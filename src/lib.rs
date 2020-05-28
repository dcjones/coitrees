
use std::time::Instant; // for debugging

// avx2 stuff
use std::arch::x86_64::{
    _lzcnt_u32,
    __m128i, __m256i, __m256,
    _mm_extract_epi32, _mm_max_epi32, _mm256_extracti128_si256,
    _mm256_set_epi32, _mm256_set1_epi32, _mm256_extract_epi32,
    _mm256_max_epi32, _mm256_cmpgt_epi32, _mm256_cmpeq_epi32, _mm256_or_si256,
    _mm256_xor_si256, _mm256_movemask_epi8 };

// number of __m256i chunks go in a leaf
const LEAF_SIZE: usize = 4;


// smallest integer >= x/y
fn cld(x: usize, y: usize) -> usize {
    return (x + (y-1)) / y;
}


// wrapper for avx integer type
struct I32x8 {
    x: __m256i,
}

impl I32x8 {
    // number of entries
    fn lanes() -> usize {
        return 8;
    }

    // new vector filled with a value
    fn fill(x: i32) -> I32x8 {
        return I32x8{
            x: unsafe { _mm256_set1_epi32(x) }
        };
    }

    // pack a slice into a vector
    fn pack(xs: [i32; 8]) -> I32x8 {
        return I32x8 {
            x: unsafe {
                _mm256_set_epi32(
                    xs[0], xs[1], xs[2], xs[3], xs[4], xs[5], xs[6], xs[7])
            }
        };
    }

    // build a vector by calling f(0), ..., f(7)
    fn build(f: impl Fn(usize) -> i32) -> I32x8 {
        let x0 = f(0);
        let x1 = f(1);
        let x2 = f(2);
        let x3 = f(3);
        let x4 = f(4);
        let x5 = f(5);
        let x6 = f(6);
        let x7 = f(7);
        return I32x8 {
            x: unsafe {
                _mm256_set_epi32(x0, x1, x2, x3, x4, x5, x6, x7)
            }
        };
    }

    // extract the first element of the vector
    fn first(&self) -> i32 {
        return unsafe {
            _mm256_extract_epi32(self.x, 0)
        };
    }

    // horizontal max
    fn hmax(&self) -> i32 {
        unsafe {
            let a: __m128i = _mm256_extracti128_si256(self.x, 0);
            let b: __m128i = _mm256_extracti128_si256(self.x, 1);
            let mab = _mm_max_epi32(a, b);

            // maybe I could permute and to avx max gain, but diminishing returns...
            let mab0 = _mm_extract_epi32(mab, 0);
            let mab1 = _mm_extract_epi32(mab, 1);
            let mab2 = _mm_extract_epi32(mab, 2);
            let mab3 = _mm_extract_epi32(mab, 3);

            return mab0.max(mab1.max(mab2.max(mab3)));
        }
    }

    // vertilal max
    fn vmax(&self, w: &I32x8) -> I32x8 {
        return I32x8{
            x: unsafe {
                _mm256_max_epi32(self.x, w.x)
            }
        };
    }
}


// check for overlaps between 2 pairs of 8 of intervals stored in vectors,
// and iterate over indexes that overlap
struct IntervalChunkOverlaps {
    // encoded where bit 4*i is set if interval pair i overlaps
    hits: u32,

    // number of low order entries to ignore
    excluded: u32,
    shift: u32,
}


fn interval_chunk_overlaps(
        first_a: &I32x8, last_a: &I32x8,
        first_b: &I32x8, last_b: &I32x8,
        mask: u8) -> IntervalChunkOverlaps {
    unsafe {
        // test if the intervals do not overlap
        let nothits256 = _mm256_or_si256(
            _mm256_cmpgt_epi32(first_a.x, last_b.x),
            _mm256_cmpgt_epi32(first_b.x, last_a.x));

        // flip all the bits
        let ones256 = _mm256_cmpeq_epi32(nothits256, nothits256);
        let hits256 = _mm256_xor_si256(nothits256, ones256);

        // compact
        let mut hits = _mm256_movemask_epi8(hits256) as u32;

        // only need one bit per hit
        hits &= 0b10001000100010001000100010001000;

        // exclude the unused bits
        let excluded = 8 - mask;
        hits >>= 4*excluded;

        eprintln!("hits: {:032b}", hits);

        // compact into an u32
        return IntervalChunkOverlaps{hits: hits, excluded: excluded as u32, shift: 0};
    }
}


impl Iterator for IntervalChunkOverlaps {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.hits == 0 {
            return None
        } else {
            let lz = unsafe {
                _lzcnt_u32(self.hits)
            };

            self.hits <<= lz + 1; // done with these bits
            self.shift += (lz + 1)/4;
            let i = self.shift - self.excluded;

            return Some(i as usize);
        }
    }
}


// find the maximum interval end in the slice of keys
fn key_max_last(keys: &[IntervalChunk]) -> i32 {
    let mut accum = I32x8::fill(0);
    for key in keys {
        accum = accum.vmax(&key.lasts);
    }
    return accum.hmax();
}


// simple interval type. Will likely need to revisit this.
pub struct Interval<T> where T: std::marker::Copy {
    pub first: i32,
    pub last: i32,
    pub metadata: T,
}

// a chunk of 8 intervals
struct IntervalChunk {
    firsts: I32x8,
    lasts: I32x8,
}

// B+-tree internal node type
struct Node {
    // Each internal node has 8 children,
    // with 32-bit keys, stored in a 256bit avx2 register.
    minfirst: I32x8,
    maxlast: I32x8,

    // number of children in [0, 8]
    num_children: u8,

    // true iff children are leaves
    leaf_children: bool,

    // if the child node is an internal node, gives the position
    // in the node vector. If the child is a leaf, gives the position
    // of the first interval in the key vector
    child_start: [u32; 8],

    // if child is a leaf node, this gives the number of intervals
    // otherwise is set to 0
    // TODO: maybe we can even shrink this down to like u8
    // gotta fiddle to see what the optimal size is.
    child_len: [u32; 8],
}


pub struct COITree<T> where T: std::marker::Copy {
    // all internal tree B+-tree nodes
    nodes: Vec<Node>,

    // all intervals stored in sorted order
    keys: Vec<IntervalChunk>,
    metadata: Vec<T>,
}

impl<T> COITree<T> where T: std::marker::Copy {
    pub fn new(intervals: Vec<Interval<T>>) -> COITree<T> {
        let n = intervals.len();
        // would like to avoid constants, but size_of is not implemented
        // let chunk_size = __m256i::size_of() / 4;
        let chunk_size = I32x8::lanes();
        let node_size = chunk_size; // just to clarify some code
        let num_chunks = cld(n, chunk_size);

        let now = Instant::now();

        // proxy sort intervals
        // (we not assume metadata is small, so avoid shuffling it around)
        let mut perm: Vec<u32> = (0..(n as u32)).collect();
        perm.sort_unstable_by_key(|i| intervals[*i as usize].first);
        eprintln!("sorting bed: {}s", now.elapsed().as_millis() as f64 / 1000.0);

        // copy metadata to metadata array
        let now = Instant::now();
        let mut metadata: Vec<T> = Vec::with_capacity(n);
        for i in &perm {
            metadata.push(intervals[*i as usize].metadata);
        }
        eprintln!("copying metadata: {}s", now.elapsed().as_millis() as f64 / 1000.0);

        // copy intervals to intervals array, packing into chunks
        let now = Instant::now();
        let mut keys: Vec<IntervalChunk> = Vec::with_capacity(num_chunks);
        for i in 0..(n/chunk_size) {
            let k = i * chunk_size;
            keys.push(IntervalChunk{
                firsts: I32x8::build(|i| intervals[perm[k+i] as usize].first),
                lasts: I32x8::build(|i| intervals[perm[k+i] as usize].last)
            });
        }

        // handle remainder
        let rem = n % chunk_size;
        if rem > 0 {
            let k = (n/chunk_size) * chunk_size;
            let firsts = I32x8::build(
                |i| if i < rem {intervals[perm[k+i] as usize].first} else {0});
            let lasts = I32x8::build(
                |i| if i < rem {intervals[perm[k+i] as usize].last} else {0});
            keys.push(IntervalChunk{firsts: firsts, lasts:lasts});

        }
        assert!(keys.len() == num_chunks);
        eprintln!("copying keys: {}s", now.elapsed().as_millis() as f64 / 1000.0);

        // figure the sizes of things
        let keys_per_leaf = LEAF_SIZE*chunk_size;
        let num_leaves = cld(n, keys_per_leaf);
        eprintln!("n: {}", n);
        eprintln!("keys_per_leaf: {}", keys_per_leaf);
        eprintln!("num_leaves: {}", num_leaves);
        let num_bottom_layer_nodes = cld(num_leaves, node_size);

        let size_guess = (chunk_size*num_bottom_layer_nodes) / (chunk_size - 1);
        let mut nodes: Vec<Node> = Vec::with_capacity(size_guess);

        // divvy up entrys into the bottom layer of internal nodes
        let now = Instant::now();
        for i in 0..num_bottom_layer_nodes {
            let first_key_chunk = i*node_size*LEAF_SIZE;
            assert!(first_key_chunk < num_chunks);

            let last_key_chunk = num_chunks.min((i+1)*node_size*LEAF_SIZE);

            let num_children = cld(last_key_chunk - first_key_chunk, LEAF_SIZE);
            assert!(num_children > 0);
            assert!(num_children <= node_size);

            // Ok, we want to pack keys[first_key_chunk..last_key_chunk]
            // into at most chunk_size internal nodes.
            let mut minfirst: [i32; 8] = [0; 8];
            let mut maxlast: [i32; 8] = [0; 8];
            let mut child_start: [u32; 8] = [0; 8];
            let mut child_len: [u32; 8] = [0; 8];

            for j in 0..num_children {
                let u = first_key_chunk+j*LEAF_SIZE;
                let v = num_chunks.min(first_key_chunk+(j+1)*LEAF_SIZE);

                // consider the leaf node made up of chunks u..v
                child_start[j] = u as u32;
                let firstchunk = &keys[u];
                minfirst[j] = firstchunk.firsts.first();

                child_len[j] = (n.min(v*chunk_size) - u*chunk_size) as u32;
                maxlast[j] = key_max_last(&keys[u..v]);
            }

            nodes.push(Node{
                minfirst: I32x8::pack(minfirst),
                maxlast: I32x8::pack(maxlast),
                num_children: num_children as u8,
                leaf_children: true,
                child_start: child_start,
                child_len: child_len
            });
        }
        eprintln!("building bottom layer nodes: {}s", now.elapsed().as_millis() as f64 / 1000.0);

        // make upper levels of the tree
        let mut curr_layer = 0..nodes.len();
        while curr_layer.len() > 1 {
            eprintln!("curr_layer {}..{}", curr_layer.start, curr_layer.end);
            // divvy current layer into a new layer
            let num_parents = cld(curr_layer.len(), node_size);

            for i in 0..num_parents {
                let num_children = curr_layer.len().min((i+1)*node_size) - i*node_size;
                assert!(num_children > 0);
                assert!(num_children <= node_size);

                let mut minfirst: [i32; 8] = [0; 8];
                let mut maxlast: [i32; 8] = [0; 8];
                let mut child_start: [u32; 8] = [0; 8];

                for j in 0..num_children {
                    let k = curr_layer.start + i*node_size + j;
                    child_start[j] = k as u32;
                    minfirst[j] = nodes[k].minfirst.first();
                    maxlast[j] = nodes[k].maxlast.hmax();
                }

                nodes.push(Node{
                    minfirst: I32x8::pack(minfirst),
                    maxlast: I32x8::pack(maxlast),
                    num_children: num_children as u8,
                    leaf_children: false,
                    child_start: child_start,
                    child_len: [0; 8]
                });
            }

            curr_layer = curr_layer.end..nodes.len()
        }

        // TODO: We should consider reordering nodes. Possible in veb order.
        // current buttome-up order is probably pretty bad.

        return COITree {
            nodes: nodes,
            keys: keys,
            metadata: metadata,
        }
    }

    pub fn query(&self, first: i32, last:i32) -> (usize, usize , usize) {
        let firstv = I32x8::fill(first);
        let lastv = I32x8::fill(last);

        let mut count = 0;
        let mut overlap = 0;
        let mut visited = 0;

        self.query_node(
            self.nodes.len()-1, firstv, lastv, &mut count, &mut overlap, &mut visited);

        return (count, overlap, visited);
    }

    fn query_node(
            &self, idx: usize, first: I32x8, last: I32x8,
            count: &mut usize, overlap: &mut usize, visited: &mut usize) {
        let node = &self.nodes[idx];

        // TODO: do interval comparison avx stuff

        // TODO: iterate through the hits somehow, traversing into children

        // TODO: call either query_node or query_leaf for each hit
    }

    fn query_leaf(
            &self, idx: usize, leaflen: usize, first: I32x8, last: I32x8,
            count: &mut usize, overlap: &mut usize, visited: &mut usize) {

        // TODO:
        // really need to figure out how to set up some kind of mask
        // to deal with unaligned data.
        // here we can do the interval comparison stuff again

    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_interval_chunk_overlaps() {
        // play around and see if this thing works
        let first_a = I32x8::fill(5);
        let last_a = I32x8::fill(70);

        let first_b = I32x8::pack([0, 80, 20, 30, 40, 50, 60, 70]);
        let last_b = I32x8::pack([9, 99, 29, 39, 49, 59, 69, 79]);

        for i in interval_chunk_overlaps(&first_a, &last_a, &first_b, &last_b, 7) {
            println!("i: {}", i);
        }
    }
}
