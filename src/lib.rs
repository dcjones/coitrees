
use std::time::Instant; // for debugging
use std::fmt;

// avx2 stuff
use std::arch::x86_64::{
    _lzcnt_u32,
    __m128i, __m256i, __m256,
    _mm_extract_epi32, _mm_max_epi32, _mm256_extracti128_si256,
    _mm256_set_epi32, _mm256_set1_epi32, _mm256_extract_epi32,
    _mm256_max_epi32, _mm256_cmpgt_epi32, _mm256_cmpeq_epi32, _mm256_or_si256,
    _mm256_xor_si256, _mm256_movemask_epi8 };

// number of __m256i chunks go in a leaf
const LEAF_SIZE: usize = 1;


// smallest integer >= x/y
fn cld(x: usize, y: usize) -> usize {
    return (x + (y-1)) / y;
}


// wrapper for avx integer type
#[derive(Copy, Clone)]
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
                    // xs[7], xs[6], xs[5], xs[4], xs[3], xs[2], xs[1], xs[0])
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
        // let x7 = f(0);
        // let x6 = f(1);
        // let x5 = f(2);
        // let x4 = f(3);
        // let x3 = f(4);
        // let x2 = f(5);
        // let x1 = f(6);
        // let x0 = f(7);
        return I32x8 {
            x: unsafe {
                _mm256_set_epi32(x0, x1, x2, x3, x4, x5, x6, x7)
            }
        };
    }

    // extract the first element of the vector
    fn first(&self) -> i32 {
        return unsafe {
            _mm256_extract_epi32(self.x, 7)
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


impl fmt::Display for I32x8 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        unsafe {
            return write!(
                f, "({}, {}, {}, {}, {}, {}, {}, {})",
                _mm256_extract_epi32(self.x, 7),
                _mm256_extract_epi32(self.x, 6),
                _mm256_extract_epi32(self.x, 5),
                _mm256_extract_epi32(self.x, 4),
                _mm256_extract_epi32(self.x, 3),
                _mm256_extract_epi32(self.x, 2),
                _mm256_extract_epi32(self.x, 1),
                _mm256_extract_epi32(self.x, 0));
        }
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

        // eprintln!("hits: {:032b}", hits);

        // exclude the unused bits
        let excluded = 8 - mask;
        hits >>= 4*excluded;

        // eprintln!("hits: {:032b} ({})", hits, excluded);

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
#[derive(Copy, Clone)]
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

        // copy metadata to metadata array
        let now = Instant::now();
        let mut metadata: Vec<T> = Vec::with_capacity(n);
        for i in &perm {
            metadata.push(intervals[*i as usize].metadata);
        }

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

        // figure the sizes of things
        let keys_per_leaf = LEAF_SIZE*chunk_size;
        let num_leaves = cld(n, keys_per_leaf);
        let num_bottom_layer_nodes = cld(num_leaves, node_size);

        let size_guess = (chunk_size*num_bottom_layer_nodes) / (chunk_size - 1);
        let mut nodes: Vec<Node> = Vec::with_capacity(size_guess);

        // TODO: this tree building algorithm is really suboptimal. We should
        // have 8 children in the top, but we usually don't, so those avx
        // comparisons are wasted. What's the right way to do this? Bottom down?
        // I think that must be the way to go. Should be trying to minimize the
        // number of nodes. We are always doing all 8 comparisons.

        // we could even do full trees with implicit coordinates.


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

        // make upper levels of the tree
        // eprintln!("--------------");
        let mut curr_layer = 0..nodes.len();
        while curr_layer.len() > 1 {
            // eprintln!("curr_layer.len: {}", curr_layer.len());
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

        // reorder nodes in pre-order
        let mut reorder: Vec<usize> = Vec::with_capacity(nodes.len());
        let mut reorder_rev = vec![0; nodes.len()];
        let mut stack: Vec<usize> = Vec::new();

        if !nodes.is_empty() {
            stack.push(nodes.len() - 1);
        }

        while !stack.is_empty() {
            let i = stack.pop().unwrap();
            reorder.push(i);
            reorder_rev[i] = reorder.len() - 1;
            if !nodes[i].leaf_children {
                for j in 0..nodes[i].num_children {
                    stack.push(nodes[i].child_start[j as usize] as usize);
                }
            }
        }
        assert!(reorder.len() == nodes.len());

        // fix indexes
        for node in &mut nodes {
            if !node.leaf_children {
                for j in 0..node.num_children {
                    node.child_start[j as usize] =
                        reorder_rev[node.child_start[j as usize] as usize] as u32;
                }
            }
        }

        // reorder
        let mut nodes_reorder: Vec<Node> = nodes.clone();
        for (i, node) in reorder_rev.iter().zip(nodes) {
            nodes_reorder[*i] = node;
        }

        return COITree {
            nodes: nodes_reorder,
            keys: keys,
            metadata: metadata,
        }
    }

    pub fn query(&self, first: i32, last:i32) -> (usize, usize , usize) {
        let query_first = I32x8::fill(first);
        let query_last = I32x8::fill(last);

        let mut count = 0;
        let mut overlap = 0;
        let mut visited = 0;

        self.query_node(
            0, &query_first, &query_last, &mut count, &mut overlap, &mut visited);
            // self.nodes.len()-1, &query_first, &query_last, &mut count, &mut overlap, &mut visited);

        return (count, overlap, visited);
    }

    fn query_node(
            &self, idx: usize, query_first: &I32x8, query_last: &I32x8,
            count: &mut usize, overlap: &mut usize, visited: &mut usize) {
        let node = &self.nodes[idx];

        if node.leaf_children {
            for i in interval_chunk_overlaps(
                    &query_first, &query_last,
                    &node.minfirst, &node.maxlast, node.num_children) {
                self.query_leaf(
                    node.child_start[i] as usize, node.child_len[i] as usize,
                    query_first, query_last, count, overlap, visited);
                *visited += 1;
            }
        } else {
            for i in interval_chunk_overlaps(
                    query_first, query_last,
                    &node.minfirst, &node.maxlast, node.num_children) {
                self.query_node(
                    node.child_start[i] as usize,
                    query_first, query_last, count, overlap, visited);
                *visited += 1;
            }
        }
    }

    fn query_leaf(
            &self, idx: usize, leaflen: usize, query_first: &I32x8, query_last: &I32x8,
            count: &mut usize, overlap: &mut usize, visited: &mut usize) {

        let chunk_size = I32x8::lanes();
        let num_chunks = cld(leaflen, chunk_size);
        for i in 0..num_chunks {
            let mask = leaflen.min((i+1)*chunk_size) - i*chunk_size;
            assert!(mask > 0);
            assert!(mask <= I32x8::lanes());
            for j in interval_chunk_overlaps(
                    query_first, query_last,
                    &self.keys[idx+i].firsts, &self.keys[idx+i].lasts, mask as u8) {
                *count += 1;
                // TODO: compute overlap (with avx?). It's kind of a pain
                // to extract specific values.
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inspect_pack() {
        let a = I32x8::pack([1, 2, 3, 4, 5, 6, 7, 8]);
        println!("a: {}", a);
    }

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
