
//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!

// use std::marker::Copy;
// use std::convert::{TryInto, TryFrom};
use std::fmt::Debug;
use std::collections::{Bound, BTreeMap, HashMap};
use std::option::Option;

extern crate num_traits;
use num_traits::Bounded;
// use std::ops::{AddAssign, SubAssign};

// use std::arch::x86_64::{
//     __m256i,
//     _mm256_load_si256,
// };


const DEFAULT_SPARSITY: f64 = 2.0;
// const MIN_FINAL_SEQ_LEN: usize = 16;
const MIN_FINAL_SEQ_LEN: usize = 1;


// TODO: Let's not use i64 by default. Use i32 

#[derive(Debug)]
pub struct Interval<I, T> {
    pub first: I,
    pub last: I,
    pub metadata: T
}

type PlainInterval<I> = Interval<I, ()>;


/// COITree data structure. A representation of a static set of intervals with
/// associated metadata, enabling fast overlap and coverage queries.
pub struct COITree<I, T> {

    // index into `intervals` according to query `last`
    // TODO: since this is static may be worth switching to a sorted list of
    // pairs and just doing interpolation search.
    index: BTreeMap<I, usize>,

    // intervals arranged to facilitate linear searches
    intervals: Vec<PlainInterval<I>>,

    // metadata associated with each interval in intervals
    metadata: Vec<T>
}


#[derive(Copy, Clone, Debug)]
pub struct SurplusTreeNode {
    // how many more points could lie below the y boundary before
    // the query region represented by this tree becomes sparse
    surplus: f64,

    // minimum surplus of any prefix in the subtree rooted at this node.
    // this is called "tolerance" in the Arge paper
    min_prefix_surplus: f64,

    // number of leaf nodes in the subtree not marked as dead
    live_nodes: usize
}


impl Default for SurplusTreeNode {
    fn default() -> Self {
        return SurplusTreeNode {
            surplus: f64::NAN,
            min_prefix_surplus: f64::NAN,
            live_nodes: 0,
        };
    }
}


/// Data structure used when constructing the COITree to keep track of the
/// "surplus" of each prefix which let's us find for the current y boundary
/// whether there are any sparse queries and where the maximum x value is
/// for that query.
pub struct SurplusTree {
    sparsity: f64,

    // these are stored in binary heap order
    nodes: Vec<SurplusTreeNode>,

    // map (x-sorted) leaf index to the corresponding (y-sorted) index
    // in the `intervals` vector. We additionally delete points when they
    // become dead to facilitate enumerating L_i

    // TODO: I don't think this needs to be ordered. Use a a hash table.
    leaf_to_index: HashMap<usize, usize>,

    // reverse direction: map node index to leaf node index
    index_to_leaf: Vec<usize>,

    // TODO: don't think we need this because we can just check leaf_to_node

    // for each leaf, true if the point is not dead
    // (a point gets killed if its in L_i - L'_i)
    // We could use the bitmap crate to make this a little faster
    // alive: Vec<bool>
}


impl SurplusTree where {
    fn new<I, T>(intervals: &Vec<Interval<I, T>>, sparsity: f64) -> Self
            where I: Ord + Copy {
        let n = intervals.len();

        // permutation that puts (y-sorted) nodes in x-sorted order.
        let mut xperm: Vec<usize> = (0..n).collect();
        xperm.sort_unstable_by_key(|i| intervals[*i].first);

        let mut nodes = vec![SurplusTreeNode::default(); 2*n-1];
        let mut leaf_to_index: HashMap<usize, usize> = HashMap::with_capacity(n);
        let mut index_to_leaf: Vec<usize> = vec![usize::MAX; n];

        // basically doing a manual closure here so I can do recursion
        struct TraverseContext<'a, I, T> {
            intervals: &'a Vec<Interval<I, T>>,
            nodes: &'a mut Vec<SurplusTreeNode>,
            leaf_to_index: &'a mut HashMap<usize, usize>,
            xperm: &'a Vec<usize>,
            sparsity: f64
        }

        // traverse the implicit tree initializing `nodes` and `index_to_leaf`.
        fn init_traverse<I, T>(
                ctx: &mut TraverseContext<I, T>,
                i: usize, leaf_count: &mut usize) -> usize where I: Ord {
            let n = ctx.nodes.len();
            if i >= n {
                return 0;
            }

            let left = 2*i+1;
            let right = 2*i+2;

            // this is a leaf node
            let mut subtree_leaf_count = 0;

            if left >= n && right >= n {
                subtree_leaf_count += 1;
                let idx = ctx.xperm[*leaf_count];
                ctx.leaf_to_index.insert(i, idx);
                *leaf_count += 1;
            } else {
                if left < n {
                    subtree_leaf_count += init_traverse(ctx, left, leaf_count);
                }
                if right < n {
                    subtree_leaf_count += init_traverse(ctx, right, leaf_count);
                }
            }

            ctx.nodes[i].surplus = (ctx.sparsity - 1.0) * (subtree_leaf_count as f64);
            ctx.nodes[i].min_prefix_surplus = ctx.sparsity - 1.0;
            ctx.nodes[i].live_nodes = subtree_leaf_count;

            return subtree_leaf_count;
        };

        let mut leaf_count = 0;
        init_traverse(
            &mut TraverseContext{
                intervals: intervals,
                nodes: &mut nodes,
                leaf_to_index: &mut leaf_to_index,
                xperm: &mut xperm,
                sparsity: sparsity },
            0, &mut leaf_count);

        // reverse the map
        for (k, v) in &leaf_to_index {
            index_to_leaf[*v] = *k;
        }

        return Self {
            sparsity: sparsity,
            nodes: nodes,
            leaf_to_index: leaf_to_index,
            index_to_leaf: index_to_leaf,
        };
    }

    fn len(&self) -> usize {
        return self.nodes.len();
    }

    fn has_sparse_prefix(&self) -> bool {
        return self.nodes[0].surplus < 0.0;
    }

    fn update_min_prefix_surplus(&mut self, j: usize) {
        let left = 2*j+1;
        let right = 2*j+2;
        let n = self.len();
        debug_assert!((left < n) == (right < n));

        if left < n {
            let left_mps = self.nodes[left].min_prefix_surplus;
            let left_surp = self.nodes[left].surplus;
            self.nodes[j].min_prefix_surplus = left_mps;

            let right_mps = self.nodes[right].min_prefix_surplus;
            if right_mps + left_surp < left_mps {
                self.nodes[j].min_prefix_surplus = right_mps + left_surp;
            }
        }
        // else {
            // leaf node, so min_prefix_surplus = surplus
            // if self.nodes[j].live_nodes > 0 {
            //     self.nodes[j].min_prefix_surplus = self.nodes[j].surplus;
            // } else {
                // self.nodes[j].min_prefix_surplus = f64::INFINITY;
            // }
        }
    // }

    // Called on nodes below the sweep line when they are also to the left of
    // the maximum sparse x query point (and thus no longer in S_i)
    fn set_node_dead(&mut self, i: usize) {
        let mut j = self.index_to_leaf[i];

        // disregard any prefix that ends on this node
        self.nodes[j].min_prefix_surplus = f64::INFINITY;

        loop {
            self.nodes[j].surplus += 1.0;
            self.nodes[j].live_nodes -= 1;
            self.update_min_prefix_surplus(j);

            if j > 0 {
                j = (j-1)/2; // parent node
            } else {
                break;
            }
        }
    }

    // Called on nodes when they fall below the sweep line
    fn set_node_useless(&mut self, i: usize) {
        let mut j = self.index_to_leaf[i];

        // disregard any prefix that ends on this node
        self.nodes[j].min_prefix_surplus = f64::INFINITY;

        loop {
            self.nodes[j].surplus -= self.sparsity;
            self.update_min_prefix_surplus(j);

            if j > 0 {
                j = (j-1)/2; // parent node
            } else {
                break;
            }
        }
    }


    fn num_live_nodes(&self) -> usize {
        return self.nodes[0].live_nodes;
    }


    fn last_live_leaf(&self) -> Option<usize> {
        let mut i = 0;

        if self.nodes[i].live_nodes == 0 {
            return None;
        }

        let num_nodes = self.len();
        let mut left = 2*i+1;
        let mut right = 2*i+2;

        while left < num_nodes {
            if self.nodes[right].live_nodes > 0 {
                i = right;
            } else {
                i = left;
            }

            left = 2*i+1;
            right = 2*i+2;
        }

        return Some(i);
    }


    fn next_live_leaf(&self, mut i: usize) -> Option<usize> {
        let num_nodes = self.len();

        let mut left = 2*i+1;
        let mut right = 2*i+2;

        if left < num_nodes {
            // internal node: climb down until we find a live leaf
            if self.nodes[i].live_nodes == 0 {
                return None;
            }

            while left < num_nodes {
                if self.nodes[left].live_nodes > 0 {
                    i = left;
                } else {
                    i = right;
                }

                left = 2*i+1;
                right = 2*i+2;
            }
        } else {
            // leaf node: climb up until we are someone's left child and right
            // is live

            loop {
                if i == 0 {
                    return None;
                }

                let parent = (i-1)/2;
                let parent_left = 2*parent+1;
                let parent_right = 2*parent+2;

                debug_assert!(i == parent_left || i == parent_right);
                if i == parent_left && self.nodes[parent_right].live_nodes > 0 {
                    i = parent_right;
                    break;
                }

                i = parent;
            }

            // now climb down and find a live node
            left = 2*i+1;
            right = 2*i+2;
            while left < num_nodes {
                if self.nodes[left].live_nodes > 0 {
                    i = left;
                } else {
                    i = right;
                }

                left = 2*i+1;
                right = 2*i+2;
            }
        }

        return Some(i);
    }

    // Iterate through a non-dead x-sorted prefix. Stop when we hit the leaf
    // node with the give index. The iterator returns indexes into the intervals
    // array
    fn prefix(&mut self, end_idx: usize) -> SurplusTreePrefixIterator {
        return SurplusTreePrefixIterator {
            i: self.next_live_leaf(0),
            end_idx: end_idx,
            tree: self
        }
    }


    fn map_prefix<F>(&mut self, end_idx: usize, mut visit: F)
            where F: FnMut(&mut Self, usize) {

        let mut i = 0;
        loop {
            if let Some(j) = self.next_live_leaf(i) {
                visit(self, self.leaf_to_index[&j]);
                if j == end_idx {
                    break;
                }
                i = j;
            } else {
                break;
            }
        }
    }


    // If there is a sparse query, return the index with the corresponding
    // maximum x boundary..
    // If there is no sparse query, return None.
    fn find_sparse_query_prefix(&self) -> Option<usize> {
        if self.nodes[0].min_prefix_surplus >= 0.0 {
            return None;
        } else {
            let n = self.len();
            let mut j = 0;
            let mut prefix_surplus = 0.0;
            loop {
                let left = 2*j+1;
                let right = 2*j+2;
                debug_assert!((left < n) == (right < n));

                if left < n {
                    let right_prefix_surplus =
                        prefix_surplus +
                        self.nodes[left].surplus +
                        self.nodes[right].min_prefix_surplus;

                    if right_prefix_surplus < 0.0 {
                        prefix_surplus += self.nodes[left].surplus;
                        j = right;
                    } else {
                        j = left;
                    }
                } else {
                    prefix_surplus += self.nodes[j].surplus;
                    break;
                }
            }

            assert!(prefix_surplus < 0.0);
            return Some(j);
        }
    }
}

// TODO: Why even bother making this an iterator. Why not just use a callback?

struct SurplusTreePrefixIterator<'a> {
    // current leaf node index
    i: Option<usize>,
    end_idx: usize,
    tree: &'a mut SurplusTree
}

impl<'a> Iterator for SurplusTreePrefixIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.i {
            let j = self.tree.leaf_to_index[&i];
            if i == self.end_idx {
                self.i = None;
            } else {
                self.i = self.tree.next_live_leaf(i);
            }

            return Some(j);
        } else {
            return None;
        }
    }
}


impl<I, T> COITree<I, T> {
    pub fn new(intervals: Vec<Interval<I, T>>) -> COITree<I, T>
            where I: Bounded + Ord + Copy + Debug, T: Copy {
        return Self::with_sparsity(intervals, DEFAULT_SPARSITY);
    }

    // TODO: is intervals getting copied here? Should I be using a mutable ref?
    pub fn with_sparsity(mut intervals: Vec<Interval<I, T>>, sparsity: f64) -> COITree<I, T>
            where I: Bounded + Ord + Copy + Debug, T: Copy {
        assert!(sparsity > 1.0);

        let n = intervals.len();

        if n == 0 {
            return Self{
                index: BTreeMap::new(),
                intervals: Vec::new(),
                metadata: Vec::new()
            }
        }

        let max_size: usize = (sparsity * (n as f64)).ceil() as usize;
        dbg!(max_size);

        let mut searchable_intervals: Vec<PlainInterval<I>> = Vec::with_capacity(max_size);
        let mut metadata: Vec<T> = Vec::with_capacity(max_size);
        let mut index: BTreeMap<I, usize> = BTreeMap::new();

        index.insert(I::min_value(), 0);

        intervals.sort_unstable_by_key(|i| i.last);

        let mut surplus_tree = SurplusTree::new(&intervals, sparsity);

        // for node in &surplus_tree.nodes {
        //     dbg!(node);
        // }

        // TODO: How do I write the sequences?
        // They are x-sorted, so it has to be done by walking though
        // the surplus tree, dumping all nodes that aren't marked dead.
        //
        // Output should be x-sorted, so I need to walk the surplus tree.

        let mut i = 0;
        while i < n && surplus_tree.num_live_nodes() > 0 {

            // can't place a boundary between equal values, so consume
            // everything that shares a value in one chunk
            // if i + 1 < n && intervals[i].last ==  intervals[i+1].last {
            //     surplus_tree.set_node_useless(i);
            //     i += 1;
            //     continue;
            // }

            let max_end_opt = if n - i <= MIN_FINAL_SEQ_LEN {
                i = n-1;
                surplus_tree.last_live_leaf()
            } else {
                surplus_tree.find_sparse_query_prefix()
            };

            // dbg!(i);
            // dbg!(surplus_tree.nodes[0].live_nodes);
            // dbg!(surplus_tree.nodes[0].min_prefix_surplus);

            if let Some(max_end) = max_end_opt {
                let boundary = intervals[i].last;
                // dbg!(boundary);
                // dbg!(intervals[surplus_tree.leaf_to_index[&max_end]].first);
                // dbg!(max_end);

                // for node in &surplus_tree.nodes {
                //     dbg!(node);
                // }

                let mut killed_count = 0;
                let mut l_count = 0;

                surplus_tree.map_prefix(max_end, |tree, idx| {
                    // TODO: wait, this is only in L if interval.first <= our
                    // max query.

                    let interval = &intervals[idx];
                    searchable_intervals.push(PlainInterval{
                        first: interval.first,
                        last: interval.last,
                        metadata: ()
                    });
                    metadata.push(interval.metadata);
                    l_count += 1;

                    if interval.last < boundary {
                        tree.set_node_dead(idx);
                        killed_count += 1;
                    }
                });

                // dbg!((killed_count, l_count));

                if !(i == n-1 || killed_count > l_count - killed_count) {
                    dbg!(&searchable_intervals[searchable_intervals.len()-l_count..]);
                }

                assert!(i == n-1 || killed_count > l_count - killed_count);

                if i < n-1 {
                    index.insert(boundary, searchable_intervals.len());
                }
            }

            // eprintln!("useless: {:?}", intervals[i].last);
            surplus_tree.set_node_useless(i);
            i += 1;

            // make sure the next value in distinct, otherwise we can choose
            // a boundary between equal values which can break things
            while i < n && intervals[i].last == intervals[i-1].last {
                surplus_tree.set_node_useless(i);
                i += 1;
            }
        }

        dbg!(&searchable_intervals[0..100]);
        dbg!(searchable_intervals.len());
        dbg!(&index);
        dbg!(index.len());

        // TODO: under this scheme we end copying metadata entries. That's bad
        // if metadata is large. We could consider adding an index to intervals
        // to avoid this. Or we could just store references to metadata.

        return Self{
            index: index,
            intervals: searchable_intervals,
            metadata: metadata
        };
    }


    fn find_search_start(&self, first: I) -> Option<usize> where I: Ord {
        if let Some((_, i)) = self.index.range((Bound::Unbounded, Bound::Included(first))).next_back() {
            return Some(*i);
        } else {
            return None;
        }
    }

    pub fn query_count(&self, first: I, last: I) -> usize
            where I: Bounded + Ord + Copy + Debug {

        let mut count = 0;

        // LAST is unused. WTF?

        if let Some(search_start) = self.find_search_start(first) {
            // dbg!((first, last, search_start));

            let mut last_hit_first = I::min_value();
            let mut misses = 0;
            for interval in &self.intervals[search_start..] {
                // dbg!(interval);
                if interval.first > last {
                    break;
                } else if interval.last >= first && interval.first >= last_hit_first {
                    debug_assert!(interval.first <= last);
                    count += 1;
                    last_hit_first = interval.first;
                } else {
                    misses += 1;
                }
            }

            // println!("count, misses = {}, {}", count, misses);
        }

        return count;
    }
}
