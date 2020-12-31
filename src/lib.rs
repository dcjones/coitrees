
//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!

// use std::marker::Copy;
// use std::convert::{TryInto, TryFrom};
use std::fmt::Debug;
use std::collections::{Bound, BTreeMap, HashMap};
use std::option::Option;
use std::time::Instant;
use std::cmp::min;

extern crate num_traits;
use num_traits::Bounded;
// use std::ops::{AddAssign, SubAssign};

// use std::arch::x86_64::{
//     __m256i,
//     _mm256_load_si256,
// };


const DEFAULT_SPARSITY: i64 = 10;
const MIN_FINAL_SEQ_LEN: usize = 16;


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
    intervals: Vec<Interval<I, u32>>,

    // metadata associated with each interval in intervals
    metadata: Vec<T>
}


#[derive(Copy, Clone, Debug)]
pub struct SurplusTreeNode {
    // how many more points could lie below the y boundary before
    // the query region represented by this tree becomes sparse
    surplus: i64,

    // minimum surplus of any prefix in the subtree rooted at this node.
    // this is called "tolerance" in the Arge paper
    min_prefix_surplus: i64,

    // number of leaf nodes in the subtree not marked as dead
    live_nodes: usize
}


impl Default for SurplusTreeNode {
    fn default() -> Self {
        return SurplusTreeNode {
            surplus: i64::MIN,
            min_prefix_surplus: i64::MIN,
            live_nodes: 0,
        };
    }
}


/// Data structure used when constructing the COITree to keep track of the
/// "surplus" of each prefix which let's us find for the current y boundary
/// whether there are any sparse queries and where the maximum x value is
/// for that query.
pub struct SurplusTree {
    sparsity: i64,

    // these are stored in binary heap order
    nodes: Vec<SurplusTreeNode>,

    // map (x-sorted) leaf index to the corresponding (y-sorted) index in the
    // `intervals` vector, along with its leaf number. We additionally delete
    // points when they become dead to facilitate enumerating L_i
    leaf_to_index: HashMap<usize, (usize, u32)>,

    // reverse direction: map node index to leaf node index
    index_to_leaf: Vec<usize>,
}


// // Iterates over ancestors of a given nodes up to the root, along with
// // each nodes immediate children
// struct AncestorIterator<'a> {
//     nodes: &'a mut Vec<SurplusTreeNode>,
//     root: Option<usize>
// }


// impl<'a, 'b> Iterator for AncestorIterator<'a> {
//     type Item = (&'b mut SurplusTreeNode, &'b SurplusTreeNode, &'b SurplusTreeNode);

//     fn next(&mut self) -> Option<Self::Item> {
//         if let Some(root) = self.root {
//             if root == 0 {
//                 self.root = None;
//             } else {
//                 self.root = Some((root-1)/2);
//             }

//             let left = 2*root+1;
//             let right = 2*root+2;

//             let ret: Self::Item = unsafe {
//                 (
//                     self.nodes.get_unchecked_mut(root),
//                     self.nodes.get_unchecked(left),
//                     self.nodes.get_unchecked(right)
//                 )
//             };

//             return Some(ret);

//         } else {
//             return None;
//         }
//     }
// }

// Function for traversing from one leaf node index to the next.
fn next_live_leaf(nodes: &Vec<SurplusTreeNode>, mut i: usize) -> Option<usize> {
    let num_nodes = nodes.len();

    let mut left = 2*i+1;
    let mut right = 2*i+2;

    if left < num_nodes {
        // internal node: climb down until we find a live leaf
        if nodes[i].live_nodes == 0 {
            return None;
        }

        while left < num_nodes {
            if nodes[left].live_nodes > 0 {
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
            if i == parent_left && nodes[parent_right].live_nodes > 0 {
                i = parent_right;
                break;
            }

            i = parent;
        }

        // now climb down and find a live node
        left = 2*i+1;
        right = 2*i+2;
        while left < num_nodes {
            if nodes[left].live_nodes > 0 {
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


impl SurplusTree where {
    fn new<I, T>(intervals: &Vec<Interval<I, T>>, sparsity: i64) -> Self
            where I: Ord + Copy {
        let n = intervals.len();
        assert!(sparsity > 1);

        // permutation that puts (y-sorted) nodes in x-sorted order.
        let mut xperm: Vec<usize> = (0..n).collect();

        // OPT: most expensive part of this function
        xperm.sort_unstable_by_key(|i| intervals[*i].first);
        // xperm.sort_unstable_by_key(|i| unsafe { intervals.get_unchecked(*i) }.first);

        let num_nodes = 2*n-1;
        let mut nodes = vec![SurplusTreeNode::default(); num_nodes];
        let mut leaf_to_index: HashMap<usize, (usize, u32)> = HashMap::with_capacity(n);
        let mut index_to_leaf: Vec<usize> = vec![usize::MAX; n];

        // go up the tree setting internal node values

        // this version won't work because we have to borrow multiple times
        // for (i_rev, node) in nodes.iter_mut().rev().enumerate() {
        //     let i = num_nodes - i_rev - 1;
        //     let left = 2*i+1;
        //     let right = 2*i+2;
        //     if left < num_nodes {
        //         node.live_nodes = nodes[left].live_nodes + nodes[right].live_nodes;
        //     } else {
        //         node.live_nodes = 1;
        //     }

        //     node.min_prefix_surplus = sparsity - 1.0;
        //     node.surplus = (sparsity - 1.0) * node.live_nodes as f64;
        // }

        let mut i = nodes.len() - 1;
        loop {
            let left = 2*i+1;
            let right = 2*i+2;

            // is a leaf node
            if left < nodes.len() {
                nodes[i].live_nodes = nodes[left].live_nodes + nodes[right].live_nodes;
            } else {
                nodes[i].live_nodes = 1;
            }

            nodes[i].min_prefix_surplus = sparsity - 1;
            nodes[i].surplus = (sparsity - 1) * nodes[i].live_nodes as i64;

            if i == 0 {
                break;
            } else {
                i -= 1;
            }
        }

        // iterate through leaves, building `leaf_to_index`, and initializing
        // leaf node values
        let mut i = 0;
        let mut leaf_count = 0;
        while let Some(j) = next_live_leaf(&nodes, i) {
            let idx = xperm[leaf_count];

            // OPT: second most expensive part of this function
            leaf_to_index.insert(j, (idx, leaf_count as u32 + 1));
            leaf_count += 1;
            i = j;
        }
        assert!(leaf_to_index.len() == n);

        // reverse the map
        for (k, (v, _)) in &leaf_to_index {
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


    // Climb up to the root from node `j`, setting each ancestors `min_prefix_surplus`
    // and calling `visit` on each ancestor. This needs to be as fast as possible,
    // so we do unsafe indexing.
    fn update_ancestors<F>(&mut self, j: usize, mut visit: F)
            where F: FnMut(&mut SurplusTreeNode) {

        assert!(j < self.nodes.len());
        let mut root = j;
        while root > 0 {
            root = (root-1)/2;

            let left = 2*root+1;
            let right = 2*root+2;

            let (left_min_prefix_surplus, left_surplus) = {
                let node_left = unsafe { self.nodes.get_unchecked(left) };
                (node_left.min_prefix_surplus, node_left.surplus)
            };

            let right_min_prefix_surplus = {
                let node_right = unsafe { self.nodes.get_unchecked(right) };
                node_right.min_prefix_surplus
            };

            let node_root = unsafe { self.nodes.get_unchecked_mut(root) };

            // node_root.min_prefix_surplus = left_min_prefix_surplus.min(
            //     left_surplus + right_min_prefix_surplus);

            let l = left_min_prefix_surplus;
            let r = left_surplus.saturating_add(right_min_prefix_surplus);
            node_root.min_prefix_surplus = min(l, r);

            visit(node_root);
        }
    }

    // Called on nodes below the sweep line when they are also to the left of
    // the maximum sparse x query point (and thus no longer in S_i)
    fn set_node_dead(&mut self, i: usize) {
        let j = self.index_to_leaf[i];

        {
            let node = &mut self.nodes[j];
            // disregard any prefix that ends on this node
            node.min_prefix_surplus = i64::MAX;
            node.surplus += 1;
            node.live_nodes -= 1;

        }
        self.update_ancestors(j, |node| {
            node.surplus += 1;
            node.live_nodes -= 1;
        });
    }

    // Called on nodes when they fall below the sweep line
    fn set_node_useless(&mut self, i: usize) {
        let j = self.index_to_leaf[i];

        {
            let node = &mut self.nodes[j];
            // disregard any prefix that ends on this node
            node.min_prefix_surplus = i64::MAX;
            node.surplus -= self.sparsity;
        }

        let sparsity = self.sparsity;
        self.update_ancestors(j, |node| {
            node.surplus -= sparsity;
        });
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


    fn map_prefix<F>(&mut self, end_idx: usize, mut visit: F)
            where F: FnMut(&mut Self, usize, u32) {

        let mut i = 0;
        loop {
            if let Some(j) = next_live_leaf(&self.nodes, i) {
                let idx = self.leaf_to_index[&j];
                visit(self, idx.0, idx.1);
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
        if self.nodes[0].min_prefix_surplus >= 0 {
            return None;
        } else {
            let n = self.len();
            let mut j = 0;
            let mut prefix_surplus: i64 = 0;
            loop {
                let left = 2*j+1;
                let right = 2*j+2;
                debug_assert!((left < n) == (right < n));

                // TODO: Since we are setting min_prefix_surplus to MAX
                // this causes the following addition to overflow, fucking
                // thing up.

                if left < n {
                    let right_prefix_surplus =
                        prefix_surplus
                            .saturating_add(self.nodes[left].surplus)
                            .saturating_add(self.nodes[right].min_prefix_surplus);

                    if right_prefix_surplus < 0 {
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

            // TODO: This doesn't really resolve the issue. Any point in keeping it?
            // find the furthest leaf node that shares the same `first` value
            // this is to avoid splitting up points that share the same values,
            // which leads to issues.
            // while let Some(k) = self.next_live_leaf(j) {
            //     if intervals[self.leaf_to_index[&j].0].first == intervals[self.leaf_to_index[&k].0].first {
            //         j = k;
            //     } else {
            //         break;
            //     }
            // }

            assert!(prefix_surplus < 0);
            return Some(j);
        }
    }
}


impl<I, T> COITree<I, T> {
    pub fn new(intervals: Vec<Interval<I, T>>) -> COITree<I, T>
            where I: Bounded + Ord + Copy + Debug, T: Copy {
        return Self::with_sparsity(intervals, DEFAULT_SPARSITY);
    }

    // TODO: is intervals getting copied here? Should I be using a mutable ref?
    pub fn with_sparsity(mut intervals: Vec<Interval<I, T>>, sparsity: i64) -> COITree<I, T>
            where I: Bounded + Ord + Copy + Debug, T: Copy {
        assert!(sparsity > 1);

        let n = intervals.len();

        if n == 0 {
            return Self{
                index: BTreeMap::new(),
                intervals: Vec::new(),
                metadata: Vec::new()
            }
        }

        let max_size: usize = ((sparsity as f64 / (sparsity - 1) as f64) * (n as f64)).ceil() as usize;
        dbg!(max_size);

        let mut searchable_intervals: Vec<Interval<I, u32>> = Vec::with_capacity(max_size);
        let mut metadata: Vec<T> = Vec::with_capacity(max_size);
        let mut index: BTreeMap<I, usize> = BTreeMap::new();

        index.insert(I::min_value(), 0);

        intervals.sort_unstable_by_key(|i| i.last);

        let now = Instant::now();
        let mut surplus_tree = SurplusTree::new(&intervals, sparsity);
        eprintln!("SurplusTree::new: {}", now.elapsed().as_millis() as f64 / 1000.0);

        // searchable intervals stores an extra integer to disambiguate and
        // avoid counting the same hits more than once.

        // TODO: NO THIS DOESN'T DISAMBIGUATE AT ALL
        // We need to assign these ids to leaves in the surplus tree
        // Ok, where do we assign them then?
        // Can `map_prefix` keep track of the leaf number?

        let now = Instant::now();
        let mut i = 0;
        while i < n && surplus_tree.num_live_nodes() > 0 {

            let max_end_opt = if n - i <= MIN_FINAL_SEQ_LEN {
                i = n-1;
                surplus_tree.last_live_leaf()
            } else {
                surplus_tree.find_sparse_query_prefix()
            };

            if let Some(max_end) = max_end_opt {
                let boundary = intervals[i].last;

                let mut killed_count = 0;
                let mut l_count = 0;

                surplus_tree.map_prefix(max_end, |tree, idx, leaf_num| {
                    let interval = &intervals[idx];
                    searchable_intervals.push(Interval{
                        first: interval.first,
                        last: interval.last,
                        metadata: leaf_num,
                    });

                    metadata.push(interval.metadata);
                    l_count += 1;

                    if interval.last < boundary {
                        tree.set_node_dead(idx);
                        killed_count += 1;
                    }
                });

                if i < n-1 {
                    index.insert(boundary, searchable_intervals.len());
                }
            }

            surplus_tree.set_node_useless(i);
            i += 1;

            // make sure the next value in distinct, otherwise we can choose
            // a boundary between equal values which can break things
            while i < n && intervals[i].last == intervals[i-1].last {
                surplus_tree.set_node_useless(i);
                i += 1;
            }
        }
        eprintln!("SurplusTree::new: {}", now.elapsed().as_millis() as f64 / 1000.0);

        dbg!(searchable_intervals.len());
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


    #[inline(always)]
    fn find_search_start(&self, first: I) -> Option<usize> where I: Ord {
        if let Some((_, i)) = self.index.range((Bound::Unbounded, Bound::Included(first))).next_back() {
            return Some(*i);
        } else {
            return None;
        }
    }

    pub fn query_count(&self, first: I, last: I) -> usize
            where I: Bounded + Ord + Copy + Debug {

        // let mut misses = 0;
        let mut count = 0;
        if let Some(search_start) = self.find_search_start(first) {
            let mut last_hit_id: u32 = 0;
            for interval in &self.intervals[search_start..] {
                if interval.first > last {
                    break;
                } else if interval.last >= first && interval.metadata > last_hit_id {
                    // debug_assert!(interval.first <= last);
                    count += 1;
                    last_hit_id = interval.metadata;
                }
                // } else {
                    // misses += 1;
                // }
            }
        }

        // eprintln!("hit prop: {}", count as f64 / (count + misses) as f64);

        return count;
    }
}
