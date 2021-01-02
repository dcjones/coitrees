
//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!

// use std::marker::Copy;
// use std::convert::{TryInto, TryFrom};
use std::fmt::Debug;
use std::collections::{Bound, BTreeMap};
use std::option::Option;
use std::time::Instant;
use std::cmp::{min, max};

extern crate num_traits;
use num_traits::Bounded;
// use std::ops::{AddAssign, SubAssign};

// use std::arch::x86_64::{
//     __m256i,
//     _mm256_load_si256,
// };


const DEFAULT_SPARSITY: f64 = 10.0;
// const DEFAULT_SPARSITY: f64 = 1.5;
const MIN_FINAL_SEQ_LEN: usize = 16;


// TODO: Let's not use i64 by default. Use i32 

#[derive(Debug)]
pub struct Interval<I, T> {
    pub first: I,
    pub last: I,
    pub metadata: T
}


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
pub struct SurplusTreeNode<I> {
    // how many more points could lie below the y boundary before
    // the query region represented by this tree becomes sparse
    surplus: f64,

    // minimum surplus of any prefix in the subtree rooted at this node.
    // this is called "tolerance" in the Arge paper
    min_prefix_surplus: f64,

    // number of leaf nodes in the subtree not marked as dead
    live_nodes: usize,

    // if this is a leaf node, the corresponding index into the intervals array
    // otherwise usize::MAX
    intervals_index: usize,

    // if this is is a leaf node, a running count of prior leaf nodes thats
    // used to avoid double counting intersections. If not a leaf node, u32::MAX
    leaf_count: u32,

    // Smallest first (x) value among all useful points
    min_useful_first: I,
}


impl<I> Default for SurplusTreeNode<I> where I: Bounded + Default {
    fn default() -> Self {
        return Self {
            surplus: f64::NAN,
            min_prefix_surplus: f64::NAN,
            live_nodes: 0,
            intervals_index: usize::MAX,
            leaf_count: u32::MAX,
            min_useful_first: I::max_value(),
        };
    }
}


/// Data structure used when constructing the COITree to keep track of the
/// "surplus" of each prefix which let's us find for the current y boundary
/// whether there are any sparse queries and where the maximum x value is
/// for that query.
pub struct SurplusTree<I> {
    sparsity: f64,

    // these are stored in binary heap order
    nodes: Vec<SurplusTreeNode<I>>,

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
fn next_live_leaf<I>(nodes: &Vec<SurplusTreeNode<I>>, mut i: usize) -> Option<usize> {
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


impl<I> SurplusTree<I> where I: Bounded + Copy + Debug + Default + Ord {
    fn new<T>(intervals: &Vec<Interval<I, T>>, sparsity: f64) -> Self {
        let n = intervals.len();

        // permutation that puts (y-sorted) nodes in x-sorted order.
        let now = Instant::now();
        let mut xperm: Vec<usize> = (0..n).collect();
        xperm.sort_by_key(|i| intervals[*i].first);
        eprintln!("second sort: {}", now.elapsed().as_millis() as f64 / 1000.0);

        let now = Instant::now();
        let num_nodes = 2*n-1;
        let mut nodes = vec![SurplusTreeNode::default(); num_nodes];
        let mut index_to_leaf: Vec<usize> = vec![usize::default(); n];
        eprintln!("surplus tree allocate: {}", now.elapsed().as_millis() as f64 / 1000.0);

        // go up the tree setting internal node values
        let now = Instant::now();
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

            nodes[i].min_prefix_surplus = sparsity - 1.0;
            nodes[i].surplus = (sparsity - 1.0) * nodes[i].live_nodes as f64;

            if i == 0 {
                break;
            } else {
                i -= 1;
            }
        }
        eprintln!("init surplus tree 1: {}", now.elapsed().as_millis() as f64 / 1000.0);

        // iterate through leaves, initializing leaf node values
        let now = Instant::now();
        let mut i = 0;
        let mut leaf_count = 0;
        while let Some(j) = next_live_leaf(&nodes, i) {
            let idx = xperm[leaf_count];
            nodes[j].min_useful_first = intervals[idx].first;
            nodes[j].intervals_index = idx;
            nodes[j].leaf_count = leaf_count as u32 + 1;
            leaf_count += 1;
            i = j;
        }
        dbg!(leaf_count);
        eprintln!("init surplus tree 2: {}", now.elapsed().as_millis() as f64 / 1000.0);

        // set min_useful_first for internal nodes
        for i in (0..num_nodes).rev() {
            let left = 2*i+1;
            let right = 2*i+2;
            if left < num_nodes {
                nodes[i].min_useful_first = min(
                    nodes[left].min_useful_first,
                    nodes[right].min_useful_first);
            }

            assert!(nodes[i].min_useful_first < I::max_value());
        }

        // reverse the map
        let now = Instant::now();
        for (i, node) in nodes.iter().enumerate() {
            if node.intervals_index != usize::MAX {
                index_to_leaf[node.intervals_index] = i;
            }
        }
        eprintln!("init reverse index: {}", now.elapsed().as_millis() as f64 / 1000.0);

        return Self {
            sparsity: sparsity,
            nodes: nodes,
            index_to_leaf: index_to_leaf,
        };
    }

    fn len(&self) -> usize {
        return self.nodes.len();
    }


    fn min_useful_first(&self) -> I {
        return self.nodes[0].min_useful_first;
    }


    // Climb up to the root from node `j`, setting each ancestors `min_prefix_surplus`
    // and calling `visit` on each ancestor. This needs to be as fast as possible,
    // so we do unsafe indexing.
    fn update_ancestors<F>(&mut self, j: usize, mut visit: F)
            where F: FnMut(&mut SurplusTreeNode<I>) {

        assert!(j < self.nodes.len());
        let mut root = j;
        while root > 0 {
            root = (root-1)/2;

            let left = 2*root+1;
            let right = 2*root+2;

            let (left_min_prefix_surplus, left_surplus, left_min_useful_first) = {
                let node_left = unsafe { self.nodes.get_unchecked(left) };
                (node_left.min_prefix_surplus, node_left.surplus, node_left.min_useful_first)
            };

            let (right_min_prefix_surplus, right_min_useful_first) = {
                let node_right = unsafe { self.nodes.get_unchecked(right) };
                (node_right.min_prefix_surplus, node_right.min_useful_first)
            };

            let node_root = unsafe { self.nodes.get_unchecked_mut(root) };

            let l = left_min_prefix_surplus;
            let r = left_surplus + right_min_prefix_surplus;
            node_root.min_prefix_surplus = if l < r { l } else { r };

            // dbg!((root, left_min_useful_first, right_min_useful_first));

            // OPT: This isn't necessary when marking dead points
            node_root.min_useful_first = min(
                left_min_useful_first, right_min_useful_first);

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
            node.min_prefix_surplus = f64::INFINITY;
            node.surplus += 1.0;
            node.live_nodes -= 1;

        }
        self.update_ancestors(j, |node| {
            node.surplus += 1.0;
            node.live_nodes -= 1;
        });
    }

    // Called on nodes when they fall below the sweep line
    fn set_node_useless(&mut self, i: usize) {
        let j = self.index_to_leaf[i];

        {
            let node = &mut self.nodes[j];
            // disregard any prefix that ends on this node
            node.min_prefix_surplus = f64::INFINITY;
            node.surplus -= self.sparsity;
            node.min_useful_first = I::max_value();
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
                let idx = self.nodes[j].intervals_index;
                let leaf_count = self.nodes[j].leaf_count;
                visit(self, idx, leaf_count);
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

            assert!(prefix_surplus < 0.0);
            return Some(j);
        }
    }
}


impl<I, T> COITree<I, T>
        where I: Bounded + Copy + Default + Ord + Debug,
              T: Copy {
    pub fn new(intervals: Vec<Interval<I, T>>) -> COITree<I, T> {
        return Self::with_sparsity(intervals, DEFAULT_SPARSITY);
    }

    // TODO: is intervals getting copied here? Should I be using a mutable ref?
    pub fn with_sparsity(mut intervals: Vec<Interval<I, T>>, sparsity: f64) -> COITree<I, T> {
        assert!(sparsity > 1.0);

        let n = intervals.len();

        if n == 0 {
            return Self{
                index: BTreeMap::new(),
                intervals: Vec::new(),
                metadata: Vec::new()
            }
        }

        let max_size: usize = ((sparsity / (sparsity - 1.0)) * (n as f64)).ceil() as usize;
        dbg!(max_size);

        let mut searchable_intervals: Vec<Interval<I, u32>> = Vec::with_capacity(max_size);
        let mut metadata: Vec<T> = Vec::with_capacity(max_size);
        let mut index: BTreeMap<I, usize> = BTreeMap::new();

        index.insert(I::min_value(), 0);

        let now = Instant::now();
        intervals.sort_unstable_by_key(|i| i.last);
        eprintln!("first sort: {}", now.elapsed().as_millis() as f64 / 1000.0);

        let now = Instant::now();
        let mut surplus_tree = SurplusTree::new(&intervals, sparsity);
        eprintln!("SurplusTree::new: {}", now.elapsed().as_millis() as f64 / 1000.0);

        let mut useless_count = 0;
        let mut max_useless_first = I::min_value();

        let now = Instant::now();
        let mut i = 0;
        while i < n && surplus_tree.num_live_nodes() > 0 {
            // eprintln!("USELESS: {} {:?}, {:?}", useless_count, max_useless_first, surplus_tree.nodes[0].min_useful_first);
            if max_useless_first <= surplus_tree.min_useful_first() && useless_count > 0 {
                let boundary = intervals[i].last;
                // dbg!(boundary);
                // eprintln!("dumping {} useless points", useless_count);

                // i-1 should be the last interval to fall below the boundary
                let max_end = surplus_tree.index_to_leaf[i];

                // TODO: Are we sure this should purge all useless points?

                surplus_tree.map_prefix(max_end, |tree, idx, leaf_num| {
                    let interval = &intervals[idx];
                    if interval.last < boundary {
                        assert!(interval.last < boundary);
                        searchable_intervals.push(Interval{
                            first: interval.first,
                            last: interval.last,
                            metadata: leaf_num,
                        });

                        metadata.push(interval.metadata);

                        tree.set_node_dead(idx);
                        useless_count -= 1;
                    }
                });

                // dbg!(max_useless_first);
                // dbg!(surplus_tree.num_live_nodes());
                // eprintln!("post useless_count = {}", useless_count);
                assert!(useless_count == 0);
                index.insert(boundary, searchable_intervals.len());
                continue;
            }

            let max_end_opt = if n - i <= MIN_FINAL_SEQ_LEN {
                i = n-1;
                surplus_tree.last_live_leaf()
            } else {
                surplus_tree.find_sparse_query_prefix()
            };

            if let Some(max_end) = max_end_opt {
                // eprintln!("OUTPUTING Li");
                let boundary = intervals[i].last;

                let mut killed_count = 0;
                let mut l_count = 0;

                // TODO: We can skip pushing intervals above the boundary
                // when the x-values of everything below the boundary
                // is <= the x-values of everything above the boundary.

                // I think we are on to something. With this scheme we are
                // oppourtunistically dropping suffixes. What we want to do though
                // is dump prefixes whenever we can, even if there isn't a strictly
                // sparse query.
                //
                // So we need some kind of scheme to keep track of the maximum
                // x value among useless points, and the minimum x value among
                // useful points.
                //
                // This runs the risk of making the index too dense though,
                // which also be counterproductive.
                //
                // I guess what we want is to dump any prefix that's separable
                // from its suffix and above some size, to avoid just indexing
                // every single value, which would be suboptimal.


                // So how do we keep track? Every time we mark a node as useless
                // update the maximum useless x. Also keep track of the number
                // of useless nodes. Easy.

                // How de we determine the minimum x value among useful nodes?
                // I think we have to use the SurplusTree to do this. Each node
                // can keep track of the minimum x value among useful points
                // in its subtree. When we mark a node as useless, we climb up
                // the tree updating these.
                //
                // Expensive, and involves more tree climbing, but I think that's
                // the only way.

                // TODO: Another possible consequence of this scheme is that we
                // could mark entire subtrees as dead in one go.

                let mut max_prefix_first = I::min_value();
                let mut min_suffice_first = I::max_value();
                surplus_tree.map_prefix(max_end, |tree, idx, leaf_num| {
                    let interval = &intervals[idx];
                    if interval.last < boundary {
                        max_prefix_first = max(max_prefix_first, interval.first);
                    } else {
                        min_suffice_first = min(min_suffice_first, interval.first);
                    }
                });
                // dbg!((max_prefix_first, min_suffice_first));
                let supress_suffix = i < n - 1 && max_prefix_first <= min_suffice_first;
                // dbg!(max_prefix_first <= min_suffice_first);
                // dbg!(boundary);
                // dbg!(supress_suffix);

                surplus_tree.map_prefix(max_end, |tree, idx, leaf_num| {
                    let interval = &intervals[idx];
                    if interval.last < boundary || !supress_suffix {
                        searchable_intervals.push(Interval{
                            first: interval.first,
                            last: interval.last,
                            metadata: leaf_num,
                        });

                        metadata.push(interval.metadata);
                        l_count += 1;

                        if interval.last < boundary {
                            tree.set_node_dead(idx);
                            useless_count -= 1;
                            killed_count += 1;
                        }
                    }
                });

                if i < n-1 {
                    index.insert(boundary, searchable_intervals.len());
                }
            }

            // mark this interval and all following intervals that share
            // a last value as useless.
            loop {
                surplus_tree.set_node_useless(i);
                useless_count += 1;
                max_useless_first = max(max_useless_first, intervals[i].first);
                i += 1;

                if i >= n || intervals[i].last != intervals[i-1].last {
                    break;
                }
            }
        }
        eprintln!("main construction loop: {}", now.elapsed().as_millis() as f64 / 1000.0);

        dbg!(searchable_intervals.len());
        dbg!(index.len());

        let index_density = index.len() as f64 /  searchable_intervals.len() as f64;
        dbg!(index_density);

        // dbg!(&searchable_intervals);
        // dbg!(&index);

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

    // TODO: disabling inlining here to make profiling easier
    #[inline(never)]
    pub fn query_count(&self, first: I, last: I) -> usize
            where I: Bounded + Ord + Copy + Debug {

        // eprintln!("query: ({:?}, {:?})", first, last);
        let mut misses1 = 0;
        let mut misses2 = 0;
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
                    // eprintln!("hit: ({:?}, {:?})", interval.first, interval.last);
                } else {
                    if interval.metadata > last_hit_id {
                        // eprintln!("type 1 miss: ({:?}, {:?})", interval.first, interval.last);
                        misses1 += 1;
                    } else {
                        // eprintln!("type 2 miss: ({:?}, {:?})", interval.first, interval.last);
                        misses2 += 1;
                    }
                }
            }
        }

        // let hit_prop = count as f64 / (count + misses1 + misses2) as f64;
        // if hit_prop < 0.5 {
        //     eprintln!("query: ({:?}, {:?})", first, last);
            // println!("misses: {}, {}",
            //     misses1 as f64 / (count + misses1 + misses2) as f64,
            //     misses2 as f64 / (count + misses1 + misses2) as f64);
        // }

        return count;
    }
}
