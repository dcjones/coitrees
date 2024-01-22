//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!
//! The data structure used a fairly standard interval tree, but with nodes stored
//! in van Emde Boas layout, which improves average cache locality, and thus
//! query performance. The downside is that building the tree is more expensive
//! so a relatively large number of queries needs to made for it to break even.
//!
//! The data structure `COITree` is constructed with an array of `IntervalNode`
//! structs which store integer, end-inclusive intervals along with associated
//! metadata. The tree can be queried directly for coverage or overlaps, or
//! through the intermediary `SortedQuerent` which keeps track of some state
//! to accelerate overlaping queries.

use super::interval::{GenericInterval, IntWithMax, Interval, IntervalTree, SortedQuerent};
use std::cmp::Ordering;
use std::fmt::Debug;
use std::iter::IntoIterator;

// Small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees. Performance isn't
// super sensitive, but is worse with a very small or very large number.
const SIMPLE_SUBTREE_CUTOFF: usize = 64;

// Very dense subtrees in which we probably intersect most of the intervals
// are more efficient to query linearly. When the expected proportion of hits
// is a above this number it becomes a simple subtree.
const SIMPLE_SUBTREE_DENSITY_CUTOFF: f32 = 0.2;

/// Internal interval node representation used by BasicCOITree
#[derive(Clone)]
pub struct IntervalNode<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    // subtree interval
    subtree_last: i32,

    // interval
    pub first: i32,
    pub last: i32,

    // when this is the root of a simple subtree, left == right is the size
    // of the subtree, otherwise they are left, right child pointers.
    left: I,
    right: I,

    pub metadata: T,
}

impl<T, I> IntervalNode<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    pub fn new(first: i32, last: i32, metadata: T) -> IntervalNode<T, I> {
        Self {
            subtree_last: last,
            first,
            last,
            left: I::MAX,
            right: I::MAX,
            metadata,
        }
    }

    fn from_interval<V>(interval: &V) -> IntervalNode<T, I>
    where
        V: GenericInterval<T>,
    {
        return Self::new(
            interval.first(),
            interval.last(),
            interval.metadata().clone(),
        );
    }

    /// Length spanned by the interval. (Interval are end-inclusive.)
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> i32 {
        (self.last - self.first + 1).max(0)
    }
}

// IntervalNodes are themselves a type of annotated interval
impl<T, I> GenericInterval<T> for IntervalNode<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    fn first(&self) -> i32 {
        self.first
    }

    fn last(&self) -> i32 {
        self.last
    }

    fn metadata(&self) -> &T {
        &self.metadata
    }
}

#[test]
fn test_interval_len() {
    fn make_interval(first: i32, last: i32) -> IntervalNode<(), u32> {
        IntervalNode::new(first, last, ())
    }

    assert_eq!(make_interval(1, -1).len(), 0);
    assert_eq!(make_interval(1, 0).len(), 0);
    assert_eq!(make_interval(1, 1).len(), 1);
    assert_eq!(make_interval(1, 2).len(), 2);
}

/// COITree data structure. A representation of a static set of intervals with
/// associated metadata, enabling fast overlap and coverage queries.
///
/// The index type `I` is a typically `usize`, but can be `u32` or `u16`.
/// It's slightly more efficient to use a smalled index type, assuming there are
/// fewer than I::MAX-1 intervals to store.
#[derive(Clone)]
pub struct BasicCOITree<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    nodes: Vec<IntervalNode<T, I>>,
    root_idx: usize,
    height: usize,
}

impl<'a, T, I> BasicCOITree<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    // Exactly the same as IntervalTree::query, but provides a reference with the
    // lifetime of the tree (not something we can guarantee with other implementations).
    fn query_static<F>(&'a self, first: i32, last: i32, mut visit: F)
    where
        F: FnMut(&'a IntervalNode<T, I>),
    {
        if !self.is_empty() {
            query_recursion(&self.nodes, self.root_idx, first, last, &mut visit);
        }
    }
}

impl<'a, T, I> IntervalTree<'a> for BasicCOITree<T, I>
where
    T: Clone + 'a,
    I: IntWithMax + 'a,
{
    type Metadata = T;
    type Index = I;
    type Item = IntervalNode<T, I>;
    type Iter = BasicCOITreeIterator<'a, T, I>;

    fn new<'c, U, V>(intervals: U) -> BasicCOITree<T, I>
    where
        U: IntoIterator<Item = &'c V>,
        V: GenericInterval<T> + 'c,
    {
        let nodes: Vec<IntervalNode<T, I>> = intervals
            .into_iter()
            .map(|node| IntervalNode::from_interval(node))
            .collect();
        if nodes.len() >= (I::MAX).to_usize() {
            panic!("BasicCOITree construction failed: more intervals than index type can enumerate")
        }

        let (nodes, root_idx, height) = veb_order(nodes);

        BasicCOITree {
            nodes,
            root_idx,
            height,
        }
    }

    /// Number of intervals in the set.
    fn len(&self) -> usize {
        self.nodes.len()
    }

    /// True iff the set is empty.
    fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Find intervals in the set overlaping the query `[first, last]` and call `visit` on every overlapping node
    fn query<F>(&'a self, first: i32, last: i32, mut visit: F)
    where
        F: FnMut(&IntervalNode<T, I>),
    {
        if !self.is_empty() {
            query_recursion(&self.nodes, self.root_idx, first, last, &mut visit);
        }
    }

    /// Count the number of intervals in the set overlapping the query `[first, last]`.
    fn query_count(&self, first: i32, last: i32) -> usize {
        if !self.is_empty() {
            query_recursion_count(&self.nodes, self.root_idx, first, last)
        } else {
            0
        }
    }

    /// Return a pair `(count, cov)`, where `count` gives the number of intervals
    /// in the set overlapping the query, and `cov` the number of positions in the query
    /// interval covered by at least one interval in the set.
    fn coverage(&self, first: i32, last: i32) -> (usize, usize) {
        assert!(last >= first);

        if self.is_empty() {
            return (0, 0);
        }

        let (mut uncov_len, last_cov, count) =
            coverage_recursion(&self.nodes, self.root_idx, first, last, first - 1);

        if last_cov < last {
            uncov_len += last - last_cov;
        }

        let cov = ((last - first + 1) as usize) - (uncov_len as usize);

        (count, cov)
    }

    /// Iterate through the interval set in sorted order by interval start position.
    fn iter(&'a self) -> BasicCOITreeIterator<'a, T, I> {
        let mut i = self.root_idx;
        let mut stack: Vec<usize> = Vec::with_capacity(self.height);
        while i < self.nodes.len()
            && self.nodes[i].left != I::MAX
            && self.nodes[i].left != self.nodes[i].right
        {
            stack.push(i);
            i = self.nodes[i].left.to_usize();
        }

        BasicCOITreeIterator {
            nodes: &self.nodes,
            i,
            count: 0,
            stack,
        }
    }
}

impl<'a, T, I> IntoIterator for &'a BasicCOITree<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    type Item = Interval<&'a T>;
    type IntoIter = BasicCOITreeIterator<'a, T, I>;

    fn into_iter(self) -> BasicCOITreeIterator<'a, T, I> {
        self.iter()
    }
}

/// Iterate through nodes in a tree in sorted order by interval start position.
pub struct BasicCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    nodes: &'a Vec<IntervalNode<T, I>>,
    i: usize,     // current node
    count: usize, // number generated so far
    stack: Vec<usize>,
}

impl<'a, T, I> Iterator for BasicCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    type Item = Interval<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.nodes.len() {
            return None;
        }

        let node = &self.nodes[self.i];

        if node.left == node.right {
            // simple node
            if node.left.to_usize() > 1 {
                self.i += 1;
            } else if let Some(i) = self.stack.pop() {
                self.i = i;
            } else {
                self.i = usize::MAX;
            }
        } else if node.right == I::MAX {
            if let Some(i) = self.stack.pop() {
                self.i = i;
            } else {
                self.i = usize::MAX;
            }
        } else {
            let mut i = node.right.to_usize();

            while self.nodes[i].left != I::MAX && self.nodes[i].left != self.nodes[i].right {
                self.stack.push(i);
                i = self.nodes[i].left.to_usize();
            }

            self.i = i;
        }

        self.count += 1;
        Some(Interval::new(node.first, node.last, &node.metadata))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.nodes.len() - self.count;
        (len, Some(len))
    }
}

impl<'a, T, I> ExactSizeIterator for BasicCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    fn len(&self) -> usize {
        self.nodes.len() - self.count
    }
}

// Recursively count overlaps between the tree specified by `nodes` and a
// query interval specified by `first`, `last`.
fn query_recursion<'a, T, I, F>(
    nodes: &'a [IntervalNode<T, I>],
    root_idx: usize,
    first: i32,
    last: i32,
    visit: &mut F,
) where
    T: Clone,
    I: IntWithMax,
    F: FnMut(&'a IntervalNode<T, I>),
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.first {
                break;
            } else if first <= node.last {
                visit(node);
            }
        }
    } else {
        if overlaps(node.first, node.last, first, last) {
            visit(node);
        }

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                query_recursion(nodes, left, first, last, visit);
            }
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.first, nodes[right].subtree_last, first, last) {
                query_recursion(nodes, right, first, last, visit);
            }
        }
    }
}

// query_recursion but just count number of overlaps
fn query_recursion_count<T, I>(
    nodes: &[IntervalNode<T, I>],
    root_idx: usize,
    first: i32,
    last: i32,
) -> usize
where
    T: Clone,
    I: IntWithMax,
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        let mut count = 0;
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.first {
                break;
            } else if first <= node.last {
                count += 1;
            }
        }
        count
    } else {
        let mut count = 0;
        if overlaps(node.first, node.last, first, last) {
            count += 1;
        }

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                count += query_recursion_count(nodes, left, first, last);
            }
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.first, nodes[right].subtree_last, first, last) {
                count += query_recursion_count(nodes, right, first, last);
            }
        }

        count
    }
}

fn coverage_recursion<T, I>(
    nodes: &[IntervalNode<T, I>],
    root_idx: usize,
    first: i32,
    last: i32,
    mut last_cov: i32,
) -> (i32, i32, usize)
where
    T: Clone,
    I: IntWithMax,
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        let mut count = 0;
        let mut uncov_len = 0;
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if overlaps(node.first, node.last, first, last) {
                if node.first > last_cov {
                    uncov_len += node.first - (last_cov + 1);
                }
                last_cov = last_cov.max(node.last);
                count += 1;
            }
        }
        (uncov_len, last_cov, count)
    } else {
        let mut uncov_len = 0;
        let mut count = 0;

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                let (left_uncov_len, left_last_cov, left_count) =
                    coverage_recursion(nodes, left, first, last, last_cov);
                last_cov = left_last_cov;
                uncov_len += left_uncov_len;
                count += left_count;
            }
        }

        if overlaps(node.first, node.last, first, last) {
            if node.first > last_cov {
                uncov_len += node.first - (last_cov + 1);
            }
            last_cov = last_cov.max(node.last);
            count += 1;
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.first, nodes[right].subtree_last, first, last) {
                let (right_uncov_len, right_last_cov, right_count) =
                    coverage_recursion(nodes, right, first, last, last_cov);
                last_cov = right_last_cov;
                uncov_len += right_uncov_len;
                count += right_count;
            }
        }

        (uncov_len, last_cov, count)
    }
}

/// An alternative query strategy that can be much faster when queries are performed
/// in sorted order and overlap.
///
/// Unilke `COITree::query`, some state is retained between queries.
/// `SortedQuerent` tracks that state. If queries are not sorted or don't
/// overlap, this strategy still works, but is slightly slower than
/// `COITree::query`.
pub struct BasicSortedQuerent<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    tree: &'a BasicCOITree<T, I>,
    prev_first: i32,
    prev_last: i32,
    overlapping_intervals: Vec<&'a IntervalNode<T, I>>,
}

impl<'a, T, I> SortedQuerent<'a> for BasicSortedQuerent<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    type Metadata = T;
    type Index = I;
    type Item = IntervalNode<T, I>;
    type Iter = BasicCOITreeIterator<'a, T, I>;
    type Tree = BasicCOITree<T, I>;

    /// Construct a new `SortedQuerent` to perform a sequence. queries.
    fn new(tree: &'a BasicCOITree<T, I>) -> BasicSortedQuerent<'a, T, I> {
        let overlapping_intervals: Vec<&IntervalNode<T, I>> = Vec::new();
        BasicSortedQuerent {
            tree,
            prev_first: -1,
            prev_last: -1,
            overlapping_intervals,
        }
    }

    /// Find intervals in the underlying `COITree` that overlap the query
    /// `[first, last]` and call `visit` on each. Works equivalently to
    /// `COITrees::query` but queries that overlap prior queries will potentially
    /// be faster.
    fn query<F>(&mut self, first: i32, last: i32, mut visit: F)
    where
        F: FnMut(&'a IntervalNode<T, I>),
    {
        if self.tree.is_empty() {
            return;
        }

        // not overlaping or preceding
        if first < self.prev_first || first > self.prev_last {
            // no overlap with previous query. have to resort to regular query strategy
            self.overlapping_intervals.clear();
            self.tree
                .query_static(first, last, |node| self.overlapping_intervals.push(node));
        } else {
            // successor query, exploit the overlap

            // delete previously overlapping intervals with end in [prev_first, first-1]
            if self.prev_first < first {
                let mut i = 0;
                while i < self.overlapping_intervals.len() {
                    if self.overlapping_intervals[i].last < first {
                        self.overlapping_intervals.swap_remove(i);
                    } else {
                        i += 1;
                    }
                }
            }

            // delete previously overlapping intervals with start in [last+1, prev_end]
            if self.prev_last > last {
                let mut i = 0;
                while i < self.overlapping_intervals.len() {
                    if self.overlapping_intervals[i].first > last {
                        self.overlapping_intervals.swap_remove(i);
                    } else {
                        i += 1;
                    }
                }
            }

            // add any interval that start in [prev_last+1, last]
            if self.prev_last < last {
                sorted_querent_query_firsts(
                    &self.tree.nodes,
                    self.tree.root_idx,
                    self.prev_last + 1,
                    last,
                    &mut self.overlapping_intervals,
                );
            }
        }

        // call visit on everything
        for overlapping_interval in &self.overlapping_intervals {
            visit(overlapping_interval);
        }

        self.prev_first = first;
        self.prev_last = last;
    }
}

// find any intervals in the tree with their first value in [first, last]
fn sorted_querent_query_firsts<'a, T, I>(
    nodes: &'a [IntervalNode<T, I>],
    root_idx: usize,
    first: i32,
    last: i32,
    overlapping_intervals: &mut Vec<&'a IntervalNode<T, I>>,
) where
    T: Clone,
    I: IntWithMax,
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.first {
                break;
            } else if first <= node.first {
                overlapping_intervals.push(node);
            }
        }
    } else {
        if first <= node.first && node.first <= last {
            overlapping_intervals.push(node);
        }

        if node.left < I::MAX && first <= node.first {
            let left = node.left.to_usize();
            sorted_querent_query_firsts(nodes, left, first, last, overlapping_intervals);
        }

        if node.right < I::MAX && last >= node.first {
            let right = node.right.to_usize();
            sorted_querent_query_firsts(nodes, right, first, last, overlapping_intervals);
        }
    }
}

// True iff the two intervals overlap.
#[inline(always)]
fn overlaps(first_a: i32, last_a: i32, first_b: i32, last_b: i32) -> bool {
    first_a <= last_b && last_a >= first_b
}

// Used by `traverse` to keep record tree metadata.
#[derive(Clone, Debug)]
struct TraversalInfo<I>
where
    I: IntWithMax,
{
    depth: u32,
    inorder: I,  // in-order visit number
    preorder: I, // pre-order visit number
    subtree_size: I,
    parent: I,
    expected_hit_proportion: f32,
}

impl<I> Default for TraversalInfo<I>
where
    I: IntWithMax,
{
    fn default() -> Self {
        TraversalInfo {
            depth: 0,
            inorder: I::default(),
            preorder: I::default(),
            subtree_size: I::one(),
            parent: I::MAX,
            expected_hit_proportion: 0.0,
        }
    }
}

// dfs traversal of an implicit bst computing dfs number, node depth, subtree
// size, and left and right pointers.
fn traverse<T, I>(nodes: &mut [IntervalNode<T, I>]) -> Vec<TraversalInfo<I>>
where
    T: Clone,
    I: IntWithMax,
{
    let n = nodes.len();
    let mut info = vec![TraversalInfo::default(); n];
    let mut inorder = 0;
    let mut preorder = 0;
    traverse_recursion(
        nodes,
        &mut info,
        0,
        n,
        0,
        I::MAX,
        &mut inorder,
        &mut preorder,
    );

    info
}

// The recursive part of the `traverse` function.
fn traverse_recursion<T, I>(
    nodes: &mut [IntervalNode<T, I>],
    info: &mut [TraversalInfo<I>],
    start: usize,
    end: usize,
    depth: u32,
    parent: I,
    inorder: &mut usize,
    preorder: &mut usize,
) -> (I, i32, f32)
where
    T: Clone,
    I: IntWithMax,
{
    if start >= end {
        return (I::MAX, i32::MAX, f32::NAN);
    }

    let root_idx = start + (end - start) / 2;
    let subtree_size = end - start;

    info[root_idx].depth = depth;
    info[root_idx].preorder = I::from_usize(*preorder);
    info[root_idx].parent = parent;
    *preorder += 1;

    let mut subtree_first = nodes[root_idx].first;

    let mut left_expected_hits = 0.0;
    let mut left_subtree_span = 0;

    if root_idx > start {
        let (left, left_subtree_first, left_expected_hits_) = traverse_recursion(
            nodes,
            info,
            start,
            root_idx,
            depth + 1,
            I::from_usize(root_idx),
            inorder,
            preorder,
        );

        left_expected_hits = left_expected_hits_;
        left_subtree_span = nodes[left.to_usize()].subtree_last - left_subtree_first + 1;

        subtree_first = left_subtree_first;
        if nodes[left.to_usize()].subtree_last > nodes[root_idx].subtree_last {
            nodes[root_idx].subtree_last = nodes[left.to_usize()].subtree_last;
        }
        nodes[root_idx].left = left;
    }

    info[root_idx].inorder = I::from_usize(*inorder);
    *inorder += 1;

    let mut right_expected_hits = 0.0;
    let mut right_subtree_span = 0;

    if root_idx + 1 < end {
        let (right, right_subtree_first, right_expected_hits_) = traverse_recursion(
            nodes,
            info,
            root_idx + 1,
            end,
            depth + 1,
            I::from_usize(root_idx),
            inorder,
            preorder,
        );

        right_expected_hits = right_expected_hits_;
        right_subtree_span = nodes[right.to_usize()].subtree_last - right_subtree_first + 1;

        if nodes[right.to_usize()].subtree_last > nodes[root_idx].subtree_last {
            nodes[root_idx].subtree_last = nodes[right.to_usize()].subtree_last;
        }
        nodes[root_idx].right = right;
    }

    info[root_idx].subtree_size = I::from_usize(subtree_size);
    let subtree_span = nodes[root_idx].subtree_last - subtree_first + 1;

    debug_assert!(left_subtree_span <= subtree_span);
    debug_assert!(right_subtree_span <= subtree_span);

    let expected_hits = ((nodes[root_idx].last - nodes[root_idx].first + 1) as f32
        + (left_subtree_span as f32) * left_expected_hits
        + (right_subtree_span as f32) * right_expected_hits)
        / subtree_span as f32;

    info[root_idx].expected_hit_proportion = expected_hits / subtree_size as f32;

    (I::from_usize(root_idx), subtree_first, expected_hits)
}

// norder partition by depth on pivot into three parts, like so
//      [ bottom left ][ top ][ bottom right ]
// where bottom left and right are the bottom subtrees with positioned to
// the left and right of the root node
fn stable_ternary_tree_partition<I>(
    input: &[I],
    output: &mut [I],
    partition: &mut [i8],
    info: &[TraversalInfo<I>],
    pivot_depth: u32,
    pivot_dfs: I,
    start: usize,
    end: usize,
) -> (usize, usize)
where
    I: IntWithMax,
{
    let n = end - start;

    // determine which partition each index goes in
    let mut bottom_left_size = 0;
    let mut top_size = 0;
    let mut bottom_right_size = 0;

    for (i, p) in input[start..end].iter().zip(&mut partition[start..end]) {
        let info_j = &info[i.to_usize()];
        if info_j.depth <= pivot_depth {
            *p = 0;
            top_size += 1;
        } else if info_j.inorder < pivot_dfs {
            *p = -1;
            bottom_left_size += 1;
        } else {
            *p = 1;
            bottom_right_size += 1;
        }
    }

    debug_assert!(bottom_left_size + top_size + bottom_right_size == n);

    // do the partition
    let mut bl = start;
    let mut t = bl + bottom_left_size;
    let mut br = t + top_size;
    for (i, p) in input[start..end].iter().zip(&partition[start..end]) {
        match p.cmp(&0) {
            Ordering::Less => {
                output[bl] = *i;
                bl += 1;
            }
            Ordering::Equal => {
                output[t] = *i;
                t += 1;
            }
            Ordering::Greater => {
                output[br] = *i;
                br += 1;
            }
        }
    }
    debug_assert!(br == end);

    (bl, t)
}

// put nodes in van Emde Boas order
fn veb_order<T, I>(mut nodes: Vec<IntervalNode<T, I>>) -> (Vec<IntervalNode<T, I>>, usize, usize)
where
    T: Clone,
    I: IntWithMax,
{
    // let now = Instant::now();
    let mut veb_nodes = nodes.clone();
    let n = veb_nodes.len();

    if veb_nodes.is_empty() {
        return (veb_nodes, 0, 0);
    }

    nodes.sort_unstable_by_key(|node| (node.first, node.last));
    // eprintln!("sorting nodes: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    let info = traverse(&mut nodes);
    // eprintln!("traversing: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    let mut max_depth = 0;
    for info_i in &info {
        if info_i.depth > max_depth {
            max_depth = info_i.depth;
        }
    }

    let idxs: &mut [I] = &mut vec![I::default(); n];

    idxs.iter_mut()
        .enumerate()
        .for_each(|(i, idx)| *idx = I::from_usize(i));

    let tmp: &mut [I] = &mut vec![I::default(); n];

    // put in dfs order
    for i in &*idxs {
        tmp[info[i.to_usize()].preorder.to_usize()] = *i;
    }
    let (idxs, tmp) = (tmp, idxs);

    // space used to by stable_ternary_tree_partition
    let partition: &mut [i8] = &mut vec![0; n];

    // let now = Instant::now();
    let root_idx = veb_order_recursion(
        idxs, tmp, partition, &info, &mut nodes, 0, n, false, 0, max_depth,
    );
    // eprintln!("computing order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[j.to_usize()] = I::from_usize(i);
    }

    // put nodes in vEB order
    for (i_, mut node) in revidx.iter().zip(nodes) {
        let i = i_.to_usize();
        if node.left != node.right {
            if node.left < I::MAX {
                node.left = revidx[node.left.to_usize()];
            }

            if node.right < I::MAX {
                node.right = revidx[node.right.to_usize()];
            }
        }
        veb_nodes[i.to_usize()] = node;
    }

    let root_idx = revidx[root_idx.to_usize()].to_usize();

    // eprintln!("ordering: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    debug_assert!(compute_tree_size(&veb_nodes, root_idx) == n);

    (veb_nodes, root_idx, max_depth as usize)
}

// Traverse the tree and return the size, used as a basic sanity check.
fn compute_tree_size<T, I>(nodes: &[IntervalNode<T, I>], root_idx: usize) -> usize
where
    T: Clone,
    I: IntWithMax,
{
    let mut subtree_size = 1;

    let node = &nodes[root_idx];
    if node.left == node.right {
        subtree_size = nodes[root_idx].right.to_usize();
    } else {
        if node.left < I::MAX {
            let left = node.left.to_usize();
            subtree_size += compute_tree_size(nodes, left);
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            subtree_size += compute_tree_size(nodes, right);
        }
    }

    subtree_size
}

// recursively reorder indexes to put it in vEB order. Called by `veb_order`
//   idxs: current permutation
//   tmp: temporary space of equal length to idxs
//   partition: space used to assist `stable_ternary_tree_partition`.
//   nodes: the interval nodes (in sorted order)
//   start, end: slice within idxs to be reordered
//   childless: true if this slice is a proper subtree and has no children below it
//   parity: true if idxs and tmp are swapped and need to be copied back,
//   min_depth, max_depth: depth extreme of the start..end slice
//
fn veb_order_recursion<T, I>(
    idxs: &mut [I],
    tmp: &mut [I],
    partition: &mut [i8],
    info: &[TraversalInfo<I>],
    nodes: &mut [IntervalNode<T, I>],
    start: usize,
    end: usize,
    parity: bool,
    min_depth: u32,
    max_depth: u32,
) -> I
where
    T: Clone,
    I: IntWithMax,
{
    let n = (start..end).len();

    assert!(n > 0);

    let childless = info[idxs[start].to_usize()].subtree_size.to_usize() == n;

    // small subtrees are put into sorted order and just searched through
    // linearly. There is a little trickiness to this because we have to
    // update the parent's child pointers and some other fields.
    if childless
        && ((info[idxs[start].to_usize()].subtree_size.to_usize() <= SIMPLE_SUBTREE_CUTOFF)
            || (info[idxs[start].to_usize()].expected_hit_proportion
                >= SIMPLE_SUBTREE_DENSITY_CUTOFF))
    {
        debug_assert!(n == info[idxs[start].to_usize()].subtree_size.to_usize());

        let old_root = idxs[start];

        idxs[start..end].sort_unstable_by_key(|i| info[i.to_usize()].inorder);
        let subtree_size = info[old_root.to_usize()].subtree_size;
        nodes[idxs[start].to_usize()].subtree_last = nodes[old_root.to_usize()].subtree_last;

        // all children nodes record the size of the remaining list
        // let mut subtree_i_size = subtree_size - i;
        let mut subtree_i_size = subtree_size;
        for idx in &idxs[start..end] {
            nodes[idx.to_usize()].left = subtree_i_size;
            nodes[idx.to_usize()].right = subtree_i_size;
            subtree_i_size -= I::one();
        }

        let parent = info[old_root.to_usize()].parent;
        if parent < I::MAX {
            if nodes[parent.to_usize()].left == old_root {
                nodes[parent.to_usize()].left = idxs[start];
            } else {
                debug_assert!(nodes[parent.to_usize()].right == old_root);
                nodes[parent.to_usize()].right = idxs[start];
            }
        }

        if parity {
            tmp[start..end].copy_from_slice(&idxs[start..end]);
        }
        return idxs[start];
    }

    // very small trees are already in order
    if n == 1 {
        if parity {
            tmp[start] = idxs[start];
        }
        return idxs[start];
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let pivot_dfs = info[idxs[start].to_usize()].inorder;

    let (top_start, bottom_right_start) = stable_ternary_tree_partition(
        idxs,
        tmp,
        partition,
        info,
        pivot_depth,
        pivot_dfs,
        start,
        end,
    );

    // tmp is not partitioned, so swap pointers
    let (tmp, idxs) = (idxs, tmp);

    // recurse on top subtree
    let top_root_idx = veb_order_recursion(
        idxs,
        tmp,
        partition,
        info,
        nodes,
        top_start,
        bottom_right_start,
        !parity,
        min_depth,
        pivot_depth,
    );

    // find on recurse on subtrees in the bottom left partition and bottom right partition
    for (part_start, part_end) in &[(start, top_start), (bottom_right_start, end)] {
        let bottom_subtree_depth = pivot_depth + 1;
        let mut i = *part_start;
        while i < *part_end {
            debug_assert!(info[idxs[i].to_usize()].depth == bottom_subtree_depth);

            let mut subtree_max_depth = info[idxs[i].to_usize()].depth;
            let mut j = *part_end;
            for (u, v) in (i + 1..*part_end).zip(&idxs[i + 1..*part_end]) {
                let depth = info[v.to_usize()].depth;
                if depth == bottom_subtree_depth {
                    j = u;
                    break;
                } else if depth > subtree_max_depth {
                    subtree_max_depth = depth;
                }
            }

            veb_order_recursion(
                idxs,
                tmp,
                partition,
                info,
                nodes,
                i,
                j,
                !parity,
                bottom_subtree_depth,
                subtree_max_depth,
            );
            i = j;
        }
    }

    top_root_idx
}
