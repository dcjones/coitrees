
use std::marker::Copy;
use std::convert::{TryInto, TryFrom};
use std::fmt::Debug;
use std::ops::{AddAssign, SubAssign};

// small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees. Performance isn't
// super sensitive, but is worse with a very small or very large number.
const SIMPLE_SUBTREE_CUTOFF: usize = 32;


// Set up trait for possible index types.
pub trait IntWithMax: TryInto<usize> + TryFrom<usize> + Copy + Default + PartialEq + Ord + AddAssign + SubAssign {
    const MAX: Self;

    // typically the branch here should be optimized out, because we are
    // converting, e.g. a u32 to a usize on 64-bit system.
    #[inline(always)]
    fn to_usize(self) -> usize {
        match self.try_into() {
            Ok(x) => return x,
            Err(_) => panic!("index conversion to usize failed")
        }
    }

    #[inline(always)]
    fn from_usize(x: usize) -> Self {
        match x.try_into() {
            Ok(y) => return y,
            Err(_) => panic!("index conversion from usize failed")
        }
    }

    fn one() -> Self {
        return Self::from_usize(1);
    }
}

impl IntWithMax for usize {
    const MAX: usize = usize::MAX;
}

impl IntWithMax for u32 {
    const MAX: u32 = u32::MAX;
}


impl IntWithMax for u16 {
    const MAX: u16 = u16::MAX;
}


// interval node structure forming the tree
// Nodes can be "simple" meaning they just give a span of sorted intervals
// rather than a left and right child.
#[derive(Copy, Clone, Debug)]
pub struct IntervalNode<T, I> where T: Copy, I: IntWithMax
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


impl<T, I> IntervalNode<T, I> where T: Copy, I: IntWithMax {
    pub fn new(first: i32, last: i32, metadata: T) -> IntervalNode<T, I> {
        return IntervalNode{
            subtree_last: last,
            first: first,
            last: last,
            left: I::MAX,
            right: I::MAX,
            metadata: metadata
        };
    }

    pub fn len(&self) -> i32 {
        return (self.last - self.first + 1).max(0);
    }
}


#[test]
fn test_interval_len() {
    fn make_interval(first: i32, last: i32) -> IntervalNode<(), u32> {
        return IntervalNode::new(first, last, ());
    }

    assert_eq!(make_interval(1, -1).len(), 0);
    assert_eq!(make_interval(1, 0).len(), 0);
    assert_eq!(make_interval(1, 1).len(), 1);
    assert_eq!(make_interval(1, 2).len(), 2);
}


pub struct COITree<T, I>  where T: Copy, I: IntWithMax {
    nodes: Vec<IntervalNode<T, I>>,
    root_idx: usize,
    height: usize
}


impl<T, I> COITree<T, I> where T: Copy, I: IntWithMax {
    pub fn new(nodes: Vec<IntervalNode<T, I>>) -> COITree<T, I> {
        if nodes.len() >= (I::MAX).to_usize() {
            panic!("COITree construction failed: more intervals that index type and enumerate")
        }

        let (nodes, root_idx, height) = veb_order(nodes);
        return COITree { nodes: nodes, root_idx: root_idx, height: height };
    }

    pub fn len(&self) -> usize {
        return self.nodes.len();
    }

    pub fn is_empty(&self) -> bool {
        return self.nodes.is_empty();
    }

    // find overlaps and call `visit` on every overlapping node
    pub fn query<'a, F>(&'a self, first: i32, last: i32, mut visit: F) where F: FnMut(&'a IntervalNode<T, I>) {
        if !self.is_empty() {
            query_recursion(&self.nodes, self.root_idx, first, last, &mut visit);
        }
    }

    // Count the number of overlaps. This can be done with `query`, but this
    // is slightly faster in cases of a large number of overlaps.
    pub fn query_count(&self, first: i32, last: i32) -> usize  {
        if !self.is_empty() {
            return query_recursion_count(&self.nodes, self.root_idx, first, last);
        } else {
            return 0;
        }
    }

    // Return the proportion of the query covered by intervals in the tree.
    pub fn coverage(&self, first: i32, last: i32) -> (usize, usize) {
        assert!(last >= first);

        if self.is_empty() {
            return (0, 0);
        }

        let (mut uncov_len, last_cov, count) = coverage_recursion(
            &self.nodes, self.root_idx, first, last, first - 1);

        if last_cov < last {
            uncov_len += last - last_cov;
        }

        let cov = ((last - first + 1) as usize) - (uncov_len as usize);

        return (count, cov);
    }

    pub fn iter(&self) -> COITreeIterator<T, I> {
        let mut i = self.root_idx;
        let mut stack: Vec<usize> = Vec::with_capacity(self.height);
        while i < self.nodes.len() &&
                self.nodes[i].left != I::MAX &&
                self.nodes[i].left != self.nodes[i].right {
            stack.push(i);
            i = self.nodes[i].left.to_usize();
        }

        return COITreeIterator{
            nodes: &self.nodes,
            i: i,
            count: 0,
            stack: stack,
        };
    }
}


impl<'a, T, I> IntoIterator for &'a COITree<T, I> where T: Copy, I: IntWithMax {
    type Item = &'a IntervalNode<T, I>;
    type IntoIter = COITreeIterator<'a, T, I>;

    fn into_iter(self) -> COITreeIterator<'a, T, I> {
        return self.iter();
    }
}


// Iterate through nodes in a tree in sorted order
pub struct COITreeIterator<'a, T, I> where T: Copy, I: IntWithMax {
    nodes: &'a Vec<IntervalNode<T, I>>,
    i: usize, // current node
    count: usize, // number generated so far
    stack: Vec<usize>,
}


impl<'a, T, I> Iterator for COITreeIterator<'a, T, I> where T: Copy, I: IntWithMax {
    type Item = &'a IntervalNode<T, I>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.nodes.len() {
            return None;
        }

        let node = &self.nodes[self.i];

        if node.left == node.right { // simple node
            if node.left.to_usize() > 1 {
                self.i += 1;
            } else {
                if let Some(i) = self.stack.pop() {
                    self.i = i;
                } else {
                    self.i = usize::max_value();
                }
            }
        } else {
            if node.right == I::MAX {
                if let Some(i) = self.stack.pop() {
                    self.i = i;
                } else {
                    self.i = usize::max_value();
                }
            } else {
                let mut i = node.right.to_usize();

                while self.nodes[i].left != I::MAX
                        && self.nodes[i].left != self.nodes[i].right {
                    self.stack.push(i);
                    i = self.nodes[i].left.to_usize();
                }

                self.i = i;
            }
        }

        self.count += 1;
        return Some(node);
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.nodes.len() - self.count;
        return (len, Some(len));
    }
}


impl<'a, T, I> ExactSizeIterator for COITreeIterator<'a, T, I> where T: Copy, I: IntWithMax {
    fn len(&self) -> usize {
        return self.nodes.len() - self.count;
    }
}


// Recursively count overlaps between the tree specified by `nodes` and a
// query interval specified by `first`, `last`.
fn query_recursion<'a, T, I, F>(
        nodes: &'a [IntervalNode<T, I>], root_idx: usize, first: i32, last: i32,
        visit: &mut F) where T: Copy, I: IntWithMax, F: FnMut(&'a IntervalNode<T, I>) {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
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
        nodes: &[IntervalNode<T, I>], root_idx: usize, first: i32, last: i32) -> usize
            where T: Copy, I: IntWithMax {

    let node = nodes[root_idx];

    if node.left == node.right { // simple subtree
        let mut count = 0;
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.first {
                break;
            } else if first <= node.last {
                count += 1;
            }
        }
        return count;
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

        return count;
    }
}


fn coverage_recursion<T, I>(
        nodes: &[IntervalNode<T, I>], root_idx: usize,
        first: i32, last: i32, mut last_cov: i32) -> (i32, i32, usize)
        where T: Copy, I: IntWithMax {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
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
        return (uncov_len, last_cov, count);
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

        return (uncov_len, last_cov, count);
    }
}


// Used to perform dense sorted queries more efficiently by leveraging
// what was learned from the previous query to speed up the current one, if
// the current query is a nearby successor to the previous.
pub struct SortedQuerent<'a, T, I> where T: Copy, I: IntWithMax {
    tree: &'a COITree<T, I>,
    prev_first: i32,
    prev_last: i32,
    overlapping_intervals: Vec<&'a IntervalNode<T, I>>,
}


impl<'a, T, I> SortedQuerent<'a, T, I> where T: Copy, I: IntWithMax {
    pub fn new(tree: &'a COITree<T, I>) -> SortedQuerent<'a, T, I> {
        let overlapping_intervals: Vec<&IntervalNode<T, I>> = Vec::new();
        return SortedQuerent {
            tree: tree,
            prev_first: -1,
            prev_last: -1,
            overlapping_intervals: overlapping_intervals,
        };
    }

    pub fn query<F>(&mut self, first: i32, last: i32, mut visit: F) where F: FnMut(&IntervalNode<T, I>) {

        if self.tree.is_empty() {
            return;
        }

        // not overlaping or preceding
        if first < self.prev_first || first > self.prev_last {
            // no overlap with previous query. have to resort to regular query strategy
            self.overlapping_intervals.clear();
            self.tree.query(first, last, |node| self.overlapping_intervals.push(node));
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
                    &self.tree.nodes, self.tree.root_idx,
                    self.prev_last+1, last, &mut self.overlapping_intervals);
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
        nodes: &'a [IntervalNode<T, I>], root_idx: usize, first: i32, last: i32,
        overlapping_intervals: &mut Vec<&'a IntervalNode<T, I>>) where T: Copy, I: IntWithMax {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.first {
                break;
            } else if first <= node.first {
                overlapping_intervals.push(node);
            }
        }
    }
    else {
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
    return first_a <= last_b && last_a >= first_b;
}


// Used by `traverse` to keep record tree metadata.
#[derive(Copy, Clone, Debug, Default)]
struct TraversalInfo<I> where I: IntWithMax {
    depth: u32,
    inorder: I, // in-order visit number
    preorder: I, // pre-order visit number
    subtree_size: I,
    left: I,
    right: I,
    parent: I,
    simple: bool // set by veb_order_recursion
}


// dfs traversal of an implicit bst computing dfs number, node depth, subtree
// size, and left and right pointers.
fn traverse<T, I>(nodes: &mut [IntervalNode<T, I>]) -> Vec<TraversalInfo<I>>
        where T: Copy, I: IntWithMax {
    let n = nodes.len();
    let mut info = vec![TraversalInfo::default(); n];
    let mut inorder = 0;
    let mut preorder = 0;
    traverse_recursion(nodes, &mut info, 0, n, 0, I::MAX, &mut inorder, &mut preorder);

    return info;
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
        preorder: &mut usize) -> I
        where T: Copy, I: IntWithMax {

    if start >= end {
        return I::MAX;
    }

    let root_idx = start + (end - start) / 2;
    let left;
    let right;
    let subtree_size = end - start;

    info[root_idx].depth = depth;
    info[root_idx].preorder = I::from_usize(*preorder);
    info[root_idx].parent = parent;
    *preorder += 1;

    if root_idx > start {
        left = traverse_recursion(
                nodes, info, start, root_idx, depth+1,
                I::from_usize(root_idx), inorder, preorder);
        if nodes[left.to_usize()].subtree_last > nodes[root_idx].subtree_last {
            nodes[root_idx].subtree_last = nodes[left.to_usize()].subtree_last;
        }
    } else {
        left = I::MAX;
    }

    info[root_idx].inorder = I::from_usize(*inorder);
    *inorder += 1;

    if root_idx + 1 < end {
        right = traverse_recursion(
            nodes, info, root_idx+1, end, depth+1,
            I::from_usize(root_idx), inorder, preorder);
        if nodes[right.to_usize()].subtree_last > nodes[root_idx].subtree_last {
            nodes[root_idx].subtree_last = nodes[right.to_usize()].subtree_last;
        }
    } else {
        right = I::MAX;
    }

    info[root_idx].subtree_size = I::from_usize(subtree_size);
    info[root_idx].left = left;
    info[root_idx].right = right;

    return I::from_usize(root_idx);
}


// norder partition by depth on pivot into three parts, like so
//      [ bottom left ][ top ][ bottom right ]
// where bottom left and right are the bottom subtrees with positioned to
// the left and right of the root node
fn stable_ternary_tree_partition<I>(
        input: &[I], output: &mut [I], partition: &mut [i8],
        info: &[TraversalInfo<I>], pivot_depth: u32, pivot_dfs: I,
        start: usize, end: usize) -> (usize, usize) where I: IntWithMax {

    let n = end - start;

    // determine which partition each index goes in
    let mut bottom_left_size = 0;
    let mut top_size = 0;
    let mut bottom_right_size = 0;
    for (i, j) in input[start..end].iter().enumerate() {
        let info_j = info[j.to_usize()];
        let p: i8;
        if info_j.depth <= pivot_depth {
            p = 0;
            top_size += 1;
        } else if info_j.inorder < pivot_dfs {
            p = -1;
            bottom_left_size += 1;
        } else {
            p = 1;
            bottom_right_size += 1;
        }
        partition[start+i] = p;
    }
    assert!(bottom_left_size + top_size + bottom_right_size == n);

    // do the partition
    let mut bl = start;
    let mut t = bl + bottom_left_size;
    let mut br = t + top_size;
    for (i, p) in input[start..end].iter().zip(&partition[start..end]) {
        if *p < 0 {
            output[bl] = *i;
            bl += 1;
        } else if *p == 0 {
            output[t] = *i;
            t += 1;
        } else {
            output[br] = *i;
            br += 1;
        }
    }
    assert!(br == end);

    return (bl, t);
}


// put nodes in van Emde Boas order
fn veb_order<T, I>(mut nodes: Vec<IntervalNode<T, I>>) -> (Vec<IntervalNode<T, I>>, usize, usize)
        where T: Copy, I: IntWithMax {

    // let now = Instant::now();
    // nodes.sort_unstable_by_key(|node| node.first);
    let mut veb_nodes = nodes.clone();
    let n = veb_nodes.len();

    if veb_nodes.is_empty() {
        return (veb_nodes, 0, 0);
    }

    let mut nodes_presorted = true;
    for i in 1..n {
        if nodes[i].first < nodes[i-1].first {
            nodes_presorted = false;
            break;
        }
    }
    if !nodes_presorted {
        radix_sort_nodes(&mut nodes, &mut veb_nodes);
    }
    // eprintln!("sorting nodes: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    let mut info = traverse(&mut nodes);
    // eprintln!("traversing: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    let mut max_depth = 0;
    for info_i in &info {
        if info_i.depth > max_depth {
            max_depth = info_i.depth;
        }
    }

    let idxs: &mut [I] = &mut vec![I::default(); n];
    for i in 0..n {
        idxs[i] = I::from_usize(i);
    }
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
        idxs, tmp, partition, &mut info, &mut nodes, 0, n, true, false, 0, max_depth);
    // eprintln!("computing order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[j.to_usize()] = I::from_usize(i);
    }

    // put nodes in vEB order in a temporary vector
    for i in 0..n {
        veb_nodes[i] = nodes[idxs[i.to_usize()].to_usize()];

        if info[idxs[i].to_usize()].simple {
            veb_nodes[i].left = info[idxs[i].to_usize()].subtree_size;
            veb_nodes[i].right = veb_nodes[i].left;
        } else {
            // update left and right pointers
            let left = info[idxs[i].to_usize()].left;
            veb_nodes[i].left = if left < I::MAX {
                revidx[left.to_usize()]
            } else {
                left
            };

            let right = info[idxs[i].to_usize()].right;
            veb_nodes[i].right = if right < I::MAX {
                revidx[right.to_usize()]
            } else {
                right
            };
        }
    }

    let root_idx = revidx[root_idx.to_usize()].to_usize();

    // eprintln!("ordering: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    assert!(compute_tree_size(&veb_nodes, root_idx) == n);

    return (veb_nodes, root_idx, max_depth as usize);
}


// Traverse the tree and return the size, used as a basic sanity check.
fn compute_tree_size<T, I>(nodes: &[IntervalNode<T, I>], root_idx: usize) -> usize
        where T: Copy, I: IntWithMax {

    let mut subtree_size = 1;

    let node = nodes[root_idx];
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

    return subtree_size;
}


// recursively reorder indexes to put it in vEB order. Called by `veb_order`
//   idxs: current permutation
//   tmp: temporary space of equal length to idxs
//   nodes: the interval nodes (in sorted order)
//   start, end: slice within idxs to be reordered
//   childless: true if this slice is a proper subtree and has no children below it
//   parity: true if idxs and tmp are swapped and need to be copied back,
//   min_depth, max_depth: depth extreme of the start..end slice
//
fn veb_order_recursion<T, I>(
        idxs: &mut [I], tmp: &mut [I], partition: &mut [i8],
        info: &mut [TraversalInfo<I>],
        nodes: &mut [IntervalNode<T, I>],
        start: usize, end: usize,
        childless: bool,
        parity: bool,
        min_depth: u32, max_depth: u32) -> I where T: Copy, I: IntWithMax {
    let n = (start..end).len();

    // small subtrees are put into sorted order and just searched through
    // linearly. There is a little trickiness to this because we have to
    // update the parent's child pointers and some other fields.
    if childless && info[idxs[start].to_usize()].subtree_size.to_usize() <= SIMPLE_SUBTREE_CUTOFF {
        assert!(n == info[idxs[start].to_usize()].subtree_size.to_usize());

        let old_root = idxs[start];

        idxs[start..end].sort_unstable_by_key(|i| info[i.to_usize()].inorder);
        info[idxs[start].to_usize()].simple = true;
        let subtree_size = info[old_root.to_usize()].subtree_size;
        info[idxs[start].to_usize()].subtree_size = subtree_size;
        nodes[idxs[start].to_usize()].subtree_last = nodes[old_root.to_usize()].subtree_last;

        // all children nodes record the size of the remaining list
        // let mut subtree_i_size = subtree_size - i;
        let mut subtree_i_size = subtree_size;
        for i in 0..info[idxs[start].to_usize()].subtree_size.to_usize() {
            info[idxs[start+i].to_usize()].subtree_size = subtree_i_size;
            info[idxs[start+i].to_usize()].simple = true;
            subtree_i_size -= I::one();
        }

        let parent = info[old_root.to_usize()].parent;
        if parent < I::MAX {
            if info[parent.to_usize()].left == old_root {
                info[parent.to_usize()].left = idxs[start];
            } else {
                assert!(info[parent.to_usize()].right == old_root);
                info[parent.to_usize()].right = idxs[start];
            }
        }

        if parity {
            tmp[start..end].copy_from_slice(&idxs[start..end]);
        }
        return idxs[start];
    }

    // very small trees are already in order
    if n <= 3 {
        if parity {
            tmp[start..end].copy_from_slice(&idxs[start..end]);
        }
        return idxs[start];
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let pivot_dfs = info[idxs[start].to_usize()].inorder;

    let (top_start, bottom_right_start) =
        stable_ternary_tree_partition(
            idxs, tmp, partition, info, pivot_depth, pivot_dfs, start, end);

    // tmp is not partitioned, so swap pointers
    let (tmp, idxs) = (idxs, tmp);

    // recurse on top subtree
    let top_root_idx = veb_order_recursion(
        idxs, tmp, partition, info, nodes, top_start, bottom_right_start, false,
        !parity, min_depth, pivot_depth);

    // find on recurse on subtrees in the bottom left partition and bottom right partition
    for (part_start, part_end) in &[(start, top_start), (bottom_right_start, end)] {
        let bottom_subtree_depth = pivot_depth + 1;
        let mut i = *part_start;
        while i < *part_end {
            assert!(info[idxs[i].to_usize()].depth == bottom_subtree_depth);
            let mut j = i+1;
            let mut subtree_max_depth = info[idxs[i].to_usize()].depth;
            while j < *part_end && info[idxs[j].to_usize()].depth != bottom_subtree_depth {
                assert!(info[idxs[j].to_usize()].depth > bottom_subtree_depth);
                if info[idxs[j].to_usize()].depth > subtree_max_depth {
                    subtree_max_depth = info[idxs[j].to_usize()].depth;
                }
                j += 1;
            }
            veb_order_recursion(
                idxs, tmp, partition, info, nodes, i, j, childless, !parity,
                bottom_subtree_depth, subtree_max_depth);
            i = j;
        }
    }

    return top_root_idx;
}


// is already sorted.
// Simple two pass radix sort of 32bit integers (16bits at a time) to sort nodes
// on start position. tmp is temporary space for the first pass of equal length
// to nodes.
fn radix_sort_nodes<T, I>(nodes: &mut [IntervalNode<T, I>], tmp: &mut [IntervalNode<T, I>])
        where T: Copy, I: IntWithMax {

    let mut max_first = 0;
    for node in &*nodes {
        max_first = max_first.max(node.first);
    }

    let mut count = 0;

    const R: usize = 16;
    const K: usize = 0xffff+1;
    const MASK: i32 = 0xffff;

    let mut shift: usize = 0;
    // let mut radix_counts: [u32; K] = [0; K];
    let mut radix_counts: Vec<u32> = vec![0; K];

    let mut from = nodes;
    let mut to   = tmp;
    while count < 32/R {
        for radix_count in &mut radix_counts {
            *radix_count = 0;
        }

        for node in &*from {
            radix_counts[((node.first >> shift) & MASK) as usize] += 1;
        }

        // make counts cumulative
        for i in 1..K {
            radix_counts[i] += radix_counts[i-1];
        }

        // change counts to offsets
        for i in 0..K-1 {
            radix_counts[K-1-i] = radix_counts[K-2-i];
        }
        radix_counts[0] = 0;

        for node in &*from {
            let radix = ((node.first >> shift) & MASK) as usize;
            to[radix_counts[radix] as usize] = *node;
            radix_counts[radix] += 1;
        }

        count += 1;
        shift += R;

        let swap_tmp = from;
        from = to;
        to = swap_tmp;
    }

    // only needed if we do an odd number of iterations
    // if count % 2 == 1 {
    //     to.copy_from_slice(from);
    // }
}

