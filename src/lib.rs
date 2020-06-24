
use std::marker::Copy;

// small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees.
const SIMPLE_SUBTREE_CUTOFF: u32 = 32;

// interval node structure forming the tree
// Nodes can be "simple" meaning they just give a span of sorted intervals
// rather than a left and right child.
#[derive(Copy, Clone, Debug)]
pub struct IntervalNode<T> where T: Copy
 {
    // subtree interval
    subtree_last: i32,

    // interval
    pub first: i32,
    pub last: i32,

    // when this is the root of a simple subtree, left == right is the size
    // of the subtree, otherwise they are left, right child pointers.
    left: u32,
    right: u32,

    pub metadata: T,
}


impl<T> IntervalNode<T> where T: Copy {
    pub fn new(first: i32, last: i32, metadata: T) -> IntervalNode<T> {
        return IntervalNode{
            subtree_last: last,
            first: first,
            last: last,
            left: u32::max_value(),
            right: u32::max_value(),
            metadata: metadata
        };
    }

    pub fn len(&self) -> i32 {
        return (self.last - self.first + 1).max(0);
    }
}


#[test]
fn test_interval_len() {
    assert_eq!(IntervalNode::new(1, -1, ()).len(), 0);
    assert_eq!(IntervalNode::new(1, 0, ()).len(), 0);
    assert_eq!(IntervalNode::new(1, 1, ()).len(), 1);
    assert_eq!(IntervalNode::new(1, 2, ()).len(), 2);
}


pub struct COITree<T>  where T: Copy {
    nodes: Vec<IntervalNode<T>>,
    root_idx: usize
}


impl<T> COITree<T> where T: Copy {
    pub fn new(nodes: Vec<IntervalNode<T>>) -> COITree<T> {
        let (nodes, root_idx) = veb_order(nodes);
        return COITree { nodes: nodes, root_idx: root_idx };
    }

    pub fn len(&self) -> usize {
        return self.nodes.len();
    }

    pub fn is_empty(&self) -> bool {
        return self.nodes.is_empty();
    }

    // find overlaps and call `visit` on every overlapping node
    pub fn query<'a, F>(&'a self, first: i32, last: i32, mut visit: F) where F: FnMut(&'a IntervalNode<T>) {
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
}


// Recursively count overlaps between the tree specified by `nodes` and a
// query interval specified by `first`, `last`.
fn query_recursion<'a, T, F>(
        nodes: &'a [IntervalNode<T>], root_idx: usize, first: i32, last: i32,
        visit: &mut F) where T: Copy, F: FnMut(&'a IntervalNode<T>) {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
        for node in &nodes[root_idx..root_idx + node.right as usize] {
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

        let left = node.left as usize;
        if left < u32::max_value() as usize {
            if nodes[left].subtree_last >= first {
                query_recursion(nodes, left, first, last, visit);
            }
        }

        let right = node.right as usize;
        if right < u32::max_value() as usize {
            if overlaps(node.first, nodes[right].subtree_last, first, last) {
                query_recursion(nodes, right, first, last, visit);
            }
        }
    }
}


// query_recursion but just count number of overlaps
fn query_recursion_count<T>(
        nodes: &[IntervalNode<T>], root_idx: usize, first: i32, last: i32) -> usize
            where T: Copy {

    let node = nodes[root_idx];

    if node.left == node.right { // simple subtree
        let mut count = 0;
        for node in &nodes[root_idx..root_idx + node.right as usize] {
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

        let left = node.left as usize;
        if left < u32::max_value() as usize {
            if nodes[left].subtree_last >= first {
                count += query_recursion_count(nodes, left, first, last);
            }
        }

        let right = node.right as usize;
        if right < u32::max_value() as usize {
            if overlaps(node.first, nodes[right].subtree_last, first, last) {
                count += query_recursion_count(nodes, right, first, last);
            }
        }

        return count;
    }
}


fn coverage_recursion<T>(
        nodes: &[IntervalNode<T>], root_idx: usize,
        first: i32, last: i32, mut last_cov: i32) -> (i32, i32, usize)
        where T: Copy {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
        let mut count = 0;
        let mut uncov_len = 0;
        for node in &nodes[root_idx..root_idx + node.right as usize] {
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

        let left = node.left as usize;
        if left < u32::max_value() as usize {
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

        let right = node.right as usize;
        if right < u32::max_value() as usize {
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
pub struct SortedQuerent<'a, T> where T: Copy {
    tree: &'a COITree<T>,
    prev_first: i32,
    prev_last: i32,
    overlapping_intervals: Vec<&'a IntervalNode<T>>,
}


impl<'a, T> SortedQuerent<'a, T> where T: Copy {
    pub fn new(tree: &'a COITree<T>) -> SortedQuerent<'a, T> {
        let overlapping_intervals: Vec<&IntervalNode<T>> = Vec::new();
        return SortedQuerent {
            tree: tree,
            prev_first: -1,
            prev_last: -1,
            overlapping_intervals: overlapping_intervals,
        };
    }

    pub fn query<F>(&mut self, first: i32, last: i32, mut visit: F) where F: FnMut(&IntervalNode<T>) {

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
fn sorted_querent_query_firsts<'a, T>(
        nodes: &'a [IntervalNode<T>], root_idx: usize, first: i32, last: i32,
        overlapping_intervals: &mut Vec<&'a IntervalNode<T>>) where T: Copy {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
        for node in &nodes[root_idx..root_idx + node.right as usize] {
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

        let left = node.left as usize;
        if left < u32::max_value() as usize && first <= node.first {
            sorted_querent_query_firsts(nodes, left, first, last, overlapping_intervals);
        }

        let right = node.right as usize;
        if right < u32::max_value() as usize && last >= node.first {
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
struct TraversalInfo {
    depth: u32,
    inorder: u32, // in-order visit number
    preorder: u32, // pre-order visit number
    subtree_size: u32,
    left: u32,
    right: u32,
    parent: u32,
    simple: bool // set by veb_order_recursion
}


// dfs traversal of an implicit bst computing dfs number, node depth, subtree
// size, and left and right pointers.
fn traverse<T>(nodes: &mut [IntervalNode<T>]) -> Vec<TraversalInfo>
        where T: Copy {
    let n = nodes.len();
    let mut info = vec![TraversalInfo::default(); n];
    let mut inorder = 0;
    let mut preorder = 0;
    traverse_recursion(nodes, &mut info, 0, n, 0, u32::MAX, &mut inorder, &mut preorder);

    return info;
}


// The recursive part of the `traverse` function.
fn traverse_recursion<T>(
        nodes: &mut [IntervalNode<T>],
        info: &mut [TraversalInfo],
        start: usize,
        end: usize,
        depth: u32,
        parent: u32,
        inorder: &mut u32,
        preorder: &mut u32) -> u32
        where T: Copy {

    if start >= end {
        return u32::max_value();
    }

    let root_idx = start + (end - start) / 2;
    let mut left = u32::max_value();
    let mut right = u32::max_value();
    let subtree_size = end - start;

    nodes[root_idx].subtree_last = nodes[root_idx].last;
    info[root_idx].depth = depth;
    info[root_idx].preorder = *preorder;
    info[root_idx].parent = parent;
    *preorder += 1;

    if root_idx > start {
        let left_ = traverse_recursion(
                nodes, info, start, root_idx, depth+1,
                root_idx as u32, inorder, preorder);
        left = left_;
        nodes[root_idx].subtree_last = nodes[root_idx].subtree_last.max(
            nodes[left as usize].subtree_last);
    }

    info[root_idx].inorder = *inorder;
    *inorder += 1;

    if root_idx + 1 < end {
        let right_ = traverse_recursion(
            nodes, info, root_idx+1, end, depth+1,
            root_idx as u32, inorder, preorder);
        right = right_;
        nodes[root_idx].subtree_last = nodes[root_idx].subtree_last.max(
            nodes[right as usize].subtree_last);
    }

    info[root_idx].subtree_size = subtree_size as u32;
    info[root_idx].left = left as u32;
    info[root_idx].right = right as u32;

    return root_idx as u32
}


// norder partition by depth on pivot into three parts, like so
//      [ bottom left ][ top ][ bottom right ]
// where bottom left and right are the bottom subtrees with positioned to
// the left and right of the root node
fn stable_ternary_tree_partition(
        input: &[usize], output: &mut [usize],
        info: &[TraversalInfo], pivot_depth: u32, pivot_dfs: u32,
        start: usize, end: usize) -> (usize, usize) {

    // bottom left
    let mut bl = start;
    for i in &input[start..end] {
        if info[*i].depth > pivot_depth && info[*i].inorder < pivot_dfs {
            output[bl] = *i;
            bl += 1;
        }
    }

    // top
    let mut t = bl;
    for i in &input[start..end] {
        if info[*i].depth <= pivot_depth {
            output[t] = *i;
            t += 1;
        }
    }

    // bottom right
    let mut br = t;
    for i in &input[start..end] {
        if info[*i].depth > pivot_depth && info[*i].inorder > pivot_dfs {
            output[br] = *i;
            br += 1;
        }
    }

    assert!(br == end);

    return (bl, t);
}


// put nodes in van Emde Boas order
fn veb_order<T>(mut nodes: Vec<IntervalNode<T>>) -> (Vec<IntervalNode<T>>, usize)
        where T: Copy {

    // let now = Instant::now();
    // nodes.sort_unstable_by_key(|node| node.first);
    let mut veb_nodes = nodes.clone();
    let n = veb_nodes.len();

    if veb_nodes.is_empty() {
        return (veb_nodes, 0);
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

    let idxs: &mut [usize] = &mut ((0..n).collect::<Vec<usize>>());
    let tmp: &mut [usize] = &mut vec![0; n];

    // put in dfs order
    for i in &*idxs {
        tmp[info[*i].preorder as usize] = *i;
    }
    let (idxs, tmp) = (tmp, idxs);

    // let now = Instant::now();
    let root_idx = veb_order_recursion(
        idxs, tmp, &mut info, &mut nodes, 0, n, true, false, 0, max_depth);
    // eprintln!("computing order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[*j] = i;
    }

    // put nodes in vEB order in a temporary vector
    for i in 0..n {
        veb_nodes[i] = nodes[idxs[i]];

        if info[idxs[i]].simple {
            veb_nodes[i].left = info[idxs[i]].subtree_size as u32;
            veb_nodes[i].right = veb_nodes[i].left;
        } else {
            // update left and right pointers
            let left = info[idxs[i]].left as u32;
            veb_nodes[i].left = if left < u32::max_value() {
                revidx[left as usize] as u32
            } else {
                left
            };

            let right = info[idxs[i]].right as u32;
            veb_nodes[i].right = if right < u32::max_value() {
                revidx[right as usize] as u32
            } else {
                right
            };
        }
    }

    let root_idx = revidx[root_idx];

    // eprintln!("ordering: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    assert!(compute_tree_size(&veb_nodes, root_idx) == n);

    return (veb_nodes, root_idx);
}


// Traverse the tree and return the size, used as a basic sanity check.
fn compute_tree_size<T>(nodes: &[IntervalNode<T>], root_idx: usize) -> usize
        where T: Copy {

    let mut subtree_size = 1;

    if nodes[root_idx].left == nodes[root_idx].right {
        subtree_size = nodes[root_idx].right as usize;
    } else {
        let left = nodes[root_idx].left as usize;
        if left < u32::max_value() as usize {
            subtree_size += compute_tree_size(nodes, left);
        }

        let right = nodes[root_idx].right as usize;
        if right < u32::max_value() as usize {
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
fn veb_order_recursion<T>(
        idxs: &mut [usize], tmp: &mut [usize],
        info: &mut [TraversalInfo],
        nodes: &mut [IntervalNode<T>],
        start: usize, end: usize,
        childless: bool,
        parity: bool,
        min_depth: u32, max_depth: u32) -> usize where T: Copy {
    let n = (start..end).len();

    // small subtrees are put into sorted order and just searched through
    // linearly. There is a little trickiness to this because we have to
    // update the parent's child pointers and some other fields.
    if childless && info[idxs[start]].subtree_size <= SIMPLE_SUBTREE_CUTOFF {
        assert!(n == info[idxs[start]].subtree_size as usize);

        let old_root = idxs[start];

        idxs[start..end].sort_unstable_by_key(|i| info[*i].inorder);
        info[idxs[start]].simple = true;
        info[idxs[start]].subtree_size = info[old_root].subtree_size;
        nodes[idxs[start]].subtree_last = nodes[old_root].subtree_last;

        let parent = info[old_root as usize].parent;
        if parent < u32::MAX {
            if info[parent as usize].left == old_root as u32 {
                info[parent as usize].left = idxs[start] as u32;
            } else {
                assert!(info[parent as usize].right == old_root as u32);
                info[parent as usize].right = idxs[start] as u32;
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
    let pivot_dfs = info[idxs[start]].inorder;

    let (top_start, bottom_right_start) =
        stable_ternary_tree_partition(
            idxs, tmp, info, pivot_depth, pivot_dfs, start, end);

    // tmp is not partitioned, so swap pointers
    let (tmp, idxs) = (idxs, tmp);

    // recurse on top subtree
    let top_root_idx = veb_order_recursion(
        idxs, tmp, info, nodes, top_start, bottom_right_start, false,
        !parity, min_depth, pivot_depth);

    // find on recurse on subtrees in the bottom left partition and bottom right partition
    for (part_start, part_end) in &[(start, top_start), (bottom_right_start, end)] {
        let bottom_subtree_depth = pivot_depth + 1;
        let mut i = *part_start;
        while i < *part_end {
            assert!(info[idxs[i]].depth == bottom_subtree_depth);
            let mut j = i+1;
            let mut subtree_max_depth = info[idxs[i]].depth;
            while j < *part_end && info[idxs[j]].depth != bottom_subtree_depth {
                assert!(info[idxs[j]].depth > bottom_subtree_depth);
                if info[idxs[j]].depth > subtree_max_depth {
                    subtree_max_depth = info[idxs[j]].depth;
                }
                j += 1;
            }
            veb_order_recursion(
                idxs, tmp, info, nodes, i, j, childless, !parity,
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
fn radix_sort_nodes<T>(nodes: &mut [IntervalNode<T>], tmp: &mut [IntervalNode<T>])
        where T: Copy {

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

