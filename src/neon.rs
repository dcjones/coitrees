use std::arch::aarch64::*;

use super::interval::{GenericInterval, IntWithMax, Interval, IntervalTree, SortedQuerent};
use std::cmp::{max, Ordering};
use std::fmt::Debug;
use std::marker::Copy;
use std::mem::transmute;

#[allow(non_camel_case_types)]
type i32x4 = int32x4_t;

const LANE_SIZE: usize = 4;

// Small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees. Performance isn't
// super sensitive, but is worse with a very small or very large number.
const SIMPLE_SUBTREE_CUTOFF: usize = 8;

// Very dense subtrees in which we probably intersect most of the intervals
// are more efficient to query linearly. When the expected proportion of hits
// is a above this number it becomes a simple subtree.
const SIMPLE_SUBTREE_DENSITY_CUTOFF: f32 = 0.2;

/// Node in the interval tree. Each node holds a chunk of 8 intervals.
#[derive(Clone)]
struct IntervalNode<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    // subtree interval
    subtree_last: i32,

    firsts: i32x4,
    lasts: i32x4,
    metadata: [T; LANE_SIZE],

    // when this is the root of a simple subtree, left == right is the size
    // of the subtree, otherwise they are left, right child pointers.
    left: I,
    right: I,
}

impl<T, I> IntervalNode<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    fn new(intervals: &[Interval<T>; LANE_SIZE]) -> Self
    where
        T: Copy,
    {
        assert!(!intervals.is_empty() && intervals.len() <= LANE_SIZE);

        unsafe {
            let max_last = intervals
                .iter()
                .max_by_key(|x| x.last)
                .map(|x| x.last)
                .unwrap();

            // the firsts and lasts are adjust by 1 here because AVX has an
            // instruction for > but not >=.
            // WARN: may not need to minus 1 <03-01-23, Yangyang Li yangyang.li@northwestern.edu>
            let firsts = vld1q_s32(
                [
                    intervals[0].first - 1,
                    intervals[1].first - 1,
                    intervals[2].first - 1,
                    intervals[3].first - 1,
                ]
                .as_ptr(),
            );

            let lasts = vld1q_s32(
                [
                    intervals[0].last + 1,
                    intervals[1].last + 1,
                    intervals[2].last + 1,
                    intervals[3].last + 1,
                ]
                .as_ptr(),
            );

            // NOTE: the order is reversed from first and lasts data <Yangyang Li yangyang.li@northwestern.edu>
            let metadata: [T; LANE_SIZE] = [
                intervals[0].metadata,
                intervals[1].metadata,
                intervals[2].metadata,
                intervals[3].metadata,
            ];

            Self {
                subtree_last: max_last,
                firsts,
                lasts,
                metadata,
                left: I::MAX,
                right: I::MAX,
            }
        }
    }

    fn first(&self, i: usize) -> i32 {
        assert!(i < LANE_SIZE);
        (unsafe { transmute::<&i32x4, &[i32; LANE_SIZE]>(&self.firsts) }[i])
    }

    fn last(&self, i: usize) -> i32 {
        assert!(i < LANE_SIZE);
        (unsafe { transmute::<&i32x4, &[i32; LANE_SIZE]>(&self.lasts) }[i])
    }

    fn min_first(&self) -> i32 {
        unsafe { vgetq_lane_s32(self.firsts, 0) + 1 }
    }

    fn max_last(&self) -> i32 {
        // faster ways to do this, but doesn't really matter
        unsafe {
            max(
                max(vgetq_lane_s32(self.lasts, 0), vgetq_lane_s32(self.lasts, 1)),
                max(vgetq_lane_s32(self.lasts, 2), vgetq_lane_s32(self.lasts, 3)),
            ) - 1
        }
    }

    #[inline(always)]
    fn query_count_chunk(&self, query_first: i32x4, query_last: i32x4) -> usize {
        unsafe {
            let cmp1 = vcgtq_s32(query_last, self.firsts);
            let cmp2 = vcgtq_s32(self.lasts, query_first);
            let cmp = vandq_u32(cmp1, cmp2);
            count_bits(cmp) / 32
        }
    }

    fn query_chunk_firsts<'a, F>(&'a self, query_first: i32x4, query_last: i32x4, mut visit: F)
    where
        F: FnMut(i32, i32, &'a T),
    {
        let cmp: u64 = unsafe {
            let cmp1 = vcgtq_s32(query_last, self.firsts);
            let cmp2 = vcgtq_s32(self.firsts, query_first);
            let cmp = vmovn_u32(vandq_u32(cmp1, cmp2));
            transmute(cmp)
        };

        let masks = [0xffff, 0xffff << 16, 0xffff << 32, 0xffff << 48];

        if cmp & masks[0] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 0) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 0) - 1 },
                &self.metadata[0],
            );
        }

        if cmp & masks[1] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 1) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 1) - 1 },
                &self.metadata[1],
            );
        }

        if cmp & masks[2] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 2) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 2) - 1 },
                &self.metadata[2],
            );
        }

        if cmp & masks[3] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 3) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 3) - 1 },
                &self.metadata[3],
            );
        }
    }

    fn query_chunk_metadata<'a, F>(&'a self, query_first: i32x4, query_last: i32x4, mut visit: F)
    where
        F: FnMut(i32, i32, &'a T),
    {
        let cmp: u64 = unsafe {
            let cmp1 = vcgtq_s32(query_last, self.firsts);
            let cmp2 = vcgtq_s32(self.lasts, query_first);
            let cmp = vmovn_u32(vandq_u32(cmp1, cmp2));
            transmute(cmp)
        };

        let masks = [0xffff, 0xffff << 16, 0xffff << 32, 0xffff << 48];

        if cmp & masks[0] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 0) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 0) - 1 },
                &self.metadata[0],
            );
        }
        if cmp & masks[1] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 1) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 1) - 1 },
                &self.metadata[1],
            );
        }
        if cmp & masks[2] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 2) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 2) - 1 },
                &self.metadata[2],
            );
        }
        if cmp & masks[3] != 0 {
            visit(
                unsafe { vgetq_lane_s32(self.firsts, 3) + 1 },
                unsafe { vgetq_lane_s32(self.lasts, 3) - 1 },
                &self.metadata[3],
            );
        }
    }

    fn query_chunk<F>(&self, query_first: i32x4, query_last: i32x4, mut visit: F)
    where
        F: FnMut(i32, i32),
    {
        let cmp: u64 = unsafe {
            let cmp1 = vcgtq_s32(query_last, self.firsts);
            let cmp2 = vcgtq_s32(self.lasts, query_first);
            let cmp = vmovn_u32(vandq_u32(cmp1, cmp2));
            transmute(cmp)
        };

        let masks = [0xffff, 0xffff << 16, 0xffff << 32, 0xffff << 48];

        if cmp & masks[0] != 0 {
            visit(unsafe { vgetq_lane_s32(self.firsts, 0) + 1 }, unsafe {
                vgetq_lane_s32(self.lasts, 0) - 1
            });
        }
        if cmp & masks[1] != 0 {
            visit(unsafe { vgetq_lane_s32(self.firsts, 1) + 1 }, unsafe {
                vgetq_lane_s32(self.lasts, 1) - 1
            });
        }
        if cmp & masks[2] != 0 {
            visit(unsafe { vgetq_lane_s32(self.firsts, 2) + 1 }, unsafe {
                vgetq_lane_s32(self.lasts, 2) - 1
            });
        }
        if cmp & masks[3] != 0 {
            visit(unsafe { vgetq_lane_s32(self.firsts, 3) + 1 }, unsafe {
                vgetq_lane_s32(self.lasts, 3) - 1
            });
        }
    }
}

fn count_bits(bits: uint32x4_t) -> usize {
    unsafe {
        let t2 = vreinterpretq_u8_u32(bits);
        let t3 = vcntq_u8(t2);
        let sum = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(t3)));
        vgetq_lane_u64(sum, 0) as usize + vgetq_lane_u64(sum, 1) as usize
    }
}

/// NeonCOITree data structure. A representation of a static set of intervals with
/// associated metadata, enabling fast overlap and coverage queries.
///
/// The index type `I` is a typically `usize`, but can be `u32` or `u16`.
/// It's slightly more efficient to use a smalled index type, assuming there are
/// fewer than I::MAX-1 intervals to store.
#[derive(Clone)]
pub struct NeonCOITree<T, I>
where
    T: Clone,
    I: IntWithMax,
{
    nodes: Vec<IntervalNode<T, I>>,
    len: usize,
    root_idx: usize,
    height: usize,
}

impl<T, I> NeonCOITree<T, I>
where
    T: Default + Copy + Clone,
    I: IntWithMax,
{
    fn chunk_intervals(intervals: Vec<Interval<T>>) -> Vec<IntervalNode<T, I>>
    where
        T: Copy,
    {
        let n = intervals.len();

        let num_chunks = (n / LANE_SIZE) + (n % LANE_SIZE != 0) as usize;

        let mut nodes: Vec<IntervalNode<T, I>> = Vec::with_capacity(num_chunks);

        // chunk initializer
        let mut chunk_init: [Interval<T>; LANE_SIZE] = [Interval {
            first: i32::MAX,
            last: i32::MIN,
            metadata: T::default(),
        }; LANE_SIZE];

        for chunk in intervals.chunks(LANE_SIZE) {
            for (j, interval) in chunk.iter().enumerate() {
                chunk_init[j] = *interval;
            }

            chunk_init.iter_mut().skip(chunk.len()).for_each(|item| {
                *item = Interval {
                    first: i32::MAX,
                    last: i32::MIN,
                    metadata: T::default(),
                };
            });

            nodes.push(IntervalNode::new(&chunk_init));
        }

        nodes
    }
}

impl<'a, T, I> IntervalTree<'a> for NeonCOITree<T, I>
where
    T: Default + Copy + Clone + 'a,
    I: IntWithMax + 'a,
{
    type Metadata = T;
    type Index = I;
    type Item = Interval<&'a T>;
    type Iter = NeonCOITreeIterator<'a, T, I>;

    fn new<'b, U, V>(intervals: U) -> NeonCOITree<T, I>
    where
        U: IntoIterator<Item = &'b V>,
        V: GenericInterval<T> + 'b,
    {
        let mut intervals: Vec<_> = intervals
            .into_iter()
            .map(|interval| {
                Interval::new(
                    interval.first(),
                    interval.last(),
                    interval.metadata().clone(),
                )
            })
            .collect();

        if intervals.len() >= (I::MAX).to_usize() {
            panic!("NeonCOITree construction failed: more intervals than index type can enumerate")
        }

        let n = intervals.len();
        intervals.sort_unstable_by_key(|interval| (interval.first, interval.last));
        let nodes = Self::chunk_intervals(intervals);

        let (nodes, root_idx, height) = veb_order(nodes);
        NeonCOITree {
            nodes,
            len: n,
            root_idx,
            height,
        }
    }

    /// Number of intervals in the set.
    fn len(&self) -> usize {
        self.len
    }

    /// True iff the set is empty.
    fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    // /// Find intervals in the set overlaping the query `[first, last]` and call `visit` on every overlapping node
    fn query<F>(&'a self, first: i32, last: i32, mut visit: F)
    where
        F: FnMut(&Interval<&'a T>),
    {
        let (firstv, lastv) = unsafe { (vdupq_n_s32(first), vdupq_n_s32(last)) };

        if !self.is_empty() {
            query_recursion(
                &self.nodes,
                self.root_idx,
                first,
                last,
                firstv,
                lastv,
                &mut visit,
            );
        }
    }

    /// Count the number of intervals in the set overlapping the query `[first, last]`.
    fn query_count(&self, first: i32, last: i32) -> usize {
        let (firstv, lastv) = unsafe { (vdupq_n_s32(first), vdupq_n_s32(last)) };

        if !self.is_empty() {
            query_recursion_count(&self.nodes, self.root_idx, first, last, firstv, lastv)
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

        let (firstv, lastv) = unsafe { (vdupq_n_s32(first), vdupq_n_s32(last)) };

        let (mut uncov_len, last_cov, count) = coverage_recursion(
            &self.nodes,
            self.root_idx,
            first,
            last,
            firstv,
            lastv,
            first - 1,
        );

        if last_cov < last {
            uncov_len += last - last_cov;
        }

        let cov = ((last - first + 1) as usize) - (uncov_len as usize);

        (count, cov)
    }

    /// Iterate through the interval set in sorted order by interval start position.
    fn iter(&self) -> NeonCOITreeIterator<T, I> {
        let mut i = self.root_idx;
        let mut stack: Vec<usize> = Vec::with_capacity(self.height);
        while i < self.nodes.len()
            && self.nodes[i].left != I::MAX
            && self.nodes[i].left != self.nodes[i].right
        {
            stack.push(i);
            i = self.nodes[i].left.to_usize();
        }

        NeonCOITreeIterator {
            nodes: &self.nodes,
            len: self.len,
            i,
            j: 0,
            count: 0,
            stack,
        }
    }
}

impl<'a, T, I> IntoIterator for &'a NeonCOITree<T, I>
where
    T: Default + Copy + Clone,
    I: IntWithMax,
{
    type Item = Interval<&'a T>;
    type IntoIter = NeonCOITreeIterator<'a, T, I>;

    fn into_iter(self) -> NeonCOITreeIterator<'a, T, I> {
        return self.iter();
    }
}

/// Iterate through nodes in a tree in sorted order by interval start position.
pub struct NeonCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    nodes: &'a Vec<IntervalNode<T, I>>,
    len: usize,   // total number of intervals
    i: usize,     // current node
    j: usize,     // offset into the current chunk
    count: usize, // number generated so far
    stack: Vec<usize>,
}

impl<'a, T, I> Iterator for NeonCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    type Item = Interval<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i * LANE_SIZE + self.j >= self.len {
            return None;
        }

        let node = &self.nodes[self.i];
        if self.j < LANE_SIZE {
            let ret = Some(Interval {
                first: node.first(self.j),
                last: node.last(self.j),
                metadata: &node.metadata[self.j],
            });
            self.count += 1;
            self.j += 1;
            return ret;
        }

        // find the next node
        self.j = 0;
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

        let node = &self.nodes[self.i];
        self.count += 1;
        Some(Interval {
            first: node.first(self.j),
            last: node.last(self.j),
            metadata: &node.metadata[self.j],
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len_remaining = self.len - self.count;
        (len_remaining, Some(len_remaining))
    }
}

impl<'a, T, I> ExactSizeIterator for NeonCOITreeIterator<'a, T, I>
where
    T: Clone,
    I: IntWithMax,
{
    fn len(&self) -> usize {
        self.len - self.count
    }
}

// // Recursively count overlaps between the tree specified by `nodes` and a
// // query interval specified by `first`, `last`.
fn query_recursion<'a, T, I, F>(
    nodes: &'a [IntervalNode<T, I>],
    root_idx: usize,
    first: i32,
    last: i32,
    firstv: i32x4,
    lastv: i32x4,
    visit: &mut F,
) where
    T: Clone,
    I: IntWithMax,
    F: FnMut(&Interval<&'a T>),
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.min_first() {
                break;
            }

            node.query_chunk_metadata(firstv, lastv, |first_hit, last_hit, metadata| {
                visit(&Interval {
                    first: first_hit,
                    last: last_hit,
                    metadata,
                });
            });
        }
    } else {
        node.query_chunk_metadata(firstv, lastv, |first_hit, last_hit, metadata| {
            visit(&Interval {
                first: first_hit,
                last: last_hit,
                metadata,
            });
        });

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                query_recursion(nodes, left, first, last, firstv, lastv, visit);
            }
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.min_first(), nodes[right].subtree_last, first, last) {
                query_recursion(nodes, right, first, last, firstv, lastv, visit);
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
    firstv: i32x4,
    lastv: i32x4,
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
            if last < node.min_first() {
                break;
            } else {
                count += node.query_count_chunk(firstv, lastv);
            }
        }
        count
    } else {
        let mut count = 0;
        count += node.query_count_chunk(firstv, lastv);

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                count += query_recursion_count(nodes, left, first, last, firstv, lastv);
            }
        }

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.min_first(), nodes[right].subtree_last, first, last) {
                count += query_recursion_count(nodes, right, first, last, firstv, lastv);
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
    firstv: i32x4,
    lastv: i32x4,
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
            if last < node.min_first() {
                break;
            }

            node.query_chunk(firstv, lastv, |first_hit, last_hit| {
                if first_hit > last_cov {
                    uncov_len += first_hit - (last_cov + 1);
                }
                last_cov = last_cov.max(last_hit);
                count += 1;
            });
        }
        (uncov_len, last_cov, count)
    } else {
        let mut uncov_len = 0;
        let mut count = 0;

        if node.left < I::MAX {
            let left = node.left.to_usize();
            if nodes[left].subtree_last >= first {
                let (left_uncov_len, left_last_cov, left_count) =
                    coverage_recursion(nodes, left, first, last, firstv, lastv, last_cov);
                last_cov = left_last_cov;
                uncov_len += left_uncov_len;
                count += left_count;
            }
        }

        node.query_chunk(firstv, lastv, |first_hit, last_hit| {
            if first_hit > last_cov {
                uncov_len += first_hit - (last_cov + 1);
            }
            last_cov = last_cov.max(last_hit);
            count += 1;
        });

        if node.right < I::MAX {
            let right = node.right.to_usize();
            if overlaps(node.min_first(), nodes[right].subtree_last, first, last) {
                let (right_uncov_len, right_last_cov, right_count) =
                    coverage_recursion(nodes, right, first, last, firstv, lastv, last_cov);
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
pub struct NeonSortedQuerent<'a, T, I>
where
    T: Default + Clone,
    I: IntWithMax,
{
    tree: &'a NeonCOITree<T, I>,
    prev_first: i32,
    prev_last: i32,
    overlapping_intervals: Vec<Interval<&'a T>>,
}

impl<'a, T, I> SortedQuerent<'a> for NeonSortedQuerent<'a, T, I>
where
    T: Default + Copy + Clone,
    I: IntWithMax,
{
    type Metadata = T;
    type Index = I;
    type Item = Interval<&'a T>;
    type Iter = NeonCOITreeIterator<'a, T, I>;
    type Tree = NeonCOITree<T, I>;

    /// Construct a new `SortedQuerent` to perform a sequence. queries.
    fn new(tree: &'a NeonCOITree<T, I>) -> NeonSortedQuerent<'a, T, I> {
        let overlapping_intervals: Vec<Interval<&'a T>> = Vec::new();
        NeonSortedQuerent {
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
        F: FnMut(&Interval<&'a T>),
    {
        if self.tree.is_empty() {
            return;
        }

        // not overlaping or preceding
        if first < self.prev_first || first > self.prev_last {
            // no overlap with previous query. have to resort to regular query strategy
            self.overlapping_intervals.clear();
            self.tree
                .query(first, last, |node| self.overlapping_intervals.push(*node));
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
                let qa = self.prev_last + 1 - 2; // -2 accounts for the adjustment made in the chunk
                let qb = last;

                let (qav, qbv) = unsafe { (vdupq_n_s32(qa), vdupq_n_s32(qb)) };

                sorted_querent_query_firsts(
                    &self.tree.nodes,
                    self.tree.root_idx,
                    qa,
                    qb,
                    qav,
                    qbv,
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
    firstv: i32x4,
    lastv: i32x4,
    overlapping_intervals: &mut Vec<Interval<&'a T>>,
) where
    T: Clone,
    I: IntWithMax,
{
    let node = &nodes[root_idx];

    if node.left == node.right {
        // simple subtree
        for node in &nodes[root_idx..root_idx + node.right.to_usize()] {
            if last < node.min_first() {
                break;
            }

            node.query_chunk_firsts(firstv, lastv, |first_hit, last_hit, metadata| {
                overlapping_intervals.push(Interval {
                    first: first_hit,
                    last: last_hit,
                    metadata,
                });
            })
        }
    } else {
        node.query_chunk_firsts(firstv, lastv, |first_hit, last_hit, metadata| {
            overlapping_intervals.push(Interval {
                first: first_hit,
                last: last_hit,
                metadata,
            });
        });

        if node.left < I::MAX && first <= node.min_first() {
            let left = node.left.to_usize();
            sorted_querent_query_firsts(
                nodes,
                left,
                first,
                last,
                firstv,
                lastv,
                overlapping_intervals,
            );
        }

        if node.right < I::MAX && last >= node.min_first() {
            let right = node.right.to_usize();
            sorted_querent_query_firsts(
                nodes,
                right,
                first,
                last,
                firstv,
                lastv,
                overlapping_intervals,
            );
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

    let mut subtree_first = nodes[root_idx].min_first();

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

    let expected_hits = ((nodes[root_idx].max_last() - nodes[root_idx].min_first() + 1) as f32
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
    (0..n).for_each(|i| idxs[i] = I::from_usize(i));

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metadata() {
        let intervals = (0..10)
            .map(|i| Interval::new(i, i + 1, i as usize))
            .collect::<Vec<Interval<usize>>>();

        let tree = NeonCOITree::<usize, u32>::new(&intervals);

        let mut metadata = Vec::new();
        tree.query(0, 3, |node| metadata.push(node.clone()));

        println!("{:?}", metadata);
    }
}
