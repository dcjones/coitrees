
use std::cmp::{max, Ord};
use std::default::Default;


// Basic interval type
pub struct Interval<T, M> {
    // inclusive interval
    pub first: T,
    pub last: T,

    // accompanying metadata
    pub metadata: M,
}


#[derive(Copy, Clone)]
struct IndexNode<T> {
    first: T,
    subtree_last: T,

    // child pointers. Whene left == right this is a leaf node
    // and left==right gives the index in data of the corresponding interval.
    // Where left!=right, this is an internal node with pointers into index.
    left: u32,
    right: u32,

    // index of first data element in this subtree
    subtree_max_data_idx: u32,
}


impl<T> IndexNode<T>
        where T: Default {
    fn new() -> IndexNode<T> {
        return IndexNode{
            first: T::default(),
            subtree_last: T::default(),
            left: u32::max_value(),
            right: u32::max_value(),
            subtree_max_data_idx: u32::max_value() };
    }
}


pub struct COITree<T, M> {

    // binary search tree in van Emde Boas layout
    index: Vec<IndexNode<T>>,

    // interval in sorted order (by first)
    data: Vec<Interval<T,M>>,
}


impl<T, M> COITree<T, M>
        where T: Ord + Default + Clone + Copy {
    pub fn new(mut data: Vec<Interval<T, M>>) -> COITree<T, M> {
        // sort intervals if needed
        let mut data_is_presorted = true;
        for i in 1..data.len() {
            if data[i].first < data[i-1].first {
                data_is_presorted = false;
                break;
            }
        }

        if !data_is_presorted {
            data.sort_unstable_by_key(|interval| interval.first);
        }

        return COITree {
            index: index_intervals(&data),
            data: data
        };
    }
}


// build a vEB layout index for `data`. Assumed `data` is sorted on `first`.
fn index_intervals<T, M>(data: &[Interval<T, M>]) -> Vec<IndexNode<T>>
        where T: Ord + Default + Copy {

    let (index, depths) = traverse(data);

    // reorder index vector
    let max_depth = *depths.iter().max().unwrap();

    // a map from vEB order to the current dfs order of the index
    let idxmap: &mut [usize] = &mut ((0..index.len()).collect::<Vec<usize>>());

    // temporary space used when constructing idxmap
    let tmp: &mut [usize] = &mut vec![0; index.len()];

    veb_order_recursion(
        idxmap, tmp, depths.as_slice(), 0, index.len(), false, 0, max_depth);

    // build dfs -> veb index map
    let revidxmap = tmp;
    for (i, j) in idxmap.iter().enumerate() {
        revidxmap[*j] = i;
    }

    // put index in VEB order
    let mut veb_index = index.clone();
    for i in 0..index.len() {
        veb_index[i] = index[idxmap[i]];

        if veb_index[i].left != veb_index[i].right {
            veb_index[i].left = revidxmap[veb_index[i].left as usize] as u32;
            veb_index[i].right = revidxmap[veb_index[i].right as usize] as u32;
        }
    }

    return veb_index;
}


// traverse data in sorted order build up index
fn traverse<T, M>(data: &[Interval<T, M>])
        -> (Vec<IndexNode<T>>, Vec<u32>) where T: Ord + Default + Copy {

    let mut depths = Vec::with_capacity(2*data.len()-1);
    let mut index = Vec::with_capacity(2*data.len()-1);

    traverse_recursion(data, &mut index, &mut depths, 0, data.len(), 0);

    return (index, depths);
}


fn traverse_recursion<T, M>(
        data: &[Interval<T, M>], index: &mut Vec<IndexNode<T>>,
        depths: &mut Vec<u32>, start: usize, end: usize,
        depth: u32) -> u32 where T: Ord + Default + Copy {

    depths.push(depth);

    // leaf node
    assert!(start < end);
    if start + 1 == end {
        index.push(IndexNode{
            first: data[start].first,
            subtree_last: data[start].last,
            left: start as u32,
            right: start as u32,
            subtree_max_data_idx: start as u32,
        });

        return (index.len() - 1) as u32;
    }

    // internal node
    let mid = start + (end - start) / 2;
    index.push(IndexNode{
        first: data[mid].first,
        subtree_last: T::default(),
        left: u32::max_value(),
        right: u32::max_value(),
        subtree_max_data_idx: u32::max_value(),
    });
    let root_idx = index.len() - 1;

    // traverse left
    index[root_idx].left =
        traverse_recursion(data, index, depths, start, mid, depth + 1);
    index[root_idx].subtree_last =
        index[index[root_idx].left as usize].subtree_last;
    index[root_idx].subtree_max_data_idx =
        index[index[root_idx].left as usize].subtree_max_data_idx;

    // traverse right
    index[root_idx].right =
        traverse_recursion(data, index, depths, mid, end, depth + 1);
    index[root_idx].subtree_last = max(
        index[root_idx].subtree_last,
        index[index[root_idx].right as usize].subtree_last);
    index[root_idx].subtree_max_data_idx = max(
        index[root_idx].subtree_max_data_idx,
        index[index[root_idx].right as usize].subtree_max_data_idx);

    return root_idx as u32;
}


fn veb_order_recursion(
        idxmap: &mut [usize],
        tmp: &mut [usize],
        depths: &[u32],
        start: usize, end: usize, parity: bool,
        min_depth: u32, max_depth: u32) {

    let n = end - start;

    // for trees <= 3 dfs order = veb order, don't need to do anything
    if n <= 3 {
        if parity {
            tmp[start..end].copy_from_slice(&idxmap[start..end]);
        }
        return;
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let top_size = stable_partition_by_depth(
        idxmap, tmp, depths, pivot_depth, start, end);

    // tmp is now partitioned by depth, so swap pointers
    let (tmp, idxmap) = (idxmap, tmp);

    // recurse on top subtree
    veb_order_recursion(
        idxmap, tmp, depths, start, start + top_size,
        !parity, min_depth, pivot_depth);

    // find and recurse on bottom subtrees
    let bottom_subtree_depth = pivot_depth + 1;
    let mut i = start + top_size;
    while i < end {
        assert!(depths[idxmap[i]] == bottom_subtree_depth);
        let mut j = i+1;
        let mut subtree_max_depth = depths[idxmap[i]];
        while j < end && depths[idxmap[j]] != bottom_subtree_depth {
            assert!(depths[idxmap[j]]> bottom_subtree_depth);
            subtree_max_depth = max(subtree_max_depth, depths[idxmap[j]]);
            j += 1;
        }

        veb_order_recursion(
            idxmap, tmp, depths, i, j, !parity,
            bottom_subtree_depth, subtree_max_depth);
        i = j;
    }
}


// partition (in the quicksort sense) indexes according do the corresponding depths
// while retaining relative order.
fn stable_partition_by_depth(
        input: &[usize], output: &mut [usize],
        depths: &[u32], pivot: u32,
        start: usize, end: usize) -> usize {
    let mut l = start;
    for i in start..end {
        if depths[input[i]] <= pivot {
            output[l] = input[i];
            l += 1;
        }
    }

    let mut r = l;
    for i in start..end {
        if depths[input[i]] > pivot {
            output[r] = input[i];
            r += 1;
        }
    }

    return l - start;
}

impl<T, M> COITree<T, M> where T: Ord + Copy {
    pub fn count_overlaps(&self, first: T, last: T) -> usize
            where T: Ord + Copy {
        let mut count = 0;

        count_overlaps_recursion(
            self, 0, first, last, 0, &mut count);

        return count;
    }
}


fn count_overlaps_recursion<T, M>(
        tree: &COITree<T, M>, root_idx: usize,
        first: T, last: T, mut min_data_idx: usize,
        count: &mut usize) -> usize
        where T: Ord + Copy {
    let node = &tree.index[root_idx];

    // leaf node
    if node.left == node.right {
        let mut i = node.left as usize;
        let data = &tree.data;
        while i < data.len() &&
                overlaps(data[i].first, data[i].last, first, last) {
            i += 1;
        }
        *count += i - node.left as usize;
        min_data_idx = i + 1;
        return min_data_idx;
    }

    // internal node
    let left_node = &tree.index[node.left as usize];
    if left_node.subtree_last >= first &&
            (left_node.subtree_max_data_idx as usize) >= min_data_idx {
        min_data_idx = count_overlaps_recursion(
            tree, node.left as usize, first, last, min_data_idx, count);
    }

    let right_node = &tree.index[node.right as usize];
    if overlaps(node.first, right_node.subtree_last, first, last) &&
            (right_node.subtree_max_data_idx as usize) >= min_data_idx {
        min_data_idx = count_overlaps_recursion(
            tree, node.right as usize, first, last, min_data_idx, count);
    }

    return min_data_idx;
}


// True iff the two intervals overlap.
fn overlaps<T>(first_a: T, last_a: T, first_b: T, last_b: T) -> bool
        where T: Ord {
    return first_a <= last_b && last_a >= first_b;
}