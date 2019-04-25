
extern crate csv; // for parsing BED
use std::io;
use std::error::Error;
use std::cmp::{min, max};

type GenericError = Box<Error>;


#[derive(Copy, Clone, Debug)]
struct IntervalNode<T> where T: std::marker::Copy
 {
    first: i32,
    last: i32,

    subtree_first: i32,
    subtree_last: i32,

    left: Option<usize>,
    right: Option<usize>,

    metadata: T,
}


struct COITree<T>  where T: std::marker::Copy {
    nodes: Vec<IntervalNode<T>>
}


impl<T> COITree<T> where T: std::marker::Copy {
    pub fn new(mut nodes: Vec<IntervalNode<T>>) -> COITree<T> {
        veb_order(&mut nodes);
        return COITree { nodes: nodes };
    }
}



#[derive(Copy, Clone, Debug)]
struct TraversalInfo {
    depth: usize,
    dfs: usize,
    left: Option<usize>,
    right: Option<usize>,
}


fn traverse(n: usize) -> Vec<TraversalInfo> {
    let mut info = vec![TraversalInfo{depth: 0, dfs: 0, left: None, right: None}; n];
    let mut dfs = 0;
    let n = info.len();
    traverse_recursion(&mut info, 0, n, 0, &mut dfs);

    return info;
}


// dfs traversal of an implicit bst computing dfs number and depth.
fn traverse_recursion(
        info: &mut [TraversalInfo],
        start: usize,
        end: usize,
        depth: usize,
        dfs: &mut usize) -> Option<usize> {
    if info.is_empty() {
        return None;
    }

    let root_idx = start + (end - start) / 2;
    info[root_idx].depth = depth;
    info[root_idx].dfs = *dfs;
    *dfs += 1;
    let mut left = None;
    let mut right = None;

    if root_idx > start {
        left = traverse_recursion(info, start, root_idx, depth+1, dfs);
    }

    if root_idx + 1 < end {
        right = traverse_recursion(info, root_idx+1, end, depth+1, dfs);
    }

    info[root_idx].left = left;
    info[root_idx].right = right;

    return Some(root_idx);
}

// partition (in the quicksort sense) indexes according do the corresponding depths
// while retaining relative order.
fn stable_partition_by_depth(
        idxs: &mut [usize], tmp: &mut [usize],
        info: &[TraversalInfo], pivot: usize,
        start: usize, end: usize) -> usize {
    // count elements <= pivot
    let mut left_size = 0;
    for i in start..end {
        if info[idxs[i]].depth <= pivot {
            left_size += 1;
        }
    }

    let mut l = start;
    let mut r = start + left_size;
    for i in start..end {
        if info[idxs[i]].depth <= pivot {
            tmp[l] = idxs[i];
            l += 1
        } else {
            tmp[r] = idxs[i];
            r += 1;
        }
    }

    idxs[start..end].copy_from_slice(&tmp[start..end]);
    return left_size;
}


// recursively reorder indexes to put it in vEB order.
fn veb_order_recursion(
        idxs: &mut [usize], tmp: &mut [usize],
        info: &[TraversalInfo],
        start: usize, end: usize,
        min_depth: usize, max_depth: usize) {
    let n = (start..end).len();
    if min_depth == max_depth {
        // TODO: probably could relax this and return if size is <= 3
        assert!(n == 1);
        return;
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let top_size = stable_partition_by_depth(
        idxs, tmp, info, pivot_depth, start, end);

    // recurse on top subtree
    veb_order_recursion(
        idxs, tmp, info, start, start + top_size, min_depth, pivot_depth);

    // find and recurse on bottom subtrees
    let bottom_subtree_depth = pivot_depth + 1;
    let mut i = start + top_size;
    while i < end {
        assert!(info[idxs[i]].depth == bottom_subtree_depth);
        let mut j = i+1;
        let mut subtree_max_depth = info[idxs[i]].depth;
        while j < end && info[idxs[j]].depth != bottom_subtree_depth {
            assert!(info[idxs[j]].depth > bottom_subtree_depth);
            if info[idxs[j]].depth > subtree_max_depth {
                subtree_max_depth = info[idxs[j]].depth;
            }
            j += 1;
        }
        veb_order_recursion(
            idxs, tmp, info, i, j, bottom_subtree_depth, subtree_max_depth);
        i = j;
    }
}


// traverse the tree filling the subtree_first and subtree_last fields.
fn compute_subtree_sizes<T>(nodes: &mut [IntervalNode<T>], root_idx: usize)
        where T: std::marker::Copy {
    let mut subtree_first = nodes[root_idx].first;
    let mut subtree_last = nodes[root_idx].last;

    if let Some(left) = nodes[root_idx].left {
        compute_subtree_sizes(nodes, left);
        subtree_first = min(subtree_first, nodes[left].subtree_first);
        subtree_last  = max(subtree_last, nodes[left].subtree_last);
    }

    if let Some(right) = nodes[root_idx].right {
        compute_subtree_sizes(nodes, right);
        subtree_first = min(subtree_first, nodes[right].subtree_first);
        subtree_last  = max(subtree_last, nodes[right].subtree_last);
    }

    nodes[root_idx].subtree_first = subtree_first;
    nodes[root_idx].subtree_last = subtree_last;
}


// put nodes in van Emde Boas order
fn veb_order<T>(nodes: &mut [IntervalNode<T>])
        where T: std::marker::Copy {
    println!("sorting nodes");
    nodes.sort_unstable_by_key(|node| node.first);
    println!("done.");

    println!("traverse");
    let info = traverse(nodes.len());
    println!("done.");

    let mut max_depth = 0;
    for info_i in &info {
        if info_i.depth > max_depth {
            max_depth = info_i.depth;
        }
    }

    let mut idxs: Vec<usize> = (0..nodes.len()).collect();
    let mut tmp: Vec<usize> = vec![0; nodes.len()];

    println!("sort on dfs");
    idxs.sort_by_key(|i| info[*i].dfs);
    println!("done");

    println!("veb_order_recursion");
    veb_order_recursion(&mut idxs, &mut tmp, &info, 0, nodes.len(), 0, max_depth);
    println!("done");

    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let mut revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[*j] = i;
    }

    // but nodes in vEB order in a temporary vector
    let mut veb_nodes = vec![nodes[0]; nodes.len()];
    for i in 0..nodes.len() {
        veb_nodes[i] = nodes[idxs[i]];

        // update left and right pointers
        veb_nodes[i].left = if let Some(left) = info[idxs[i]].left {
            Some(revidx[left])
        } else {
            None
        };

        veb_nodes[i].right = if let Some(right) = info[idxs[i]].right {
            Some(revidx[right])
        } else {
            None
        };
    }

    // copy reordered nodes back to the original vector
    nodes.copy_from_slice(&mut veb_nodes);

    compute_subtree_sizes(nodes, 0);
}


fn overlaps(first_a: i32, last_a: i32, first_b: i32, last_b: i32) -> bool {
    return first_a <= last_b && last_a >= first_b;
}


fn query_recursion(
        nodes: &[IntervalNode<()>], root_idx: usize, first: i32, last: i32,
        count: &mut usize, overlap: &mut usize) {
    if overlaps(nodes[root_idx].first, nodes[root_idx].last, first, last) {
        *count += 1;
        // *overlap +=
        //     (min(nodes[root_idx].last, last) -
        //     max(nodes[root_idx].first, first)) as usize;
    }

    if let Some(left) = nodes[root_idx].left {
        if overlaps(
                nodes[left].subtree_first, nodes[left].subtree_last,
                first, last) {
            query_recursion(nodes, left, first, last, count, overlap);
        }
    }

    if let Some(right) = nodes[root_idx].right {
        if overlaps(
                nodes[right].subtree_first, nodes[right].subtree_last,
                first, last) {
            query_recursion(nodes, right, first, last, count, overlap);
        }
    }
}

// super simple query which prints every overlap
fn query(tree: &COITree<()>, first: i32, last: i32) -> (usize, usize) {
    let mut count = 0;
    let mut overlap = 0;
    query_recursion(&tree.nodes, 0, first, last, &mut count, &mut overlap);
    return (count, overlap);
}


fn read_bed_file(path: &str) -> Result<COITree<()>, GenericError> {
    let mut nodes: Vec<IntervalNode<()>> = Vec::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(path)?;

    println!("reading bed");
    for result in rdr.records() {
        let record = result?;
        if record.len() < 3 {
            return Err(GenericError::from(
                io::Error::new(io::ErrorKind::Other, "Bed line too short")));
        }

        // TODO: handle sequence: (by building multiple trees in a hashmap)

        let first: i32 = record[1].parse().unwrap();
        let mut last: i32 = record[2].parse().unwrap();
        last -= 1; // bed is end-exclusive
        nodes.push(IntervalNode{
            first: first, last: last,
            subtree_first: first, subtree_last: last,
            left: None, right: None, metadata: ()});
    }
    println!("done.");

    return Ok(COITree::new(nodes));
}


fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let tree = read_bed_file(filename_a)?;

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(filename_b)?;

    for result in rdr.records() {
        let record = result?;
        if record.len() < 3 {
            return Err(GenericError::from(
                io::Error::new(io::ErrorKind::Other, "Bed line too short")));
        }

        let first: i32 = record[1].parse().unwrap();
        let mut last: i32 = record[2].parse().unwrap();
        last -= 1; // bed is end-exclusive

        let count_overlap = query(&tree, first, last);
        println!(
            "{}\t{}\t{}\t{}\t{}",
            &record[0], &record[1], &record[2], count_overlap.0, count_overlap.1);
    }

    return Ok(());
}


fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        println!("Must specify file name.");
        std::process::exit(1);
    }

    let result = query_bed_files(&args[1], &args[2]);
    if let Err(err) = result {
        println!("error: {}", err)
    }
}

