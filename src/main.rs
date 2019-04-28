
use std::cmp::max;
use std::error::Error;
use std::str;
use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader};

extern crate fnv;
use fnv::FnvHashMap;

extern crate libc;

type GenericError = Box<Error>;


// small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees.
const SIMPLE_SUBTREE_CUTOFF: usize = 32;


#[derive(Copy, Clone, Debug)]
struct IntervalNode<T> where T: std::marker::Copy
 {
    // subtree interval
    subtree_last: i32,

    // interval
    first: i32,
    last: i32,

    // when this is the root of a simple subtree, left == right is the size
    // of the subtree, otherwise they are left, right childe pointers.
    left: u32,
    right: u32,

    metadata: T,
}


struct COITree<T>  where T: std::marker::Copy {
    nodes: Vec<IntervalNode<T>>
}


impl<T> COITree<T> where T: std::marker::Copy {
    pub fn new(nodes: Vec<IntervalNode<T>>) -> COITree<T> {
        return COITree { nodes: veb_order(nodes) };
    }
}


// Used by `traverse` to keep record tree metadata.
#[derive(Copy, Clone, Debug)]
struct TraversalInfo {
    depth: usize,
    dfs: usize,
    subtree_size: usize,
    left: usize,
    right: usize,
    simple: bool // set by veb_order_recursion
}


// dfs traversal of an implicit bst computing dfs number, node depth, subtree
// size, and left and right pointers.
fn traverse<T>(nodes: &mut [IntervalNode<T>]) -> Vec<TraversalInfo>
        where T: std::marker::Copy {
    let n = nodes.len();
    let mut info = vec![TraversalInfo{
        depth: 0, dfs: 0, subtree_size: 0, left: 0, right: 0, simple: false}; n];
    let mut dfs = 0;
    traverse_recursion(nodes, &mut info, 0, n, 0, &mut dfs);

    return info;
}


// The recursive part of the `traverse` function.
fn traverse_recursion<T>(
        nodes: &mut [IntervalNode<T>],
        info: &mut [TraversalInfo],
        start: usize,
        end: usize,
        depth: usize,
        dfs: &mut usize) -> (usize, usize)
        where T: std::marker::Copy {

    if start >= end {
        return (usize::max_value(), 0);
    }

    let root_idx = start + (end - start) / 2;
    info[root_idx].depth = depth;
    info[root_idx].dfs = *dfs;
    nodes[root_idx].subtree_last = nodes[root_idx].last;
    *dfs += 1;
    let mut left = usize::max_value();
    let mut right = usize::max_value();
    let mut subtree_size = 1;

    if root_idx > start {
        let (left_, left_subtree_size) =
            traverse_recursion(nodes, info, start, root_idx, depth+1, dfs);
        left = left_;
        subtree_size += left_subtree_size;
        nodes[root_idx].subtree_last = max(nodes[root_idx].subtree_last, nodes[left].subtree_last);
    }

    if root_idx + 1 < end {
        let (right_, right_subtree_size) =
            traverse_recursion(nodes, info, root_idx+1, end, depth+1, dfs);
        right = right_;
        subtree_size += right_subtree_size;
        nodes[root_idx].subtree_last = max(nodes[root_idx].subtree_last, nodes[right].subtree_last);
    }

    info[root_idx].subtree_size = subtree_size;
    info[root_idx].left = left;
    info[root_idx].right = right;

    return (root_idx, subtree_size)
}

// partition (in the quicksort sense) indexes according do the corresponding depths
// while retaining relative order.
fn stable_partition_by_depth(
        input: &[usize], output: &mut [usize],
        info: &[TraversalInfo], pivot: usize,
        start: usize, end: usize) -> usize {

    let mut l = start;
    for i in start..end {
        if info[input[i]].depth <= pivot {
            output[l] = input[i];
            l += 1;
        }
    }

    let mut r = l;
    for i in start..end {
        if info[input[i]].depth > pivot {
            output[r] = input[i];
            r += 1;
        }
    }

    return l - start;
}


// recursively reorder indexes to put it in vEB order. Called by `veb_order`
//   idxs: current permutation
//   tmp: temporary space of equal length to idxs
//   nodes: the interval nodes (in sorted order)
//   start, end: slice within idxs to be reordered
//   childless: true if this slice is a proper subtree and has no children below it
//   parity: true if idxs and tmp are swapped and need to be copied back,
//   min_depth, max_depth: depth extreme of the start..end slice
fn veb_order_recursion<T>(
        idxs: &mut [usize], tmp: &mut [usize],
        info: &mut [TraversalInfo],
        nodes: &[IntervalNode<T>],
        start: usize, end: usize,
        childless: bool,
        parity: bool,
        min_depth: usize, max_depth: usize)
        where T: std::marker::Copy {
    let n = (start..end).len();

    // small subtrees are just put in sorted order
    if childless && info[idxs[start]].subtree_size <= SIMPLE_SUBTREE_CUTOFF {
        assert!(n == info[idxs[start]].subtree_size);
        info[idxs[start]].simple = true;

        if parity {
            tmp[start..end].copy_from_slice(&idxs[start..end]);
        }
        return;
    }

    // very small trees are already in order
    if n <= 3 {
        if parity {
            tmp[start..end].copy_from_slice(&idxs[start..end]);
        }
        return;
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let top_size = stable_partition_by_depth(
        idxs, tmp, info, pivot_depth, start, end);

    // tmp is now partitioned by depth, so swap pointers
    let (tmp, idxs) = (idxs, tmp);

    // recurse on top subtree
    veb_order_recursion(
        idxs, tmp, info, nodes, start,
        start + top_size, false, !parity, min_depth, pivot_depth);

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
            idxs, tmp, info, nodes, i, j, childless, !parity,
            bottom_subtree_depth, subtree_max_depth);
        i = j;
    }
}


// Traverse the tree and return the size, used as a basic sanity check.
fn compute_tree_size<T>(nodes: &[IntervalNode<T>], root_idx: usize) -> usize
        where T: std::marker::Copy {

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


// Simple two pass radix sort of 32bit integers (16bits at a time) to sort nodes
// on start position. tmp is temporary space for the first pass of equal length
// to nodes.
fn radix_sort_nodes<T>(nodes: &mut [IntervalNode<T>], tmp: &mut [IntervalNode<T>])
        where T: std::marker::Copy {

    let mut max_first = 0;
    for node in &*nodes {
        max_first = max(max_first, node.first);
    }

    let mut count = 0;
    let n = nodes.len();

    const R: usize = 16;
    const K: usize = 0xffff+1;
    const MASK: i32 = 0xffff;

    let mut shift: usize = 0;
    // let mut radix_counts: [u32; K] = [0; K];
    let mut radix_counts: Vec<u32> = vec![0; K];

    let mut from = nodes;
    let mut to   = tmp;
    while count < 32/R {
        for i in 0..K {
            radix_counts[i] = 0;
        }

        for i in 0..n {
            radix_counts[((from[i].first >> shift) & MASK) as usize] += 1;
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

        for i in 0..n {
            let radix = ((from[i].first >> shift) & MASK) as usize;
            to[radix_counts[radix] as usize] = from[i];
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



// put nodes in van Emde Boas order
fn veb_order<T>(mut nodes: Vec<IntervalNode<T>>) -> Vec<IntervalNode<T>>
        where T: std::marker::Copy {

    // let now = Instant::now();
    // nodes.sort_unstable_by_key(|node| node.first);
    let mut veb_nodes = nodes.clone();
    radix_sort_nodes(&mut nodes, &mut veb_nodes);
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

    let idxs: &mut [usize] = &mut ((0..nodes.len()).collect::<Vec<usize>>());
    let tmp: &mut [usize] = &mut vec![0; nodes.len()];

    // put in dfs order
    for i in &*idxs {
        tmp[info[*i].dfs] = *i;
    }
    let (idxs, tmp) = (tmp, idxs);

    // let now = Instant::now();
    veb_order_recursion(
        idxs, tmp, &mut info, &nodes, 0, nodes.len(), true, false, 0, max_depth);
    // eprintln!("computing order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    // let now = Instant::now();
    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[*j] = i;
    }

    // but nodes in vEB order in a temporary vector
    for i in 0..nodes.len() {
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

    // eprintln!("ordering: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    assert!(compute_tree_size(&veb_nodes, 0) == veb_nodes.len());

    return veb_nodes;
}


// True iff the two intervals overlap.
fn overlaps(first_a: i32, last_a: i32, first_b: i32, last_b: i32) -> bool {
    return first_a <= last_b && last_a >= first_b;
}


// Recursively count overlaps between the tree specified by `nodes` and a
// query interval specified by `first`, `last`.
fn query_recursion(
        nodes: &[IntervalNode<()>], root_idx: usize, first: i32, last: i32,
        count: &mut usize, overlap: &mut usize, visited: &mut usize) {

    let node = &nodes[root_idx];

    if node.left == node.right { // simple subtree
        for k in root_idx..root_idx + node.right as usize {
            let node = &nodes[k];
            if overlaps(node.first, node.last, first, last) {
                *count += 1;
                *visited += 1;
            }
        }
    } else {
        *visited += 1;
        if overlaps(node.first, node.last, first, last) {
            *count += 1;
        }

        let left = node.left as usize;
        if left < u32::max_value() as usize {
            if nodes[left].subtree_last >= first {
                query_recursion(nodes, left, first, last, count, overlap, visited);
            }
        }

        let right = node.right as usize;
        if right < u32::max_value() as usize {
            if overlaps(
                    node.first, nodes[right].subtree_last,
                    first, last) {
                query_recursion(nodes, right, first, last, count, overlap, visited);
            }
        }
    }
}

// super simple query which prints every overlap
fn query(tree: &COITree<()>, first: i32, last: i32) -> (usize, usize, usize) {
    let mut count = 0;
    let mut overlap = 0;
    let mut visited = 0;
    query_recursion(
        &tree.nodes, 0, first, last, &mut count, &mut overlap, &mut visited);
    return (count, overlap, visited);
}


fn parse_bed_line(line: &[u8]) -> (&str, &str, &str, i32, i32) {
    let mut p = 0;
    while p < line.len()-1 && line[p] != b'\t' {
        p += 1;
    }
    let seqname = unsafe {
        str::from_utf8_unchecked(&line[..p])
    };
    p += 1;
    let p0 = p;

    while p < line.len()-1 && line[p] != b'\t' {
        p += 1;
    }
    let first_str = unsafe {
        str::from_utf8_unchecked(&line[p0..p])
    };
    let first = first_str.parse::<i32>().unwrap();
    p += 1;
    let p0 = p;

    while p < line.len()-1 && line[p] != b'\t' {
        p += 1;
    }
    let last_str = unsafe {
        str::from_utf8_unchecked(&line[p0..p])
    };
    let last = last_str.parse::<i32>().unwrap() - 1;

    return (seqname, first_str, last_str, first, last);
}


// Read a bed file into a COITree
fn read_bed_file(path: &str) -> Result<FnvHashMap<String, COITree<()>>, GenericError> {
    let mut nodes = FnvHashMap::<String, Vec<IntervalNode<()>>>::default();

    let now = Instant::now();

    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line_count = 0;
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, _, _, first, last) =
            parse_bed_line(&line);

        let node_arr = if let Some(node_arr) = nodes.get_mut(seqname) {
            node_arr
        } else {
            nodes.entry(seqname.to_string()).or_insert(Vec::new())
        };

        node_arr.push(IntervalNode{
            first: first, last: last,
            subtree_last: last,
            left: u32::max_value(),
            right: u32::max_value(),
            metadata: ()});

        line_count += 1;
        line.clear();
    }

    eprintln!("reading bed: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("lines: {}", line_count);
    eprintln!("sequences: {}", nodes.len());

    let now = Instant::now();
    let mut trees = FnvHashMap::<String, COITree<()>>::default();
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(seqname_nodes));
    }
    eprintln!("veb_order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    return Ok(trees);
}


fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let tree = read_bed_file(filename_a)?;

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();

    let mut total_visits = 0;
    let mut total_count = 0;
    let now = Instant::now();

    // let stdout = io::stdout();
    // let mut out = stdout.lock();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, _first_str, _last_str, first, last) =
            parse_bed_line(&line);

        let mut count = 0;
        let mut overlap = 0;
        let mut visited = 0;
        if let Some(seqname_tree) = tree.get(seqname) {

            let count_overlap = query(&seqname_tree, first, last);
            count = count_overlap.0;
            overlap = count_overlap.1;
            visited = count_overlap.2;
        }

        // out.write(&line[..line.len()-1])?;
        // writeln!(out, "\t{}", count)?;

        // unfortunately printing in c is quite a bit faster than rust
        unsafe {
            let linelen = line.len();
            line[linelen-1] = b'\0';
            libc::printf(
                b"%s\t%u\n\0".as_ptr() as *const i8,
                line.as_ptr() as *const i8,
                count as u32);
        }

        total_visits += visited;
        total_count += count;

        line.clear();
    }

    eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("total overlaps: {}", total_count);
    eprintln!("total visits: {}", total_visits);

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

