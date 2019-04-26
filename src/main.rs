
use std::cmp::{min, max};
use std::error::Error;
use std::io;
use std::str;
use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

extern crate fnv;
use fnv::FnvHashMap;

type GenericError = Box<Error>;


// small subtrees at the bottom of the tree are stored in sorted order
// This gives the upper bound on the size of such subtrees.
const SIMPLE_SUBTREE_CUTOFF: usize = 16;


#[derive(Copy, Clone, Debug)]
struct IntervalNode<T> where T: std::marker::Copy
 {
    // subtree interval
    subtree_first: i32,
    subtree_last: i32,

    // interval
    first: i32,
    last: i32,

    // when `simple` is false, these are pointers to children
    // when `simple` is true, left has no meaning, and right
    // is the number of nodes in the sequence.
    left: u32,
    right: u32,

    simple: bool,

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
    subtree_size: usize,
    left: usize,
    right: usize,
    simple: bool // set by veb_order_recursion
}


fn traverse(n: usize) -> Vec<TraversalInfo> {
    let mut info = vec![TraversalInfo{
        depth: 0, dfs: 0, subtree_size: 0, left: 0, right: 0, simple: false}; n];
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
        dfs: &mut usize) -> (usize, usize) {

    if start >= end {
        return (usize::max_value(), 0);
    }

    let root_idx = start + (end - start) / 2;
    info[root_idx].depth = depth;
    info[root_idx].dfs = *dfs;
    *dfs += 1;
    let mut left = usize::max_value();
    let mut right = usize::max_value();
    let mut subtree_size = 1;

    if root_idx > start {
        let (left_, left_subtree_size) =
            traverse_recursion(info, start, root_idx, depth+1, dfs);
        left = left_;
        subtree_size += left_subtree_size;
    }

    if root_idx + 1 < end {
        let (right_, right_subtree_size) =
            traverse_recursion(info, root_idx+1, end, depth+1, dfs);
        right = right_;
        subtree_size += right_subtree_size;
    }

    info[root_idx].subtree_size = subtree_size;
    info[root_idx].left = left;
    info[root_idx].right = right;

    return (root_idx, subtree_size)
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
fn veb_order_recursion<T>(
        idxs: &mut [usize], tmp: &mut [usize],
        info: &mut [TraversalInfo],
        nodes: &[IntervalNode<T>],
        start: usize, end: usize,
        childless: bool, // true if no nodes below the subtree indicated start..end
        min_depth: usize, max_depth: usize)
        where T: std::marker::Copy {
    let n = (start..end).len();
    if min_depth == max_depth {
        // TODO: probably could relax this and return if size is <= 3
        assert!(n == 1);
        return;
    }

    // eprintln!("{} {} {} {} {}",
    //     start, end,
    //     childless,
    //     info[idxs[start]].depth,
    //     info[idxs[start]].subtree_size);

    // small subtrees are just put in sorted order
    if childless && info[idxs[start]].subtree_size <= SIMPLE_SUBTREE_CUTOFF {
        assert!(n == info[idxs[start]].subtree_size);
        info[idxs[start]].simple = true;
        return;
    }

    let pivot_depth = min_depth + (max_depth - min_depth) / 2;
    let top_size = stable_partition_by_depth(
        idxs, tmp, info, pivot_depth, start, end);

    // recurse on top subtree
    veb_order_recursion(
        idxs, tmp, info, nodes, start,
        start + top_size, false, min_depth, pivot_depth);

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
            idxs, tmp, info, nodes, i, j, childless,
            bottom_subtree_depth, subtree_max_depth);
        i = j;
    }
}


// traverse the tree filling the subtree_first and subtree_last fields.
fn compute_subtree_sizes<T>(nodes: &mut [IntervalNode<T>], root_idx: usize)
        where T: std::marker::Copy {

    let mut subtree_first = nodes[root_idx].first;
    let mut subtree_last = nodes[root_idx].last;

    if nodes[root_idx].simple {
        for k in 1..nodes[root_idx].right as usize {
            subtree_first = min(subtree_first, nodes[root_idx+k].first);
            subtree_last = max(subtree_last, nodes[root_idx+k].last);
        }
    } else {
        let left = nodes[root_idx].left as usize;
        if left < u32::max_value() as usize {
            compute_subtree_sizes(nodes, left);
            subtree_first = min(subtree_first, nodes[left].subtree_first);
            subtree_last  = max(subtree_last, nodes[left].subtree_last);
        }

        let right = nodes[root_idx].right as usize;
        if right < u32::max_value() as usize {
            compute_subtree_sizes(nodes, right);
            subtree_first = min(subtree_first, nodes[right].subtree_first);
            subtree_last  = max(subtree_last, nodes[right].subtree_last);
        }
    }

    nodes[root_idx].subtree_first = subtree_first;
    nodes[root_idx].subtree_last = subtree_last;
}


fn compute_tree_size<T>(nodes: &[IntervalNode<T>], root_idx: usize) -> usize
        where T: std::marker::Copy {

    let mut subtree_size = 1;

    if nodes[root_idx].simple {
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


// put nodes in van Emde Boas order
fn veb_order<T>(nodes: &mut [IntervalNode<T>])
        where T: std::marker::Copy {

    // it seems to not matter all that much how this is sorted
    nodes.sort_unstable_by_key(|node| node.first);
    // nodes.sort_unstable_by_key(|node| node.last);
    // nodes.sort_unstable_by_key(|node| (node.first, node.last));
    // nodes.sort_unstable_by_key(|node| (node.first, -node.last));
    // nodes.sort_unstable_by_key(|node| (node.last + node.first)/2);

    let mut info = traverse(nodes.len());

    let mut max_depth = 0;
    for info_i in &info {
        if info_i.depth > max_depth {
            max_depth = info_i.depth;
        }
    }

    let mut idxs: Vec<usize> = (0..nodes.len()).collect();
    let mut tmp: Vec<usize> = vec![0; nodes.len()];

    for i in &idxs {
        tmp[info[*i].dfs] = *i;
    }
    idxs.copy_from_slice(&mut tmp);

    veb_order_recursion(
        &mut idxs, &mut tmp, &mut info, &nodes, 0, nodes.len(), true, 0, max_depth);

    // idxs is now a vEB -> sorted order map. Build the reverse here.
    let mut revidx = tmp;
    for (i, j) in idxs.iter().enumerate() {
        revidx[*j] = i;
    }

    // but nodes in vEB order in a temporary vector
    let mut veb_nodes = vec![nodes[0]; nodes.len()];
    for i in 0..nodes.len() {
        veb_nodes[i] = nodes[idxs[i]];

        if info[idxs[i]].subtree_size <= SIMPLE_SUBTREE_CUTOFF {
            veb_nodes[i].simple = true;
            veb_nodes[i].left = u32::max_value();
            veb_nodes[i].right = info[idxs[i]].subtree_size as u32;
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

    // copy reordered nodes back to the original vector
    nodes.copy_from_slice(&mut veb_nodes);
    compute_subtree_sizes(nodes, 0);
}


fn overlaps(first_a: i32, last_a: i32, first_b: i32, last_b: i32) -> bool {
    return first_a <= last_b && last_a >= first_b;
}


fn query_recursion(
        nodes: &[IntervalNode<()>], root_idx: usize, first: i32, last: i32,
        count: &mut usize, overlap: &mut usize, visited: &mut usize) {
    // println!("{} {:?} {:?} {} {}",
        // root_idx, nodes[root_idx].left, nodes[root_idx].right,
        // nodes[root_idx].first, nodes[root_idx].last);

    let node = &nodes[root_idx];

    *visited += 1;
    if overlaps(node.first, node.last, first, last) {
        *count += 1;
    }

    if node.simple {
        for k in 1..node.right as usize {
            let node = &nodes[root_idx+k];
            if overlaps(node.first, node.last, first, last) {
                *count += 1;
            }
        }
    } else {
        let left = node.left as usize;
        if left < u32::max_value() as usize {
            if overlaps(
                    nodes[left].subtree_first, nodes[left].subtree_last,
                    first, last) {
                query_recursion(nodes, left, first, last, count, overlap, visited);
            }
        }

        let right = node.right as usize;
        if right < u32::max_value() as usize {
            if overlaps(
                    nodes[right].subtree_first, nodes[right].subtree_last,
                    first, last) {
                query_recursion(nodes, right, first, last, count, overlap, visited);
            }
        }
    }

    // alternative check, using just subtree_last
    // if first <= nodes[root_idx].subtree_last {
    //     if last >= nodes[root_idx].first {
    //         if first <= nodes[root_idx].last {
    //             println!("hit!");
    //             *count += 1;
    //         }

    //         if let Some(right) = nodes[root_idx].right {
    //             query_recursion(nodes, right, first, last, count, overlap);
    //         }
    //     }

    //     if let Some(left) = nodes[root_idx].left {
    //         query_recursion(nodes, left, first, last, count, overlap);
    //     }
    // }
}

// super simple query which prints every overlap
fn query(tree: &COITree<()>, first: i32, last: i32) -> (usize, usize, usize) {
    // println!("QUERY");
    let mut count = 0;
    let mut overlap = 0;
    let mut visited = 0;
    query_recursion(
        &tree.nodes, 0, first, last, &mut count, &mut overlap, &mut visited);
    return (count, overlap, visited);
}


// This is about the same as the recursive version but relies on a fixed
// size stack space.
// fn inlined_query(
//         tree: &COITree<()>,
//         first: i32, last: i32) -> (usize, usize, usize) {
//     let mut count = 0;
//     let mut overlap = 0;
//     let mut visited = 0;
//     let nodes = &tree.nodes;
//     let stack: &mut [u32]  = &mut [0; 64];
//     let mut p: isize = 0;

//     while p >= 0 {
//         let root_idx = stack[p as usize] as usize;
//         p -= 1;

//         visited += 1;
//         if overlaps(nodes[root_idx].first, nodes[root_idx].last, first, last) {
//             count += 1;
//         }

//         let left = nodes[root_idx].left as usize;
//         if left < u32::max_value() as usize {
//             if overlaps(
//                     nodes[left].subtree_first, nodes[left].subtree_last,
//                     first, last) {
//                 p += 1;
//                 stack[p as usize] = left as u32;
//             }
//         }

//         let right = nodes[root_idx].right as usize;
//         if right < u32::max_value() as usize {
//             if overlaps(
//                     nodes[right].subtree_first, nodes[right].subtree_last,
//                     first, last) {
//                 p += 1;
//                 stack[p as usize] = right as u32;
//             }
//         }
//     }

//     return (count, overlap, visited);
// }


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
            subtree_first: first,
            subtree_last: last,
            left: u32::max_value(),
            right: u32::max_value(),
            simple: false,
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

    let stdout = io::stdout();
    let mut out = stdout.lock();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first_str, last_str, first, last) =
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

        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            seqname, first_str, last_str, count)?;

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

