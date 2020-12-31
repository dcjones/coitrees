
use std::error::Error;
use std::str;
use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ffi::CString;
use coitrees::{COITree, Interval};

extern crate fnv;
use fnv::FnvHashMap;

extern crate clap;
use clap::{Arg, App};

extern crate libc;

type GenericError = Box<dyn Error>;


// Parse a i32 with no checking whatsoever. (e.g. non-number characters will just)
fn i32_from_bytes_uncheckd(s: &[u8]) -> i32 {
    if s.is_empty() {
        return 0;
    } else if s[0] == b'-' {
        return -s[1..].iter().fold(0, |a, b| a*10 + (b & 0x0f) as i32);
    } else {
        return s.iter().fold(0, |a, b| a*10 + (b & 0x0f) as i32);
    }
}


fn parse_bed_line(line: &[u8]) -> (&str, i32, i32) {
    let n = line.len() - 1;
    let mut p = 0;
    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let seqname = unsafe {
        str::from_utf8_unchecked(&line[..p])
    };
    p += 1;
    let p0 = p;

    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let first = i32_from_bytes_uncheckd(&line[p0..p]);
    p += 1;
    let p0 = p;

    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let last = i32_from_bytes_uncheckd(&line[p0..p]) - 1;

    return (seqname, first, last);
}


// Read a bed file into a COITree
fn read_bed_file(path: &str) -> Result<FnvHashMap<String, COITree<i32, ()>>, GenericError> {
    let mut nodes = FnvHashMap::<String, Vec<Interval<i32, ()>>>::default();

    let now = Instant::now();

    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line_count = 0;
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) =
            parse_bed_line(&line);

        let node_arr = if let Some(node_arr) = nodes.get_mut(seqname) {
            node_arr
        } else {
            nodes.entry(seqname.to_string()).or_insert(Vec::new())
        };

        node_arr.push(Interval {
            first: first, last: last, metadata: () });

        line_count += 1;
        line.clear();
    }

    eprintln!("reading bed: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("lines: {}", line_count);
    eprintln!("sequences: {}", nodes.len());

    let now = Instant::now();
    let mut trees = FnvHashMap::<String, COITree<i32, ()>>::default();
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(seqname_nodes));
    }
    eprintln!("veb_order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    return Ok(trees);
}


fn read_bed_file_numbered(path: &str) -> Result<FnvHashMap<String, COITree<i32, usize>>, GenericError> {
    let mut nodes = FnvHashMap::<String, Vec<Interval<i32, usize>>>::default();

    let now = Instant::now();

    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line_count = 0;
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) =
            parse_bed_line(&line);

        let node_arr = if let Some(node_arr) = nodes.get_mut(seqname) {
            node_arr
        } else {
            nodes.entry(seqname.to_string()).or_insert(Vec::new())
        };

        node_arr.push(Interval {
            first: first, last: last, metadata: node_arr.len() });

        line_count += 1;
        line.clear();
    }

    eprintln!("reading bed: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("lines: {}", line_count);
    eprintln!("sequences: {}", nodes.len());

    let now = Instant::now();
    let mut trees = FnvHashMap::<String, COITree<i32, usize>>::default();
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

    let mut total_count: usize = 0;
    let now = Instant::now();

    // let stdout = io::stdout();
    // let mut out = stdout.lock();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) =
            parse_bed_line(&line);

        let mut count: usize = 0;

        if let Some(seqname_tree) = tree.get(seqname) {
            // seqname_tree.query(first, last, |_| count += 1);
            count = seqname_tree.query_count(first, last);
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

        total_count += count;

        line.clear();
    }

    eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("total overlaps: {}", total_count);

    return Ok(());
}


// fn query_bed_files_tvt(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
//     let a_trees = read_bed_file(filename_a)?;
//     let b_trees = read_bed_file_numbered(filename_b)?;

//     let mut a_querents = FnvHashMap::<String, SortedQuerent<(), u32>>::default();
//     for (seqname, a_tree) in &a_trees {
//         a_querents.insert(seqname.clone(), SortedQuerent::new(a_tree));
//     }

//     let mut total_count = 0;
//     for (seqname, b_tree) in &b_trees {
//         if let Some(a_querent) = a_querents.get_mut(seqname) {
//             let c_seqname = CString::new(seqname.clone()).expect("CString::new failed");

//             for b_node in b_tree {
//                 let mut count = 0;
//                 a_querent.query(b_node.first, b_node.last, |_| count += 1);
//                 total_count += count;

//                 unsafe {
//                     libc::printf(
//                         b"%s\t%d\t%d\t%u\n\0".as_ptr() as *const i8,
//                         c_seqname.as_bytes_with_nul().as_ptr() as *const i8,
//                         b_node.first,
//                         b_node.last,
//                         count as u32);
//                 }
//             }
//         }
//     }

//     eprintln!("total overlaps: {}", total_count);

//     return Ok(());
// }


// fn query_bed_files_coverage(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
//     let tree = read_bed_file(filename_a)?;

//     let file = File::open(filename_b)?;
//     let mut rdr = BufReader::new(file);
//     let mut line = Vec::new();

//     let mut total_count: usize = 0;
//     let now = Instant::now();

//     // let stdout = io::stdout();
//     // let mut out = stdout.lock();

//     while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
//         let (seqname, first, last) =
//             parse_bed_line(&line);

//         let mut cov: usize = 0;
//         let mut count: usize = 0;

//         if let Some(seqname_tree) = tree.get(seqname) {
//             let countcov = seqname_tree.coverage(first, last);
//             count = countcov.0;
//             cov = countcov.1;
//         }

//         // out.write(&line[..line.len()-1])?;
//         // writeln!(out, "\t{}", count)?;

//         // unfortunately printing in c is quite a bit faster than rust
//         unsafe {
//             let linelen = line.len();
//             line[linelen-1] = b'\0';
//             libc::printf(
//                 b"%s\t%u\t%u\n\0".as_ptr() as *const i8,
//                 line.as_ptr() as *const i8,
//                 count as u32,
//                 cov);
//         }

//         total_count += count;

//         line.clear();
//     }

//     eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
//     eprintln!("total overlaps: {}", total_count);

//     return Ok(());
// }


// fn query_bed_files_with_sorted_querent(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
//     let trees = read_bed_file(filename_a)?;

//     let file = File::open(filename_b)?;
//     let mut rdr = BufReader::new(file);
//     let mut line = Vec::new();

//     let mut total_count: usize = 0;
//     let now = Instant::now();

//     let mut querents = FnvHashMap::<String, SortedQuerent<(), u32>>::default();
//     for (seqname, tree) in &trees {
//         querents.insert(seqname.clone(), SortedQuerent::new(tree));
//     }

//     while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
//         let (seqname, first, last) =
//             parse_bed_line(&line);

//         let mut count: usize = 0;
//         if let Some(querent) = querents.get_mut(seqname) {
//             querent.query(first, last, |_| count += 1);
//         }

//         // unfortunately printing in c is quite a bit faster than rust
//         unsafe {
//             let linelen = line.len();
//             line[linelen-1] = b'\0';
//             libc::printf(
//                 b"%s\t%u\n\0".as_ptr() as *const i8,
//                 line.as_ptr() as *const i8,
//                 count as u32);
//         }

//         total_count += count;

//         line.clear();
//     }

//     eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
//     eprintln!("total overlaps: {}", total_count);

//     return Ok(());
// }


fn main() {
    let matches = App::new("coitrees")
        .about("Find overlaps between two groups of intervals")
        .arg(Arg::new("input1")
            .about("intervals to index")
            .value_name("intervals.bed")
            .required(true)
            .index(1))
        .arg(Arg::new("input2")
            .about("query intervals")
            .value_name("queries.bed")
            .required(true)
            .index(2))
        .arg(Arg::new("use_sorted_querent")
            .long("--sorted")
            .short('s')
            .about("use alternative search strategy that's faster if queries are sorted and tend to overlap"))
        .arg(Arg::new("tree_vs_tree")
            .long("--tree-vs-tree")
            .short('t')
            .about("load both interval sets into memory instead of streaming queries"))
        .arg(Arg::new("coverage")
            .long("--coverage")
            .short('c')
            .about("compute proportion of queries covered"))
        .get_matches();

    let input1 = matches.value_of("input1").unwrap();
    let input2 = matches.value_of("input2").unwrap();

    let result;
    if matches.is_present("coverage") {
        panic!("coverage not implemented")
        // result = query_bed_files_coverage(input1, input2);
    } else if matches.is_present("use_sorted_querent") {
        panic!("sorted querent not implemented")
        // result = query_bed_files_with_sorted_querent(input1, input2);
    } else if matches.is_present("tree_vs_tree") {
        panic!("tree vs tree not implemented")
        // result = query_bed_files_tvt(input1, input2);
    } else {
        result = query_bed_files(input1, input2);
    }

    if let Err(err) = result {
        println!("error: {}", err)
    }
}

