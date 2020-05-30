
use std::error::Error;
use std::str;
use std::time::Instant;
use std::fs::File;
use std::io::{BufRead, BufReader};
use coitrees::{COITree, IntervalNode, SortedQuerent};

extern crate fnv;
use fnv::FnvHashMap;

extern crate clap;
use clap::{Arg, App};

extern crate libc;

type GenericError = Box<dyn Error>;


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

        node_arr.push(IntervalNode::new(first, last, ()));

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

    let mut total_count: usize = 0;
    let now = Instant::now();

    // let stdout = io::stdout();
    // let mut out = stdout.lock();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, _first_str, _last_str, first, last) =
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


fn query_bed_files_with_sorted_querent(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let trees = read_bed_file(filename_a)?;

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();

    let mut total_count: usize = 0;
    let now = Instant::now();

    let mut querents = FnvHashMap::<String, SortedQuerent<()>>::default();
    for (seqname, tree) in &trees {
        querents.insert(seqname.clone(), SortedQuerent::new(tree));
    }

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, _first_str, _last_str, first, last) =
            parse_bed_line(&line);

        let mut count: usize = 0;
        if let Some(querent) = querents.get_mut(seqname) {
            querent.query(first, last, |_| count += 1);
        }

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


fn main() {
    let matches = App::new("coitrees")
        .about("Find overlaps between two groups of intervals")
        .arg(Arg::with_name("input1")
            .about("intervals to index")
            .value_name("intervals.bed")
            .required(true)
            .index(1))
        .arg(Arg::with_name("input2")
            .about("query intervals")
            .value_name("queries.bed")
            .required(true)
            .index(2))
        .arg(Arg::with_name("use_sorted_querent")
            .long("--sorted")
            .short('s')
            .about("use alternative search strategy that's faster if queries are sorted and tend to overlap"))
        .get_matches();

    let input1 = matches.value_of("input1").unwrap();
    let input2 = matches.value_of("input2").unwrap();

    if matches.is_present("use_sorted_querent") {
        let result = query_bed_files_with_sorted_querent(input1, input2);
        if let Err(err) = result {
            println!("error: {}", err)
        }
    } else {
        let result = query_bed_files(input1, input2);
        if let Err(err) = result {
            println!("error: {}", err)
        }
    }
}

