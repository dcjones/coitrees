
# COITrees: Cache Oblivious Interval Trees

COITrees implements a data structure for very fast overlap queries of a
static set of integer intervals, with genomic intervals in mind.

Borrowing from [cgranges](https://github.com/lh3/cgranges), this data
structure stores intervals in contiguous memory, but improves query
performance by storing the nodes in in-order [van Emde Boas
layout](http://erikdemaine.org/papers/FOCS2000b/paper.pdf) layout. Computing
the layout requires some extra time and memory, but improves average cache
locality for queries of the tree. If the interval set is relatively large,
and a sufficiently large number of queries are performed, it tends to out-perform
other data structures.

The `SortedQuerent` type implements an alternative query strategy that keeps track
of the results of the previous query. When a query overlaps the previous one,
the results from that previous query can be reused to dramatically accelerate
the current one. (In the benchmarks, this is the `--sorted` option.)

Some operations can further be sped up using SIMD instructions. Two COITree
variants are implemented to exloit AVX2 instructions on x86-64 cpus
(`AVXCOITree`), and Neon instructions on ARM cpus (`NeonCOITree`). The `COITree`
type is oppurtunistically defined to one of these types if the right instruction
set is detected. Typically it's necessary to compile with the environment
variable `RUSTFLAGS="-Ctarget-cpu=native"` set for this to work. The fallback
implemntation (`BasicCOITree`) supports any platform rust compiles to and
remains highly efficient.

# Trying Out

This is primary a library for use in other programs, but for benchmarking
purposes it includes a program for intersecting BED files.

To try out, just clone this repo and run:
```shell
cargo run --release --example bed-intersect -- test1.bed test2.bed > intersections.bed
```

# Benchmarks

`A` is 2,755,864 intervals from Ensembl's human genome annotations, `B` is
62,159,484 intervals from some RNA-Seq alignments, and `B'` is the first 2
million lines of `B`.

## Intervals in sorted order

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees AVX                        |      10.1s |   **3.6s** | **0.6s** |      2.3s  |
| coitrees AVX (`--sorted`)           |       6.5s |       4.9s | **0.6s** |  **0.6s**  |
| coitrees                            |      12.0s |       5.9s |     0.8s |      7.7s  |
| coitrees (`--sorted`)               |   **6.0s** |       6.0s | **0.6s** |  **0.6s**  |
| cgranges (`bedcov-cr -c`)           |      34.9s |       6.6s |     2.0s |     16.0s  |
| AIList                              |      13.2s |      10.1s |     1.1s |     24.3s  |
| CITree                              |      20.6s |      14.0s |     1.5s |     45.6s  |
| NCList                              |      23.3s |      19.9s |     2.0s |     41.7s  |
| AITree                              |      22.2s |      26.4s |     2.2s |     60.2s  |
| `bedtools coverage -counts -sorted` |     257.5s |     295.6s |    71.6s |   2130.9s  |
| `bedtools coverage -counts`         |     322.4s |     378.5s |    75.0s |   3595.9s  |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees AVX                        |      16.0s |   **4.8s** | **1.0s** | **12.0s**  |
| coitrees                            |  **14.2s** |       6.8s | **1.0s** |     14.6s  |
| cgranges                            |      37.3s |       7.8s |     2.2s |     29.4s  |
| CITree                              |      23.2s |      26.5s |     2.0s |    161.1s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A   | B' vs B'  |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees AVX                        |  **20.8s** |   **7.5s** | **1.5s** |  **3.3s** |
| coitrees                            |      24.2s |      10.1s |     1.9s |      9.3s |
| cgranges (`bedcov-cr -c`)           |      57.6s |      11.1s |     3.4s |     18.4s |
| AIList                              |      33.4s |      21.5s |     2.6s |     24.9s |
| CITree                              |      38.2s |      19.3s |     2.8s |     45.8s |
| NCList                              |      40.0s |      27.9s |     3.6s |     43.4s |
| AITree                              |     228.7s |     147.0s |    16.3s |    883.4s |
| `bedtools coverage -counts`         |    1160.4s |     849.6s |   104.5s |   9254.6s |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees AVX                        |  **31.0s** |   **9.5s** | **2.0s** | **13.9s**  |
| coitrees                            |      31.8s |      10.5s |     2.2s |     16.5s  |
| cgranges                            |      57.3s |      12.7s |     3.5s |     30.9s  |
| CITree                              |      47.5s |      33.3s |     3.7s |    170.8s  |


All benchmarks run on a ryzen 5950x.

# Discussion

These benchmarks are somewhat realistic in that they use real data, but are
not entirely apples-to-apples because they all involve parsing and writing
BED files. Most of the programs (including the one implemented in coitrees)
have incomplete BED parsers, and some use other shortcuts like assuming a
fixed set of chromosomes with specific naming schemes.

`bedtools` carries the disadvantage of being an actually useful tool, rather
than implemented being implemented entirely for the purpose of winning benchmark
games. It seems clear it could be a lot faster, but there no doubt some cost can
be chalked up to featurefulness, completeness, and safety.

If you have a BED intersection program you suspect may be faster (or just
interesting), please let me know and I'll try to benchmark it.
