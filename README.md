
# COITrees: Cache Oblivious Interval Trees

COITrees implements a data structure for very fast overlap queries of a
static set of integer intervals, with genomic intervals in mind.

Borrowing from [cgranges](https://github.com/lh3/cgranges), this data
structure stores intervals in contiguous memory, but improves query
performance by storing the nodes in in-order [van Emde Boas
layout](http://erikdemaine.org/papers/FOCS2000b/paper.pdf). Computing
the layout requires some extra time and memory, but improves average cache
locality for queries of the tree. If the interval set is relatively large,
and a sufficiently large number of queries are performed, it tends to out-perform
other data structures.

The `SortedQuerent` type implements an alternative query strategy that keeps track
of the results of the previous query. When a query overlaps the previous one,
the results from that previous query can be reused to dramatically accelerate
the current one. (In the benchmarks, this is the `--sorted` option.)

Some operations can further be sped up using SIMD instructions. Two COITree
variants are implemented to exploit AVX2 instructions on x86-64 cpus
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
| coitrees AVX                        |      11.8s |   **3.7s** |      0.7 |      5.3s  |
| coitrees AVX (`--sorted`)           |       6.4s |       4.2s | **0.6s** |  **0.5s**  |
| coitrees                            |      11.4s |       5.2s |     0.8s |      8.3s  |
| coitrees (`--sorted`)               |   **5.8s** |       5.4s | **0.6s** |  **0.5s**  |
| cgranges (`bedcov-cr -c`)           |      35.4s |       6.6s |     2.0s |     17.6s  |
| AIList                              |      13.8s |      10.1s |     1.1s |     18.4s  |
| CITree                              |      20.1s |      13.5s |     1.6s |     45.7s  |
| NCList                              |      22.5s |      16.8s |     1.9s |     39.8s  |
| AITree                              |      23.8s |      26.3s |     2.1s |     63.4s  |
| `bedtools coverage -counts -sorted` |     257.5s |     295.6s |    71.6s |   2130.9s  |
| `bedtools coverage -counts`         |     322.4s |     378.5s |    75.0s |   3595.9s  |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees AVX                        |      18.2s |   **4.8s** |     1.1s |      16.0s |
| coitrees                            |  **14.6s** |       5.7s | **1.0s** |  **12.0s** |
| cgranges                            |      38.4s |       8.1s |     2.2s |     31.0s  |
| CITree                              |      23.2s |      25.6s |     2.0s |    160.4s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A   | B' vs B'  |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees AVX                        |  **23.9s** |   **7.2s** | **1.6s** |  **6.1s** |
| coitrees                            |      24.2s |       8.9s |     1.9s |      9.4s |
| cgranges (`bedcov-cr -c`)           |      55.7s |      11.1s |     3.3s |     19.6s |
| AIList                              |      31.2s |      18.2s |     2.3s |     19.3s |
| CITree                              |      39.4s |      19.0s |     2.9s |     47.1s |
| NCList                              |      42.7s |      23.8s |     3.4s |     44.0s |
| AITree                              |     225.3s |     134.8s |    14.7s |    921.6s |
| `bedtools coverage -counts`         |    1160.4s |     849.6s |   104.5s |   9254.6s |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees AVX                        |      34.3s |   **8.8s** | **2.2s** |     16.3s  |
| coitrees                            |  **29.6s** |       9.7s |     2.3s | **13.1s**  |
| cgranges                            |      57.6s |      12.5s |     3.6s |     32.6s  |
| CITree                              |      50.0s |      32.5s |     3.8s |    170.4s  |


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
