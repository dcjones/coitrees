
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
| coitrees (`--sorted`)               |  **11.9s** |      11.2s | **1.2s** |  **1.1s**  |
| coitrees                            |      21.6s |   **9.9s** |     1.5s |     11.8s  |
| cgranges (`bedcov-cr -c`)           |      65.0s |      11.8s |     3.8s |     30.3s  |
| AIList                              |      25.0s |      15.0s |     1.9s |     30.1s  |
| CITree                              |      35.0s |      26.2s |     2.7s |     83.4s  |
| NCList                              |      41.8s |      29.2s |     3.5s |     56.2s  |
| AITree                              |      53.2s |      44.6s |     3.9s |    180.4s  |
| `bedtools coverage -counts -sorted` |     868.3s |     371.9s |   315.9s |   4865.8s  |
| `bedtools coverage -counts`         |     977.0s |     626.2s |   329.4s |   4690.8s  |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees                            |  **29.5s** |  **11.2s** | **1.9s** | **22.6s**  |
| cgranges                            |      69.2s |      14.1s |     4.0s |     54.0s  |
| CITree                              |      43.9s |      46.7s |     3.7s |    272.8s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A   | B' vs B'  |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees                            |  **47.3s** |  **16.0s** | **3.6s** | **18.2s** |
| cgranges (`bedcov-cr -c`)           |     102.0s |      19.8s |     6.2s |     36.7s |
| AIList                              |      60.0s |      31.3s |     4.4s |     32.6s |
| CITree                              |      76.3s |      35.0s |     5.3s |     80.3s |
| NCList                              |      75.2s |      40.8s |     6.0s |     62.6s |
| AITree                              |     332.8s |     253.2s |    23.0s |   1463.4s |
| `bedtools coverage -counts`         |    2040.1s |    1721.7s |   379.9s |  11098.0s |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees                            |  **61.0s** |  **17.2s** | **4.3s** | **27.5s**  |
| cgranges                            |     109.1s |      22.4s |     6.6s |     62.8s  |
| CITree                              |      92.6s |      56.5s |     6.6s |    280.0s  |

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
