
# COITrees: Cache Oblivious Interval Trees

This is a prototype data structure inspired by
[cgranges](https://github.com/lh3/cgranges) and my previous work on
[IntevalTrees.jl](https://github.com/BioJulia/IntervalTrees.jl). It's a very
fast data structure for storing and querying collections of intervals.

Like cgranges this data structure is an augmented binary search tree stored
in contiguous memory, but following from IntervalTrees.jl (which uses
B+-trees) it tries to improve performance by reducing cache misses.

It does this by storing the tree in [van Emde Boas
layout](http://erikdemaine.org/papers/FOCS2000b/paper.pdf). That adds some
overhead in computing the layout and storing child pointers, but makes
searching through large trees significantly faster.

The one additional optimization trick is that small subtrees at the bottom
are stored in no particular order and searched through linearly.

# Using

This is just a prototype implemented as a program that takes two bed files and
prints the number of intervals in the first file that overlap each interval in the second.

To try out, just clone this repo and run:
```shell
cargo run --release -- test1.bed test2.bed > intersections.bed
```

# Benchmarks

`A` is 2,755,864 intervals from Ensembl's human genome annotations, `B` is
62,159,484 intervals from some RNA-Seq alignments, and `B'` is the first 2
million lines of `B`.

## Intervals in sorted order

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees (`--sorted`)               |  **13.6s** |      12.2s | **1.3s** |  **1.2s**  |
| coitrees                            |      22.8s |  **10.8s** |     1.6s |     15.6s  |
| cgranges (`bedcov-cr -c`)           |      64.8s |      12.0s |     3.7s |     30.3s  |
| AIList                              |      25.2s |      15.2s |     1.9s |     29.9s  |
| CITree                              |      34.4s |      25.0s |     2.7s |     73.3s  |
| NCList                              |      41.6s |      28.8s |     3.5s |     55.8s  |
| AITree                              |      54.9s |      44.6s |     3.9s |    184.8s  |
| `bedtools coverage -counts -sorted` |     868.3s |     371.9s |   315.9s |   4865.8s  |
| `bedtools coverage -counts`         |     977.0s |     626.2s |   329.4s |   4690.8s  |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees                            |  **29.7s** |  **11.7s** | **2.0s** | **18.8s**  |
| cgranges                            |      69.9s |      14.1s |     4.0s |     54.8s  |
| CITree                              |      50.2s |      47.1s |     3.8s |    268.6s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees                            |  **48.0s** |  **16.4s** | **3.6s** | **22.9s** |
| cgranges (`bedcov-cr -c`)           |     101.7s |      20.2s |     6.2s |     36.8s |
| AIList                              |      59.5s |      31.5s |     4.3s |     32.3s |
| CITree                              |      79.8s |      34.3s |     5.6s |    114.0s |
| NCList                              |      77.2s |      41.6s |     6.1s |     63.7s |
| AITree                              |     341.1s |     259.5s |    23.6s |   1583.9s |
| `bedtools coverage -counts`         |    2040.1s |    1721.7s |   379.9s |  11098.0s |

### With coverage

|                                     |     A vs B |     B vs A |  A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | ---------: |
| coitrees                            |  **62.0s** |  **17.4s** | **4.3s** | **26.7s**  |
| cgranges                            |     106.1s |      22.3s |     6.6s |     61.7s  |
| CITree                              |      91.8s |      55.4s |     6.6s |    312.7s  |

# Discussion

To put things in context, AIList, AITree, and NCList all sleazily use
hardcoded chromosome names to avoid hashing names, and skip mitochondrial
intervals for some reason, which excludes 126 intervals in `A` and 3 million
in `B`. coitrees and cgranges don't make these assumptions but are optimized
to this particular benchmark. bedtools on the other hand is an actual useful
tool with a lot of features.

It should also be noted that many of these data structures are so fast that
precisely how the BED files are read and written is a significant factor.
Even when piping output to `/dev/null`, as was done in these benchmarks,
reading/writing can be the bulk of the runtime. I had to optimize this in
coitrees by calling out to `printf` in an `unsafe` block.

A fair benchmark would have each of these methods implemented as a library
with a standard interface, then a benchmark program would link and call
each. If someone wanted to implement that, I'd probably help.

In [cgranges](https://github.com/lh3/cgranges) there are also benchmarks that
compute covered fraction, not just intersection counts. I haven't implemented
this yet, but I suspect cgranges will have an advantage in keeping intervals
in sorted order.

Another good benchmark would be to store some metadata with each interval,
count intersections that pass some filter based on that metadata. That way no
counting shortcuts could be made and every intersection would need to be
inspected. Not all these methods support storing arbitrary data with the
interval, so that too would take some work.
