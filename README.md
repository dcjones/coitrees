
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
| coitrees                            |  **22.4s** |      11.7s | **1.6s** |     29.7s  |
| cgranges (`bedcov-cr -c`)           |      66.3s |  **10.9s** |     3.6s | **28.4s**  |
| AIList                              |      24.9s |      15.5s |     1.9s |     37.2s  |
| CITree                              |      33.9s |      23.7s |     2.6s |    102.6s  |
| NCList                              |      41.1s |      28.7s |     3.3s |     57.0s  |
| AITree                              |      41.3s |      43.4s |     3.5s |     94.8s  |
| `bedtools coverage -counts`         |     807.6s |    1080.5s |   311.4s |   4820.0s  |
| `bedtools coverage -counts -sorted` |     584.7s |     660.7s |   304.8s |   2846.8s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees                            |  **51.8s** |  **18.5s** | **3.6s** |     32.8s |
| cgranges (`bedcov-cr -c`)           |      99.0s |      19.4s |     6.3s |     35.4s |
| AIList                              |      58.6s |      31.2s |     4.2s | **30.8s** |
| CITree                              |      77.8s |      32.4s |     5.3s |     82.0s |
| NCList                              |      73.4s |      40.0s |     5.7s |     59.7s |
| AITree                              |     376.7s |     314.3s |    26.3s |   1575.3s |
| `bedtools coverage -counts`         |   1809.2s  |    1511.5s |   358.0s |  13056.0s |

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
