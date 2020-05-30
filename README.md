
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
| coitrees (`--sorted`)               |  **13.0s** |      11.2s | **1.2s** |  **1.1s**  |
| coitrees                            |      22.8s |  **10.7s** |     1.5s |     19.2s  |
| cgranges (`bedcov-cr -c`)           |      65.7s |      11.8s |     3.7s |     29.6s  |
| AIList                              |      25.1s |      14.9s |     2.0s |     29.4s  |
| CITree                              |      34.0s |      25.5s |     2.7s |    105.6s  |
| NCList                              |      28.7s |      41.4s |     3.5s |     55.7s  |
| AITree                              |      54.3s |      45.0s |     3.9s |    185.7s  |
| `bedtools coverage -counts -sorted` |     868.3s |     371.9s |   315.9s |   4865.8s  |
| `bedtools coverage -counts`         |     977.0s |     626.2s |   329.4s |   4690.8s  |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees                            |  **48.6s** |  **16.1s** | **3.6s** | **25.4s** |
| cgranges (`bedcov-cr -c`)           |     100.3s |      20.4s |     6.3s |     35.6s |
| AIList                              |      58.5s |      31.2s |     4.3s |     31.1s |
| CITree                              |      74.3s |      34.9s |     5.3s |     82.2s |
| NCList                              |      40.2s |      74.0s |     5.9s |     61.8s |
| AITree                              |     324.0s |     242.0s |    22.7s |   1473.0s |
| `bedtools coverage -counts`         |    2040.1s |    1721.7s |   379.9s |  11098.0s |

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
