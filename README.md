
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
| coitrees                            |  **23.3s** |      14.7s | **1.8s** | **30.0s**  |
| cgranges (`bedcov-cr -c`)           |      69.6s |  **12.1s** |   3.8s   |  33.9s     |
| AIList                              |      26.8s |      15.3s |   2.0s   |  36.9s     |
| CITree                              |      34.4s |      25.0s |   2.7s   |  67.6s     |
| NCList                              |      40.2s |      30.2s |   3.5s   |  55.4s     |
| AITree                              |      39.9s |      37.0s |   3.0s   |  95.1s     |
| `bedtools coverage -counts`         |     807.6s |    1080.5s |   311.4s |  4820.0s   |
| `bedtools coverage -counts -sorted` |     584.7s |     660.7s |   304.8s |  2846.8s   |

## Intervals in randomized order

|                                     |     A vs B |     B vs A | A vs A  | B' vs B'   |
| ----------------------------------- | ---------: | ---------: | -------: | --------: |
| coitrees                            |  **51.7s** |  **18.5s** | **3.6s** | **32.9s** |
| cgranges (`bedcov-cr -c`)           |     104.0s |      19.9s |     5.9s |     34.6s |
| AIList                              |      59.4s |      32.6s |     4.3s |     33.3s |
| CITree                              |      76.5s |      31.9s |     5.1s |     73.0s |
| NCList                              |      74.1s |      42.1s |     5.8s |     64.9s |
| AITree                              |     366.2s |     306.2s |    25.9s |   1577.0s |
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
