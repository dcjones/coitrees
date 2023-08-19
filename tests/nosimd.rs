mod tests {
    use coitrees::*;
    use rand::{thread_rng, Rng};

    mod iteration {
        use super::*;
        use std::collections::HashSet;

        fn random_interval(
            min_first: i32,
            max_last: i32,
            min_len: i32,
            max_len: i32,
        ) -> (i32, i32) {
            let mut rng = thread_rng();
            let len = rng.gen_range(min_len..max_len + 1);
            let start = rng.gen_range(min_first..max_last - len + 1);
            (start, start + len - 1)
        }

        fn check_iteration(n: usize) {
            let mut b: Vec<IntervalNode<usize, usize>> = Vec::with_capacity(n);
            let min_first = 0;
            let max_last = 10000000;
            let min_len = 1;
            let max_len = 10000;
            for i in 0..n {
                let (first, last) = random_interval(min_first, max_last, min_len, max_len);
                b.push(IntervalNode::new(first, last, i));
            }

            let a: COITree<usize, usize> = COITree::new(&b);

            // check that intervals are sorted and that every value is generated
            let mut last_first = i32::min_value();
            let mut seen: HashSet<usize> = HashSet::new();
            for node in &a {
                assert!(last_first <= node.first);
                last_first = node.first;
                seen.insert(*node.metadata);
            }
            assert_eq!(seen.len(), n);
        }

        #[test]
        fn check_iteration_empty() {
            check_iteration(0);
        }

        #[test]
        fn check_iteration_medium() {
            check_iteration(10000);
        }
    }

    mod query {
        use super::*;

        // True iff the two intervals overlap.
        #[inline(always)]
        fn overlaps(first_a: i32, last_a: i32, first_b: i32, last_b: i32) -> bool {
            first_a <= last_b && last_a >= first_b
        }

        // Find overlapping intervals by simply checking every single one.
        // We test against this algorithm which we assume to be correct.
        fn brute_force_query<T, I, F>(
            intervals: &[IntervalNode<T, I>],
            query_first: i32,
            query_last: i32,
            mut visit: F,
        ) where
            T: Copy,
            I: IntWithMax,
            F: FnMut(&IntervalNode<T, I>),
        {
            for interval in intervals {
                if overlaps(interval.first, interval.last, query_first, query_last) {
                    visit(interval);
                }
            }
        }

        // Brute coverage calculation. `intervals` must be sorted.
        fn brute_force_coverage<T, I>(
            intervals: &[IntervalNode<T, I>],
            query_first: i32,
            query_last: i32,
        ) -> (usize, usize)
        where
            T: Copy,
            I: IntWithMax,
        {
            let mut last_cov = query_first - 1;
            let mut uncov_len = 0;
            let mut count = 0;
            for interval in intervals {
                if overlaps(interval.first, interval.last, query_first, query_last) {
                    if interval.first > last_cov {
                        uncov_len += interval.first - (last_cov + 1);
                    }
                    last_cov = last_cov.max(interval.last);
                    count += 1;
                }
            }
            if last_cov < query_last {
                uncov_len += query_last - last_cov;
            }

            let cov = ((query_last - query_first + 1) as usize) - (uncov_len as usize);
            (count, cov)
        }

        // Run queries against both a COITree and by brute force and check that
        // they get the same results.
        fn check_queries<I>(
            a: &COITree<u32, I>,
            b: &[IntervalNode<u32, I>],
            queries: &mut [(i32, i32)],
        ) where
            I: IntWithMax,
        {
            let mut a_hits: Vec<u32> = Vec::new();
            let mut b_hits: Vec<u32> = Vec::new();

            for (query_first, query_last) in queries {
                a_hits.clear();
                b_hits.clear();

                a.query(*query_first, *query_last, |node| {
                    a_hits.push(node.metadata.clone())
                });

                brute_force_query(b, *query_first, *query_last, |node| {
                    b_hits.push(node.metadata)
                });

                a_hits.sort();
                b_hits.sort();

                assert_eq!(a_hits, b_hits);
            }
        }

        fn check_coverage<I>(
            a: &COITree<u32, I>,
            b: &[IntervalNode<u32, I>],
            queries: &mut [(i32, i32)],
        ) where
            I: IntWithMax,
        {
            for (query_first, query_last) in queries {
                let (a_count, a_cover) = a.coverage(*query_first, *query_last);
                let (b_count, b_cover) = brute_force_coverage(b, *query_first, *query_last);

                assert_eq!(a_cover, b_cover);
                assert_eq!(a_count, b_count);
            }
        }

        fn check_count_queries<I>(
            a: &COITree<u32, I>,
            b: &[IntervalNode<u32, I>],
            queries: &mut [(i32, i32)],
        ) where
            I: IntWithMax,
        {
            for (query_first, query_last) in queries {
                let a_cnt = a.query_count(*query_first, *query_last);

                let mut b_cnt = 0;
                brute_force_query(b, *query_first, *query_last, |_| {
                    b_cnt += 1;
                });

                assert_eq!(a_cnt, b_cnt);
            }
        }

        // check SortedQuerent queries against brute force
        fn check_sorted_querent_queries<I>(
            a: &COITree<u32, I>,
            b: &[IntervalNode<u32, I>],
            queries: &mut [(i32, i32)],
        ) where
            I: IntWithMax,
        {
            let mut a_hits: Vec<u32> = Vec::new();
            let mut b_hits: Vec<u32> = Vec::new();

            let mut qa = COITreeSortedQuerent::new(a);

            queries.sort();

            for (query_first, query_last) in queries {
                a_hits.clear();
                b_hits.clear();

                qa.query(*query_first, *query_last, |node| {
                    a_hits.push(node.metadata.clone())
                });

                brute_force_query(b, *query_first, *query_last, |node| {
                    b_hits.push(node.metadata)
                });

                a_hits.sort();
                b_hits.sort();

                assert_eq!(a_hits, b_hits);
            }
        }

        // check that SortedQuerent still works even when queries are unsorted
        fn check_sorted_querent_unsorted_queries<I>(
            a: &COITree<u32, I>,
            b: &[IntervalNode<u32, I>],
            queries: &mut [(i32, i32)],
        ) where
            I: IntWithMax,
        {
            let mut a_hits: Vec<u32> = Vec::new();
            let mut b_hits: Vec<u32> = Vec::new();

            let mut qa = COITreeSortedQuerent::new(a);

            for (query_first, query_last) in queries {
                a_hits.clear();
                b_hits.clear();

                qa.query(*query_first, *query_last, |node| {
                    a_hits.push(node.metadata.clone())
                });

                brute_force_query(b, *query_first, *query_last, |node| {
                    b_hits.push(node.metadata)
                });

                a_hits.sort();
                b_hits.sort();

                assert_eq!(a_hits, b_hits);
            }
        }

        fn random_interval(
            min_first: i32,
            max_last: i32,
            min_len: i32,
            max_len: i32,
        ) -> (i32, i32) {
            let mut rng = thread_rng();
            let len = rng.gen_range(min_len..max_len + 1);
            let start = rng.gen_range(min_first..max_last - len + 1);
            (start, start + len - 1)
        }

        fn check_random_queries<I, F>(
            n: usize,
            num_queries: usize,
            max_last: i32,
            min_len: i32,
            max_len: i32,
            query_min_len: i32,
            query_max_len: i32,
            check: F,
        ) where
            I: IntWithMax,
            F: Fn(&COITree<u32, I>, &[IntervalNode<u32, I>], &mut [(i32, i32)]),
        {
            let min_first = 0;

            let mut b: Vec<IntervalNode<u32, I>> = Vec::with_capacity(n);
            for i in 0..n {
                let (first, last) = random_interval(min_first, max_last, min_len, max_len);
                b.push(IntervalNode::new(first, last, i as u32));
            }
            b.sort_unstable_by_key(|node| node.first);

            let a = COITree::new(&b);

            let mut queries: Vec<(i32, i32)> = Vec::with_capacity(num_queries);
            for _ in 0..num_queries {
                queries.push(random_interval(
                    min_first,
                    max_last,
                    query_min_len,
                    query_max_len,
                ));
            }

            check(&a, &b, &mut queries);
        }

        fn check_random_queries_default<I, F>(n: usize, num_queries: usize, check: F)
        where
            I: IntWithMax,
            F: Fn(&COITree<u32, I>, &[IntervalNode<u32, I>], &mut [(i32, i32)]),
        {
            let max_last = 1000000;
            let min_len = 20;
            let max_len = 2000;

            check_random_queries::<I, F>(
                n,
                num_queries,
                max_last,
                min_len,
                max_len,
                min_len,
                max_len,
                check,
            );
        }

        const CHECKS: [fn(&COITree<u32, usize>, &[IntervalNode<u32, usize>], &mut [(i32, i32)]);
            4] = [
            check_queries,
            check_count_queries,
            check_sorted_querent_queries,
            check_sorted_querent_unsorted_queries,
        ];

        #[test]
        fn query_empty_tree() {
            for check in &CHECKS {
                check_random_queries_default(0, 1000, check);
            }
            check_random_queries_default::<usize, _>(0, 1000, check_coverage);
        }

        #[test]
        fn query_small_trees() {
            for n in 1..16 {
                for check in &CHECKS {
                    check_random_queries_default(n, 1000, check);
                }
                check_random_queries_default::<usize, _>(n, 1000, check_coverage);
            }
        }

        #[test]
        fn query_medium_tree() {
            for check in &CHECKS {
                check_random_queries_default(10000, 1000, check);
            }
            check_random_queries_default::<usize, _>(10000, 1000, check_coverage);
        }

        #[test]
        fn query_singeton_intervals() {
            for check in &CHECKS {
                check_random_queries(10000, 1000, 1000, 1, 1, 1, 1, check);
                check_random_queries(10000, 1000, 1000, 1, 1, 10, 100, check);
            }
            check_random_queries::<usize, _>(10000, 1000, 1000, 1, 1, 1, 1, check_coverage);
            check_random_queries::<usize, _>(10000, 1000, 1000, 1, 1, 10, 100, check_coverage);
        }

        #[test]
        fn query_empty_intervals() {
            for check in &CHECKS {
                check_random_queries(10000, 1000, 1000, 0, 0, 0, 0, check);
                check_random_queries(10000, 1000, 1000, 0, 0, 10, 100, check);
            }
        }

        const CHECKS_U16: [fn(&COITree<u32, u16>, &[IntervalNode<u32, u16>], &mut [(i32, i32)]);
            4] = [
            check_queries,
            check_count_queries,
            check_sorted_querent_queries,
            check_sorted_querent_unsorted_queries,
        ];

        #[test]
        fn test_largest_tree() {
            // make sure nothing breaks when have the largest tree that
            // can be represented by an index type.

            let n = (u16::MAX - 1) as usize;
            for check in &CHECKS_U16 {
                check_random_queries_default::<u16, _>(n, 1000, check);
            }
            check_random_queries_default::<u16, _>(n, 1000, check_coverage);
        }

        #[test]
        #[should_panic]
        fn index_type_overflow() {
            // test that tree construction panics if there are more nodes than can be
            // indexed.

            let n = u16::MAX as usize;
            let mut b: Vec<IntervalNode<u32, u16>> = Vec::with_capacity(n);
            for i in 0..n {
                let (first, last) = random_interval(0, 10000, 10, 100);
                b.push(IntervalNode::new(first, last, i as u32));
            }

            let _tree: COITree<u32, u16> = COITree::new(&b);
        }
    }
}
