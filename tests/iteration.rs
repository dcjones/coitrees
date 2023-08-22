mod tests {
    use coitrees::{COITree, Interval, IntervalTree};

    extern crate rand;
    use rand::{thread_rng, Rng};
    use std::collections::HashSet;

    fn random_interval(min_first: i32, max_last: i32, min_len: i32, max_len: i32) -> (i32, i32) {
        let mut rng = thread_rng();
        let len = rng.gen_range(min_len..(max_len + 1));
        let start = rng.gen_range(min_first..(max_last - len + 1));

        (start, start + len - 1)
    }
    fn check_iteration(n: usize) {
        let min_first = 0;
        let max_last = 10000000;

        let min_len = 1;
        let max_len = 10000;

        let b = (0..n)
            .map(|i| {
                let (first, last) = random_interval(min_first, max_last, min_len, max_len);
                Interval {
                    first,
                    last,
                    metadata: i,
                }
            })
            .collect::<Vec<Interval<usize>>>();

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
