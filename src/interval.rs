use std::iter::IntoIterator;
use std::ops::{AddAssign, SubAssign};

pub trait GenericInterval<T>
where
    T: Clone,
{
    fn first(&self) -> i32;
    fn last(&self) -> i32;
    fn metadata(&self) -> &T;

    fn len(&self) -> i32 {
        0.max(self.last() - self.first() + 1)
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// An interval with associated metadata.
///
/// Intervals in `COITree` are treated as end-inclusive.
///
/// Metadata can be an arbitrary type `T`, but because nodes are stored in contiguous
/// memory, it may be better to store large metadata outside the node and
/// use a pointer or reference for the metadata.
///
/// # Examples
/// ```
/// use coitrees::Interval;
/// use coitrees::GenericInterval;
///
/// #[derive(Clone)]
/// struct MyMetadata {
///     chrom: String,
///     posstrand: bool
/// }
///
/// let some_interval = Interval{
///     first: 10, last: 24000,
///     metadata: MyMetadata{chrom: String::from("chr1"), posstrand: false}};
///
/// assert_eq!(some_interval.len(), 23991);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Interval<T>
where
    T: Clone,
{
    pub first: i32,
    pub last: i32,
    pub metadata: T,
}

impl<T> Interval<T>
where
    T: Clone,
{
    pub fn new(first: i32, last: i32, metadata: T) -> Interval<T> {
        Self {
            first,
            last,
            metadata,
        }
    }
}

impl<T> GenericInterval<T> for Interval<T>
where
    T: Clone,
{
    fn first(&self) -> i32 {
        self.first
    }

    fn last(&self) -> i32 {
        self.last
    }

    fn metadata(&self) -> &T {
        &self.metadata
    }
}

impl<'a, T> GenericInterval<T> for Interval<&'a T>
where
    T: Clone,
{
    fn first(&self) -> i32 {
        self.first
    }

    fn last(&self) -> i32 {
        self.last
    }

    fn metadata(&self) -> &T {
        self.metadata
    }
}

#[test]
fn test_interval_len() {
    fn make_interval(first: i32, last: i32) -> Interval<()> {
        Interval {
            first,
            last,
            metadata: (),
        }
    }

    assert_eq!(make_interval(1, -1).len(), 0);
    assert_eq!(make_interval(1, 0).len(), 0);
    assert_eq!(make_interval(1, 1).len(), 1);
    assert_eq!(make_interval(1, 2).len(), 2);
}

/// A trait facilitating COITree index types.
pub trait IntWithMax:
    TryInto<usize> + TryFrom<usize> + Copy + Default + PartialEq + Ord + AddAssign + SubAssign
{
    const MAX: Self;

    // typically the branch here should be optimized out, because we are
    // converting, e.g. a u32 to a usize on 64-bit system.
    #[inline(always)]
    fn to_usize(self) -> usize {
        match self.try_into() {
            Ok(x) => x,
            Err(_) => panic!("index conversion to usize failed"),
        }
    }

    #[inline(always)]
    fn from_usize(x: usize) -> Self {
        match x.try_into() {
            Ok(y) => y,
            Err(_) => panic!("index conversion from usize failed"),
        }
    }

    fn one() -> Self {
        Self::from_usize(1)
    }
}

impl IntWithMax for usize {
    const MAX: usize = usize::MAX;
}

impl IntWithMax for u32 {
    const MAX: u32 = u32::MAX;
}

impl IntWithMax for u16 {
    const MAX: u16 = u16::MAX;
}

/// Basic interval tree interface supported by each COITree implementation.
pub trait IntervalTree<'a> {
    type Metadata: Clone + 'a;
    type Index: IntWithMax;
    type Item: GenericInterval<Self::Metadata> + 'a;
    type Iter: Iterator<Item = Interval<&'a Self::Metadata>>;

    fn new<'b, U, V>(intervals: U) -> Self
    where
        U: IntoIterator<Item = &'b V>,
        V: GenericInterval<Self::Metadata> + 'b;

    fn len(&self) -> usize;

    fn is_empty(&self) -> bool;

    fn query<F>(&'a self, first: i32, last: i32, visit: F)
    where
        F: FnMut(&Self::Item);

    fn query_count(&self, first: i32, last: i32) -> usize;
    fn coverage(&self, first: i32, last: i32) -> (usize, usize);

    fn iter(&'a self) -> Self::Iter;
}

pub trait SortedQuerent<'a> {
    type Metadata: Clone + 'a;
    type Index: IntWithMax;
    type Item: GenericInterval<Self::Metadata> + 'a;
    type Iter: Iterator<Item = Interval<&'a Self::Metadata>>;
    type Tree: IntervalTree<
            'a,
            Metadata = Self::Metadata,
            Index = Self::Index,
            Item = Self::Item,
            Iter = Self::Iter,
        > + 'a;

    fn new(tree: &'a Self::Tree) -> Self;

    fn query<F>(&mut self, first: i32, last: i32, visit: F)
    where
        F: FnMut(&Self::Item);
}
