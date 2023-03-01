//! # COITrees
//! `coitrees` implements a fast static interval tree data structure with genomic
//! data in mind.
//!
//! The data structure used a fairly standard interval tree, but with nodes stored
//! in van Emde Boas layout, which improves average cache locality, and thus
//! query performance. The downside is that building the tree is more expensive
//! so a relatively large number of queries needs to made for it to break even.
//!
//! The data structure `COITree` is constructed with an array of `IntervalNode`
//! structs which store integer, end-inclusive intervals along with associated
//! metadata. The tree can be queried directly for coverage or overlaps, or
//! through the intermediary `SortedQuerenty` which keeps track of some state
//! to accelerate overlaping queries.

#[cfg(feature = "avx")]
pub mod avx;
#[cfg(feature = "avx")]
pub use avx::*;

#[cfg(feature = "default")]
mod default;
#[cfg(feature = "default")]
pub use default::*;

#[cfg(feature = "neon")]
pub mod neon;
#[cfg(feature = "neon")]
pub use neon::*;
