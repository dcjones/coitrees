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

#[cfg(all(target_feature = "avx2", feature = "default"))]
mod avx;

#[cfg(all(target_feature = "avx2", feature = "default"))]
pub use avx::*;

#[cfg(all(feature = "nosimd", not(feature = "default")))]
mod nosimd;
#[cfg(all(feature = "nosimd", not(feature = "default")))]
pub use nosimd::*;

#[cfg(all(target_feature = "neon", feature = "default"))]
mod neon;
#[cfg(all(target_feature = "neon", feature = "default"))]
pub use neon::*;
