//! B-rep operations, split into independently owned port units.
//!
//! Upstream: `source/SMLib/inc/SmBrep.h`, `source/SMLib/src/SmBrep.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`p5`] | B-rep editing and transformation operations, part 5 |
//!
//! Deviations: the arena record remains [`super::Brep`]. Supporting notification
//! tracking integration is documented in part 5 and `plan/modules.toml`.

pub mod p5;
