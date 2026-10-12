//! Advanced routines that act on NLib NURBS surfaces: extension, normal derivatives,
//! skinning, sweeping, Coons and Gordon surfaces.
//!
//! Upstream: `source/SMLib/src/NL_SrfAdv.cpp` (declarations in
//! `source/SMLib/inc/NL_SrfAdv.h`). The file is split over several plan units
//! (`nlib.srf_adv_p1`, `_p2`, ...); each lives in its own submodule.

pub mod p2;
