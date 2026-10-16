//! The topology intersector `SmTopologyIntersector`: creates the edges and vertices
//! corresponding to the intersections of two B-reps in each of them, and keeps the
//! relationships ("mates") between the topology that represents common geometry.
//!
//! Upstream: `source/SMLib/src/SmTopologyIntersector.cpp`,
//! `source/SMLib/inc/SmTopologyIntersector.h`. The file is split over several plan
//! units (`boolean.topology_intersector_p1`, `_p2`, `_p3`); each lives in its own
//! submodule and adds methods to the types of [`p1`].

pub mod p1;

pub use p1::{
    EdgeFaceIntersect, ExistingSsi, FaceFaceCoincidence, FaceFaceIntersect, ImprintFace,
    IntersectionResults, IntersectionTable, IntersectorBrep, Reason, SplitTwoEdges,
    TopologyIntersector, VertexFaceIntersect,
};
