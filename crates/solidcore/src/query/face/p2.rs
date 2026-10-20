//! `SmFace` queries, part 2: classification of a parameter point
//! (`SmFace::PointClassify`) and of a curve (`SmFace::CurveOnClassify`) against the
//! face boundary, the two-face parameter classification of the proto topology
//! (`sm_ClassifyParamOnFaces`), and the point containment of a proto vertex.
//!
//! Upstream: `source/SMLib/src/SmFace.cpp`, `source/SMLib/src/SmProtoTopology.cpp`,
//! `source/SMLib/inc/SmProtoTopology.h`.
//!
//! | Rust                                  | Upstream                                          |
//! |---------------------------------------|---------------------------------------------------|
//! | [`point_classify`], [`PointClassifyOptions`] | `SmFace::PointClassify`                    |
//! | [`curve_on_classify`], [`CurveOnClassifyObjects`] | `SmFace::CurveOnClassify`             |
//! | [`classify_param_on_faces`]           | `sm_ClassifyParamOnFaces` (file static)           |
//! | [`proto_vertex_distance_to`]          | `SmProtoVertex::DistanceTo`                       |
//! | [`proto_vertex_capture_distance`]     | `SmProtoVertex::GetCaptureDistance`               |
//! | [`proto_vertex_contains_point`]       | `SmProtoVertex::ContainsPoint`                    |
//! | [`ProtoTopologyManager::capture_distance`] | `SmProtoTopologyManager::GetCaptureDistance` (already in `topology::proto_topology`) |
//! | [`face_vertices`]                     | `SmFace::GetVertices(SmTArray&, ..)` (CO-PORTED, `topology.face_p4`) |
//! | `ClassObject for EntityRef`           | the `SmObject*` of a classification of B-rep topology |
//!
//! Deviations:
//! - **Object model.** The face is a record of a [`Brep`]; the functions take the
//!   B-rep and the face's id. Classification objects are [`EntityRef`]s (the face, edge
//!   or vertex). `SmFace::PointClassify`'s status is the result: `Ok(())` for
//!   `SM_SUCCESS`, [`SolidError::Upstream`] `SM_ERR` when no loop classifies the point
//!   (the classification is still updated, as upstream).
//! - **Unused parameters.** `PointClassify`'s `sSrcZoneTol3d`,
//!   `bDoBoundaryIntersections` and `bUVSpaceOkay` are `NotUsed` upstream; they are
//!   kept in [`PointClassifyOptions`] and ignored. `bUseUVClassification` only chooses
//!   between two identical `SM_ERR` returns.
//! - **Trimmed surface cache.** `SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF, surface)`
//!   has no counterpart: B-rep surfaces carry no object cache in the port, so the
//!   caller passes the face surface's [`TrimSrfCache`] (if any) in
//!   [`PointClassifyOptions::trim_cache`]. The test
//!   `GetSurface()->GetFace() == this` is "the face's surface has a face owner"
//!   (the owner of a face's surface is that face).
//! - **Out-of-bounds evaluation.** `SmTemporaryChangeValue` on the surface's
//!   `GetOutOfBoundsEnabled()` becomes a copy of the surface with out-of-bounds
//!   evaluation enabled (B-spline kinds only), used for the evaluation and the loop
//!   containment tests; the B-rep is not changed.
//! - **Loop containment.** [`crate::query::loop::contains_uv_point`] walks the loop in
//!   the B-rep, which therefore is taken mutably (`SmEdgeuse::GetOrCreateUVTrimCurve`
//!   attaches the UV trim curves the loops lack); [`point_classify`] and
//!   [`classify_param_on_faces`] take `&mut Brep`. The surface it is given is an
//!   owned copy of the face's. Its zone tolerance (`SmTol::GetZoneTol3d(pSurface)`) is
//!   the loop face's tolerance. An edge use without edge or an edge without curve is
//!   [`SolidError::InvalidInput`] (upstream dereferences them).
//! - **Ignored statuses.** `Evaluate1stDerivatives` of the surface and `EvaluatePoint`
//!   of the on-boundary edge curve: on failure the point stays `SmPoint3d()`
//!   ([`Vec3::UNDEF`]); [`SolidError::Unported`] still propagates. The `PointClassify`
//!   statuses `sm_ClassifyParamOnFaces` ignores are ignored the same way.
//!   Part-9 boundary views also ignore ordinary `ComputeUVPoint` failures and keep
//!   the initialized [`Point2::UNDEF`]; unported errors still propagate.
//! - **Null dereferences** upstream (a face, loop, edge use or vertex use without the
//!   object it points to) are [`SolidError::InvalidInput`]; a face without surface is
//!   `NER` ([`SolidError::null_pointer`]).
//! - **`sm_ClassifyParamOnFaces`** takes the faces' tolerances for
//!   `SmTol::GetZoneTol3d_ForXSectResult(pFace1, pFace2)` (`SM_USE_OLDTOL`:
//!   `SmFace::GetTolerance`).
//! - **`SmProtoVertex`** keeps no pointer to its manager; the manager is passed.
//! - **`CurveOnClassify`** stores handles of the caller's choosing in the curve
//!   classification; [`CurveOnClassifyObjects`] makes the face and surface handles and
//!   resolves curve and loop handles, and the caller's [`crate::intersect::point_class::ClassGeometry`] resolves them
//!   for `ClassifyFaceInterval` (`lFaceIndx` left at its default).
//!   `SmCurveClassification::ReSet` is [`CurveClassification::reset`]; its callees
//!   `IntersectTSBound` / `IntersectEdges` (part 9, wired through borrowed boundary
//!   views) and `SmFace::GetEdges` (`topology.face_p4`) retain their dependency
//!   boundaries. `pbOptModifiedCurve` is the result. `IntersectVertices` and
//!   `ClassifyFaceInterval` call their landed implementations through the caller's
//!   classification handles and [`crate::intersect::point_class::ClassGeometry`] callbacks.
//! - Boundary views use required handle mappings and the same [`crate::intersect::point_class::ClassGeometry`]
//!   implementation as part 5. Lazy UV and gap calculations retain their writes
//!   in the original mutable B-rep through the same owner provider.
//! - No recursion.

use crate::curve::curve::Curve;
use crate::error::{ignored, EntityRef, Result, SolidError, SM_ERR};
use crate::geom::extent2d::Extent2d;
use crate::geom::point_set::PointContainment;
use crate::geom::tree::{NodeClass, NodeData};
use crate::intersect::curve_class::p1::{src_zone_tol_3d_of, CurveClassification};
use crate::intersect::curve_class::p2::ResetOptions;
use crate::intersect::curve_class::p9::{
    BoundaryEdge, BoundaryEdgeuse, BoundaryFace, BoundaryLoop, BoundaryVertex,
    IntersectEdgesOptions,
};
use crate::intersect::point_class::{
    ClassObject, ClassObjectKind, PointClassType, PointClassification,
};
use crate::kernel::Kernel;
use crate::math::{max2, Point2, Point3, Vec3, BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQRT};
use crate::query::r#loop::{contains_uv_point, UvContainment};
use crate::surface::surface::{Surface, SurfaceKind};
use crate::surface::surface_cache::p1::SurfaceCacheContext;
use crate::surface::surface_cache::p1::{FaceContainment, SurfaceCacheOptions};
use crate::surface::trim_srf_cache::TrimSrfCache;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::face::p2 as fp2;
use crate::topology::loopuse;
use crate::topology::proto_topology::{ProtoTopologyManager, ProtoVertex};
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, LoopId, VertexId, VertexuseId};

const OP_CLASSIFY: &str = "query::face::p2::point_classify";
const OP_CURVE: &str = "query::face::p2::curve_on_classify";
const OP_VERTICES: &str = "query::face::p2::face_vertices";

impl ClassObject for EntityRef {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Face(_) => ClassObjectKind::Face,
            Self::Edge(_) => ClassObjectKind::Edge,
            Self::Vertex(_) => ClassObjectKind::Vertex,
            Self::Loop(_) => ClassObjectKind::Loop,
            Self::Region(_) => ClassObjectKind::Region,
            Self::Faceuse(_) => ClassObjectKind::Faceuse,
            Self::Loopuse(_) => ClassObjectKind::Loopuse,
            Self::Edgeuse(_) => ClassObjectKind::Edgeuse,
            Self::Vertexuse(_) => ClassObjectKind::Vertexuse,
            Self::Brep(_) | Self::Shell(_) => ClassObjectKind::Other,
        }
    }
}

fn null(operation: &'static str, reason: &'static str) -> SolidError {
    SolidError::invalid_input(operation, reason)
}

// ---- co-ported topology access ----------------------------------------------------------

/// `SmFace::GetVertices(SmTArray&, NULL, bTgtLoopsOnly, pOptTgtLoops)` (CO-PORTED from
/// unit `topology.face_p4`): the vertices of the loop uses of the upward face use
/// selected as [`eu1::faceuse_loopuses`] selects them ([`loopuse::vertices`]: the vertex
/// of a vertex loop use, else the distinct start vertices of its edge uses in CCW
/// order), appended (a vertex shared by two loops is listed twice). None without face
/// use.
///
/// # Errors
///
/// Those of [`eu1::faceuse_loopuses`] and [`loopuse::vertices`];
/// [`SolidError::InvalidInput`] for a vertex use without vertex (upstream dereferences
/// it).
pub fn face_vertices(
    b: &Brep,
    f: FaceId,
    tgt_loops_only: bool,
    tgt_loops: Option<&[LoopId]>,
) -> Result<Vec<VertexId>> {
    // init output
    let mut out = Vec::new();

    // for every upwardFaceuse->Loopuse
    let Some(fu) = b.face(f)?.faceuse else {
        return Ok(out);
    };
    for lu in eu1::faceuse_loopuses(b, fu, tgt_loops_only, tgt_loops)? {
        // accumulate vertices
        for v in loopuse::vertices(b, lu, None)? {
            out.push(v.ok_or_else(|| null(OP_VERTICES, "vertex use without vertex"))?);
        }
    }
    Ok(out)
}

/// A copy of `surface` with out-of-bounds evaluation enabled, when it is an
/// `SmBSplineSurface` (B-spline, plane, surface of revolution) without it; `None`
/// when nothing changes (`SmTemporaryChangeValue` on `GetOutOfBoundsEnabled()`).
fn with_out_of_bounds_enabled(surface: &Surface) -> Option<Surface> {
    let enabled = match surface.kind() {
        SurfaceKind::BSpline(s) => s.is_out_of_bounds_enabled(),
        SurfaceKind::Plane(p) => p.base().is_out_of_bounds_enabled(),
        SurfaceKind::SurfOfRevolution(s) => s.base().is_out_of_bounds_enabled(),
        SurfaceKind::OffsetSurface(_)
        | SurfaceKind::CurveBounded(_)
        | SurfaceKind::Step(_)
        | SurfaceKind::SrfInVolume(_) => {
            return None;
        }
    };
    if enabled {
        return None;
    }
    let mut copy = surface.clone();
    match copy.kind_mut() {
        SurfaceKind::BSpline(s) => *s.out_of_bounds_enabled_mut() = true,
        SurfaceKind::Plane(p) => *p.out_of_bounds_enabled_mut() = true,
        SurfaceKind::SurfOfRevolution(s) => *s.base_mut().out_of_bounds_enabled_mut() = true,
        SurfaceKind::OffsetSurface(_)
        | SurfaceKind::CurveBounded(_)
        | SurfaceKind::Step(_)
        | SurfaceKind::SrfInVolume(_) => {}
    }
    // the copy constructor resets the owner; the temporary change keeps it
    copy.set_owner(surface.owner());
    Some(copy)
}

// ---- SmFace::PointClassify --------------------------------------------------------------

/// Decision from the original retained trimmed-surface cache.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FacePointCacheDecision {
    /// No classified leaf answered the point; continue the normal face query.
    NoDecision,
    /// The first classified leaf is outside.
    Outside,
    /// The first classified leaf is inside, retaining its exact optional face.
    Inside(Option<EntityRef>),
}

/// Narrow mutable access to a caller's retained cache, independent of PC handles.
/// Cache construction receives the original B-rep only for this call.
pub trait FacePointCache: std::fmt::Debug {
    /// Native ready flags, inspected before the surface's owner query.
    fn ready(&self) -> bool;

    /// Build the actual tree if absent, then inspect its leaves in native order.
    ///
    /// # Errors
    /// Propagates construction, original topology and exact identity-map errors.
    fn classify_ready(
        &mut self,
        k: &Kernel,
        brep: &mut Brep,
        face: FaceId,
        uv: Point2,
        keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
        owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
    ) -> Result<FacePointCacheDecision>;
}

/// Borrow a retained cache with its own handle type; no tree or dispatch is copied.
/// Both mappings name the same original operand and preserve complete PC state.
pub struct MappedFacePointCache<'a, C> {
    cache: &'a mut TrimSrfCache<C>,
    to_entity: &'a dyn Fn(C) -> Result<EntityRef>,
    from_entity: &'a dyn Fn(EntityRef) -> Result<C>,
}

impl<C> std::fmt::Debug for MappedFacePointCache<'_, C> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("MappedFacePointCache")
            .finish_non_exhaustive()
    }
}

impl<'a, C> MappedFacePointCache<'a, C> {
    /// Lend the original cache and exact bidirectional topology identity mappings.
    pub fn new(
        cache: &'a mut TrimSrfCache<C>,
        to_entity: &'a dyn Fn(C) -> Result<EntityRef>,
        from_entity: &'a dyn Fn(EntityRef) -> Result<C>,
    ) -> Self {
        Self {
            cache,
            to_entity,
            from_entity,
        }
    }
}

struct FaceCacheContext<'a, 'b, C> {
    context: super::p1::FacePointContext<'a>,
    face: FaceId,
    to_entity: &'b dyn Fn(C) -> Result<EntityRef>,
    from_entity: &'b dyn Fn(EntityRef) -> Result<C>,
}

fn cache_topology_object(entity: EntityRef) -> Result<crate::topology::brep_cache::CacheObject> {
    use crate::topology::brep_cache::CacheObject;
    match entity {
        EntityRef::Face(face) => Ok(CacheObject::Face(face)),
        EntityRef::Edge(edge) => Ok(CacheObject::Edge(edge)),
        EntityRef::Vertex(vertex) => Ok(CacheObject::Vertex(vertex)),
        _ => Err(SolidError::invalid_input(
            OP_CLASSIFY,
            "unsupported trim-cache topology handle",
        )),
    }
}

fn cache_topology_entity(object: crate::topology::brep_cache::CacheObject) -> Result<EntityRef> {
    use crate::topology::brep_cache::CacheObject;
    match object {
        CacheObject::Face(face) => Ok(EntityRef::Face(face)),
        CacheObject::Edge(edge) => Ok(EntityRef::Edge(edge)),
        CacheObject::Vertex(vertex) => Ok(EntityRef::Vertex(vertex)),
        _ => Err(SolidError::invalid_input(
            OP_CLASSIFY,
            "non-topology trim-cache classification",
        )),
    }
}

impl<C: Copy + PartialEq> SurfaceCacheContext<C> for FaceCacheContext<'_, '_, C> {
    fn surface_face(&self, surface: &Surface, _: Option<C>) -> Result<Option<C>> {
        use crate::topology::brep_cache::CacheObject;
        self.context
            .surface_face(surface, Some(CacheObject::Surface(self.face)))?
            .map(|object| (self.from_entity)(cache_topology_entity(object)?))
            .transpose()
    }
    fn surface_has_direct_face_owner(&self, surface: &Surface, _: Option<C>) -> Result<bool> {
        use crate::topology::brep_cache::CacheObject;
        self.context
            .surface_has_direct_face_owner(surface, Some(CacheObject::Surface(self.face)))
    }
    fn surface_zone_tol_3d(&self, k: &Kernel, surface: &Surface, _: Option<C>) -> Result<f64> {
        use crate::topology::brep_cache::CacheObject;
        self.context
            .surface_zone_tol_3d(k, surface, Some(CacheObject::Surface(self.face)))
    }
    fn surface_owner_zone_tol_3d(
        &self,
        k: &Kernel,
        surface: &Surface,
        _: Option<C>,
    ) -> Result<f64> {
        use crate::topology::brep_cache::CacheObject;
        self.context
            .surface_owner_zone_tol_3d(k, surface, Some(CacheObject::Surface(self.face)))
    }
    fn object_zone_tol_3d(&self, k: &Kernel, object: Option<C>) -> Result<f64> {
        let object = object
            .map(|object| cache_topology_object((self.to_entity)(object)?))
            .transpose()?;
        self.context.object_zone_tol_3d(k, object)
    }
    fn face_brep_editing_enabled(&self, face: C) -> Result<bool> {
        self.context
            .face_brep_editing_enabled(cache_topology_object((self.to_entity)(face)?)?)
    }
    fn set_face_brep_editing_enabled(&mut self, face: C, editing: bool) -> Result<()> {
        self.context
            .set_face_brep_editing_enabled(cache_topology_object((self.to_entity)(face)?)?, editing)
    }
    fn face_create_uv_trim_curves(&mut self, k: &Kernel, face: C) -> Result<()> {
        self.context
            .face_create_uv_trim_curves(k, cache_topology_object((self.to_entity)(face)?)?)
    }
    fn face_point_classify(
        &mut self,
        k: &Kernel,
        face: C,
        uv: Point2,
        tolerance: f64,
        do_boundaries: bool,
        use_uv: bool,
        classification: &mut PointClassification<C>,
    ) -> Result<()> {
        let mut local = classification
            .try_map_objects(|object| cache_topology_object((self.to_entity)(object)?))?;
        let result = self.context.face_point_classify(
            k,
            cache_topology_object((self.to_entity)(face)?)?,
            uv,
            tolerance,
            do_boundaries,
            use_uv,
            &mut local,
        );
        *classification =
            local.try_map_objects(|object| (self.from_entity)(cache_topology_entity(object)?))?;
        result
    }
}

fn classify_retained_cache<C: Copy + PartialEq>(
    cache: &mut TrimSrfCache<C>,
    k: &Kernel,
    b: &mut Brep,
    face: FaceId,
    uv: Point2,
    keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
    owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
    to_entity: &dyn Fn(C) -> Result<EntityRef>,
    from_entity: &dyn Fn(EntityRef) -> Result<C>,
) -> Result<FacePointCacheDecision> {
    if cache.base().tree().is_none() {
        let mut context = FaceCacheContext {
            context: super::p1::FacePointContext::new(b, keys, owner_of),
            face,
            to_entity,
            from_entity,
        };
        cache.build_tree_with_context(k, &mut context)?;
    }
    let nodes = cache
        .base()
        .find_uv_nodes_built(k, &Extent2d::from_point(uv), false)?;
    let tree = cache
        .base()
        .tree()
        .ok_or_else(|| null(OP_CLASSIFY, "trimmed surface cache without tree"))?;
    let mut leaf_nodes = Vec::new();
    for node in nodes {
        if tree.node(node)?.child1.is_none() {
            leaf_nodes.push(node);
        }
    }
    for node in leaf_nodes {
        let NodeData::Surface(aux) = &tree.node(node)?.data else {
            return Err(null(OP_CLASSIFY, "leaf node without surface data"));
        };
        if aux.node_class == NodeClass::Outside {
            return Ok(FacePointCacheDecision::Outside);
        }
        if aux.node_class == NodeClass::Inside {
            return Ok(FacePointCacheDecision::Inside(
                aux.face.map(to_entity).transpose()?,
            ));
        }
    }
    Ok(FacePointCacheDecision::NoDecision)
}

impl FacePointCache for TrimSrfCache<EntityRef> {
    fn ready(&self) -> bool {
        self.base().have_t_surface_cache
            && self.base().face_containment_done == FaceContainment::Done
            && !self.base().face_was_modified
    }
    fn classify_ready(
        &mut self,
        k: &Kernel,
        b: &mut Brep,
        face: FaceId,
        uv: Point2,
        keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
        owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
    ) -> Result<FacePointCacheDecision> {
        classify_retained_cache(self, k, b, face, uv, keys, owner_of, &Ok, &Ok)
    }
}

impl<C: Copy + PartialEq> FacePointCache for MappedFacePointCache<'_, C> {
    fn ready(&self) -> bool {
        self.cache.base().have_t_surface_cache
            && self.cache.base().face_containment_done == FaceContainment::Done
            && !self.cache.base().face_was_modified
    }
    fn classify_ready(
        &mut self,
        k: &Kernel,
        b: &mut Brep,
        face: FaceId,
        uv: Point2,
        keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
        owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
    ) -> Result<FacePointCacheDecision> {
        classify_retained_cache(
            self.cache,
            k,
            b,
            face,
            uv,
            keys,
            owner_of,
            self.to_entity,
            self.from_entity,
        )
    }
}

/// The optional arguments of [`point_classify`]; [`Default`] gives upstream's
/// defaults for the defaulted ones (`bTgtLoopsOnly = FALSE`, `pOptTgtLoops = NULL`,
/// `bUVSpaceOkay = TRUE`) and `0` / `false` for the others.
#[derive(Debug)]
#[allow(clippy::struct_excessive_bools)] // upstream's SmBoolean arguments
pub struct PointClassifyOptions<'a> {
    /// `sSrcZoneTol3d` (`NotUsed` upstream).
    pub src_zone_tol_3d: f64,
    /// `bDoBoundaryIntersections` (`NotUsed` upstream).
    pub do_boundary_intersections: bool,
    /// `bUseUVClassification` (no longer used: both of its branches return
    /// `SM_ERR`).
    pub use_uv_classification: bool,
    /// `bTgtLoopsOnly`: only classify against [`PointClassifyOptions::tgt_loops`]
    /// (default `false`: all loops of the face).
    pub tgt_loops_only: bool,
    /// `pOptTgtLoops`: with `tgt_loops_only`, the loops to classify against (`None` or
    /// empty: the outer loop).
    pub tgt_loops: Option<&'a [LoopId]>,
    /// `bUVSpaceOkay` (`NotUsed` upstream; upstream default `TRUE`).
    pub uv_space_okay: bool,
    /// The face surface's trimmed surface cache
    /// (`SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF, GetSurface())`), `None` when the
    /// surface has none.
    pub trim_cache: Option<&'a mut dyn FacePointCache>,
}

impl Default for PointClassifyOptions<'_> {
    fn default() -> Self {
        Self {
            src_zone_tol_3d: 0.0,
            do_boundary_intersections: false,
            use_uv_classification: false,
            tgt_loops_only: false,
            tgt_loops: None,
            uv_space_okay: true,
            trim_cache: None,
        }
    }
}

/// `SmFace::PointClassify`: classifies the parameter point `uv` of face `f` into `pc`:
/// [`PointClassType::Face`] with the face when inside, [`PointClassType::Edge`] /
/// [`PointClassType::Vertex`] with the edge (and its parameter) or vertex when on the
/// boundary, [`PointClassType::Unknown`] when outside (a point outside the face's UV
/// domain, outside the outer loop or inside an inner loop). The face surface's
/// trimmed surface cache answers first when it is up to date; then the point's
/// surface position is tested against every face vertex (within the larger of the
/// face, vertex and `pc` tolerances), then the point against every loop
/// ([`contains_uv_point`]). With `tgt_loops_only`, only the selected loops (and their
/// vertices) are used. `pc` keeps its other members; its tolerance
/// (`GetTolerance`) enters the vertex and loop tolerances.
///
/// # Errors
///
/// [`SolidError::Upstream`] `SM_ERR` when the point is neither outside nor inside nor
/// on the checked loops (a loop that cannot classify it, or only inner loops checked
/// without the point being outside all of them); `SM_ERR_NULL_POINTER` for a face
/// without surface; those of the trimmed surface cache search, of [`face_vertices`]
/// and [`contains_uv_point`]; [`SolidError::InvalidInput`] for an
/// on-boundary edge without curve.
/// This plain entry point retains the arena-identity convention. Call
/// [`point_classify_with_context`] for registered operand namespaces and owners.
pub fn point_classify(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    uv: Point2,
    pc: &mut PointClassification<EntityRef>,
    opts: PointClassifyOptions<'_>,
) -> Result<()> {
    let mut owners = std::collections::BTreeMap::new();
    if opts.trim_cache.is_some() {
        for (edge, record) in b.edges().iter() {
            if let Some(owner) = record.curve.as_ref().and_then(Curve::owner) {
                owners.insert(
                    crate::topology::brep_cache::arena_identity(
                        crate::topology::brep_cache::CacheObject::Curve(edge),
                    ),
                    owner,
                );
            }
        }
        for (face, record) in b.faces().iter() {
            if let Some(owner) = record.surface.as_ref().and_then(Surface::owner) {
                owners.insert(
                    crate::topology::brep_cache::arena_identity(
                        crate::topology::brep_cache::CacheObject::Surface(face),
                    ),
                    owner,
                );
            }
        }
    }
    point_classify_with_context(
        k,
        b,
        f,
        uv,
        pc,
        opts,
        &crate::topology::brep_cache::arena_identity,
        &|object| owners.get(&object).copied(),
    )
}

/// [`point_classify`] over the caller's original retained cache, whose stored
/// handles are translated only for topology queries and cached output fields.
/// Construction uses the exact caller keys and owner links on the original B-rep.
/// No cache copy is constructed; lazy tree writes and subdivision dispatch stay
/// in the same cache on success and on error.
///
/// # Errors
/// Those of [`point_classify`] and the operand-scoped cache handle translator.
pub fn point_classify_with_context(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    uv: Point2,
    pc: &mut PointClassification<EntityRef>,
    opts: PointClassifyOptions<'_>,
    keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
    owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
) -> Result<()> {
    let PointClassifyOptions {
        tgt_loops_only,
        tgt_loops,
        use_uv_classification,
        trim_cache,
        ..
    } = opts;
    let face = b.face(f)?;
    let this = EntityRef::Face(f);
    let (face_tolerance, face_uv_domain) = (face.tolerance, face.uv_domain);

    // init output - SM_PC_UNKNOWN == 'Outside' Face
    pc.set_class_object(PointClassType::Unknown, None);

    // No Work: TargetPoint outside of Face->UVDomain
    let scaled_zero = EFF_ZERO_SQRT * (1.0 + face_uv_domain.size().length());
    if !face_uv_domain.contains_point(uv, scaled_zero) {
        // TargetPoint is outside of Face
        return Ok(());
    }

    let has_surface = face.surface.is_some();
    if !tgt_loops_only && has_surface {
        if let Some(cache) = trim_cache {
            if cache.ready()
                && b.face(f)?.surface.as_ref().map_or(Ok(false), |surface| {
                    super::p1::surface_face_with_context(b, surface, keys, owner_of)
                        .map(|owner| owner == Some(f))
                })?
            {
                match cache.classify_ready(k, b, f, uv, keys, owner_of)? {
                    FacePointCacheDecision::NoDecision => {}
                    FacePointCacheDecision::Outside => {
                        pc.set_class_object(PointClassType::Unknown, None);
                        return Ok(());
                    }
                    FacePointCacheDecision::Inside(face) => {
                        pc.set_class_object(PointClassType::Face, face);
                        pc.set_gap_3d(0.0);
                        pc.set_pre_snap_param(BIG_DOUBLE);
                        pc.set_uv_param(uv);
                        return Ok(());
                    }
                }
            }
        }
    }

    // arrive here when PointClassify by TrimSrfCache failed or when bTgtLoopsOnly

    // init
    pc.set_class_object(PointClassType::Unknown, None);

    // We must enable out-of-bounds evaluation, otherwise a point outside the surface
    // domain will be classified on a boundary.
    let surface = b
        .face(f)?
        .surface
        .as_ref()
        .ok_or(SolidError::null_pointer(OP_CLASSIFY))?;
    // (an owned surface: the loop classification edits the B-rep while it walks)
    let surface: Surface = if let Some(copy) = with_out_of_bounds_enabled(surface) {
        copy
    } else {
        // (the copy constructor resets the owner, which the face's surface keeps)
        let mut copy = surface.clone();
        copy.set_owner(surface.owner());
        copy
    };
    let surface = &surface;

    // Eval Surface at TargetPoint (status ignored)
    let surf_pt = ignored(surface.evaluate_1st_derivatives(uv, false, false))?
        .map_or(Vec3::UNDEF, |(p, _, _)| p);

    {
        // begin Scope for TargetPoint/Vertex intersection
        let verts = face_vertices(b, f, tgt_loops_only, tgt_loops)?;

        // see if the target point is on a vertex location
        for v in verts {
            let vert = b.vertex(v)?;
            let mut vert_tol = max2(face_tolerance, vert.tolerance);
            vert_tol = max2(vert_tol, pc.tolerance());

            if vert.point.closer_than(vert_tol, surf_pt) {
                pc.set_class_object(PointClassType::Vertex, Some(EntityRef::Vertex(v)));

                // Need to calculate the deviation.
                pc.set_deviation(vert.point.distance(surf_pt));
                pc.set_src_zone_tol_3d(vert_tol);
                pc.set_pre_snap_param(BIG_DOUBLE);
                return Ok(());
            }
        } // end iter every face vertex
    }

    // arrive here when TargetPoint is not on a vertex,
    // next: check TargetPoint Loop containment.
    let mut type_inner = PointContainment::Unknown;
    let mut type_outer = PointContainment::Unknown;
    let mut checked_outer_loop = false;
    let outer_loop = fp2::outer_loops(b, f)?.first().copied();
    let loops = fp2::loops(b, f, tgt_loops_only, tgt_loops)?;

    // no work - no loops to check
    if loops.is_empty() {
        return Ok(());
    }

    // Get tolerance
    let tol_3d = max2(face_tolerance, pc.tolerance());

    // For every loop: check containment
    for l in loops {
        let is_outer = Some(l) == outer_loop;

        // check TargetPoint containment
        let UvContainment {
            containment,
            on_edge,
            edge_param,
            ..
        } = contains_uv_point(k, b, l, Some(surface), uv, tol_3d)?;
        // when checking an outerloop
        if is_outer {
            type_outer = containment;

            // remember outer loop has been checked
            checked_outer_loop = true;

            // Exit Case: when outside the outer loop - outside the face
            if type_outer == PointContainment::Outside {
                return Ok(());
            }
            // Exit Case: when on outer loop - on a Face->Edge - Get Deviation
            else if type_outer == PointContainment::OnBoundary {
                return set_on_edge(b, pc, on_edge, edge_param, surf_pt, tol_3d);
            }
            // arrive here when TargetPoint is 'in' OuterLoop
        } else {
            // when checking an InnerLoop branch
            type_inner = containment;

            // Exit Case: Inside an InnerLoop = Outside Face
            if type_inner == PointContainment::Inside {
                return Ok(());
            }
            // Exit Case: On an InnerLoop bound = On a Face->Edge - Get Deviation
            else if type_inner == PointContainment::OnBoundary {
                return set_on_edge(b, pc, on_edge, edge_param, surf_pt, tol_3d);
            }
            // arrive here when TargetPoint is 'Out' of this InnerLoop
        }
    } // end iter every Loop to Check

    // When TargetPoint is 'in' Face
    if type_outer == PointContainment::Inside
        || (tgt_loops_only && !checked_outer_loop && type_inner == PointContainment::Outside)
    {
        // inside Face
        pc.set_class_object(PointClassType::Face, Some(this));
        pc.set_pre_snap_param(BIG_DOUBLE);
        pc.set_deviation(0.0);
        pc.set_src_zone_tol_3d(tol_3d);
        pc.set_uv_param(uv);
        return Ok(());
    }

    // PointClassify() fails for BadDatabase loops that cross seam boundaries
    let _ = use_uv_classification;
    Err(SolidError::Upstream {
        operation: OP_CLASSIFY,
        code: SM_ERR,
        entity: Some(this),
    })
}

/// The on-boundary exit of [`point_classify`]: the edge, its parameter, the gap from
/// the edge point to the surface point and the loop tolerance.
fn set_on_edge(
    b: &Brep,
    pc: &mut PointClassification<EntityRef>,
    on_edge: Option<EdgeId>,
    edge_param: f64,
    surf_pt: Point3,
    tol_3d: f64,
) -> Result<()> {
    pc.set_class_object(PointClassType::Edge, on_edge.map(EntityRef::Edge));
    pc.set_t_param(edge_param, false);

    // Need to calculate the deviation.
    let on_edge = on_edge.ok_or_else(|| null(OP_CLASSIFY, "on-boundary without edge"))?;
    let curve = b
        .edge(on_edge)?
        .curve
        .as_ref()
        .ok_or_else(|| null(OP_CLASSIFY, "edge without curve"))?;
    let edge_pt = ignored(curve.evaluate_point(edge_param))?.unwrap_or(Vec3::UNDEF);
    pc.set_deviation(edge_pt.distance(surf_pt));
    pc.set_src_zone_tol_3d(tol_3d);
    pc.set_pre_snap_param(BIG_DOUBLE);
    Ok(())
}

// ---- sm_ClassifyParamOnFaces ----------------------------------------------------------------

/// `sm_ClassifyParamOnFaces`: which of two faces the parameter point `uv` lies on: `1`
/// when it is in `domain1` only or classifies in `face1`, `2` when it is in `domain2`
/// only or classifies in `face2`, `0` when it is on the boundary of both or classifies
/// cleanly to neither. `src_zone_tol_3d` above `SM_EFF_ZERO` is the classification
/// tolerance, else the larger face tolerance. The statuses of the face classifications
/// are ignored.
///
/// # Errors
///
/// [`SolidError::Unported`] from the classifications; their
/// [`SolidError::InvalidInput`] and [`SolidError::NotFound`] (malformed B-rep).
pub fn classify_param_on_faces(
    k: &Kernel,
    b: &mut Brep,
    uv: Point2,
    src_zone_tol_3d: f64,
    face1: FaceId,
    domain1: &Extent2d,
    face2: FaceId,
    domain2: &Extent2d,
) -> Result<i32> {
    classify_param_on_faces_with_context(
        k,
        b,
        uv,
        src_zone_tol_3d,
        face1,
        domain1,
        face2,
        domain2,
        &crate::topology::brep_cache::arena_identity,
        &|_| None,
    )
}

/// Classify a shared UV parameter with the caller's cache identity and owner chain.
///
/// # Errors
/// Returns the original face query error when either classification fails.
pub fn classify_param_on_faces_with_context(
    k: &Kernel,
    b: &mut Brep,
    uv: Point2,
    src_zone_tol_3d: f64,
    face1: FaceId,
    domain1: &Extent2d,
    face2: FaceId,
    domain2: &Extent2d,
    keys: &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
    owner_of: &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
) -> Result<i32> {
    // quick BBox containment check
    let in_face1 = domain1.contains_point(uv, EFF_ZERO);
    let in_face2 = domain2.contains_point(uv, EFF_ZERO);

    // low work case - TgtUVParam is in only one TgtFace Domain BBox
    if in_face1 && !in_face2 {
        return Ok(1);
    } else if !in_face1 && in_face2 {
        return Ok(2);
    }

    // ZoneTol3d. When input ZoneTol3d is zero - set ZoneTol3d = Max(pFace1->ZoneTol3d,
    // pFace2->ZoneTol3d)
    let src_tol = if src_zone_tol_3d > EFF_ZERO {
        src_zone_tol_3d
    } else {
        max2(fp2::tolerance(b, face1)?, fp2::tolerance(b, face2)?)
    };

    // classify sFaceUVParam against pFace1
    let mut classify = |face: FaceId| -> Result<PointClassType> {
        let mut pc = PointClassification::new(src_tol);
        let opts = PointClassifyOptions {
            src_zone_tol_3d: src_tol,
            do_boundary_intersections: true,
            use_uv_classification: true,
            uv_space_okay: true,
            ..PointClassifyOptions::default()
        };
        match point_classify_with_context(k, b, face, uv, &mut pc, opts, keys, owner_of) {
            Err(
                e @ (SolidError::Unported { .. }
                | SolidError::InvalidInput { .. }
                | SolidError::NotFound { .. }),
            ) => Err(e),
            _ => Ok(pc.point_class()),
        }
    };
    let pc1 = classify(face1)?;
    if pc1 == PointClassType::Face {
        return Ok(1);
    }

    // classify sFaceUVParam against pFace2
    let pc2 = classify(face2)?;
    if pc2 == PointClassType::Face {
        return Ok(2);
    }

    // when UVParam is on both FaceBndrys its probably on the XSectCrve between the
    // faces - no action required
    if pc1 == PointClassType::Edge && pc2 == PointClassType::Edge {
        return Ok(0);
    }

    // arrive here when the UVPoint did not classify cleanly to one face or the other
    Ok(0)
}

// ---- SmProtoVertex --------------------------------------------------------------------------

/// `SmProtoVertex::DistanceTo`: the distance from the proto vertex's position to
/// `position`.
#[must_use]
pub fn proto_vertex_distance_to(pv: &ProtoVertex, position: Point3) -> f64 {
    pv.position().distance(position)
}

/// `SmProtoVertex::GetCaptureDistance`: the capture distance of the owning manager.
#[must_use]
pub fn proto_vertex_capture_distance<O: ClassObject>(mgr: &ProtoTopologyManager<O>) -> f64 {
    mgr.capture_distance()
}

/// `SmProtoVertex::ContainsPoint`: whether `position` is within the larger of the
/// proto vertex's own tolerance (`m_dTol`) and the manager's capture distance of its
/// position.
#[must_use]
pub fn proto_vertex_contains_point<O: ClassObject>(
    mgr: &ProtoTopologyManager<O>,
    pv: &ProtoVertex,
    position: Point3,
) -> bool {
    // m_dTol >= Local object Tols. Takes precedence over CaptureDistance [B673]
    pv.position().closer_than(
        max2(pv.raw_tolerance(), proto_vertex_capture_distance(mgr)),
        position,
    )
}

// ---- SmFace::CurveOnClassify ------------------------------------------------------------

/// What [`curve_on_classify`] needs of the caller's classification handles.
pub trait CurveOnClassifyObjects<O: ClassObject>:
    crate::intersect::curve_class::p7::FaceIntervalQueries<O>
{
    /// The handle of face `f` (`this`).
    fn face(&self, f: FaceId) -> O;
    /// The classification handle of vertex `v`.
    fn vertex(&self, v: VertexId) -> O;
    /// The classification handle of edge `e`.
    fn edge(&self, e: EdgeId) -> O;
    /// The identity of vertex use `vu`, also used by trim-cache lists.
    fn vertexuse(&self, vu: VertexuseId) -> O;
    /// The identity of edge use `eu`, also used by trim-cache lists.
    fn edgeuse(&self, eu: EdgeuseId) -> O;
    /// The handle of edge `e`'s curve, resolvable through [`crate::intersect::point_class::ClassGeometry::curve`].
    fn edge_curve(&self, e: EdgeId) -> O;
    /// The handle of the trim curve stored on `eu`, which may be the selected
    /// edgeuse's mate. Resolvable through [`crate::intersect::point_class::ClassGeometry::curve`].
    fn trim_curve(&self, eu: EdgeuseId) -> O;
    /// The handle of face `f`'s surface, `None` without one.
    fn face_surface(&self, f: FaceId) -> Option<O>;
    /// The loop behind a target loop handle.
    fn loop_id(&self, tgt_loop: O) -> Option<LoopId>;
    /// Existing face trim cache, if supplied by the caller's cache storage.
    fn trim_cache(&self, _f: FaceId) -> Option<&TrimSrfCache<O>> {
        None
    }
}

/// Resolve stored edge data and its optional gap on this face. The snapshot
/// retains lazy gap-cache updates in the caller's original mutable B-rep.
fn boundary_edge<O: ClassObject>(
    k: &Kernel,
    objects: &mut impl CurveOnClassifyObjects<O>,
    face: FaceId,
    edge: EdgeId,
    has_big_gaps: bool,
) -> Result<BoundaryEdge<O>> {
    let face_gap = if has_big_gaps {
        match crate::topology::edge::edgeuse_of_face(
            objects.classification_brep()?,
            edge,
            Some(face),
        )? {
            Some(eu) => eu1::max_edge_face_gap(k, objects.classification_brep_mut()?, eu, false)?
                .map(|g| g.gap_3d),
            None => None,
        }
    } else {
        None
    };
    let record = objects.classification_brep()?.edge(edge)?;
    if record.curve.is_none() {
        return Err(SolidError::null_pointer(OP_CURVE));
    }
    Ok(BoundaryEdge {
        edge: objects.edge(edge),
        curve: objects.edge_curve(edge),
        interval: record.interval,
        zone_tol_3d: record.tolerance,
        face_gap,
    })
}

/// Build real ordered topology views; only handle identity is caller-specific.
/// UV computation uses the existing vertexuse port on the original mutable arena.
fn boundary_loops<O: ClassObject>(
    k: &Kernel,
    objects: &mut impl CurveOnClassifyObjects<O>,
    face: FaceId,
) -> Result<Vec<BoundaryLoop<O>>> {
    let fu = fp2::upward_faceuse(objects.classification_brep()?, face)?
        .ok_or(SolidError::null_pointer(OP_CURVE))?;
    let cached = objects.trim_cache(face).is_some_and(|cache| {
        let base = cache.base();
        !base.face_was_modified
            && base.have_t_surface_cache
            && base.face_containment_done != FaceContainment::NotDone
    });
    let mut loops = Vec::new();
    let loopuses = objects.classification_brep()?.faceuse(fu)?.loopuses.clone();
    for lu in loopuses {
        let mut vertices = Vec::new();
        for vu in loopuse::vertexuses(objects.classification_brep()?, lu, None)? {
            let vu = vu.ok_or(SolidError::null_pointer(OP_CURVE))?;
            let vertex = objects
                .classification_brep()?
                .vertexuse(vu)?
                .vertex
                .ok_or(SolidError::null_pointer(OP_CURVE))?;
            let record = objects.classification_brep()?.vertex(vertex)?;
            let (point, zone_tol_3d) = (record.point, record.tolerance);
            vertices.push(BoundaryVertex {
                vertexuse: objects.vertexuse(vu),
                vertex: objects.vertex(vertex),
                // Cached traversal supplies its own UV only for visited vertices.
                uv: if cached {
                    Point2::UNDEF
                } else {
                    ignored(crate::topology::vertexuse::compute_uv_point(
                        k,
                        objects.classification_brep_mut()?,
                        vu,
                        true,
                    ))?
                    .unwrap_or(Point2::UNDEF)
                },
                point,
                zone_tol_3d,
            });
        }
        let mut edges = Vec::new();
        for eu in loopuse::edgeuses(objects.classification_brep()?, lu, None, None)? {
            let edge = objects
                .classification_brep()?
                .edgeuse(eu)?
                .edge
                .ok_or(SolidError::null_pointer(OP_CURVE))?;
            let (trim, owner) = eu1::uv_trim_curve(objects.classification_brep()?, eu)?;
            let trim_curve = if trim.is_some() {
                owner.map(|owner| objects.trim_curve(owner))
            } else {
                None
            };
            edges.push(BoundaryEdgeuse {
                edgeuse: objects.edgeuse(eu),
                edge: boundary_edge(k, objects, face, edge, false)?,
                trim_curve,
            });
        }
        loops.push(BoundaryLoop { vertices, edges });
    }
    Ok(loops)
}

/// `SmFace::CurveOnClassify`: splits the curve of `cc` into intervals at its
/// intersections with the face boundary (the trim curves of a 2D curve, else the face's
/// vertices and edges, redone once when the intersection changed the curve) and
/// classifies each interval in, out or on the face (by the interval end
/// classifications, or by point classification with `do_point_classify`). Returns
/// `pbOptModifiedCurve`: whether the curve was modified and the classification rebuilt.
///
/// # Errors
///
/// [`SolidError::Unported`]: the boundary intersections and interval classification (see
/// the module Deviations); `SM_ERR` when the trim boundary intersection fails;
/// [`SolidError::InvalidInput`] for a classification without (resolvable) curve or an
/// unresolvable target loop.
#[allow(clippy::fn_params_excessive_bools)] // upstream's SmBoolean arguments
pub fn curve_on_classify<O: ClassObject>(
    k: &Kernel,
    f: FaceId,
    objects: &mut impl CurveOnClassifyObjects<O>,
    do_point_classify: bool,
    cc: &mut CurveClassification<O>,
    has_big_gaps: bool,
    do_check_3d: bool,
    uv_space_okay: bool,
) -> Result<bool> {
    // init output
    reset_classification(k, cc)?;
    cc.set_face(Some(objects.face(f)));
    let mut modified = false;

    // for two iterations to allow the surf/surf/surf intersection heuristic a chance to
    // update curve and edge shapes.
    for ii in 0..2 {
        // skip 2nd classification when its not needed else set up for DoOver
        if ii > 0 {
            if cc.is_curve_current() {
                continue;
            }
            modified = true;
            reset_classification(k, cc)?;
        }

        let curve = classification_curve(objects, cc)?;
        if curve.dim() == 2 {
            // 2D Curve branch - work in surface's uv space directly on the Face's UVcurves
            cc.set_surface1(objects.face_surface(f));
            match intersect_ts_bound(k, objects, cc, f, do_check_3d) {
                Ok(()) => {}
                Err(e) if e.is_unported() => return Err(e),
                Err(_) => return Err(SolidError::sm_err(OP_CURVE)),
            }
        } else {
            // 3D Curve Branch
            // 1st intersect face vertices with the curve
            let tgt_loops = cc
                .tgt_loops()
                .iter()
                .map(|&o| {
                    objects
                        .loop_id(o)
                        .ok_or_else(|| null(OP_CURVE, "unresolvable target loop"))
                })
                .collect::<Result<Vec<_>>>()?;
            let vertices = face_vertices(
                objects.classification_brep()?,
                f,
                cc.tgt_loops_only(),
                Some(&tgt_loops),
            )?;
            let surface = objects
                .classification_brep()?
                .face(f)?
                .surface
                .as_ref()
                .ok_or_else(|| null(OP_CURVE, "face without surface"))?;
            let vertices = vertices
                .iter()
                .map(|&v| objects.vertex(v))
                .collect::<Vec<_>>();
            cc.intersect_vertices(k, objects, &vertices, surface, has_big_gaps, true)?;

            // 2nd intersect face edges with the curve
            let edges = face_edges(
                objects.classification_brep()?,
                f,
                cc.tgt_loops_only(),
                Some(&tgt_loops),
            )?;
            intersect_edges(k, objects, cc, &edges, f, has_big_gaps)?;
        }
    }

    // 2. classify each interval relative to the face
    let src_zone_tol_3d = src_zone_tol_3d_of(k, Some(cc));
    let face_handle = objects.face(f);
    for ii in 0..cc.size() {
        let degenerate = classification_curve(objects, cc)?.is_degenerate(EFF_ZERO, None)?;
        let ivl = cc.at_mut(ii)?;

        // set Face_BackPtr and SrcZoneTol3d values in interval's PointClassifications
        for p in [&mut ivl.start, &mut ivl.mid, &mut ivl.end] {
            p.set_face(Some(face_handle));
            p.set_src_zone_tol_3d(src_zone_tol_3d);
        }

        // when curve is degenerate
        if degenerate {
            ivl.mid = ivl.start;
            ivl.end = ivl.start;
        }

        // low work - go to next interval
        if !do_point_classify
            && ivl.start.point_class() == PointClassType::Unknown
            && ivl.end.point_class() == PointClassType::Unknown
        {
            continue;
        }

        // classify the interval - avoiding point classification if possible
        if ivl.mid.point_class() == PointClassType::Unknown {
            cc.classify_face_interval(
                k,
                objects,
                ii,
                do_point_classify,
                face_handle,
                None,
                uv_space_okay,
                None,
            )?;
        }
    }
    Ok(modified)
}

/// The curve of `cc` (`GetCurve()`, which upstream dereferences).
fn classification_curve<'o, O: ClassObject>(
    objects: &'o impl CurveOnClassifyObjects<O>,
    cc: &CurveClassification<O>,
) -> Result<&'o Curve> {
    cc.curve()
        .and_then(|c| objects.curve(c))
        .ok_or_else(|| null(OP_CURVE, "classification without curve"))
}

/// `SmCurveClassification::ReSet()` with its defaults.
fn reset_classification<O: ClassObject>(k: &Kernel, cc: &mut CurveClassification<O>) -> Result<()> {
    cc.reset(k, &ResetOptions::default())
}

/// `SmCurveClassification::IntersectTSBound(this, bDoCheck3d)`.
fn intersect_ts_bound<O: ClassObject>(
    k: &Kernel,
    objects: &mut impl CurveOnClassifyObjects<O>,
    cc: &mut CurveClassification<O>,
    face: FaceId,
    do_check_3d: bool,
) -> Result<()> {
    let surface = objects
        .classification_brep()?
        .face(face)?
        .surface
        .as_ref()
        .ok_or(SolidError::null_pointer(OP_CURVE))?
        .clone_for_evaluation()?;
    let loops = boundary_loops(k, objects, face)?;
    let fresh;
    let cache = if let Some(cache) = objects.trim_cache(face) {
        cache
    } else {
        fresh = TrimSrfCache::new(surface.clone(), SurfaceCacheOptions::default());
        &fresh
    };
    cc.intersect_ts_bound(
        k,
        objects,
        &BoundaryFace {
            surface: &surface,
            loops: &loops,
            cache,
        },
        do_check_3d,
    )
}

/// `SmFace::GetEdges(SmTArray&, NULL, bTgtLoopsOnly, pOptTgtLoops)`.
fn face_edges(
    _b: &Brep,
    _f: FaceId,
    _tgt_loops_only: bool,
    _tgt_loops: Option<&[LoopId]>,
) -> Result<Vec<EdgeId>> {
    // UNPORTED(topology.face_p4): SmFace::GetEdges
    Err(SolidError::unported("SmFace::GetEdges", "topology.face_p4"))
}

/// `SmCurveClassification::IntersectEdges(edges, *this, bHasBigGaps)`.
fn intersect_edges<O: ClassObject>(
    k: &Kernel,
    objects: &mut impl CurveOnClassifyObjects<O>,
    cc: &mut CurveClassification<O>,
    edges: &[EdgeId],
    face: FaceId,
    has_big_gaps: bool,
) -> Result<()> {
    let surface = objects
        .classification_brep()?
        .face(face)?
        .surface
        .as_ref()
        .ok_or(SolidError::null_pointer(OP_CURVE))?
        .clone_for_evaluation()?;
    let mut views = Vec::new();
    for &edge in edges {
        views.push(boundary_edge(k, objects, face, edge, has_big_gaps)?);
    }
    cc.intersect_edges(
        k,
        objects,
        &views,
        &surface,
        IntersectEdgesOptions {
            has_big_gaps,
            ..IntersectEdgesOptions::default()
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::curve::Curve;
    use crate::curve::line::{BoundedLineOptions, Line};
    use crate::geom::extent1d::Extent1d;
    use crate::geom::vector2d::Vec2;
    use crate::intersect::point_class::ClassGeometry;
    use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
    use crate::surface::plane::Plane;
    use crate::topology::{
        Edge, Edgeuse, EdgeuseId, EdgeuseOwner, Face, Faceuse, FaceuseId, Loop, Loopuse, LoopuseId,
        LoopuseStart, Vertex, Vertexuse, VertexuseId, VertexuseOwner,
    };

    /// A unit-square planar face over `[0, 1]^2` with UV trim lines.
    fn square(b: &mut Brep) -> FaceId {
        let s = b.add_shell().unwrap();
        let o = Point3::new(0.0, 0.0, 0.0);
        let dom = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let plane = Plane::new(o, Vec3::X, Vec3::Y, Vec2::new(1.0, 1.0), &dom);
        let f = b
            .add_face_record(Face {
                surface: Some(Surface::from(plane)),
                faceuse: Some(FaceuseId::from_index(0)),
                uv_domain: dom,
                tolerance: 1e-6,
                ..Face::new(s)
            })
            .unwrap();
        let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        for &(x, y) in &corners {
            b.add_vertex_record(Vertex {
                point: Point3::new(x, y, 0.0),
                tolerance: 1e-6,
                ..Vertex::default()
            })
            .unwrap();
        }
        b.add_faceuse_record(Faceuse {
            shell: Some(s),
            face: Some(f),
            mate: None,
            orientation: OrientType::Same,
            loopuses: vec![LoopuseId::from_index(0)],
            sibling_links: None,
        })
        .unwrap();
        b.add_loopuse_record(Loopuse {
            faceuse: Some(FaceuseId::from_index(0)),
            loop_: Some(LoopId::from_index(0)),
            mate: None,
            start: LoopuseStart::Edgeuse(Some(EdgeuseId::from_index(0))),
            ..Loopuse::default()
        })
        .unwrap();
        b.add_loop_record(Loop {
            face: f,
            edges: Vec::new(),
            loopuse: Some(LoopuseId::from_index(0)),
        })
        .unwrap();
        for i in 0..4u32 {
            let (x0, y0) = corners[i as usize];
            let (x1, y1) = corners[((i + 1) % 4) as usize];
            let a = Point3::new(x0, y0, 0.0);
            let c = Point3::new(x1, y1, 0.0);
            let unit = Extent1d::from_min_max(0.0, 1.0);
            let line = |a: Point3, c: Point3| {
                Curve::from(Line::bounded(
                    a,
                    c - a,
                    &unit,
                    &BoundedLineOptions::default(),
                ))
            };
            b.add_edge_record(Edge {
                vertices: [VertexId::from_index(i), VertexId::from_index((i + 1) % 4)],
                edgeuses: vec![EdgeuseId::from_index(i)],
                interval: Extent1d::from_min_max(0.0, 1.0),
                tolerance: 1e-6,
                curve: Some(line(a, c)),
            })
            .unwrap();
            let uv_line = Curve::from(Line::bounded(
                a,
                c - a,
                &unit,
                &BoundedLineOptions {
                    dim: 2,
                    ..BoundedLineOptions::default()
                },
            ));
            b.add_edgeuse_record(Edgeuse {
                edge: Some(EdgeId::from_index(i)),
                owner: EdgeuseOwner::Loopuse(Some(LoopuseId::from_index(0))),
                ccw: Some(EdgeuseId::from_index((i + 1) % 4)),
                cw: Some(EdgeuseId::from_index((i + 3) % 4)),
                vertexuse: Some(VertexuseId::from_index(i)),
                orientation: OrientType::Same,
                uv_trim_curve: Some(uv_line),
                ..Edgeuse::default()
            })
            .unwrap();
            b.add_vertexuse_record(Vertexuse {
                vertex: Some(VertexId::from_index(i)),
                owner: VertexuseOwner::Edgeuse(Some(EdgeuseId::from_index(i))),
                ..Vertexuse::default()
            })
            .unwrap();
        }
        f
    }

    fn classify(
        b: &mut Brep,
        f: FaceId,
        uv: Point2,
    ) -> (Result<()>, PointClassification<EntityRef>) {
        let k = Kernel::default();
        let mut pc = PointClassification::new(1e-6);
        let r = point_classify(&k, b, f, uv, &mut pc, PointClassifyOptions::default());
        (r, pc)
    }

    #[test]
    fn outside_the_domain_is_unknown_without_work() {
        let mut b = Brep::new();
        let f = square(&mut b);
        let (r, pc) = classify(&mut b, f, Point2::new(2.0, 0.5));
        assert!(r.is_ok());
        assert_eq!(pc.point_class(), PointClassType::Unknown);
        assert_eq!(pc.object(), None);
    }

    #[test]
    fn vertices_are_found_first() {
        let mut b = Brep::new();
        let f = square(&mut b);
        assert_eq!(
            face_vertices(&b, f, false, None).unwrap(),
            (0..4).map(VertexId::from_index).collect::<Vec<_>>()
        );
        let (r, pc) = classify(&mut b, f, Point2::new(1.0, 1.0));
        assert!(r.is_ok());
        assert_eq!(pc.point_class(), PointClassType::Vertex);
        assert_eq!(
            pc.object(),
            Some(EntityRef::Vertex(VertexId::from_index(2)))
        );
        assert_eq!(pc.pre_snap_param(), BIG_DOUBLE);
    }

    #[test]
    fn interior_point_is_on_the_face() {
        let mut b = Brep::new();
        let f = square(&mut b);
        let (r, pc) = classify(&mut b, f, Point2::new(0.25, 0.5));
        r.unwrap();
        assert_eq!(pc.point_class(), PointClassType::Face);
        assert_eq!(pc.object(), Some(EntityRef::Face(f)));
        assert_eq!(pc.uv_param(), Point2::new(0.25, 0.5));
        assert_eq!(pc.src_zone_tol_3d(), 1e-6);
    }

    #[test]
    fn classify_param_on_faces_by_domain() {
        let mut b = Brep::new();
        let f = square(&mut b);
        let k = Kernel::default();
        let d1 = Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0);
        let d2 = Extent2d::from_bounds(2.0, 2.0, 3.0, 3.0);
        let uv = Point2::new(0.5, 0.5);
        assert_eq!(
            classify_param_on_faces(&k, &mut b, uv, 0.0, f, &d1, f, &d2).unwrap(),
            1
        );
        assert_eq!(
            classify_param_on_faces(&k, &mut b, uv, 0.0, f, &d2, f, &d1).unwrap(),
            2
        );
        assert_eq!(
            classify_param_on_faces(&k, &mut b, uv, 0.0, f, &d1, f, &d1).unwrap(),
            1
        );
    }

    #[test]
    fn missing_surface_is_a_null_pointer() {
        let mut b = Brep::new();
        let f = square(&mut b);
        b.face_mut(f).unwrap().surface = None;
        let (r, _) = classify(&mut b, f, Point2::new(0.5, 0.5));
        assert!(matches!(
            r,
            Err(SolidError::Upstream {
                code: crate::error::SM_ERR_NULL_POINTER,
                ..
            })
        ));
    }

    #[test]
    fn proto_vertex_containment_uses_the_larger_tolerance() {
        let e = Brep::new();
        let mut mgr: ProtoTopologyManager<EntityRef> = ProtoTopologyManager::new(
            &Kernel::default(),
            (EntityRef::Face(FaceId::from_index(0)), &e),
            (EntityRef::Face(FaceId::from_index(1)), &e),
            EntityRef::Face(FaceId::from_index(2)),
            1e-3,
            0.5,
            0.0,
        )
        .unwrap();
        let id = mgr
            .new_proto_vertex(Point3::new(0.0, 0.0, 0.0), 0.1)
            .unwrap();
        let pv = mgr.proto_vertex(id).unwrap().clone();
        assert!(proto_vertex_contains_point(
            &mgr,
            &pv,
            Point3::new(0.4, 0.0, 0.0)
        ));
        assert!(!proto_vertex_contains_point(
            &mgr,
            &pv,
            Point3::new(0.6, 0.0, 0.0)
        ));
        assert_eq!(
            proto_vertex_distance_to(&pv, Point3::new(3.0, 4.0, 0.0)),
            5.0
        );
        mgr.set_capture_distance(0.0);
        assert_eq!(proto_vertex_capture_distance(&mgr), 0.0);
        assert!(!proto_vertex_contains_point(
            &mgr,
            &pv,
            Point3::new(0.4, 0.0, 0.0)
        ));
    }
    #[derive(Clone, Copy, PartialEq)]
    enum Handle {
        Topology(EntityRef),
        Source,
        EdgeCurve(EdgeId),
        Trim(EdgeuseId),
        Surface(FaceId),
    }
    impl ClassObject for Handle {
        fn class_kind(&self) -> ClassObjectKind {
            match self {
                Self::Topology(entity) => entity.class_kind(),
                Self::Source | Self::EdgeCurve(_) | Self::Trim(_) => ClassObjectKind::Curve,
                Self::Surface(_) => ClassObjectKind::Surface,
            }
        }
    }
    struct Handles<'a> {
        b: &'a mut Brep,
        source: Curve,
        cache: Option<TrimSrfCache<Handle>>,
    }
    impl ClassGeometry<Handle> for Handles<'_> {
        fn curve(&self, h: Handle) -> Option<&Curve> {
            match h {
                Handle::Source => Some(&self.source),
                Handle::EdgeCurve(e) => self.b.edge(e).ok()?.curve.as_ref(),
                Handle::Trim(eu) => eu1::uv_trim_curve(self.b, eu).ok()?.0,
                _ => None,
            }
        }
        fn surface(&self, h: Handle) -> Option<&Surface> {
            if let Handle::Surface(f) = h {
                self.b.face(f).ok()?.surface.as_ref()
            } else {
                None
            }
        }
        fn volume(&self, _: Handle) -> Option<&crate::surface::volume::Volume> {
            None
        }
        fn point(&self, _: Handle) -> Option<Point3> {
            None
        }
        fn topology_zone_tol_3d(&self, h: Handle) -> Option<f64> {
            match h {
                Handle::Topology(EntityRef::Vertex(v)) => {
                    self.b.vertex(v).ok().map(|v| v.tolerance)
                }
                Handle::Topology(EntityRef::Edge(e)) => self.b.edge(e).ok().map(|e| e.tolerance),
                _ => None,
            }
        }
    }
    impl crate::intersect::curve_class::p7::FaceQueryOwner<Handle> for Handles<'_> {
        fn classification_brep(&self) -> Result<&Brep> {
            Ok(self.b)
        }
        fn classification_brep_mut(&mut self) -> Result<&mut Brep> {
            Ok(self.b)
        }
        fn classification_entity(&self, _: Handle, object: Handle) -> Result<EntityRef> {
            match object {
                Handle::Topology(entity) => Ok(entity),
                _ => Err(null(OP_CURVE, "not a fixture topology handle")),
            }
        }
        fn classification_handle(&self, _: Handle, entity: EntityRef) -> Result<Handle> {
            Ok(Handle::Topology(entity))
        }
        fn query_on_face<T>(
            &mut self,
            kernel: &Kernel,
            face: Handle,
            query: impl FnOnce(
                &Kernel,
                &mut Brep,
                FaceId,
                &dyn Fn(crate::topology::brep_cache::CacheObject) -> crate::geom::object::ObjectRef,
                &dyn Fn(crate::geom::object::ObjectRef) -> Option<crate::geom::object::ObjectRef>,
                Option<&mut dyn FacePointCache>,
            ) -> Result<T>,
        ) -> Result<T> {
            let Handle::Topology(EntityRef::Face(face)) = face else {
                return Err(null(OP_CURVE, "not a fixture face"));
            };
            let decode = |handle| match handle {
                Handle::Topology(entity) => Ok(entity),
                _ => Err(null(OP_CURVE, "not a fixture cache topology handle")),
            };
            let encode = |entity| Ok(Handle::Topology(entity));
            let owner_of = |_| None;
            if let Some(cache) = self.cache.as_mut() {
                let mut view = MappedFacePointCache::new(cache, &decode, &encode);
                query(
                    kernel,
                    self.b,
                    face,
                    &crate::topology::brep_cache::arena_identity,
                    &owner_of,
                    Some(&mut view),
                )
            } else {
                query(
                    kernel,
                    self.b,
                    face,
                    &crate::topology::brep_cache::arena_identity,
                    &owner_of,
                    None,
                )
            }
        }
    }

    impl CurveOnClassifyObjects<Handle> for Handles<'_> {
        fn face(&self, f: FaceId) -> Handle {
            Handle::Topology(EntityRef::Face(f))
        }
        fn vertex(&self, v: VertexId) -> Handle {
            Handle::Topology(EntityRef::Vertex(v))
        }
        fn edge(&self, e: EdgeId) -> Handle {
            Handle::Topology(EntityRef::Edge(e))
        }
        fn vertexuse(&self, vu: VertexuseId) -> Handle {
            Handle::Topology(EntityRef::Vertexuse(vu))
        }
        fn edgeuse(&self, eu: EdgeuseId) -> Handle {
            Handle::Topology(EntityRef::Edgeuse(eu))
        }
        fn edge_curve(&self, e: EdgeId) -> Handle {
            Handle::EdgeCurve(e)
        }
        fn trim_curve(&self, eu: EdgeuseId) -> Handle {
            Handle::Trim(eu)
        }
        fn face_surface(&self, f: FaceId) -> Option<Handle> {
            Some(Handle::Surface(f))
        }
        fn loop_id(&self, _: Handle) -> Option<LoopId> {
            None
        }
        fn trim_cache(&self, _: FaceId) -> Option<&TrimSrfCache<Handle>> {
            self.cache.as_ref()
        }
    }
    #[test]
    fn nonempty_boundary_adapters_reach_both_part9_methods() {
        let k = Kernel::default();
        let mut b = Brep::new();
        let f = square(&mut b);
        // Edge-loop UV computation retains its genuine NormalizedEvaluate
        // dependency. An isolated vertex loop can use the existing projection.
        b.vertexuse_mut(VertexuseId::from_index(0)).unwrap().owner =
            VertexuseOwner::Loopuse(Some(LoopuseId::from_index(0)));
        b.loopuse_mut(LoopuseId::from_index(0)).unwrap().start =
            LoopuseStart::Vertexuse(Some(VertexuseId::from_index(0)));
        let mut objects = Handles {
            b: &mut b,
            cache: None,
            source: Curve::from(Line::from_points(
                Point3::new(0.2, 1.5, 0.),
                Point3::new(0.8, 1.5, 0.),
                2,
            )),
        };
        let loops = boundary_loops(&k, &mut objects, f).unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].vertices.len(), 1);
        assert!(loops[0].edges.is_empty());
        assert!(
            loops[0].vertices[0].vertex
                == Handle::Topology(EntityRef::Vertex(VertexId::from_index(0)))
        );
        let mut cc = CurveClassification::with_curve_interval(
            &k,
            Handle::Source,
            Extent1d::from_min_max(0., 1.),
            1e-6,
        );
        intersect_ts_bound(&k, &mut objects, &mut cc, f, false).unwrap();
        assert_eq!(cc.size(), 1);
        let mut objects = Handles {
            b: &mut b,
            cache: None,
            source: Curve::from(Line::from_points(
                Point3::new(0.2, 0.5, 2.),
                Point3::new(0.8, 0.5, 2.),
                3,
            )),
        };
        intersect_edges(
            &k,
            &mut objects,
            &mut cc,
            &[EdgeId::from_index(0)],
            f,
            false,
        )
        .unwrap();
        assert_eq!(cc.size(), 1);
    }
    #[test]
    fn cached_boundary_views_resolve_real_edge_and_trim_handles() {
        let k = Kernel::default();
        let mut b = Brep::new();
        let f = square(&mut b);
        let mut cache = TrimSrfCache::new(
            b.face(f).unwrap().surface.clone().unwrap(),
            SurfaceCacheOptions::default(),
        );
        cache.base_mut().have_t_surface_cache = true;
        cache.base_mut().face_was_modified = false;
        cache.base_mut().face_containment_done = FaceContainment::Done;
        let mut objects = Handles {
            b: &mut b,
            cache: Some(cache),
            source: Curve::from(Line::from_points(Vec3::ZERO, Vec3::X, 2)),
        };
        let loops = boundary_loops(&k, &mut objects, f).unwrap();
        assert_eq!(loops[0].vertices.len(), 4);
        assert_eq!(loops[0].edges.len(), 4);
        for eu in &loops[0].edges {
            assert_eq!(objects.curve(eu.edge.curve).unwrap().dim(), 3);
            assert_eq!(objects.curve(eu.trim_curve.unwrap()).unwrap().dim(), 2);
        }
    }
}

#[cfg(test)]
mod retained_point_cache_tests {
    use super::*;
    use crate::geom::attribute::ObjectKey;
    use crate::geom::object::ObjectRef;
    use crate::geom::tree::NodeId;
    use crate::surface::plane::Plane;
    use crate::surface::surface_cache::p1::SurfaceCacheOptions;
    use crate::topology::brep_cache::{arena_identity, CacheObject};

    #[derive(Debug, Clone, Copy, PartialEq)]
    enum CacheHandle {
        Surface,
        Topology(EntityRef),
    }

    fn caller_key(object: CacheObject) -> ObjectRef {
        let mut key = arena_identity(object);
        key.key = ObjectKey(key.key.0 + 123_456);
        key
    }

    fn to_entity(handle: CacheHandle) -> Result<EntityRef> {
        match handle {
            CacheHandle::Topology(entity) => Ok(entity),
            CacheHandle::Surface => Err(SolidError::invalid_input(
                OP_CLASSIFY,
                "a surface handle is not a topology entity",
            )),
        }
    }

    fn ready_without_tree() -> (Brep, FaceId, TrimSrfCache<CacheHandle>) {
        let mut b = Brep::new();
        let shell = b.add_shell().unwrap();
        let face = b.add_face(shell).unwrap();
        let domain = Extent2d::from_bounds(-1.0, -1.0, 1.0, 1.0);
        let mut surface = Surface::from(Plane::new(
            Point3::ZERO,
            Vec3::X,
            Vec3::Y,
            crate::math::Vec2::new(1.0, 1.0),
            &domain,
        ));
        surface.set_owner(Some(caller_key(CacheObject::Face(face))));
        let cache_surface = surface.clone_for_evaluation().unwrap();
        b.face_mut(face).unwrap().surface = Some(surface);
        b.face_mut(face).unwrap().uv_domain = domain;
        let mut cache = TrimSrfCache::new(cache_surface, SurfaceCacheOptions::default());
        cache
            .base_mut()
            .set_surface_handle(Some(CacheHandle::Surface));
        cache.base_mut().have_t_surface_cache = true;
        cache.base_mut().face_containment_done = FaceContainment::Done;
        (b, face, cache)
    }

    #[test]
    fn ready_missing_tree_builds_in_place_and_restores_native_error_scopes() {
        let (mut b, face, mut cache) = ready_without_tree();
        let kernel = Kernel::default();
        let mut pc = PointClassification::new(0.125);
        pc.set_gap_3d(17.0);
        let encode = |entity| Ok(CacheHandle::Topology(entity));
        let mut mapped = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
        let result = point_classify_with_context(
            &kernel,
            &mut b,
            face,
            Point2::new(0.0, 0.0),
            &mut pc,
            PointClassifyOptions {
                trim_cache: Some(&mut mapped),
                ..PointClassifyOptions::default()
            },
            &caller_key,
            &|_| None,
        );
        assert!(matches!(
            result,
            Err(SolidError::Unported {
                unit: "topology.face_p4",
                ..
            })
        ));
        assert!(cache.base().tree().is_some());
        assert_eq!(cache.base().face_containment_done, FaceContainment::Done);
        assert!(cache.base().point_test_enabled);
        assert!(!b.editing_enabled());
        assert_eq!(
            b.face(face).unwrap().surface.as_ref().unwrap().owner(),
            Some(caller_key(CacheObject::Face(face)))
        );
        assert_eq!(pc.point_class(), PointClassType::Unknown);
        let mut expected = PointClassification::new(0.125);
        expected.set_gap_3d(17.0);
        assert_eq!(pc, expected);
        assert_eq!(cache.base().surface_handle(), Some(CacheHandle::Surface));
    }

    #[test]
    fn external_surface_owner_chain_enters_the_real_ready_cache_gate() {
        let (mut b, face, _) = ready_without_tree();
        let intermediate = ObjectRef {
            type_code: b.face(face).unwrap().surface.as_ref().unwrap().type_code(),
            key: ObjectKey(987_654),
        };
        let surface = b.face_mut(face).unwrap().surface.as_mut().unwrap();
        surface.set_owner(Some(intermediate));
        assert_eq!(surface.face().unwrap(), None);
        let mut cache = TrimSrfCache::new(
            surface.clone_for_evaluation().unwrap(),
            SurfaceCacheOptions::default(),
        );
        cache
            .base_mut()
            .set_surface_handle(Some(CacheHandle::Surface));
        cache.base_mut().have_t_surface_cache = true;
        cache.base_mut().face_containment_done = FaceContainment::Done;
        let owner_of =
            |owner| (owner == intermediate).then_some(caller_key(CacheObject::Face(face)));
        let encode = |entity| Ok(CacheHandle::Topology(entity));
        let mut mapped = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
        let mut pc = PointClassification::new(0.125);
        let result = point_classify_with_context(
            &Kernel::default(),
            &mut b,
            face,
            Point2::new(0.0, 0.0),
            &mut pc,
            PointClassifyOptions {
                trim_cache: Some(&mut mapped),
                ..PointClassifyOptions::default()
            },
            &caller_key,
            &owner_of,
        );
        // A skipped gate would leave the actual cache tree absent. The real
        // build retains its tree before the genuine topology.face_p4 dependency.
        assert!(matches!(
            result,
            Err(SolidError::Unported {
                unit: "topology.face_p4",
                ..
            })
        ));
        assert!(cache.base().tree().is_some());
        assert_eq!(cache.base().face_containment_done, FaceContainment::Done);
        assert!(!b.editing_enabled());
    }

    #[test]
    fn existing_cache_maps_only_leaf_faces_and_preserves_inside_null_face() {
        let (mut b, face, mut cache) = ready_without_tree();
        b.set_editing_enabled(true);
        let kernel = Kernel::default();
        let encode = |entity| Ok(CacheHandle::Topology(entity));
        {
            let mut mapped = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
            assert_eq!(
                mapped
                    .classify_ready(
                        &kernel,
                        &mut b,
                        face,
                        Point2::new(0.0, 0.0),
                        &caller_key,
                        &|_| None
                    )
                    .unwrap(),
                FacePointCacheDecision::NoDecision
            );
        }
        assert!(cache.base().tree().is_some());
        let node_count = cache.base().tree().unwrap().nodes().len();
        for stored_face in [Some(CacheHandle::Topology(EntityRef::Face(face))), None] {
            let tree = cache.base_mut().tree_mut().unwrap();
            for index in 0..tree.nodes().len() {
                let NodeData::Surface(aux) = &mut tree.node_mut(NodeId(index)).unwrap().data else {
                    panic!("surface tree")
                };
                aux.node_class = NodeClass::Inside;
                aux.face = stored_face;
            }
            cache.base_mut().have_t_surface_cache = true;
            cache.base_mut().face_containment_done = FaceContainment::Done;
            cache.base_mut().face_was_modified = false;
            let mut mapped = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
            let mut pc = PointClassification::new(0.125);
            point_classify_with_context(
                &kernel,
                &mut b,
                face,
                Point2::new(0.0, 0.0),
                &mut pc,
                PointClassifyOptions {
                    trim_cache: Some(&mut mapped),
                    ..PointClassifyOptions::default()
                },
                &caller_key,
                &|_| None,
            )
            .unwrap();
            assert_eq!(pc.point_class(), PointClassType::Face);
            assert_eq!(
                pc.object(),
                stored_face.map(|handle| to_entity(handle).unwrap())
            );
            assert_eq!(cache.base().tree().unwrap().nodes().len(), node_count);
            assert_eq!(cache.base().surface_handle(), Some(CacheHandle::Surface));
        }
        assert!(b.editing_enabled());
    }

    #[test]
    fn ready_cache_of_a_different_surface_owner_is_not_consumed() {
        for external_chain in [false, true] {
            let (mut b, face, mut cache) = ready_without_tree();
            let kernel = Kernel::default();
            let encode = |entity| Ok(CacheHandle::Topology(entity));
            b.set_editing_enabled(true);
            {
                let mut view = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
                assert_eq!(
                    view.classify_ready(
                        &kernel,
                        &mut b,
                        face,
                        Point2::new(0.0, 0.0),
                        &caller_key,
                        &|_| None
                    )
                    .unwrap(),
                    FacePointCacheDecision::NoDecision
                );
            }
            // A real ready tree would immediately classify the target as inside
            // if the erroneous any-face-owner gate consumed this cache.
            let tree = cache.base_mut().tree_mut().unwrap();
            for index in 0..tree.nodes().len() {
                let NodeData::Surface(aux) = &mut tree.node_mut(NodeId(index)).unwrap().data else {
                    panic!("surface tree")
                };
                aux.node_class = NodeClass::Inside;
                aux.face = Some(CacheHandle::Topology(EntityRef::Face(face)));
            }
            let node_count = cache.base().tree().unwrap().nodes().len();
            cache.base_mut().have_t_surface_cache = true;
            cache.base_mut().face_containment_done = FaceContainment::Done;
            cache.base_mut().face_was_modified = false;
            let shell = b.add_shell().unwrap();
            let other = b.add_face(shell).unwrap();
            let intermediate = ObjectRef {
                type_code: b.face(face).unwrap().surface.as_ref().unwrap().type_code(),
                key: ObjectKey(777_123),
            };
            let real_owner = caller_key(CacheObject::Face(other));
            b.face_mut(face)
                .unwrap()
                .surface
                .as_mut()
                .unwrap()
                .set_owner(Some(if external_chain {
                    intermediate
                } else {
                    real_owner
                }));
            let owner_of = |owner| (external_chain && owner == intermediate).then_some(real_owner);
            let mut view = MappedFacePointCache::new(&mut cache, &to_entity, &encode);
            let mut pc = PointClassification::new(0.125);
            point_classify_with_context(
                &kernel,
                &mut b,
                face,
                Point2::new(0.0, 0.0),
                &mut pc,
                PointClassifyOptions {
                    trim_cache: Some(&mut view),
                    ..PointClassifyOptions::default()
                },
                &caller_key,
                &owner_of,
            )
            .unwrap();
            assert_eq!(pc.point_class(), PointClassType::Unknown);
            assert_eq!(pc.object(), None);
            assert_eq!(cache.base().tree().unwrap().nodes().len(), node_count);
            assert!(cache.base().tree().unwrap().nodes().iter().all(|node| {
                !matches!(&node.data, NodeData::Surface(aux) if aux.node_class != NodeClass::Inside)
            }));
        }
    }
}
