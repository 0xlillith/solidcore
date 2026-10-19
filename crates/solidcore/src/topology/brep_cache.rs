//! Spatial trees of B-rep vertices, curves and surfaces.
//!
//! Upstream: `source/SMLib/inc/SmBrepCache.h`, `source/SMLib/src/SmBrepCache.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`BrepCache::new`] | `SmBrepCache::SmBrepCache` |
//! | [`BrepCache::build_trees`] | `SmBrepCache::BuildTrees` |
//! | [`BrepCache::resize_trees`] | `SmBrepCache::ReSizeTrees` |
//! | [`BrepCache::vertex_tree`], [`BrepCache::curve_tree`], [`BrepCache::surface_tree`] | `GetVertexTree`, `GetCurveTree`, `GetSurfaceTree` |
//! | [`make_brep_cache`] | `SmSAGObject::CacheMakeOrValidate` (new-cache path) |
//! | [`BrepObjects`] | `SmTopologySolver` object/cache access |
//! | [`TYPE`], [`is_kind_of`] | `GetType`, `IsKindOf` |
//! | `Drop` (implicit) | `SmBrepCache::~SmBrepCache` |
//!
//! Deviations: trees own typed object handles instead of pointers. The source is
//! supplied when building, so the cache does not borrow its geometry permanently;
//! the caller retains source identity and invalidates the cache after edits. Rust
//! drops trees automatically. Exact B-rep lists remain explicit
//! `topology.brep_p1` back-edges; arena allocation order is not list order.
//! Rebuilding an initialized cache and resizing an
//! uninitialized cache return `InvalidInput` (upstream asserts). Ignored statuses
//! propagate Unported; an ordinary ignored box failure keeps the previous box,
//! as upstream. No new recursion; spatial-tree traversal uses geom.tree's bounds.
//! B-rep and face boxes call `query::brep` and `query::face::p1`; face queries
//! (`topology.face_p4`) remain explicit back-edges. The polygon path is ported.

use crate::error::{Result, SolidError};
use crate::geom::a_object::{BREP_TYPE, POLY_BREP_TYPE};
use crate::geom::extent3d::Extent3d;
use crate::geom::object;
use crate::geom::tree::{NodeGeomType, Tree, TreeOwner};
use crate::kernel::Kernel;
use crate::math::EFF_ZERO;
use crate::nurbs::b_spline_curve::BoundingBoxRequest;
use crate::query::poly::{PolyBrep, PolyEdgeId, PolyFaceId, PolyVertexId};
use crate::topology::{edge, Brep, EdgeId, FaceId, VertexId};

const OP: &str = "topology::brep_cache";
/// `SmBrepCache_TYPE`.
pub const TYPE: i64 = 16400;
/// `SmBrepCache::GetTypeString`.
pub const TYPE_STRING: &str = "SmBrepCache_TYPE";
/// `SmBrepCache::GetClassString`.
pub const CLASS_STRING: &str = "SmBrepCache";

/// `SmBrepCache::IsKindOf` delegates directly to `SmObject`, skipping `SmCacheObj`.
pub const fn is_kind_of(t: i64) -> bool {
    t == TYPE || object::is_kind_of(t)
}

/// Geometry held in a B-rep cache; curve and surface ids name their owning records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CacheObject {
    /// An owner outside the adapter's exact topology records.
    Other(object::ObjectRef),
    /// The source B-rep of a topology solver.
    Brep,
    /// An edge record (distinct from its curve).
    Edge(EdgeId),
    /// A face record (distinct from its surface).
    Face(FaceId),
    /// An exact vertex.
    Vertex(VertexId),
    /// A curve belonging to an edge.
    Curve(EdgeId),
    /// A surface belonging to a face.
    Surface(FaceId),
    /// A polygon vertex.
    PolyVertex(PolyVertexId),
    /// A polygon edge.
    PolyEdge(PolyEdgeId),
    /// A polygon face.
    PolyFace(PolyFaceId),
}

/// The source object of `BuildTrees`.
#[derive(Clone, Copy)]
pub enum CacheSource<'a> {
    /// An exact B-rep.
    Brep(&'a Brep),
    /// A polygonal B-rep.
    PolyBrep(&'a PolyBrep),
}

impl CacheSource<'_> {
    fn type_code(self) -> i64 {
        match self {
            Self::Brep(_) => BREP_TYPE,
            Self::PolyBrep(_) => POLY_BREP_TYPE,
        }
    }
}

/// `SmBrepCache`: three optional spatial trees. Object identities can be mapped
/// to a solver adapter's handles during construction.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepCache<O = CacheObject> {
    target_type: i64,
    vertices: Option<Tree<O>>,
    curves: Option<Tree<O>>,
    surfaces: Option<Tree<O>>,
}

impl<O> BrepCache<O> {
    /// Constructs an empty cache for an `SmBrep` or `SmPolyBrep` type code.
    /// Other codes return upstream's `SM_ERR`.
    pub fn new(target_type: i64) -> Result<Self> {
        if target_type != BREP_TYPE && target_type != POLY_BREP_TYPE {
            return Err(SolidError::sm_err(OP));
        }
        Ok(Self {
            target_type,
            vertices: None,
            curves: None,
            surfaces: None,
        })
    }

    /// `GetVertexTree`: absent until `BuildTrees` succeeds.
    pub const fn vertex_tree(&self) -> Option<&Tree<O>> {
        self.vertices.as_ref()
    }
    /// `GetCurveTree`: absent until `BuildTrees` succeeds.
    pub const fn curve_tree(&self) -> Option<&Tree<O>> {
        self.curves.as_ref()
    }
    /// `GetSurfaceTree`: absent until `BuildTrees` succeeds.
    pub const fn surface_tree(&self) -> Option<&Tree<O>> {
        self.surfaces.as_ref()
    }
    /// Mutable `GetVertexTree` (the tree a B-rep notification edits).
    pub fn vertex_tree_mut(&mut self) -> Option<&mut Tree<O>> {
        self.vertices.as_mut()
    }
    /// Mutable `GetCurveTree` (the tree a B-rep notification edits).
    pub fn curve_tree_mut(&mut self) -> Option<&mut Tree<O>> {
        self.curves.as_mut()
    }
    /// Mutable `GetSurfaceTree` (the tree a B-rep notification edits).
    pub fn surface_tree_mut(&mut self) -> Option<&mut Tree<O>> {
        self.surfaces.as_mut()
    }
}

impl<O: Copy + PartialEq> BrepCache<O> {
    /// Builds all trees atomically in upstream vertex/edge/face order. `map`
    /// translates arena ids to caller-owned object handles, without changing order.
    /// `keys` identifies exact edge/face owners independently of arena numbering.
    pub fn build_trees(
        &mut self,
        k: &Kernel,
        source: CacheSource<'_>,
        mut map: impl FnMut(CacheObject) -> O,
        keys: &dyn Fn(CacheObject) -> object::ObjectRef,
    ) -> Result<()> {
        if self.vertices.is_some()
            || self.curves.is_some()
            || self.surfaces.is_some()
            || self.target_type != source.type_code()
        {
            return Err(SolidError::invalid_input(
                OP,
                "cache already built or source type mismatch",
            ));
        }
        let mut bbox = match source {
            CacheSource::Brep(b) => {
                let mut bbox = brep_bounding_box(k, b)?;
                bbox.expand_absolute(b.tolerance() * 100.0);
                let size = bbox.size();
                bbox.add_point(bbox.max() + 0.15 * size);
                bbox.add_point(bbox.min() - 0.15 * size);
                bbox
            }
            CacheSource::PolyBrep(b) => {
                let mut bbox = b
                    .calculate_bounding_box(false, None, false, None)?
                    .normal_box;
                bbox.expand_absolute(b.tolerance() * 100.0);
                bbox
            }
        };
        // Preserve the source's inside-out empty box, including its arithmetic.
        let mut vertices = Tree::new_spatial(&bbox);
        let mut curves = Tree::new_spatial(&bbox);
        let mut surfaces = Tree::new_spatial(&bbox);
        for (tree, geom) in [
            (&mut vertices, None),
            (&mut curves, Some(NodeGeomType::LineSeg)),
            (&mut surfaces, Some(NodeGeomType::Polygon)),
        ] {
            tree.set_owner_object(Some(TreeOwner::BrepCache));
            if let Some(geom) = geom {
                let top = tree.top_node().ok_or(SolidError::null_pointer(OP))?;
                tree.node_mut(top)?.geom_type = geom;
            }
        }
        match source {
            CacheSource::PolyBrep(b) => {
                for v in b.poly_vertices() {
                    bbox = Extent3d::from_point(b.vertex_point(v)?);
                    bbox.expand_absolute(b.vertex_tolerance(v)?);
                    vertices.add_to_spatial_tree(
                        &bbox,
                        Some(map(CacheObject::PolyVertex(v))),
                        NodeGeomType::Default,
                    )?;
                }
                for e in b.poly_edges(None)? {
                    bbox = b.edge_bounding_box(e)?;
                    bbox.expand_absolute(b.edge_tolerance(e)?);
                    curves.add_to_spatial_tree(
                        &bbox,
                        Some(map(CacheObject::PolyEdge(e))),
                        NodeGeomType::LineSeg,
                    )?;
                }
                for f in b.poly_faces(None)? {
                    bbox = b.face_bounding_box(f)?;
                    surfaces.add_to_spatial_tree(
                        &bbox,
                        Some(map(CacheObject::PolyFace(f))),
                        NodeGeomType::Polygon,
                    )?;
                }
            }
            CacheSource::Brep(b) => {
                for v in brep_vertices(b)? {
                    let data = b.vertex(v)?;
                    bbox = Extent3d::from_point(data.point);
                    bbox.expand_absolute(data.tolerance);
                    vertices.add_to_spatial_tree(
                        &bbox,
                        Some(map(CacheObject::Vertex(v))),
                        NodeGeomType::Default,
                    )?;
                }
                for e in brep_edges(b)? {
                    let data = b.edge(e)?;
                    let curve = data.curve.as_ref().ok_or(SolidError::null_pointer(OP))?;
                    if curve.owner() == Some(keys(CacheObject::Edge(e))) {
                        bbox = curve
                            .calculate_bounding_box(
                                k,
                                &data.interval,
                                BoundingBoxRequest {
                                    normal: true,
                                    pseudo: false,
                                    polar: false,
                                },
                                true,
                            )?
                            .normal
                            .ok_or(SolidError::null_pointer(OP))?;
                        bbox.expand_absolute(data.tolerance);
                        curves.add_to_spatial_tree(
                            &bbox,
                            Some(map(CacheObject::Curve(e))),
                            NodeGeomType::Default,
                        )?;
                    }
                }
                for f in brep_faces(b)? {
                    let data = b.face(f)?;
                    let surface = data.surface.as_ref().ok_or(SolidError::null_pointer(OP))?;
                    if surface.owner() == Some(keys(CacheObject::Face(f))) {
                        bbox = face_bounding_box(k, b, f)?;
                        bbox.expand_absolute(data.tolerance);
                        surfaces.add_to_spatial_tree(
                            &bbox,
                            Some(map(CacheObject::Surface(f))),
                            NodeGeomType::Default,
                        )?;
                    }
                }
            }
        }
        self.vertices = Some(vertices);
        self.curves = Some(curves);
        self.surfaces = Some(surfaces);
        Ok(())
    }

    /// Expands the three trees to include exact faces, edges and vertices. This
    /// does not insert objects; all entity tolerances contribute to the SAME
    /// maximum used for every tree, as upstream. Empty slices are null arrays.
    /// `keys` supplies the identities used by the geometry owner links.
    /// `owner_of` resolves intermediate generated-curve and edgeuse owners.
    pub fn resize_trees(
        &mut self,
        k: &Kernel,
        b: &Brep,
        faces: &[FaceId],
        edges: &[EdgeId],
        vertices: &[VertexId],
        keys: &dyn Fn(CacheObject) -> object::ObjectRef,
        owner_of: &dyn Fn(object::ObjectRef) -> Option<object::ObjectRef>,
    ) -> Result<()> {
        if self.vertices.is_none() || self.curves.is_none() || self.surfaces.is_none() {
            return Err(SolidError::invalid_input(OP, "cache has not been built"));
        }
        let mut surface_union = Extent3d::new();
        let mut curve_union = Extent3d::new();
        let mut vertex_union = Extent3d::new();
        let mut surface_box = Extent3d::new();
        let mut curve_box = Extent3d::new();
        let mut tol = 0.0;
        for &f in faces {
            let face = b.face(f)?;
            // The original surface's owning face, not necessarily the input face.
            let owner = face
                .surface
                .as_ref()
                .ok_or(SolidError::null_pointer(OP))?
                .face()?
                .ok_or(SolidError::null_pointer(OP))?;
            let owner = face_id(b, owner, keys)?;
            match face_bounding_box(k, b, owner) {
                Ok(bbox) => surface_box = bbox,
                Err(e) if e.is_unported() => return Err(e),
                Err(_) => {}
            }
            surface_union = surface_union.union(&surface_box);
            if face.tolerance > tol {
                tol = face.tolerance;
            }
            for e in face_edges(b, f)? {
                resize_edge(
                    k,
                    b,
                    e,
                    &mut curve_box,
                    &mut curve_union,
                    &mut vertex_union,
                    &mut tol,
                    keys,
                    owner_of,
                )?;
            }
        }
        for &e in edges {
            resize_edge(
                k,
                b,
                e,
                &mut curve_box,
                &mut curve_union,
                &mut vertex_union,
                &mut tol,
                keys,
                owner_of,
            )?;
        }
        for &v in vertices {
            let vertex = b.vertex(v)?;
            if vertex.tolerance > tol {
                tol = vertex.tolerance;
            }
            vertex_union.add_point(vertex.point);
        }
        for (tree, mut bbox, containment_tol) in [
            (&mut self.surfaces, surface_union, 0.0),
            (&mut self.curves, curve_union, 0.0),
            (&mut self.vertices, vertex_union, EFF_ZERO),
        ] {
            if !bbox.has_negative_volume() {
                bbox.expand_absolute(100.0 * tol);
                let tree = tree.as_mut().ok_or(SolidError::null_pointer(OP))?;
                let top = tree.top_node().ok_or(SolidError::null_pointer(OP))?;
                if !bbox.is_contained_by(&tree.node(top)?.bbox, containment_tol) {
                    crate::error::ignore_status(tree.add_to_spatial_tree(
                        &bbox,
                        None,
                        NodeGeomType::Default,
                    ))?;
                }
            }
        }
        Ok(())
    }
}

fn resize_edge(
    k: &Kernel,
    b: &Brep,
    e: EdgeId,
    curve_box: &mut Extent3d,
    curve_union: &mut Extent3d,
    vertex_union: &mut Extent3d,
    tol: &mut f64,
    keys: &dyn Fn(CacheObject) -> object::ObjectRef,
    owner_of: &dyn Fn(object::ObjectRef) -> Option<object::ObjectRef>,
) -> Result<()> {
    let data = b.edge(e)?;
    if data.tolerance > *tol {
        *tol = data.tolerance;
    }
    let owner = data
        .curve
        .as_ref()
        .ok_or(SolidError::null_pointer(OP))?
        .edge(owner_of)?
        .ok_or(SolidError::null_pointer(OP))?;
    let owner = edge_id(b, owner, keys)?;
    match crate::query::edge::calculate_bounding_box(
        k,
        b,
        owner,
        BoundingBoxRequest {
            normal: true,
            pseudo: false,
            polar: false,
        },
        false,
    ) {
        Ok(boxes) => {
            if let Some(bbox) = boxes.normal {
                *curve_box = bbox;
            }
        }
        Err(e) if e.is_unported() => return Err(e),
        Err(_) => {}
    }
    *curve_union = curve_union.union(curve_box);
    let first = edge::vertex(b, e)?.ok_or(SolidError::null_pointer(OP))?;
    let second = edge::other_vertex(b, e, Some(first))?.ok_or(SolidError::null_pointer(OP))?;
    for v in [first, second] {
        let data = b.vertex(v)?;
        if data.tolerance > *tol {
            *tol = data.tolerance;
        }
        vertex_union.add_point(data.point);
    }
    Ok(())
}

fn brep_bounding_box(k: &Kernel, b: &Brep) -> Result<Extent3d> {
    crate::query::brep::calculate_loose_bounding_box(k, b, false, true)
}

fn brep_vertices(b: &Brep) -> Result<Vec<VertexId>> {
    crate::topology::brep::p1::vertices(b, None)
}

fn brep_edges(b: &Brep) -> Result<Vec<EdgeId>> {
    crate::topology::brep::p1::edges(b, None)
}

fn brep_faces(b: &Brep) -> Result<Vec<FaceId>> {
    crate::topology::brep::p1::faces(b, None)
}

fn face_bounding_box(k: &Kernel, b: &Brep, f: FaceId) -> Result<Extent3d> {
    crate::query::face::p1::calculate_loose_bounding_box(k, b, f)
}

fn face_edges(b: &Brep, f: FaceId) -> Result<Vec<EdgeId>> {
    crate::topology::topology_traverser::face_edges(b, f)
}

/// Builds a new cache for `SmSAGObject::CacheMakeOrValidate`. The source geometry
/// is required, as for the curve/surface cache factories.
pub fn make_brep_cache(
    k: &Kernel,
    source: CacheSource<'_>,
    keys: &dyn Fn(CacheObject) -> object::ObjectRef,
) -> Result<BrepCache> {
    let mut cache = BrepCache::new(source.type_code())?;
    cache.build_trees(k, source, |object| object, keys)?;
    Ok(cache)
}

impl crate::geom::solution_array::SolutionObject for CacheObject {
    fn kind(&self) -> crate::geom::solution_array::ObjectKind {
        use crate::geom::solution_array::ObjectKind;
        match self {
            Self::Vertex(_) => ObjectKind::Vertex,
            Self::Curve(_) => ObjectKind::Curve,
            Self::Surface(_) => ObjectKind::Surface,
            Self::Edge(_) => ObjectKind::Edge,
            Self::Face(_) => ObjectKind::Face,
            Self::PolyVertex(_) => ObjectKind::PolyVertex,
            Self::PolyEdge(_) => ObjectKind::PolyEdge,
            Self::PolyFace(_) => ObjectKind::PolyFace,
            Self::Brep | Self::Other(_) => ObjectKind::Other,
        }
    }
}

/// `SmTol::GetZoneTol3d(const SmObject*)` for exact B-rep query handles.
/// A null source uses the system tolerance. Geometry uses its actual topology
/// owner (`Curve::edge` may follow intermediate owners), not its storage record.
///
/// # Errors
/// Missing records or geometry, an unresolved topology owner, or a handle outside
/// the supplied B-rep. Such a handle is never treated as a null source.
pub fn object_zone_tol_3d(
    kernel: &Kernel,
    brep: &Brep,
    source: Option<CacheObject>,
    keys: &dyn Fn(CacheObject) -> object::ObjectRef,
    owner_of: &dyn Fn(object::ObjectRef) -> Option<object::ObjectRef>,
) -> Result<f64> {
    use crate::geom::tol::{zone_tol_3d, SmallTopoType};
    let system = || zone_tol_3d(kernel, SmallTopoType::NormalSize);
    match source {
        None => Ok(system()),
        Some(CacheObject::Brep) => Ok(brep.tolerance()),
        Some(CacheObject::Vertex(v)) => Ok(brep.vertex(v)?.tolerance),
        Some(CacheObject::Edge(e)) => Ok(brep.edge(e)?.tolerance),
        Some(CacheObject::Face(f)) => Ok(brep.face(f)?.tolerance),
        Some(CacheObject::Curve(e)) => {
            let curve = brep
                .edge(e)?
                .curve
                .as_ref()
                .ok_or(SolidError::null_pointer(OP))?;
            match curve.edge(owner_of)? {
                Some(owner) => Ok(brep.edge(edge_id(brep, owner, keys)?)?.tolerance),
                None => Ok(system()),
            }
        }
        Some(CacheObject::Surface(f)) => {
            let surface = brep
                .face(f)?
                .surface
                .as_ref()
                .ok_or(SolidError::null_pointer(OP))?;
            match surface.face()? {
                Some(owner) => Ok(brep.face(face_id(brep, owner, keys)?)?.tolerance),
                None => Ok(system()),
            }
        }
        Some(_) => Err(SolidError::invalid_input(
            OP,
            "object is not in source B-rep",
        )),
    }
}

/// Exact B-rep geometry and its cache as seen by the topology solver. Face
/// classification may populate UV trim data through the mutable source borrow;
/// it does not invalidate the separately borrowed spatial trees. Geometry/topology
/// edits remain the caller's cache-invalidation responsibility. Both borrows must
/// name the same B-rep, and owner resolution uses the caller's identities.
pub struct BrepObjects<'a> {
    /// Source arenas.
    pub brep: &'a mut Brep,
    /// Identity of the source B-rep in the caller's shape store.
    pub brep_id: crate::topology::BrepId,
    /// Cache built from these arenas with identity handles.
    pub cache: &'a BrepCache,
    /// Identity of each edge/face record (independent of its arena index).
    pub keys: &'a dyn Fn(CacheObject) -> object::ObjectRef,
    /// Owners of intermediate generated curves and edge uses, as for `Curve::edge`.
    pub owner_of: &'a dyn Fn(object::ObjectRef) -> Option<object::ObjectRef>,
}

impl crate::intersect::topology_solver::p1::TopologyObjects<CacheObject> for BrepObjects<'_> {
    fn with_brep_trees(
        &mut self,
        _kernel: &Kernel,
        first: CacheObject,
        second: Option<CacheObject>,
        solve: &mut crate::intersect::topology_solver::p1::BrepTreeSolve<'_, CacheObject>,
    ) -> Result<crate::geom::solution_array::SolutionArray<CacheObject>> {
        if first != CacheObject::Brep {
            return Err(SolidError::invalid_input(OP, "expected B-rep handle"));
        }
        // The cache lives outside this adapter. Copy its reference, not its trees,
        // so the live source can be borrowed mutably while they remain published.
        let cache = self.cache;
        let trees = crate::intersect::topology_solver::p2::TopologyTrees {
            vertices: cache.vertex_tree().ok_or(SolidError::null_pointer(OP))?,
            curves: cache.curve_tree().ok_or(SolidError::null_pointer(OP))?,
            surfaces: cache.surface_tree().ok_or(SolidError::null_pointer(OP))?,
        };
        let second = match second {
            Some(CacheObject::Brep) => Some(trees),
            Some(_) => return Err(SolidError::invalid_input(OP, "expected B-rep handle")),
            None => None,
        };
        solve(self, trees, second)
    }

    fn object_zone_tol_3d(&self, kernel: &Kernel, object: Option<CacheObject>) -> Result<f64> {
        object_zone_tol_3d(kernel, self.brep, object, self.keys, self.owner_of)
    }

    fn face_point_classify(
        &mut self,
        kernel: &Kernel,
        face: CacheObject,
        uv: crate::math::Vec2,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        use_uv: bool,
        classification: &mut crate::intersect::point_class::PointClassification<crate::EntityRef>,
    ) -> Result<()> {
        let CacheObject::Face(face) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        crate::query::face::p2::point_classify(
            kernel,
            self.brep,
            face,
            uv,
            classification,
            crate::query::face::p2::PointClassifyOptions {
                src_zone_tol_3d: source_zone_tolerance,
                do_boundary_intersections: do_boundaries,
                use_uv_classification: use_uv,
                ..crate::query::face::p2::PointClassifyOptions::default()
            },
        )
    }

    fn face_point_3d_classify(
        &mut self,
        kernel: &Kernel,
        face: CacheObject,
        point: crate::math::Point3,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        classification: &mut crate::intersect::point_class::PointClassification<crate::EntityRef>,
    ) -> Result<()> {
        let CacheObject::Face(id) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        let Some(uv) = crate::intersect::topology_solver::p1::face_point_3d_uv_with_context(
            kernel,
            self,
            face,
            id,
            point,
            source_zone_tolerance,
            classification,
            None,
        )?
        else {
            return Ok(());
        };
        crate::query::face::p1::point_3d_classify_finish(
            kernel,
            self.brep,
            id,
            uv,
            point,
            source_zone_tolerance,
            do_boundaries,
            classification,
            crate::query::face::p1::Point3dClassifyOptions::default(),
        )
    }

    fn face_brep_editing_enabled(&self, face: CacheObject) -> Result<bool> {
        let CacheObject::Face(face) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        self.brep.face(face)?;
        Ok(self.brep.editing_enabled())
    }

    fn set_face_brep_editing_enabled(&mut self, face: CacheObject, editing: bool) -> Result<()> {
        let CacheObject::Face(face) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        self.brep.face(face)?;
        self.brep.set_editing_enabled(editing);
        Ok(())
    }

    fn face_create_uv_trim_curves(&mut self, kernel: &Kernel, face: CacheObject) -> Result<()> {
        let CacheObject::Face(face) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        crate::topology::face::p3::create_uv_trim_curves(kernel, self.brep, face, true, None)
            .map(|_| ())
    }

    fn face_classification_object(
        &self,
        face: CacheObject,
        entity: crate::EntityRef,
    ) -> Result<CacheObject> {
        let CacheObject::Face(face) = face else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        self.brep.face(face)?;
        match entity {
            crate::EntityRef::Vertex(id) => {
                self.brep.vertex(id)?;
                Ok(CacheObject::Vertex(id))
            }
            crate::EntityRef::Edge(id) => {
                self.brep.edge(id)?;
                Ok(CacheObject::Edge(id))
            }
            crate::EntityRef::Face(id) => {
                self.brep.face(id)?;
                Ok(CacheObject::Face(id))
            }
            _ => Err(SolidError::invalid_input(
                OP,
                "classification entity has no cached topology handle",
            )),
        }
    }

    fn brep_cache(&self, _kernel: &Kernel, object: CacheObject) -> Result<&BrepCache> {
        if object != CacheObject::Brep {
            return Err(SolidError::invalid_input(OP, "expected B-rep handle"));
        }
        Ok(self.cache)
    }
    fn curve(&self, object: CacheObject) -> Result<&crate::curve::curve::Curve> {
        let CacheObject::Curve(e) = object else {
            return Err(SolidError::invalid_input(OP, "expected curve handle"));
        };
        self.brep
            .edge(e)?
            .curve
            .as_ref()
            .ok_or(SolidError::null_pointer(OP))
    }
    fn analytic_curve_handle(&self, object: CacheObject) -> Option<CacheObject> {
        Some(object)
    }
    fn vertex_point(&self, object: CacheObject) -> Result<crate::math::Point3> {
        let CacheObject::Vertex(v) = object else {
            return Err(SolidError::invalid_input(OP, "expected vertex handle"));
        };
        Ok(self.brep.vertex(v)?.point)
    }
    fn vertex_tolerance(&self, object: CacheObject) -> Result<f64> {
        let CacheObject::Vertex(v) = object else {
            return Err(SolidError::invalid_input(OP, "expected vertex handle"));
        };
        Ok(self.brep.vertex(v)?.tolerance)
    }
    fn curve_edge(&self, object: CacheObject) -> Result<Option<CacheObject>> {
        let owner = self.curve(object)?.edge(self.owner_of)?;
        owner
            .map(|o| edge_id(self.brep, o, self.keys).map(CacheObject::Edge))
            .transpose()
    }
    fn edge_curve(&self, object: CacheObject) -> Result<Option<CacheObject>> {
        let CacheObject::Edge(e) = object else {
            return Err(SolidError::invalid_input(OP, "expected edge handle"));
        };
        Ok(self
            .brep
            .edge(e)?
            .curve
            .as_ref()
            .map(|_| CacheObject::Curve(e)))
    }
    fn edge_interval(&self, object: CacheObject) -> Result<crate::geom::extent1d::Extent1d> {
        let CacheObject::Edge(e) = object else {
            return Err(SolidError::invalid_input(OP, "expected edge handle"));
        };
        Ok(self.brep.edge(e)?.interval)
    }
    fn surface_face(&self, object: CacheObject) -> Result<Option<CacheObject>> {
        let owner = self.surface(object)?.face()?;
        owner
            .map(|o| face_id(self.brep, o, self.keys).map(CacheObject::Face))
            .transpose()
    }
    fn face_surface(&self, object: CacheObject) -> Result<CacheObject> {
        let CacheObject::Face(f) = object else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        self.brep
            .face(f)?
            .surface
            .as_ref()
            .ok_or(SolidError::null_pointer(OP))?;
        Ok(CacheObject::Surface(f))
    }
    fn surface(&self, object: CacheObject) -> Result<&crate::surface::surface::Surface> {
        let CacheObject::Surface(f) = object else {
            return Err(SolidError::invalid_input(OP, "expected surface handle"));
        };
        self.brep
            .face(f)?
            .surface
            .as_ref()
            .ok_or(SolidError::null_pointer(OP))
    }
    fn face_uv_domain(&self, object: CacheObject) -> Result<crate::geom::extent2d::Extent2d> {
        let CacheObject::Face(f) = object else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        Ok(self.brep.face(f)?.uv_domain)
    }
    fn face_tolerance(&self, object: CacheObject) -> Result<f64> {
        let CacheObject::Face(f) = object else {
            return Err(SolidError::invalid_input(OP, "expected face handle"));
        };
        Ok(self.brep.face(f)?.tolerance)
    }
    fn surface_owner(&self, object: CacheObject) -> Result<Option<CacheObject>> {
        self.surface(object)?
            .owner()
            .map(|owner| {
                if owner.type_code == crate::topology::face::p2::TYPE {
                    face_id(self.brep, owner, self.keys).map(CacheObject::Face)
                } else {
                    Ok(CacheObject::Other(owner))
                }
            })
            .transpose()
    }
    fn shape_member(&self, object: CacheObject) -> object::ShapeMember {
        use crate::geom::object::ShapeMember;
        let b = self.brep_id;
        match object {
            CacheObject::Vertex(v) => ShapeMember::Vertex((b, v)),
            CacheObject::Edge(e) => ShapeMember::Edge((b, e)),
            CacheObject::Face(f) => ShapeMember::Face((b, f)),
            _ => ShapeMember::Other,
        }
    }
}

/// Conventional identities for a single B-rep: the type code and raw arena index.
/// Callers using other `ObjectKey`s pass their own mapping to build/resize instead.
pub fn arena_identity(object: CacheObject) -> object::ObjectRef {
    use crate::geom::attribute::ObjectKey;
    let (type_code, key) = match object {
        CacheObject::Edge(e) | CacheObject::Curve(e) => (edge::TYPE, u64::from(e.index())),
        CacheObject::Face(f) | CacheObject::Surface(f) => {
            (crate::topology::face::p2::TYPE, u64::from(f.index()))
        }
        CacheObject::Vertex(v) => (crate::topology::vertex::TYPE, u64::from(v.index())),
        _ => (BREP_TYPE, 0),
    };
    object::ObjectRef {
        type_code,
        key: ObjectKey(key),
    }
}

fn edge_id(
    b: &Brep,
    owner: object::ObjectRef,
    keys: &dyn Fn(CacheObject) -> object::ObjectRef,
) -> Result<EdgeId> {
    b.edges()
        .iter()
        .find(|(id, _)| keys(CacheObject::Edge(*id)) == owner)
        .map(|(id, _)| id)
        .ok_or_else(|| SolidError::invalid_input(OP, "edge owner is not in source B-rep"))
}
fn face_id(
    b: &Brep,
    owner: object::ObjectRef,
    keys: &dyn Fn(CacheObject) -> object::ObjectRef,
) -> Result<FaceId> {
    b.faces()
        .iter()
        .find(|(id, _)| keys(CacheObject::Face(*id)) == owner)
        .map(|(id, _)| id)
        .ok_or_else(|| SolidError::invalid_input(OP, "face owner is not in source B-rep"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A triangular trim strictly inside the plane's rectangular UV domain.
    fn triangle(b: &mut Brep) -> FaceId {
        use crate::curve::curve::Curve;
        use crate::curve::line::{BoundedLineOptions, Line};
        use crate::geom::extent1d::Extent1d;
        use crate::geom::extent2d::Extent2d;
        use crate::math::{Point3, Vec2, Vec3};
        use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
        use crate::surface::plane::Plane;
        use crate::surface::surface::Surface;
        use crate::topology::{
            Edge, Edgeuse, EdgeuseId, EdgeuseOwner, Face, Faceuse, FaceuseId, Loop, LoopId,
            Loopuse, LoopuseId, LoopuseStart, Vertex, Vertexuse, VertexuseId, VertexuseOwner,
        };
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
        let corners = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
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
        for i in 0..3u32 {
            let (x0, y0) = corners[i as usize];
            let (x1, y1) = corners[((i + 1) % 3) as usize];
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
                vertices: [VertexId::from_index(i), VertexId::from_index((i + 1) % 3)],
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
                ccw: Some(EdgeuseId::from_index((i + 1) % 3)),
                cw: Some(EdgeuseId::from_index((i + 2) % 3)),
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
        for id in (0..3).map(EdgeId::from_index) {
            b.edge_mut(id)
                .unwrap()
                .curve
                .as_mut()
                .unwrap()
                .set_owner(Some(arena_identity(CacheObject::Edge(id))), None);
        }
        b.face_mut(f)
            .unwrap()
            .surface
            .as_mut()
            .unwrap()
            .set_owner(Some(arena_identity(CacheObject::Face(f))));
        f
    }

    #[test]
    fn borrowed_cache_queries_classify_the_live_trimmed_face() {
        use crate::geom::global_solver::{SolutionRequested, SolverOperation};
        use crate::intersect::topology_solver::p1::TopologySolver;
        use crate::intersect::topology_solver::p2::TopologySolveOptions;
        use crate::math::{Point3, Vec3, BIG_DOUBLE};
        use crate::topology::BrepId;

        let kernel = Kernel::default();
        let mut brep = Brep::new();
        let face = triangle(&mut brep);
        let cache = make_brep_cache(&kernel, CacheSource::Brep(&brep), &arena_identity).unwrap();
        let vectors = [Vec3::new(0.0, 0.0, -1.0)];
        let options = TopologySolveOptions {
            operation: SolverOperation::Rayfire,
            requested: SolutionRequested::All,
            tolerance_3d: 1e-6,
            best_answer_so_far_sq: BIG_DOUBLE,
            opt_vectors: &vectors,
        };
        let mut objects = BrepObjects {
            brep: &mut brep,
            brep_id: BrepId::from_index(9),
            cache: &cache,
            keys: &arena_identity,
            owner_of: &|_| None,
        };
        // Both rays meet the plane inside its UV domain. Only the first meets
        // the actual triangular face. Reuse the same published spatial trees.
        for (point, expected) in [
            (Point3::new(0.25, 0.25, 1.0), 1),
            (Point3::new(0.75, 0.75, 1.0), 0),
            (Point3::new(0.25, 0.25, 1.0), 1),
        ] {
            let answers = TopologySolver::brep_point_solve(
                &kernel,
                &mut objects,
                Some(CacheObject::Brep),
                point,
                options,
            )
            .unwrap();
            assert_eq!(answers.len(), expected);
            if let Some(answer) = answers.as_slice().first() {
                assert_eq!(answer.objects[0], Some(CacheObject::Face(face)));
                assert_eq!(answer.num_variables, 2);
                assert!((answer.start.solution_value - 1.0).abs() < 1e-12);
                assert!((answer.start.param(0).unwrap() - 0.25).abs() < 1e-12);
                assert!((answer.start.param(1).unwrap() - 0.25).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn source_tolerance_follows_actual_owners_and_null_uses_system() {
        use crate::curve::curve::Curve;
        use crate::curve::line::Line;
        use crate::geom::attribute::ObjectKey;
        use crate::geom::extent2d::Extent2d;
        use crate::kernel::KernelOptions;
        use crate::math::{Point3, Vec2, Vec3};
        use crate::surface::plane::Plane;
        use crate::surface::surface::Surface;
        use crate::topology::{Edge, Face, Vertex};

        let mut options = KernelOptions::default();
        options.tolerances.zone_tol_3d = 0.125;
        let kernel = Kernel::new(&options);
        let mut brep = Brep::new();
        let vertex = brep
            .add_vertex_record(Vertex {
                tolerance: 0.03125,
                ..Vertex::default()
            })
            .unwrap();
        let storage = brep.add_edge_record(Edge::new([vertex, vertex])).unwrap();
        let owner = brep.add_edge_record(Edge::new([vertex, vertex])).unwrap();
        brep.edge_mut(storage).unwrap().tolerance = 0.0625;
        brep.edge_mut(owner).unwrap().tolerance = 0.25;
        let keys = |o| {
            let mut key = arena_identity(o);
            key.key.0 += 1_000_000;
            key
        };
        let intermediate = object::ObjectRef {
            type_code: crate::curve::curve::p2::TYPE,
            key: ObjectKey(700),
        };
        let owner_of = |o| (o == intermediate).then(|| keys(CacheObject::Edge(owner)));
        let mut curve = Curve::from(
            Line::create_line_segment(3, Point3::ZERO, Point3::new(1.0, 0.0, 0.0), None).unwrap(),
        );
        curve.set_owner(Some(intermediate), None);
        brep.edge_mut(storage).unwrap().curve = Some(curve);
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face_record(Face::new(shell)).unwrap();
        let surface_owner = brep.add_face_record(Face::new(shell)).unwrap();
        brep.face_mut(face).unwrap().tolerance = 0.5;
        brep.face_mut(surface_owner).unwrap().tolerance = 1.0;
        let mut surface = Surface::from(Plane::new(
            Point3::ZERO,
            Vec3::X,
            Vec3::Y,
            Vec2::new(1.0, 1.0),
            &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
        ));
        surface.set_owner(Some(keys(CacheObject::Face(surface_owner))));
        brep.face_mut(face).unwrap().surface = Some(surface);
        for (source, expected) in [
            (None, 0.125),
            (Some(CacheObject::Vertex(vertex)), 0.03125),
            (Some(CacheObject::Edge(storage)), 0.0625),
            (Some(CacheObject::Curve(storage)), 0.25),
            (Some(CacheObject::Surface(face)), 1.0),
        ] {
            assert_eq!(
                object_zone_tol_3d(&kernel, &brep, source, &keys, &owner_of).unwrap(),
                expected
            );
        }
        assert!(matches!(
            object_zone_tol_3d(
                &kernel,
                &brep,
                Some(CacheObject::Curve(storage)),
                &arena_identity,
                &owner_of
            ),
            Err(SolidError::InvalidInput { .. })
        ));
        brep.edge_mut(storage)
            .unwrap()
            .curve
            .as_mut()
            .unwrap()
            .set_owner(None, None);
        assert_eq!(
            object_zone_tol_3d(
                &kernel,
                &brep,
                Some(CacheObject::Curve(storage)),
                &keys,
                &owner_of,
            )
            .unwrap(),
            0.125
        );
    }

    #[test]
    fn solver_adapter_uses_caller_owner_chain_and_brep_identity() {
        use crate::curve::curve::Curve;
        use crate::curve::line::Line;
        use crate::geom::attribute::ObjectKey;
        use crate::intersect::topology_solver::p1::TopologyObjects;
        use crate::math::Point3;
        use crate::topology::{BrepId, Edge, Vertex};

        let mut b = Brep::new();
        let a = b.add_vertex_record(Vertex::default()).unwrap();
        let z = b.add_vertex_record(Vertex::default()).unwrap();
        let e = b.add_edge_record(Edge::new([a, z])).unwrap();
        let keys = |o| {
            let mut key = arena_identity(o);
            key.key.0 += 1000000000;
            key
        };
        let generated = object::ObjectRef {
            type_code: crate::curve::curve::p2::TYPE,
            key: ObjectKey(500),
        };
        let edgeuse = object::ObjectRef {
            type_code: crate::topology::edgeuse::p1::TYPE,
            key: ObjectKey(600),
        };
        let owner_of = |o| {
            if o == generated {
                Some(edgeuse)
            } else if o == edgeuse {
                Some(keys(CacheObject::Edge(e)))
            } else {
                None
            }
        };
        let mut curve = Curve::from(
            Line::create_line_segment(3, Point3::ZERO, Point3::new(1.0, 0.0, 0.0), None).unwrap(),
        );
        curve.set_owner(Some(generated), None);
        b.edge_mut(e).unwrap().curve = Some(curve);
        let cache = BrepCache::new(BREP_TYPE).unwrap();
        let brep_id = BrepId::from_index(7);
        let mut objects = BrepObjects {
            brep: &mut b,
            brep_id,
            cache: &cache,
            keys: &keys,
            owner_of: &owner_of,
        };
        assert_eq!(
            objects.curve_edge(CacheObject::Curve(e)).unwrap(),
            Some(CacheObject::Edge(e))
        );
        assert_eq!(
            objects.shape_member(CacheObject::Edge(e)),
            object::ShapeMember::Edge((brep_id, e))
        );
        objects.keys = &arena_identity;
        assert!(matches!(
            objects.curve_edge(CacheObject::Curve(e)),
            Err(SolidError::InvalidInput { .. })
        ));
    }

    #[test]
    fn initialization_contract_and_atomic_build_failure() {
        let k = Kernel::default();
        assert!(matches!(
            BrepCache::<CacheObject>::new(0),
            Err(SolidError::Upstream { code: 1001, .. })
        ));
        let mut cache = BrepCache::new(BREP_TYPE).unwrap();
        assert!(cache.vertex_tree().is_none());
        let mut b = Brep::new();
        // A real incomplete face fails its box calculation before the cache can
        // publish any trees; an empty B-rep no longer stops at a query.brep stub.
        let shell = b.add_shell().unwrap();
        b.add_face(shell).unwrap();
        assert!(matches!(
            cache.resize_trees(&k, &b, &[], &[], &[], &arena_identity, &|_| None),
            Err(SolidError::InvalidInput { .. })
        ));
        assert!(matches!(
            cache.build_trees(&k, CacheSource::Brep(&b), |o| o, &arena_identity),
            Err(SolidError::InvalidInput { .. })
        ));
        assert!(
            cache.vertex_tree().is_none()
                && cache.curve_tree().is_none()
                && cache.surface_tree().is_none()
        );
        let p = PolyBrep::new(&k, 1e-6);
        assert!(matches!(
            cache.build_trees(&k, CacheSource::PolyBrep(&p), |o| o, &arena_identity),
            Err(SolidError::InvalidInput { .. })
        ));
    }

    #[test]
    fn resize_retains_object_identity_and_rejects_dangling_ids() {
        let k = Kernel::default();
        let mut p = PolyBrep::new(&k, 1e-6);
        let v = p.add_vertex(&k, crate::math::Point3::ZERO, 1e-6).unwrap();
        let mut cache = make_brep_cache(&k, CacheSource::PolyBrep(&p), &arena_identity).unwrap();
        assert!(matches!(
            cache.build_trees(&k, CacheSource::PolyBrep(&p), |o| o, &arena_identity),
            Err(SolidError::InvalidInput { .. })
        ));
        let mut b = Brep::new();
        let extra = b
            .add_vertex_record(crate::topology::Vertex {
                point: crate::math::Point3::new(100.0, 0.0, 0.0),
                tolerance: 0.25,
                ..crate::topology::Vertex::default()
            })
            .unwrap();
        cache
            .resize_trees(&k, &b, &[], &[], &[extra], &arena_identity, &|_| None)
            .unwrap();
        let t = cache.vertex_tree().unwrap();
        let list = t.object_list().unwrap();
        assert_eq!(list.len(), 1);
        assert_eq!(
            t.list(list[0]).unwrap().object,
            Some(CacheObject::PolyVertex(v))
        );
        assert!(t.node(t.top_node().unwrap()).unwrap().bbox.max().x >= 125.0);
        assert!(matches!(
            cache.resize_trees(
                &k,
                &b,
                &[],
                &[],
                &[VertexId::from_index(u32::MAX)],
                &arena_identity,
                &|_| None,
            ),
            Err(SolidError::NotFound { .. })
        ));
    }
}
