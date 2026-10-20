//! Arena access for the existing generic topology intersector used by merge p4.
//!
//! Upstream: `source/SMLib/src/SmTopologyIntersector.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | `ArenaPair` (private) | operand access through the intersector's topology handles |
//!
//! Deviations: borrowed arenas replace pointers; generated curves have local handles.
//! This adapts the ported operations without copying their geometry algorithms.
//! No input-dependent recursion is introduced.

use crate::boolean::topology_intersector::p1 as ti;
use crate::boolean::topology_intersector::p1::IntersectorBrep;
use crate::boolean::topology_intersector::p2::IntersectorBrepP2;
use crate::curve::curve::Curve as SmCurve;
use crate::geom::a_object::BREP_TYPE;
use crate::geom::attribute::{AttributeStore, ObjectKey};
use crate::geom::context::Context;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::geom::object::{ObjectRef, ShapeMember};
use crate::geom::solution_array::{ObjectKind, SolutionArray, SolutionObject};
use crate::geom::tol::{update_object_tolerance, TopologyToUpdate};
use crate::intersect::curve_class::p1::CurveClassification;
use crate::intersect::point_class::{
    ClassGeometry, ClassObject, ClassObjectKind, PointClassification,
};
use crate::intersect::topology_solver::p1::TopologyObjects;
use crate::kernel::Kernel;
use crate::math::{Point2, Point3};
use crate::nurbs::b_spline_curve::p2::BoundingBoxRequest;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::query::face::p2::{
    curve_on_classify, point_classify, CurveOnClassifyObjects, PointClassifyOptions,
};
use crate::surface::surface::Surface;
use crate::surface::volume::Volume;
use crate::topology::brep::p1 as bp1;
use crate::topology::brep::p7;
use crate::topology::brep_cache::{BrepCache, CacheObject, CacheSource};
use crate::topology::edge as ed;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::face::p2 as fp2;
use crate::topology::face::p3 as fp3;
use crate::topology::face::p3::{FaceIntersectOptions, IntersectObject};
use crate::topology::proto_topology::p1::ProtoTopology;
use crate::topology::proto_topology::{ProtoBrep, ProtoGeometry};
use crate::topology::topology::Topology;
use crate::topology::vertex as vx;
use crate::topology::{
    Brep, BrepId, EdgeId, EdgeuseId, FaceId, FaceuseId, LoopId, LoopuseId, RegionId, ShellId,
    VertexId, VertexuseId,
};
use crate::{EntityRef, SolidError};
use std::collections::BTreeMap;
const K_VERTEX: u8 = 0;
const K_EDGE: u8 = 1;
const K_FACE: u8 = 2;
const K_EDGE_CURVE: u8 = 3;
const K_FACE_SURFACE: u8 = 4;
const K_CURVE: u8 = 5;
const K_BREP: u8 = 6;
const K_EXTRA: u8 = 7;
const K_EDGEUSE: u8 = 8;
const K_VERTEXUSE: u8 = 9;
const K_LOOP: u8 = 10;
const K_TI: u8 = 11;
const K_TRIM_CURVE: u8 = 12;
const K_NULL: u8 = 13;
const K_REGION: u8 = 14;
const K_FACEUSE: u8 = 15;
const K_LOOPUSE: u8 = 16;
const K_SHELL: u8 = 17;
/// An object handle: the B-rep (`2` for the curve stores, `-1` for the intersector),
/// the kind and the index. Ordered by operand, class and arena index. A
/// curve handle also carries its curve's natural interval (what `SmSolutionArray`
/// reads through `SmCurve::GetNaturalInterval`); it takes no part in identity.
#[derive(Debug, Clone, Copy)]
pub(super) struct H {
    b: i8,
    kind: u8,
    i: u32,
    ivl: [f64; 2],
}

impl PartialEq for H {
    fn eq(&self, o: &Self) -> bool {
        (self.b, self.kind, self.i) == (o.b, o.kind, o.i)
    }
}

impl Eq for H {}

impl PartialOrd for H {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}

impl Ord for H {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        (self.b, self.kind, self.i).cmp(&(o.b, o.kind, o.i))
    }
}

pub(super) const fn hd(b: i8, kind: u8, i: u32) -> H {
    H {
        b,
        kind,
        i,
        ivl: [f64::NAN, f64::NAN],
    }
}

/// A curve handle carrying the curve's natural interval.
fn curve_h(b: i8, kind: u8, i: u32, c: &SmCurve) -> H {
    let ivl = c.natural_interval();
    H {
        b,
        kind,
        i,
        ivl: [ivl.min(), ivl.max()],
    }
}

impl ClassObject for H {
    fn class_kind(&self) -> ClassObjectKind {
        match self.kind {
            K_REGION => ClassObjectKind::Region,
            K_FACEUSE => ClassObjectKind::Faceuse,
            K_LOOPUSE => ClassObjectKind::Loopuse,
            K_VERTEX => ClassObjectKind::Vertex,
            K_EDGE => ClassObjectKind::Edge,
            K_FACE => ClassObjectKind::Face,
            K_EDGE_CURVE | K_CURVE | K_EXTRA | K_TRIM_CURVE => ClassObjectKind::Curve,
            K_FACE_SURFACE => ClassObjectKind::Surface,
            K_EDGEUSE => ClassObjectKind::Edgeuse,
            K_VERTEXUSE => ClassObjectKind::Vertexuse,
            K_LOOP => ClassObjectKind::Loop,
            _ => ClassObjectKind::Other,
        }
    }
}

impl SolutionObject for H {
    fn kind(&self) -> ObjectKind {
        match self.kind {
            K_VERTEX => ObjectKind::Vertex,
            K_EDGE => ObjectKind::Edge,
            K_FACE => ObjectKind::Face,
            K_EDGE_CURVE | K_CURVE | K_EXTRA | K_TRIM_CURVE => ObjectKind::Curve,
            K_FACE_SURFACE => ObjectKind::Surface,
            _ => ObjectKind::Other,
        }
    }
    fn natural_interval(&self) -> Option<Extent1d> {
        match self.kind {
            K_EDGE_CURVE | K_CURVE | K_EXTRA if self.ivl.first()?.is_finite() => Some(
                Extent1d::from_min_max(*self.ivl.first()?, *self.ivl.get(1)?),
            ),
            _ => None,
        }
    }
}

pub(super) const TI_HANDLE: H = hd(-1, K_TI, 0);

pub(super) const fn brep_h(b: i8) -> H {
    hd(b, K_BREP, b as u32)
}

#[derive(Debug)]
pub(super) struct ArenaPair<'a> {
    breps: [Option<&'a mut Brep>; 2],
    curves: &'a mut Vec<SmCurve>,
    extra: &'a mut Vec<SmCurve>,
    records: &'a mut BTreeMap<ObjectKey, Topology>,
    attributes: &'a mut AttributeStore,
    owner_links: &'a BTreeMap<ObjectRef, ObjectRef>,
    geometry_checks: &'a [Option<Box<dyn super::kernel_adapter::OperandGeometryCheck>>; 2],
    notification_checks: &'a [Option<Box<dyn super::kernel_adapter::OperandGeometryCheck>>; 2],
    primary_context: &'a mut Context,
    other_context: Option<&'a mut Context>,
    namespace: u32,
    other_kernel: Option<&'a Kernel>,
    point_caches: &'a mut [Option<BrepCache>; 2],
    region_queries: &'a mut [super::kernel_adapter::RegionQuerySlot; 2],
    caches: [std::cell::OnceCell<std::sync::Arc<BrepCache<H>>>; 2],
}
fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input("boolean::merge::p4::arena_pair", reason)
}
const fn vid(h: H) -> VertexId {
    VertexId::from_index(h.i)
}
const fn eid(h: H) -> EdgeId {
    EdgeId::from_index(h.i)
}
const fn fid(h: H) -> FaceId {
    FaceId::from_index(h.i)
}
impl ArenaPair<'_> {
    fn require_notifications(&self, slot: usize) -> crate::Result<()> {
        self.notification_checks
            .get(slot)
            .and_then(Option::as_deref)
            .map_or(Ok(()), super::kernel_adapter::OperandGeometryCheck::check)
    }

    fn require_query_geometry(&self, slot: usize) -> crate::Result<()> {
        self.geometry_checks
            .get(slot)
            .and_then(Option::as_deref)
            .map_or(Ok(()), super::kernel_adapter::OperandGeometryCheck::check)
    }

    fn classification_entity(&self, slot: i8, h: H) -> crate::Result<EntityRef> {
        if h.b != slot {
            return Err(invalid("classification handle is outside this B-rep"));
        }
        let b = self.brep_r(h)?;
        match h.kind {
            K_VERTEX => b.vertex(vid(h)).map(|_| EntityRef::Vertex(vid(h))),
            K_EDGE => b.edge(eid(h)).map(|_| EntityRef::Edge(eid(h))),
            K_FACE => b.face(fid(h)).map(|_| EntityRef::Face(fid(h))),
            K_LOOP => {
                let id = LoopId::from_index(h.i);
                b.loop_(id).map(|_| EntityRef::Loop(id))
            }
            K_EDGEUSE => {
                let id = EdgeuseId::from_index(h.i);
                b.edgeuse(id).map(|_| EntityRef::Edgeuse(id))
            }
            K_VERTEXUSE => {
                let id = VertexuseId::from_index(h.i);
                b.vertexuse(id).map(|_| EntityRef::Vertexuse(id))
            }
            K_REGION => {
                let id = RegionId::from_index(h.i);
                b.region(id).map(|_| EntityRef::Region(id))
            }
            K_FACEUSE => {
                let id = FaceuseId::from_index(h.i);
                b.faceuse(id).map(|_| EntityRef::Faceuse(id))
            }
            K_LOOPUSE => {
                let id = LoopuseId::from_index(h.i);
                b.loopuse(id).map(|_| EntityRef::Loopuse(id))
            }
            K_SHELL => {
                let id = ShellId::from_index(h.i);
                b.shell(id).map(|_| EntityRef::Shell(id))
            }
            K_BREP if u32::try_from(slot).ok() == Some(h.i) => {
                Ok(EntityRef::Brep(BrepId::from_index(h.i)))
            }
            _ => Err(invalid("classification handle is not registered topology")),
        }
    }

    pub(super) fn new<'a>(
        primary: &'a mut Brep,
        other: Option<&'a mut Brep>,
        curves: &'a mut Vec<SmCurve>,
        extra: &'a mut Vec<SmCurve>,
        records: &'a mut BTreeMap<ObjectKey, Topology>,
        attributes: &'a mut AttributeStore,
        owner_links: &'a BTreeMap<ObjectRef, ObjectRef>,
        geometry_checks: &'a [Option<Box<dyn super::kernel_adapter::OperandGeometryCheck>>; 2],
        notification_checks: &'a [Option<Box<dyn super::kernel_adapter::OperandGeometryCheck>>; 2],
        other_kernel: Option<&'a Kernel>,
        point_caches: &'a mut [Option<BrepCache>; 2],
        region_queries: &'a mut [super::kernel_adapter::RegionQuerySlot; 2],
        primary_context: &'a mut Context,
        other_context: Option<&'a mut Context>,
        namespace: u32,
    ) -> ArenaPair<'a> {
        ArenaPair {
            breps: [Some(primary), other],
            curves,
            extra,
            records,
            attributes,
            owner_links,
            geometry_checks,
            notification_checks,
            primary_context,
            other_context,
            namespace,
            other_kernel,
            point_caches,
            region_queries,
            caches: [std::cell::OnceCell::new(), std::cell::OnceCell::new()],
        }
    }
    fn brep(&self, h: H) -> Option<&Brep> {
        self.breps.get(usize::try_from(h.b).ok()?)?.as_deref()
    }

    fn brep_mut(&mut self, h: H) -> Option<&mut Brep> {
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        self.breps
            .get_mut(usize::try_from(h.b).ok()?)?
            .as_deref_mut()
    }

    fn brep_r(&self, h: H) -> crate::Result<&Brep> {
        self.brep(h)
            .ok_or_else(|| invalid("handle names an absent B-rep"))
    }

    fn brep_m(&mut self, h: H) -> crate::Result<&mut Brep> {
        self.brep_mut(h)
            .ok_or_else(|| invalid("handle names an absent B-rep"))
    }

    /// The edge's curve handle (`ProtoBrep::edge_curve`'s shape).
    #[allow(clippy::option_option)]
    fn edge_curve_h(&self, e: H) -> Option<Option<H>> {
        let b = self.brep(e)?;
        let edge = b.edge(eid(e)).ok()?;
        Some(
            edge.curve
                .as_ref()
                .map(|c| curve_h(e.b, K_EDGE_CURVE, e.i, c)),
        )
    }

    /// The face's surface handle (`ProtoBrep::face_surface`'s shape).
    #[allow(clippy::option_option)]
    fn face_surface_h(&self, f: H) -> Option<Option<H>> {
        let b = self.brep(f)?;
        let face = b.face(fid(f)).ok()?;
        Some(face.surface.as_ref().map(|_| hd(f.b, K_FACE_SURFACE, f.i)))
    }

    /// Gives a curve made by the port a handle in the curve store.
    fn register_curve(&mut self, curve: SmCurve) -> crate::Result<H> {
        self.insert_curve(curve)
            .ok_or_else(|| invalid("curve store is full"))
    }

    /// [`Two::register_curve`] for each curve present.
    fn register_curves(&mut self, curves: Vec<Option<SmCurve>>) -> crate::Result<Vec<Option<H>>> {
        curves
            .into_iter()
            .map(|c| c.map(|c| self.register_curve(c)).transpose())
            .collect()
    }

    fn entity_h(b: i8, e: EntityRef) -> Option<H> {
        match e {
            EntityRef::Region(v) => Some(hd(b, K_REGION, v.index())),
            EntityRef::Shell(v) => Some(hd(b, K_SHELL, v.index())),
            EntityRef::Faceuse(v) => Some(hd(b, K_FACEUSE, v.index())),
            EntityRef::Loopuse(v) => Some(hd(b, K_LOOPUSE, v.index())),
            EntityRef::Vertexuse(v) => Some(hd(b, K_VERTEXUSE, v.index())),
            EntityRef::Vertex(v) => Some(hd(b, K_VERTEX, v.index())),
            EntityRef::Edge(x) => Some(hd(b, K_EDGE, x.index())),
            EntityRef::Face(f) => Some(hd(b, K_FACE, f.index())),
            EntityRef::Loop(l) => Some(hd(b, K_LOOP, l.index())),
            EntityRef::Edgeuse(eu) => Some(hd(b, K_EDGEUSE, eu.index())),
            EntityRef::Brep(_) => None,
        }
    }
}
impl ClassGeometry<H> for ArenaPair<'_> {
    fn curve(&self, o: H) -> Option<&SmCurve> {
        match o.kind {
            K_EDGE_CURVE => self.brep(o)?.edge(eid(o)).ok()?.curve.as_ref(),
            K_CURVE => self.curves.get(o.i as usize),
            K_TRIM_CURVE => {
                eu1::uv_trim_curve(self.brep(o)?, EdgeuseId::from_index(o.i))
                    .ok()?
                    .0
            }
            K_EXTRA => self.extra.get(o.i as usize),
            _ => None,
        }
    }
    fn surface(&self, o: H) -> Option<&Surface> {
        match o.kind {
            K_FACE_SURFACE => self.brep(o)?.face(fid(o)).ok()?.surface.as_ref(),
            _ => None,
        }
    }
    fn volume(&self, _o: H) -> Option<&Volume> {
        None
    }
    fn point(&self, _o: H) -> Option<Point3> {
        None
    }
    // Classification reads the actual vertex points and stored tolerances.
    fn classification_vertex_point(&self, o: H) -> Option<Point3> {
        if o.kind != K_VERTEX {
            return None;
        }
        vx::point(self.brep(o)?, vid(o)).ok()
    }
    fn topology_zone_tol_3d(&self, o: H) -> Option<f64> {
        match o.kind {
            K_VERTEX | K_EDGE | K_FACE => self.topology_tolerance(o),
            _ => None,
        }
    }
}

impl ProtoGeometry<H> for ArenaPair<'_> {
    fn curve_mut(&mut self, o: H) -> Option<&mut SmCurve> {
        match o.kind {
            K_EDGE_CURVE => self.brep_mut(o)?.edge_mut(eid(o)).ok()?.curve.as_mut(),
            K_CURVE => self.curves.get_mut(o.i as usize),
            K_EXTRA => self.extra.get_mut(o.i as usize),
            _ => None,
        }
    }
}

/// The face-side objects of `curve_on_classify` for one B-rep.
struct Objs<'a, 'b> {
    two: &'a mut ArenaPair<'b>,
    b: i8,
}

impl ClassGeometry<H> for Objs<'_, '_> {
    fn curve(&self, o: H) -> Option<&SmCurve> {
        ClassGeometry::curve(self.two, o)
    }
    fn surface(&self, o: H) -> Option<&Surface> {
        ClassGeometry::surface(self.two, o)
    }
    fn volume(&self, o: H) -> Option<&Volume> {
        ClassGeometry::volume(self.two, o)
    }
    fn point(&self, o: H) -> Option<Point3> {
        ClassGeometry::point(self.two, o)
    }
}

impl crate::intersect::curve_class::p7::FaceQueryOwner<H> for Objs<'_, '_> {
    fn classification_brep(&self) -> crate::Result<&Brep> {
        self.two.brep_r(hd(self.b, K_BREP, u32::from(self.b == 1)))
    }
    fn classification_brep_mut(&mut self) -> crate::Result<&mut Brep> {
        let slot = usize::try_from(self.b).map_err(|_| invalid("bad operand"))?;
        // Only lazy query state is edited. Keep published spatial trees alive.
        self.two
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))
    }
    fn classification_entity(&self, face: H, object: H) -> crate::Result<EntityRef> {
        if face.b != self.b || face.kind != K_FACE {
            return Err(invalid("not a face of this operand"));
        }
        self.two.classification_entity(self.b, object)
    }
    fn classification_handle(&self, face: H, entity: EntityRef) -> crate::Result<H> {
        if face.b != self.b || face.kind != K_FACE {
            return Err(invalid("not a face of this operand"));
        }
        ArenaPair::entity_h(self.b, entity)
            .ok_or_else(|| invalid("unresolved classification identity"))
    }
    fn query_on_face<T>(
        &mut self,
        kernel: &Kernel,
        face: H,
        query: impl FnOnce(
            &Kernel,
            &mut Brep,
            FaceId,
            &dyn Fn(CacheObject) -> ObjectRef,
            &dyn Fn(ObjectRef) -> Option<ObjectRef>,
            Option<&mut dyn crate::query::face::p2::FacePointCache>,
        ) -> crate::Result<T>,
    ) -> crate::Result<T> {
        if face.b != self.b || face.kind != K_FACE {
            return Err(invalid("not a face of this operand"));
        }
        let slot = usize::try_from(self.b).map_err(|_| invalid("bad operand"))?;
        self.two.require_query_geometry(slot)?;
        let namespace = self.two.namespace;
        let keys = |object| super::kernel_adapter::query_key(namespace, slot == 1, object);
        let links = self.two.owner_links;
        let owner_of = |owner| links.get(&owner).copied();
        let kernel = if slot == 1 {
            self.two.other_kernel.unwrap_or(kernel)
        } else {
            kernel
        };
        let brep = self
            .two
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        query(kernel, brep, fid(face), &keys, &owner_of, None)
    }
}

impl CurveOnClassifyObjects<H> for Objs<'_, '_> {
    fn face(&self, f: FaceId) -> H {
        hd(self.b, K_FACE, f.index())
    }
    fn vertex(&self, v: VertexId) -> H {
        hd(self.b, K_VERTEX, v.index())
    }
    fn edge(&self, e: EdgeId) -> H {
        hd(self.b, K_EDGE, e.index())
    }
    fn vertexuse(&self, vu: VertexuseId) -> H {
        hd(self.b, K_VERTEXUSE, vu.index())
    }
    fn edgeuse(&self, eu: EdgeuseId) -> H {
        hd(self.b, K_EDGEUSE, eu.index())
    }
    fn edge_curve(&self, e: EdgeId) -> H {
        hd(self.b, K_EDGE_CURVE, e.index())
    }
    fn trim_curve(&self, eu: EdgeuseId) -> H {
        hd(self.b, K_TRIM_CURVE, eu.index())
    }
    fn face_surface(&self, f: FaceId) -> Option<H> {
        self.two
            .face_surface_h(hd(self.b, K_FACE, f.index()))
            .flatten()
    }
    fn loop_id(&self, tgt_loop: H) -> Option<LoopId> {
        (tgt_loop.kind == K_LOOP && tgt_loop.b == self.b).then(|| LoopId::from_index(tgt_loop.i))
    }
}

impl ProtoBrep<H> for ArenaPair<'_> {
    fn brep_point_3d_classify(
        &mut self,
        k: &Kernel,
        brep: H,
        point: Point3,
        point_zone_tolerance: f64,
        boundary_intersections: bool,
        classification: &mut PointClassification<H>,
    ) -> crate::Result<()> {
        if brep.kind != K_BREP {
            return Err(invalid("not a B-rep handle"));
        }
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        let namespace = self.namespace;
        let checks = self.geometry_checks;
        let check = || {
            checks
                .get(slot)
                .and_then(Option::as_deref)
                .map_or(Ok(()), super::kernel_adapter::OperandGeometryCheck::check)
        };
        let keys = |object| super::kernel_adapter::query_key(namespace, slot == 1, object);
        let links = self.owner_links;
        let owner_of = |owner| links.get(&owner).copied();
        let handle = |entity| {
            Self::entity_h(brep.b, entity)
                .ok_or_else(|| invalid("classification entity is outside this B-rep"))
        };
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        let cache = self
            .point_caches
            .get_mut(slot)
            .ok_or_else(|| invalid("bad B-rep cache slot"))?;
        let mut objects = crate::query::brep::BrepPointObjects::with_handles(
            b,
            BrepId::from_index(namespace * 2 + u32::from(slot == 1)),
            cache,
            &keys,
            &owner_of,
            &handle,
        )
        .with_geometry_check(&check);
        crate::query::brep::point_3d_classify(
            if slot == 1 {
                self.other_kernel.unwrap_or(k)
            } else {
                k
            },
            &mut objects,
            point,
            point_zone_tolerance,
            boundary_intersections,
            classification,
        )
    }

    fn edge_interval(&self, edge: H) -> Option<Extent1d> {
        ed::interval(self.brep(edge)?, eid(edge)).ok()
    }
    fn edge_curve(&self, edge: H) -> Option<Option<H>> {
        self.edge_curve_h(edge)
    }
    fn topology_tolerance(&self, topo: H) -> Option<f64> {
        let b = self.brep(topo)?;
        match topo.kind {
            K_VERTEX => vx::tolerance(b, vid(topo)).ok(),
            K_EDGE => ed::tolerance(b, eid(topo)).ok(),
            K_FACE => fp2::tolerance(b, fid(topo)).ok(),
            _ => None,
        }
    }
    fn face_uv_domain(&self, face: H) -> Option<Extent2d> {
        fp2::uv_domain(self.brep(face)?, fid(face)).ok()
    }
    fn set_face_uv_domain(&mut self, face: H, domain: Extent2d) -> Option<()> {
        fp2::set_uv_domain(self.brep_mut(face)?, fid(face), &domain).ok()
    }
    fn face_surface(&self, face: H) -> Option<Option<H>> {
        self.face_surface_h(face)
    }
    fn edge_edgeuse_of_face(&self, edge: H, face: H) -> crate::Result<Option<H>> {
        let b = self.brep_r(edge)?;
        let face = (face.b == edge.b && face.kind == K_FACE).then(|| fid(face));
        Ok(ed::edgeuse_of_face(b, eid(edge), face)?.map(|eu| hd(edge.b, K_EDGEUSE, eu.index())))
    }
    fn edge_is_connected_to_face(&self, edge: H, face: H) -> crate::Result<bool> {
        let b = self.brep_r(edge)?;
        let face = (face.b == edge.b && face.kind == K_FACE).then(|| fid(face));
        ed::is_connected_to_face(b, eid(edge), face)
    }
    fn insert_curve(&mut self, curve: SmCurve) -> Option<H> {
        let i = u32::try_from(self.extra.len()).ok()?;
        let hh = curve_h(2, K_EXTRA, i, &curve);
        self.extra.push(curve);
        Some(hh)
    }
    fn vertex_common_edges(&self, v1: H, v2: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(v1)?;
        let other = (v2.b == v1.b && v2.kind == K_VERTEX).then(|| vid(v2));
        Ok(vx::common_edges(b, vid(v1), other, None)?
            .into_iter()
            .map(|e| hd(v1.b, K_EDGE, e.index()))
            .collect())
    }
    fn face_curve_on_classify(
        &mut self,
        k: &Kernel,
        face: H,
        do_point_classify: bool,
        cc: &mut CurveClassification<H>,
    ) -> crate::Result<bool> {
        let mut objs = Objs {
            b: face.b,
            two: self,
        };
        curve_on_classify(
            k,
            fid(face),
            &mut objs,
            do_point_classify,
            cc,
            false,
            false,
            true,
        )
    }
}

/// The B-rep walks over one operand's own entities (as `Brep`'s own impl, through the
/// handles); an unknown handle is an invalid-input error, as a null dereference upstream.
impl ProtoTopology<H> for ArenaPair<'_> {
    fn vertex_use_faces(&self, vertex: H) -> crate::Result<Vec<Option<H>>> {
        let b = self.brep_r(vertex)?;
        let mut out = Vec::new();
        for vu in vx::vertexuses(b, vid(vertex))? {
            let face = match crate::topology::vertexuse::faceuse(b, vu)? {
                Some(fu) => eu1::faceuse_face(b, fu)?,
                None => None,
            };
            out.push(face.map(|f| hd(vertex.b, K_FACE, f.index())));
        }
        Ok(out)
    }
    fn edge_use_faces(&self, edge: H) -> crate::Result<Vec<Option<H>>> {
        let b = self.brep_r(edge)?;
        let mut out = Vec::new();
        for eu in ed::edgeuses(b, eid(edge), None, None)? {
            out.push(eu1::face(b, eu)?.map(|f| hd(edge.b, K_FACE, f.index())));
        }
        Ok(out)
    }
    fn edge_start_vertex(&self, edge: H) -> crate::Result<Option<H>> {
        Ok(ed::start_vertex(self.brep_r(edge)?, eid(edge))?
            .map(|v| hd(edge.b, K_VERTEX, v.index())))
    }
    fn edge_other_vertex(&self, edge: H, vertex: Option<H>) -> crate::Result<Option<H>> {
        let b = self.brep_r(edge)?;
        let v = vertex
            .filter(|v| v.b == edge.b && v.kind == K_VERTEX)
            .map(vid);
        Ok(ed::other_vertex(b, eid(edge), v)?.map(|v| hd(edge.b, K_VERTEX, v.index())))
    }
    fn vertex_point(&self, vertex: H) -> crate::Result<Point3> {
        vx::point(self.brep_r(vertex)?, vid(vertex))
    }
    fn edge_brep(&self, edge: H) -> crate::Result<Option<ObjectKey>> {
        Ok(eu1::edge_brep(self.brep_r(edge)?, eid(edge))?.map(Brep::key))
    }
    fn edge_faces(&self, edge: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(edge)?;
        Ok(ed::faces(b, eid(edge), None)?
            .into_iter()
            .map(|f| hd(edge.b, K_FACE, f.index()))
            .collect())
    }
}

impl crate::topology::proto_topology::p3::BrepGlue<H> for ArenaPair<'_> {
    fn glue_vertices(
        &mut self,
        k: &mut Kernel,
        brep: H,
        v1: H,
        v2: H,
        opt_wire_edges_in_brep: Option<&[H]>,
    ) -> crate::Result<()> {
        if brep.kind != K_BREP
            || v1.kind != K_VERTEX
            || v2.kind != K_VERTEX
            || v1.b != brep.b
            || v2.b != brep.b
        {
            return Err(invalid("glue vertex owner mismatch"));
        }
        // Native Glue returns before inspecting the optional wire list or
        // notification inputs when both vertex identities are the same.
        if v1 == v2 {
            self.brep_r(brep)?;
            return Ok(());
        }
        let wires = opt_wire_edges_in_brep
            .map(|edges| {
                edges
                    .iter()
                    .map(|&e| {
                        if e.kind == K_EDGE && e.b == brep.b {
                            Ok(eid(e))
                        } else {
                            Err(invalid("glue wire edge owner mismatch"))
                        }
                    })
                    .collect::<crate::Result<Vec<_>>>()
            })
            .transpose()?;
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        self.require_notifications(slot)?;
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("absent B-rep"))?;
        let namespace = self.namespace;
        let key_of =
            |entity| super::kernel_adapter::notification_key(namespace, brep.b == 1, entity);
        let mut env = ed::NotifyEnv {
            records: self.records,
            attributes: self.attributes,
            context: Some(if brep.b == 1 {
                self.other_context
                    .as_deref_mut()
                    .unwrap_or(self.primary_context)
            } else {
                self.primary_context
            }),
            key_of: &key_of,
        };
        crate::topology::brep::p5::glue_with_notifications(
            k,
            b,
            Some(vid(v1)),
            Some(vid(v2)),
            wires.as_deref(),
            Some(&mut env),
        )
    }
}

impl IntersectorBrep<H> for ArenaPair<'_> {
    fn brep_arena(&self, brep: H) -> crate::Result<&Brep> {
        if brep.kind != K_BREP {
            return Err(invalid("not a B-rep handle"));
        }
        self.brep_r(brep)
    }

    fn face_bounding_box(
        &mut self,
        k: &Kernel,
        face: H,
        output: &mut Extent3d,
    ) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad face B-rep"))?;
        self.require_query_geometry(slot)?;
        let kernel = if slot == 1 {
            self.other_kernel.unwrap_or(k)
        } else {
            k
        };
        *output = crate::query::face::p1::calculate_loose_bounding_box(
            kernel,
            self.brep_r(face)?,
            fid(face),
        )?;
        Ok(())
    }

    fn face_point_3d_classify(
        &mut self,
        k: &Kernel,
        face: H,
        point: Point3,
        src_zone_tol_3d: f64,
        do_boundary_intersections: bool,
        pc: &mut PointClassification<H>,
        guess_uv: Option<Point2>,
    ) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_query_geometry(slot)?;
        let mut local = pc.try_map_objects(|h| self.classification_entity(face.b, h))?;
        let result = (|| {
            // Keep the canonical Rust guard for an absent surface. This branch
            // stops after its one domain read and PC initialisation, before any
            // cache access; valid surfaces use the live context below.
            if matches!(self.face_surface_h(face), Some(None)) {
                let b = self
                    .breps
                    .get_mut(slot)
                    .and_then(Option::as_deref_mut)
                    .ok_or_else(|| invalid("handle names an absent B-rep"))?;
                return crate::query::face::p1::point_3d_classify(
                    k,
                    b,
                    fid(face),
                    point,
                    src_zone_tol_3d,
                    do_boundary_intersections,
                    &mut local,
                    crate::query::face::p1::Point3dClassifyOptions {
                        guess_uv,
                        ..crate::query::face::p1::Point3dClassifyOptions::default()
                    },
                );
            }
            let Some(uv) = crate::intersect::topology_solver::p1::face_point_3d_uv_with_context(
                k,
                self,
                face,
                fid(face),
                point,
                src_zone_tol_3d,
                &mut local,
                guess_uv,
            )?
            else {
                return Ok(());
            };
            // Only lazy UV trim data changes. Keep both published spatial caches
            // alive; brep_m is reserved for topology/3D geometry edits.
            let b = self
                .breps
                .get_mut(slot)
                .and_then(Option::as_deref_mut)
                .ok_or_else(|| invalid("handle names an absent B-rep"))?;
            crate::query::face::p1::point_3d_classify_finish(
                k,
                b,
                fid(face),
                uv,
                point,
                src_zone_tol_3d,
                do_boundary_intersections,
                &mut local,
                crate::query::face::p1::Point3dClassifyOptions {
                    guess_uv,
                    ..crate::query::face::p1::Point3dClassifyOptions::default()
                },
            )
        })();
        // Copy every field back even when the native call returns an error:
        // its caller deliberately ignores ordinary statuses during refreshes.
        *pc = local.try_map_objects(|entity| {
            Self::entity_h(face.b, entity)
                .ok_or_else(|| invalid("classification entity is outside this B-rep"))
        })?;
        result
    }
    fn make_vertex_loop(&mut self, brep: H, face: H, point: Point3) -> crate::Result<(H, H)> {
        if brep.kind != K_BREP || face.kind != K_FACE || brep.b != face.b {
            return Err(invalid("vertex loop owner mismatch"));
        }
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_notifications(slot)?;
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("absent B-rep"))?;

        let mut env = crate::topology::edge::NotifyEnv {
            records: self.records,
            attributes: self.attributes,
            context: Some(if brep.b == 1 {
                self.other_context
                    .as_deref_mut()
                    .unwrap_or(self.primary_context)
            } else {
                self.primary_context
            }),
            key_of: &|e| super::kernel_adapter::notification_key(self.namespace, brep.b == 1, e),
        };
        let (loop_, vertex) = crate::topology::brep::p9::make_vertex_loop(
            b,
            fid(face),
            point,
            OrientType::Opposite,
            &mut env,
        )?;
        Ok((
            hd(brep.b, K_LOOP, loop_.index()),
            hd(brep.b, K_VERTEX, vertex.index()),
        ))
    }
    fn make_vertex_split_edge(
        &mut self,
        k: &mut Kernel,
        brep: H,
        edge: H,
        parameter: f64,
    ) -> crate::Result<(H, H, H)> {
        if brep.kind != K_BREP || edge.kind != K_EDGE || brep.b != edge.b {
            return Err(invalid("split edge owner mismatch"));
        }
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_notifications(slot)?;
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("absent B-rep"))?;

        let mut env = crate::topology::edge::NotifyEnv {
            records: self.records,
            attributes: self.attributes,
            context: Some(if brep.b == 1 {
                self.other_context
                    .as_deref_mut()
                    .unwrap_or(self.primary_context)
            } else {
                self.primary_context
            }),
            key_of: &|e| super::kernel_adapter::notification_key(self.namespace, brep.b == 1, e),
        };
        let out = crate::topology::brep::p9::make_vertex_split_edge(
            k,
            b,
            eid(edge),
            parameter,
            &mut env,
            &|e, kind| ObjectRef {
                key: super::kernel_adapter::query_key(
                    self.namespace,
                    brep.b == 1,
                    CacheObject::Curve(e),
                )
                .key,
                type_code: kind,
            },
        )?;
        if !out.split {
            return Err(SolidError::sm_err("boolean::merge::p4::arena_pair::split"));
        }
        let [a, z] = out
            .edges
            .ok_or_else(|| invalid("successful split without edges"))?;
        let v = out
            .vertex
            .ok_or_else(|| invalid("successful split without vertex"))?;
        Ok((
            hd(brep.b, K_EDGE, a.index()),
            hd(brep.b, K_EDGE, z.index()),
            hd(brep.b, K_VERTEX, v.index()),
        ))
    }

    fn brep_edges(&self, brep: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(brep)?;
        Ok(bp1::edges(b, None)?
            .into_iter()
            .map(|e| hd(brep.b, K_EDGE, e.index()))
            .collect())
    }
    fn brep_vertices(&self, brep: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(brep)?;
        Ok(bp1::vertices(b, None)?
            .into_iter()
            .map(|v| hd(brep.b, K_VERTEX, v.index()))
            .collect())
    }
    fn vertex_brep(&self, vertex: H) -> crate::Result<Option<H>> {
        let b = self.brep_r(vertex)?;
        Ok(vx::brep(b, vid(vertex))?.map(|_| brep_h(vertex.b)))
    }
    fn edge_brep_handle(&self, edge: H) -> crate::Result<Option<H>> {
        let b = self.brep_r(edge)?;
        Ok(ed::brep(b, eid(edge))?.map(|_| brep_h(edge.b)))
    }
    fn set_vertex_point(&mut self, _k: &Kernel, vertex: H, point: Point3) -> crate::Result<()> {
        let slot = usize::try_from(vertex.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_notifications(slot)?;
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("absent B-rep"))?;
        let this = super::kernel_adapter::notification_key(
            self.namespace,
            vertex.b == 1,
            EntityRef::Vertex(vid(vertex)),
        );
        vx::set_point(
            b,
            vid(vertex),
            point,
            self.records,
            self.attributes,
            Some(if vertex.b == 1 {
                self.other_context
                    .as_deref_mut()
                    .unwrap_or(self.primary_context)
            } else {
                self.primary_context
            }),
            this,
            &|vu| {
                super::kernel_adapter::notification_key(
                    self.namespace,
                    vertex.b == 1,
                    EntityRef::Vertexuse(vu),
                )
            },
        )
    }
    fn set_topology_tolerance(
        &mut self,
        _k: &Kernel,
        topo: H,
        tolerance: f64,
        update_only_if_larger: bool,
        cascade_to_boundaries: bool,
    ) -> crate::Result<()> {
        let b = self.brep_m(topo)?;
        match topo.kind {
            K_VERTEX => vx::set_tolerance(
                b,
                vid(topo),
                tolerance,
                update_only_if_larger,
                cascade_to_boundaries,
            ),
            K_EDGE => ed::set_tolerance(
                b,
                eid(topo),
                tolerance,
                update_only_if_larger,
                cascade_to_boundaries,
            ),
            _ => Err(invalid("SetTolerance of a non-vertex, non-edge")),
        }
    }
    fn update_object_tolerance(
        &mut self,
        k: &Kernel,
        topo: H,
        this_gap: f64,
        this_tol: f64,
        check_gaps: bool,
    ) -> crate::Result<()> {
        let target = match topo.kind {
            K_VERTEX => TopologyToUpdate::Vertex(vid(topo)),
            K_EDGE => TopologyToUpdate::Edge(eid(topo)),
            K_FACE => TopologyToUpdate::Face(fid(topo)),
            _ => return Err(invalid("UpdateObjectTolerance of a non-topology object")),
        };
        let b = self.brep_m(topo)?;
        update_object_tolerance(k, b, Some(target), this_gap, this_tol, check_gaps)
    }
    fn vertex_is_shell_vertex(&self, vertex: H) -> crate::Result<bool> {
        vx::is_shell_vertex(self.brep_r(vertex)?, vid(vertex))
    }
    fn edge_is_wire(&self, edge: H) -> crate::Result<bool> {
        ed::is_wire(self.brep_r(edge)?, eid(edge))
    }
    fn edge_is_lamina(&self, edge: H) -> crate::Result<bool> {
        ed::is_lamina(self.brep_r(edge)?, eid(edge))
    }
    fn edge_vertex(&self, edge: H) -> crate::Result<Option<H>> {
        Ok(ed::vertex(self.brep_r(edge)?, eid(edge))?.map(|v| hd(edge.b, K_VERTEX, v.index())))
    }
    fn edge_primary_edgeuse_orientation(&self, edge: H) -> crate::Result<Option<OrientType>> {
        let b = self.brep_r(edge)?;
        match eu1::primary_edgeuse(b, eid(edge))? {
            Some(eu) => Ok(Some(eu1::orientation(b, eu)?)),
            None => Ok(None),
        }
    }
    fn edge_bounding_box(&self, k: &Kernel, edge: H) -> crate::Result<Extent3d> {
        let b = self.brep_r(edge)?;
        let boxes = crate::query::edge::calculate_bounding_box(
            k,
            b,
            eid(edge),
            BoundingBoxRequest {
                normal: true,
                ..BoundingBoxRequest::default()
            },
            false,
        )?;
        boxes.normal.ok_or_else(|| invalid("no normal box"))
    }
    fn face_uv_domain_from_uv_trim_curves(
        &mut self,
        k: &Kernel,
        face: H,
    ) -> crate::Result<Extent2d> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let b = self.brep_m(face)?;
        Ok(
            crate::topology::face::p1::calculate_uv_domain_from_uv_trim_curves(
                k,
                b,
                fid(face),
                true,
                false,
                false,
            )?
            .domain,
        )
    }
    fn face_point_classify(
        &mut self,
        k: &Kernel,
        face: H,
        uv: Point2,
        pc: &mut PointClassification<H>,
    ) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad face B-rep"))?;
        self.require_query_geometry(slot)?;
        let mut local = pc.try_map_objects(|h| self.classification_entity(face.b, h))?;
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        let result = point_classify(
            k,
            b,
            fid(face),
            uv,
            &mut local,
            PointClassifyOptions::default(),
        );
        *pc = local.try_map_objects(|entity| {
            Self::entity_h(face.b, entity)
                .ok_or_else(|| invalid("classification entity is outside this B-rep"))
        })?;
        result
    }
    fn face_intersect(
        &mut self,
        k: &Kernel,
        face1: H,
        face2: H,
        approx_tol_3d: f64,
        ang_tol_deg: f64,
        opt_coincidence_flag: Option<bool>,
    ) -> crate::Result<ti::FaceIntersection<H>> {
        if face1.kind != K_FACE || face2.kind != K_FACE || face1.b != face2.b {
            return Err(invalid("FaceIntersect needs two faces of one B-rep"));
        }
        let options = FaceIntersectOptions {
            want_this_classes: true,
            want_other_classes: true,
            want_cc_curves: true,
            coincidence: opt_coincidence_flag,
            ..FaceIntersectOptions::default()
        };
        let namespace = self.namespace;
        let other_operand = face1.b == 1;
        let owner_links = self.owner_links;
        let keys = |object| super::kernel_adapter::query_key(namespace, other_operand, object);
        let slot = usize::try_from(face1.b).map_err(|_| invalid("bad B-rep"))?;
        // FaceIntersect may create lazy UV trims; spatial topology/3D geometry
        // is unchanged. Keep both published spatial caches on errors as well.
        let brep = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        let out = fp3::face_intersect_with_context(
            k,
            brep,
            fid(face1),
            fid(face2),
            approx_tol_3d,
            ang_tol_deg,
            &options,
            &keys,
            &|object| owner_links.get(&object).copied(),
        )?;
        let bi = face1.b;

        // every curve the call made gets a handle in the curve store
        let mut curves_3d = Vec::new();
        for c in out.curves_3d {
            curves_3d.push(self.register_curve(c)?);
        }
        let cc3d: Vec<Option<H>> = self.register_curves(
            out.cc_curves_3d
                .map(|v| v.into_iter().map(Some).collect())
                .unwrap_or_default(),
        )?;
        let uv1 = self.register_curves(out.cc_this_uv_curves.unwrap_or_default())?;
        let uv2 = self.register_curves(out.cc_other_uv_curves.unwrap_or_default())?;

        // the classifications name their objects by handle
        let mut map = |o: IntersectObject| -> crate::Result<H> {
            let bad = || invalid("classification object without a handle");
            match o {
                IntersectObject::Curve(i) => cc3d.get(i).copied().flatten().ok_or_else(bad),
                IntersectObject::UvCurve(i, false) => uv1.get(i).copied().flatten().ok_or_else(bad),
                IntersectObject::UvCurve(i, true) => uv2.get(i).copied().flatten().ok_or_else(bad),
                IntersectObject::Surface(f) => Ok(hd(bi, K_FACE_SURFACE, f.index())),
                IntersectObject::EdgeCurve(e) => self
                    .edge_curve_h(hd(bi, K_EDGE, e.index()))
                    .flatten()
                    .ok_or_else(bad),
                IntersectObject::Trim(eu) => Ok(hd(bi, K_TRIM_CURVE, eu.index())),
                IntersectObject::Topology(e) => match e {
                    EntityRef::Vertexuse(vu) => Ok(hd(bi, K_VERTEXUSE, vu.index())),
                    other => Self::entity_h(bi, other).ok_or_else(bad),
                },
            }
        };
        let mut cc1 = Vec::new();
        for c in out.this_classes.unwrap_or_default() {
            cc1.push(c.map_objects(&mut map)?);
        }
        let mut cc2 = Vec::new();
        for c in out.other_classes.unwrap_or_default() {
            cc2.push(c.map_objects(&mut map)?);
        }
        Ok(ti::FaceIntersection {
            curves_3d,
            cc1,
            cc2,
            cc_curves: [
                cc3d.into_iter().flatten().collect(),
                uv1.into_iter().flatten().collect(),
                uv2.into_iter().flatten().collect(),
            ],
        })
    }
    fn editing_enabled(&self, brep: H) -> crate::Result<bool> {
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        if brep.kind != K_BREP {
            return Err(invalid("not a B-rep handle"));
        }
        self.breps
            .get(slot)
            .and_then(Option::as_ref)
            .map(|b| b.editing_enabled())
            .ok_or_else(|| invalid("not a B-rep handle"))
    }
    fn brep_tolerance(&self, brep: H) -> crate::Result<f64> {
        Ok(self.brep_r(brep)?.tolerance())
    }
    fn set_editing_enabled(&mut self, brep: H, enabled: bool) -> crate::Result<bool> {
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        if brep.kind != K_BREP {
            return Err(invalid("not a B-rep handle"));
        }
        let b = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("not a B-rep handle"))?;
        let old = b.editing_enabled();
        b.set_editing_enabled(enabled);
        Ok(old)
    }
    fn topology_objects(&mut self) -> &mut dyn TopologyObjects<H> {
        self
    }
    fn brep_glue_edges_general(
        &mut self,
        k: &mut Kernel,
        brep: H,
        surviving: Option<H>,
        relative_orientation: OrientType,
        edge_gap: f64,
        do_region_nesting: bool,
        to_delete: Option<H>,
        inside_faceuses: Option<&[H]>,
    ) -> crate::Result<()> {
        if brep.kind != K_BREP {
            return Err(invalid("not a B-rep handle"));
        }
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        let edge = |h: Option<H>| -> crate::Result<Option<EdgeId>> {
            match h {
                None => Ok(None),
                Some(h) if h.kind == K_EDGE && h.b == brep.b => Ok(Some(eid(h))),
                Some(_) => Err(invalid("not an edge of this B-rep")),
            }
        };
        let (surviving, to_delete) = (edge(surviving)?, edge(to_delete)?);
        // (no face-use handles exist here: the only caller passes NULL)
        let inside: Option<&[FaceuseId]> = match inside_faceuses {
            None => None,
            Some([]) => Some(&[]),
            Some(_) => return Err(invalid("face-use handles are not modelled")),
        };
        self.require_notifications(slot)?;
        for cache in &mut self.caches {
            let _ = cache.take();
        }
        *self.point_caches = [None, None];
        let namespace = self.namespace;
        let keys = |object| super::kernel_adapter::query_key(namespace, slot == 1, object);
        let links = self.owner_links;
        let owner_of = |owner| links.get(&owner).copied();
        let breps = &mut self.breps;
        let b = breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;

        let mut notify = ed::NotifyEnv {
            records: self.records,
            attributes: self.attributes,
            context: Some(if brep.b == 1 {
                self.other_context
                    .as_deref_mut()
                    .unwrap_or(self.primary_context)
            } else {
                self.primary_context
            }),
            key_of: &|e| super::kernel_adapter::notification_key(self.namespace, brep.b == 1, e),
        };
        let queries = self
            .region_queries
            .get_mut(slot)
            .ok_or_else(|| invalid("bad region query operand"))?;
        p7::glue_edges_general(
            k,
            b,
            surviving,
            relative_orientation,
            edge_gap,
            do_region_nesting,
            to_delete,
            inside,
            &mut notify,
            queries,
            crate::intersect::curve_class::p7::FaceQueryIdentity {
                keys: &keys,
                owner_of: &owner_of,
            },
        )
    }
}

impl TopologyObjects<H> for ArenaPair<'_> {
    fn with_brep_trees(
        &mut self,
        k: &Kernel,
        first: H,
        second: Option<H>,
        solve: &mut crate::intersect::topology_solver::p1::BrepTreeSolve<'_, H>,
    ) -> crate::Result<SolutionArray<H>> {
        fn trees(
            cache: &BrepCache<H>,
        ) -> crate::Result<crate::intersect::topology_solver::p2::TopologyTrees<'_, H>> {
            Ok(crate::intersect::topology_solver::p2::TopologyTrees {
                vertices: cache
                    .vertex_tree()
                    .ok_or_else(|| invalid("missing vertex tree"))?,
                curves: cache
                    .curve_tree()
                    .ok_or_else(|| invalid("missing curve tree"))?,
                surfaces: cache
                    .surface_tree()
                    .ok_or_else(|| invalid("missing surface tree"))?,
            })
        }
        let cache = |world: &Self, handle: H| -> crate::Result<std::sync::Arc<BrepCache<H>>> {
            world.brep_cache(k, handle)?;
            let slot = usize::try_from(handle.b).map_err(|_| invalid("bad B-rep"))?;
            world
                .caches
                .get(slot)
                .and_then(std::cell::OnceCell::get)
                .cloned()
                .ok_or_else(|| invalid("B-rep cache not stored"))
        };
        let first_cache = cache(self, first)?;
        let first_trees = trees(&first_cache)?;
        let second_cache = second.map(|handle| cache(self, handle)).transpose()?;
        let second_trees = second_cache.as_deref().map(trees).transpose()?;
        // Owning only Arc handles keeps the same caches published while the
        // classifiers edit UV trim state in the original B-reps.
        solve(self, first_trees, second_trees)
    }
    fn object_zone_tol_3d(&self, k: &Kernel, object: Option<H>) -> crate::Result<f64> {
        let system =
            || crate::geom::tol::zone_tol_3d(k, crate::geom::tol::SmallTopoType::NormalSize);
        let Some(mut handle) = object else {
            return Ok(system());
        };
        if handle.kind == K_CURVE {
            let Some(owner) = self.curve_edge(handle)? else {
                return Ok(system());
            };
            handle = owner;
        }
        let mapped = match handle.kind {
            K_VERTEX => CacheObject::Vertex(vid(handle)),
            K_EDGE => CacheObject::Edge(eid(handle)),
            K_FACE => CacheObject::Face(fid(handle)),
            K_EDGE_CURVE => CacheObject::Curve(eid(handle)),
            K_FACE_SURFACE => CacheObject::Surface(fid(handle)),
            K_BREP => CacheObject::Brep,
            _ => return Err(invalid("unsupported topology tolerance handle")),
        };
        let slot = usize::try_from(handle.b).map_err(|_| invalid("bad B-rep"))?;
        crate::topology::brep_cache::object_zone_tol_3d(
            k,
            self.brep_r(handle)?,
            Some(mapped),
            &|object| super::kernel_adapter::query_key(self.namespace, slot == 1, object),
            &|owner| self.owner_links.get(&owner).copied(),
        )
    }
    fn face_point_classify(
        &mut self,
        k: &Kernel,
        face: H,
        uv: Point2,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        use_uv: bool,
        classification: &mut PointClassification<EntityRef>,
    ) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_query_geometry(slot)?;
        // Creating a missing UV trim does not invalidate the 3D query trees.
        let brep = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        point_classify(
            k,
            brep,
            fid(face),
            uv,
            classification,
            PointClassifyOptions {
                src_zone_tol_3d: source_zone_tolerance,
                do_boundary_intersections: do_boundaries,
                use_uv_classification: use_uv,
                ..PointClassifyOptions::default()
            },
        )
    }
    fn face_point_3d_classify(
        &mut self,
        k: &Kernel,
        face: H,
        point: Point3,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        classification: &mut PointClassification<EntityRef>,
    ) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad B-rep"))?;
        self.require_query_geometry(slot)?;
        let Some(uv) = crate::intersect::topology_solver::p1::face_point_3d_uv_with_context(
            k,
            self,
            face,
            fid(face),
            point,
            source_zone_tolerance,
            classification,
            None,
        )?
        else {
            return Ok(());
        };
        let brep = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        crate::query::face::p1::point_3d_classify_finish(
            k,
            brep,
            fid(face),
            uv,
            point,
            source_zone_tolerance,
            do_boundaries,
            classification,
            crate::query::face::p1::Point3dClassifyOptions::default(),
        )
    }

    fn face_brep_editing_enabled(&self, face: H) -> crate::Result<bool> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let brep = self.brep_r(face)?;
        brep.face(fid(face))?;
        Ok(brep.editing_enabled())
    }

    fn set_face_brep_editing_enabled(&mut self, face: H, editing: bool) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad B-rep"))?;
        // A cache's temporary editing scope changes neither topology nor its
        // published 3D trees. Do not use brep_mut, which invalidates those trees.
        let brep = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        brep.face(fid(face))?;
        brep.set_editing_enabled(editing);
        Ok(())
    }

    fn face_create_uv_trim_curves(&mut self, k: &Kernel, face: H) -> crate::Result<()> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        let slot = usize::try_from(face.b).map_err(|_| invalid("bad B-rep"))?;
        let brep = self
            .breps
            .get_mut(slot)
            .and_then(Option::as_deref_mut)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        crate::topology::face::p3::create_uv_trim_curves(k, brep, fid(face), true, None).map(|_| ())
    }

    fn face_classification_object(&self, face: H, entity: EntityRef) -> crate::Result<H> {
        if face.kind != K_FACE {
            return Err(invalid("not a face handle"));
        }
        self.brep_r(face)?.face(fid(face))?;
        let object = Self::entity_h(face.b, entity)
            .ok_or_else(|| invalid("classification entity is outside this B-rep"))?;
        self.classification_entity(face.b, object)?;
        Ok(object)
    }

    /// The B-rep's cache (`SmBrep::GetBrepCache`), built from the B-rep on first use (the
    /// C++ topology solver builds it on demand) with this adapter's handles. The build
    /// starts with `SmBrep::CalculateBoundingBox` (unit `query.brep`), so until that
    /// lands the solves needing a cache defer there.
    fn brep_cache(&self, k: &Kernel, brep: H) -> crate::Result<&BrepCache<H>> {
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        let cell = self.caches.get(slot).ok_or_else(|| invalid("bad B-rep"))?;
        if let Some(cache) = cell.get() {
            return Ok(cache.as_ref());
        }
        let b = self
            .breps
            .get(slot)
            .and_then(Option::as_deref)
            .ok_or_else(|| invalid("handle names an absent B-rep"))?;
        self.require_query_geometry(slot)?;
        let bi = brep.b;
        let mut cache = BrepCache::new(BREP_TYPE)?;
        cache.build_trees(
            if slot == 1 {
                self.other_kernel.unwrap_or(k)
            } else {
                k
            },
            CacheSource::Brep(b),
            |o| match o {
                CacheObject::Vertex(v) => hd(bi, K_VERTEX, v.index()),
                CacheObject::Edge(e) => hd(bi, K_EDGE, e.index()),
                CacheObject::Face(f) => hd(bi, K_FACE, f.index()),
                CacheObject::Curve(e) => b
                    .edge(e)
                    .ok()
                    .and_then(|x| x.curve.as_ref())
                    .map_or(hd(bi, K_EDGE_CURVE, e.index()), |c| {
                        curve_h(bi, K_EDGE_CURVE, e.index(), c)
                    }),
                CacheObject::Surface(f) => hd(bi, K_FACE_SURFACE, f.index()),
                CacheObject::Brep => brep_h(bi),
                _ => hd(bi, K_NULL, 0),
            },
            &|object| super::kernel_adapter::query_key(self.namespace, slot == 1, object),
        )?;
        let _ = cell.set(std::sync::Arc::new(cache));
        cell.get()
            .map(std::sync::Arc::as_ref)
            .ok_or_else(|| invalid("B-rep cache not stored"))
    }
    fn curve(&self, object: H) -> crate::Result<&SmCurve> {
        if object.b == 0 || object.b == 1 {
            self.require_query_geometry(
                usize::try_from(object.b).map_err(|_| invalid("bad curve owner"))?,
            )?;
        }
        ClassGeometry::curve(self, object).ok_or_else(|| invalid("not a curve handle"))
    }
    fn analytic_curve_handle(&self, _curve: H) -> Option<H> {
        None
    }
    fn vertex_point(&self, vertex: H) -> crate::Result<Point3> {
        vx::point(self.brep_r(vertex)?, vid(vertex))
    }
    fn vertex_tolerance(&self, vertex: H) -> crate::Result<f64> {
        vx::tolerance(self.brep_r(vertex)?, vid(vertex))
    }
    fn curve_edge(&self, curve: H) -> crate::Result<Option<H>> {
        self.require_query_geometry(
            usize::try_from(curve.b).map_err(|_| invalid("bad curve B-rep"))?,
        )?;
        let geometry =
            ClassGeometry::curve(self, curve).ok_or_else(|| invalid("not a curve handle"))?;
        let Some(owner) = geometry.edge(|o| self.owner_links.get(&o).copied())? else {
            return Ok(None);
        };
        let b = self.brep_r(curve)?;
        b.edges()
            .ids()
            .find(|&id| {
                super::kernel_adapter::query_key(
                    self.namespace,
                    curve.b == 1,
                    CacheObject::Edge(id),
                ) == owner
            })
            .map(|id| Some(hd(curve.b, K_EDGE, id.index())))
            .ok_or_else(|| invalid("curve owner is outside this B-rep"))
    }
    fn edge_curve(&self, edge: H) -> crate::Result<Option<H>> {
        self.edge_curve_h(edge)
            .ok_or_else(|| invalid("not an edge handle"))
    }
    fn edge_interval(&self, edge: H) -> crate::Result<Extent1d> {
        ed::interval(self.brep_r(edge)?, eid(edge))
    }
    fn surface_face(&self, surface: H) -> crate::Result<Option<H>> {
        self.require_query_geometry(
            usize::try_from(surface.b).map_err(|_| invalid("bad surface B-rep"))?,
        )?;
        let geometry =
            ClassGeometry::surface(self, surface).ok_or_else(|| invalid("not a surface handle"))?;
        let Some(owner) = geometry.face()? else {
            return Ok(None);
        };
        let b = self.brep_r(surface)?;
        b.faces()
            .ids()
            .find(|&id| {
                super::kernel_adapter::query_key(
                    self.namespace,
                    surface.b == 1,
                    CacheObject::Face(id),
                ) == owner
            })
            .map(|id| Some(hd(surface.b, K_FACE, id.index())))
            .ok_or_else(|| invalid("surface owner is outside this B-rep"))
    }
    fn face_surface(&self, face: H) -> crate::Result<H> {
        self.face_surface_h(face)
            .flatten()
            .ok_or_else(|| SolidError::null_pointer("boolean::merge::p4::arena_pair::face_surface"))
    }
    fn surface(&self, surface: H) -> crate::Result<&Surface> {
        self.require_query_geometry(
            usize::try_from(surface.b).map_err(|_| invalid("bad surface owner"))?,
        )?;
        ClassGeometry::surface(self, surface).ok_or_else(|| invalid("not a surface handle"))
    }
    fn face_uv_domain(&self, face: H) -> crate::Result<Extent2d> {
        fp2::uv_domain(self.brep_r(face)?, fid(face))
    }
    fn face_tolerance(&self, face: H) -> crate::Result<f64> {
        fp2::tolerance(self.brep_r(face)?, fid(face))
    }
    fn surface_owner(&self, surface: H) -> crate::Result<Option<H>> {
        self.surface_face(surface)
    }
    fn shape_member(&self, _object: H) -> ShapeMember {
        ShapeMember::Other
    }
}

impl IntersectorBrepP2<H> for ArenaPair<'_> {
    fn live_brep(&self, brep: H) -> crate::Result<&Brep> {
        self.brep_arena(brep)
    }

    fn check_edge_curve(&self, brep: H, edge: EdgeId) -> crate::Result<()> {
        self.brep_arena(brep)?.edge(edge)?;
        let slot = usize::try_from(brep.b).map_err(|_| invalid("bad B-rep"))?;
        self.geometry_checks
            .get(slot)
            .and_then(Option::as_deref)
            .map_or(Ok(()), |check| check.check_edge_curve(edge))
    }

    fn brep_faces(
        &self,
        brep: H,
        attribute: Option<&crate::topology::edge::AttributeFilter<'_>>,
    ) -> crate::Result<Vec<H>> {
        bp1::faces(self.brep_arena(brep)?, attribute)?
            .into_iter()
            .map(|face| Ok(hd(brep.b, K_FACE, face.index())))
            .collect()
    }

    fn vertex_edges(&self, vertex: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(vertex)?;
        Ok(vx::edges(b, vid(vertex), None)?
            .into_iter()
            .map(|e| hd(vertex.b, K_EDGE, e.index()))
            .collect())
    }
    fn vertex_faces(&self, vertex: H) -> crate::Result<Vec<H>> {
        let b = self.brep_r(vertex)?;
        Ok(vx::faces(b, vid(vertex))?
            .into_iter()
            .map(|f| hd(vertex.b, K_FACE, f.index()))
            .collect())
    }
}

pub(super) fn encode(entity: EntityRef, other: bool) -> H {
    let (kind, i) = match entity {
        EntityRef::Face(f) => (K_FACE, f.index()),
        EntityRef::Edge(e) => (K_EDGE, e.index()),
        EntityRef::Vertex(v) => (K_VERTEX, v.index()),
        EntityRef::Region(v) => (K_REGION, v.index()),
        EntityRef::Shell(v) => (K_SHELL, v.index()),
        EntityRef::Faceuse(v) => (K_FACEUSE, v.index()),
        EntityRef::Loopuse(v) => (K_LOOPUSE, v.index()),
        EntityRef::Vertexuse(v) => (K_VERTEXUSE, v.index()),
        EntityRef::Loop(v) => (K_LOOP, v.index()),
        EntityRef::Edgeuse(v) => (K_EDGEUSE, v.index()),
        EntityRef::Brep(_) => (K_BREP, u32::from(other)),
    };
    hd(i8::from(other), kind, i)
}
pub(super) fn decode(h: H) -> EntityRef {
    match h.kind {
        K_FACE => EntityRef::Face(fid(h)),
        K_EDGE => EntityRef::Edge(eid(h)),
        K_VERTEX => EntityRef::Vertex(vid(h)),
        K_REGION => EntityRef::Region(RegionId::from_index(h.i)),
        K_SHELL => EntityRef::Shell(ShellId::from_index(h.i)),
        K_FACEUSE => EntityRef::Faceuse(FaceuseId::from_index(h.i)),
        K_LOOPUSE => EntityRef::Loopuse(LoopuseId::from_index(h.i)),
        K_VERTEXUSE => EntityRef::Vertexuse(VertexuseId::from_index(h.i)),
        K_LOOP => EntityRef::Loop(LoopId::from_index(h.i)),
        K_EDGEUSE => EntityRef::Edgeuse(EdgeuseId::from_index(h.i)),
        _ => EntityRef::Brep(BrepId::from_index(h.i)),
    }
}

pub(super) fn surface_handle(other: bool, face: FaceId) -> H {
    hd(i8::from(other), K_FACE_SURFACE, face.index())
}

// A nonempty subset whose external identities match no surface must still
// exclude every actual surface. Operand 2 has no B-rep surfaces in this adapter.
pub(super) fn unmatched_surface() -> H {
    hd(2, K_FACE_SURFACE, 0)
}

pub(super) fn deleted_vertex(h: H) -> crate::Result<VertexId> {
    if h.b == 0 && h.kind == K_VERTEX {
        Ok(vid(h))
    } else {
        Err(invalid("deleted vertex is not a primary vertex"))
    }
}
pub(super) fn deleted_edge(h: H) -> crate::Result<EdgeId> {
    if h.b == 0 && h.kind == K_EDGE {
        Ok(eid(h))
    } else {
        Err(invalid("deleted edge is not a primary edge"))
    }
}

#[cfg(test)]
mod face_intersection_cache_tests {
    use super::*;

    #[test]
    fn early_face_intersection_error_retains_both_published_spatial_caches() {
        let mut primary = Brep::new();
        let shell = primary.add_shell().unwrap();
        let face = primary.add_face(shell).unwrap();
        let mut other = Brep::new();
        let mut curves = Vec::new();
        let mut extra = Vec::new();
        let mut records = BTreeMap::new();
        let mut attributes = AttributeStore::default();
        let owners = BTreeMap::new();
        let geometry = [None, None];
        let notifications = [None, None];
        let mut points = [
            Some(BrepCache::new(BREP_TYPE).unwrap()),
            Some(BrepCache::new(BREP_TYPE).unwrap()),
        ];
        let mut region_queries = [
            super::super::kernel_adapter::RegionQuerySlot::default(),
            super::super::kernel_adapter::RegionQuerySlot::default(),
        ];
        let mut context = Context::default();
        let mut pair = ArenaPair::new(
            &mut primary,
            Some(&mut other),
            &mut curves,
            &mut extra,
            &mut records,
            &mut attributes,
            &owners,
            &geometry,
            &notifications,
            None,
            &mut points,
            &mut region_queries,
            &mut context,
            None,
            42,
        );
        let caches = [
            std::sync::Arc::new(BrepCache::<H>::new(BREP_TYPE).unwrap()),
            std::sync::Arc::new(BrepCache::<H>::new(BREP_TYPE).unwrap()),
        ];
        for (slot, cache) in pair.caches.iter().zip(&caches) {
            slot.set(std::sync::Arc::clone(cache)).unwrap();
        }
        let handle = hd(0, K_FACE, face.index());
        assert!(IntersectorBrep::face_intersect(
            &mut pair,
            &Kernel::default(),
            handle,
            handle,
            1e-5,
            5.0,
            None,
        )
        .is_err());
        for (slot, cache) in pair.caches.iter().zip(&caches) {
            assert!(std::sync::Arc::ptr_eq(slot.get().unwrap(), cache));
        }
        assert!(pair.point_caches.iter().all(Option::is_some));
    }
}
