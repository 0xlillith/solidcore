//! Construction of a shape from live topology and its caller-owned marks.
//!
//! Upstream: `source/SMLib/src/SmShape.cpp`, `SmShape::SmShape`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`Shape::new_with_objects`] | `SmShape::SmShape` |
//! | [`ShapeObjects`] | The caller's topology, context and geometry owner links |
//!
//! Deviations: topology and geometry are borrowed through typed ids and explicit
//! owner lookups instead of pointers. The caller supplies the persistent mark
//! records and kernel. Geometry tree entries retain the edge or face that stores
//! the native curve or surface. Unavailable dependencies return `Unported`;
//! ordinary ignored native statuses retain the same initialized output state.

use super::{EdgeKey, FaceKey, ObjectRef, Shape, ShapeMember, VertexKey};
use crate::error::{ignore_status, Result, SolidError};
use crate::geom::context::{with_mark_lock, MarkType, ALL_MARKS};
use crate::geom::extent3d::Extent3d;
use crate::geom::tree::{NodeGeomType, Tree};
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3};
use crate::nurbs::b_spline_curve::BoundingBoxRequest;
use crate::topology::brep::p1;
use crate::topology::topology::Topology;
use crate::topology::{Brep, BrepId};

const OP: &str = "geom::object::Shape::new_with_objects";

/// Live geometry, owner resolution and mark records used by shape construction.
///
/// Every returned B-rep and topology record must belong to the kernel passed to
/// [`Shape::new_with_objects`]. Owner keys are resolved by the application; their
/// numeric values are never treated as arena ids. Marks are the same persistent
/// `SmTopology` records used by the caller's other traversals and notifications.
/// The constructor neither creates missing records nor copies source geometry.
pub trait ShapeObjects {
    /// The B-rep registered under `id` in the caller's store.
    fn brep(&self, id: BrepId) -> Result<&Brep>;
    /// The edge referred to by the actual `SmCurve::GetOwner()` link.
    fn edge_owner(&self, owner: ObjectRef) -> Result<EdgeKey>;
    /// The face referred to by the actual surface owner link.
    fn face_owner(&self, owner: ObjectRef) -> Result<FaceKey>;
    /// The object's persistent topology base record.
    fn topology(&self, member: ShapeMember) -> Result<&Topology>;
    /// Mutable access to that same record, for `SmTopology::Mark`.
    fn topology_mut(&mut self, member: ShapeMember) -> Result<&mut Topology>;
}

fn required<T>(value: Option<T>, reason: &'static str) -> Result<T> {
    value.ok_or_else(|| SolidError::invalid_input(OP, reason))
}

fn normal_request() -> BoundingBoxRequest {
    BoundingBoxRequest {
        normal: true,
        pseudo: false,
        polar: false,
    }
}

fn add_unique<T: PartialEq>(items: &mut Vec<T>, item: T) {
    if !items.contains(&item) {
        items.push(item);
    }
}

/// Bounding-box status is ignored by the native constructor. Retain its initialized
/// output on ordinary failure, but never hide an unavailable dependency.
fn brep_box(k: &Kernel, b: &Brep) -> Result<Extent3d> {
    let mut bbox = Extent3d::new();
    ignore_status(
        crate::query::brep::calculate_loose_bounding_box(k, b, false, true)
            .map(|value| bbox = value),
    )?;
    Ok(bbox)
}

fn curve_owner(objects: &impl ShapeObjects, edge: EdgeKey) -> Result<EdgeKey> {
    let data = objects.brep(edge.0)?.edge(edge.1)?;
    let curve = required(data.curve.as_ref(), "edge without curve")?;
    objects.edge_owner(required(curve.owner(), "edge curve without owner")?)
}

fn surface_owner(objects: &impl ShapeObjects, face: FaceKey) -> Result<FaceKey> {
    let data = objects.brep(face.0)?.face(face.1)?;
    let surface = required(data.surface.as_ref(), "face without surface")?;
    objects.face_owner(required(surface.owner(), "face surface without owner")?)
}

fn edge_box(
    k: &Kernel,
    objects: &impl ShapeObjects,
    edge: EdgeKey,
    out: &mut Extent3d,
) -> Result<()> {
    // SmEdge::CalculateBoundingBox initializes pNormalBox even if a child fails.
    *out = Extent3d::new();
    ignore_status(
        crate::query::edge::calculate_bounding_box(
            k,
            objects.brep(edge.0)?,
            edge.1,
            normal_request(),
            false,
        )
        .and_then(|boxes| {
            *out = required(boxes.normal, "missing requested edge box")?;
            Ok(())
        }),
    )
}

fn face_box(
    k: &Kernel,
    objects: &impl ShapeObjects,
    face: FaceKey,
    out: &mut Extent3d,
) -> Result<()> {
    // The loose face call writes only on success; one precompute loop reuses out.
    ignore_status(
        crate::query::face::p1::calculate_loose_bounding_box(k, objects.brep(face.0)?, face.1)
            .map(|value| *out = value),
    )
}

fn curve_box(k: &Kernel, objects: &impl ShapeObjects, edge: EdgeKey) -> Result<Extent3d> {
    let owner = curve_owner(objects, edge)?;
    let interval = objects.brep(owner.0)?.edge(owner.1)?.interval;
    let data = objects.brep(edge.0)?.edge(edge.1)?;
    let curve = required(data.curve.as_ref(), "edge without curve")?;
    let mut bbox = Extent3d::new();
    ignore_status(
        curve
            .calculate_bounding_box(k, &interval, normal_request(), true)
            .and_then(|boxes| {
                bbox = required(boxes.normal, "missing requested curve box")?;
                Ok(())
            }),
    )?;
    Ok(bbox)
}

fn face_edges(objects: &impl ShapeObjects, face: FaceKey) -> Result<Vec<EdgeKey>> {
    crate::topology::topology_traverser::face_edges(objects.brep(face.0)?, face.1)
        .map(|edges| edges.into_iter().map(|edge| (face.0, edge)).collect())
}

fn face_vertices(objects: &impl ShapeObjects, face: FaceKey) -> Result<Vec<VertexKey>> {
    crate::query::face::p2::face_vertices(objects.brep(face.0)?, face.1, false, None).map(
        |vertices| {
            vertices
                .into_iter()
                .map(|vertex| (face.0, vertex))
                .collect()
        },
    )
}

fn point(objects: &impl ShapeObjects, vertex: VertexKey) -> Result<Point3> {
    Ok(objects.brep(vertex.0)?.vertex(vertex.1)?.point)
}

impl Shape {
    /// Builds `SmShape` from live B-reps and their actual owner and mark records.
    ///
    /// The caller supplies the context shared by the input topology. Whole B-reps
    /// use the effective native lists from `topology::brep::p1`; explicit lists
    /// retain their order and add each object once. The same context mark is
    /// advanced, locked and written to the caller's records as in the native
    /// constructor, and its lock is released even when a dependency returns an
    /// error. Tree entries use the original curves and surfaces through their
    /// owning edge/face keys.
    ///
    /// # Errors
    /// Unavailable callees remain `Unported`. The co-ported native face edge and
    /// vertex getters are reused. Missing live objects, owners or mark records
    /// are explicit context errors. Native ordinary box and tree
    /// errors are ignored, retaining their initialized or partial state.
    pub fn new_with_objects(
        kernel: &mut Kernel,
        objects: &mut impl ShapeObjects,
        breps: &[BrepId],
        faces: &[FaceKey],
        edges: &[EdgeKey],
        vertices: &[VertexKey],
    ) -> Result<Self> {
        let mut bbox_union = Extent3d::new();
        let mut has_topology = false;
        let mut listed_faces = Vec::new();
        let mut listed_edges = Vec::new();
        let mut listed_vertices = Vec::new();
        for &id in breps {
            let b = objects.brep(id)?;
            has_topology = b.infinite_region().is_some();
            bbox_union = bbox_union.union(&brep_box(kernel, b)?);
        }
        for &face in faces {
            has_topology = true;
            add_unique(&mut listed_faces, face);
            let mut bbox = Extent3d::new();
            let owner = surface_owner(objects, face)?;
            face_box(kernel, objects, owner, &mut bbox)?;
            bbox_union = bbox_union.union(&bbox);
            for edge in face_edges(objects, face)? {
                let owner = curve_owner(objects, edge)?;
                edge_box(kernel, objects, owner, &mut bbox)?;
                bbox_union = bbox_union.union(&bbox);
            }
        }
        for &edge in edges {
            has_topology = true;
            add_unique(&mut listed_edges, edge);
            let mut bbox = Extent3d::new();
            let owner = curve_owner(objects, edge)?;
            edge_box(kernel, objects, owner, &mut bbox)?;
            bbox_union = bbox_union.union(&bbox);
        }
        for &vertex in vertices {
            has_topology = true;
            add_unique(&mut listed_vertices, vertex);
            bbox_union = bbox_union.union(&Extent3d::from_point(point(objects, vertex)?));
        }
        if !has_topology {
            bbox_union.set_min_max(Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0))?;
        }
        let mut shape = Self {
            vertex_tree: Tree::new_spatial(&bbox_union),
            curve_tree: Tree::new_spatial(&bbox_union),
            surface_tree: Tree::new_spatial(&bbox_union),
            faces: listed_faces,
            edges: listed_edges,
            vertices: listed_vertices,
        };
        if !breps.is_empty() || !faces.is_empty() || !edges.is_empty() || !vertices.is_empty() {
            with_mark_lock(kernel, ALL_MARKS, |k, mark| {
                shape.load(k, objects, mark, breps, faces, edges, vertices)
            })??;
        }
        Ok(shape)
    }

    fn load_vertex(
        &mut self,
        k: &Kernel,
        objects: &mut impl ShapeObjects,
        mark: MarkType,
        vertex: VertexKey,
        list: bool,
    ) -> Result<()> {
        let member = ShapeMember::Vertex(vertex);
        if objects.topology(member)?.is_marked(&k.context, mark) {
            return Ok(());
        }
        if list {
            add_unique(&mut self.vertices, vertex);
        }
        objects.topology_mut(member)?.mark(Some(&k.context), mark);
        let bbox = Extent3d::from_point(point(objects, vertex)?);
        ignore_status(self.vertex_tree.add_to_spatial_tree(
            &bbox,
            Some(vertex),
            NodeGeomType::Default,
        ))
    }

    fn load_edge(
        &mut self,
        k: &Kernel,
        objects: &mut impl ShapeObjects,
        mark: MarkType,
        edge: EdgeKey,
        list: bool,
    ) -> Result<bool> {
        let member = ShapeMember::Edge(edge);
        if objects.topology(member)?.is_marked(&k.context, mark) {
            return Ok(false);
        }
        if list {
            add_unique(&mut self.edges, edge);
        }
        objects.topology_mut(member)?.mark(Some(&k.context), mark);
        let bbox = curve_box(k, objects, edge)?;
        ignore_status(self.curve_tree.add_to_spatial_tree(
            &bbox,
            Some(edge),
            NodeGeomType::Default,
        ))?;
        Ok(true)
    }

    fn load_face(
        &mut self,
        k: &Kernel,
        objects: &mut impl ShapeObjects,
        mark: MarkType,
        face: FaceKey,
        list: bool,
    ) -> Result<bool> {
        let member = ShapeMember::Face(face);
        if objects.topology(member)?.is_marked(&k.context, mark) {
            return Ok(false);
        }
        if list {
            add_unique(&mut self.faces, face);
        }
        objects.topology_mut(member)?.mark(Some(&k.context), mark);
        let owner = surface_owner(objects, face)?;
        let mut bbox = Extent3d::new();
        face_box(k, objects, owner, &mut bbox)?;
        ignore_status(self.surface_tree.add_to_spatial_tree(
            &bbox,
            Some(face),
            NodeGeomType::Default,
        ))?;
        Ok(true)
    }

    #[allow(clippy::too_many_arguments)] // the constructor's ordered input lists
    fn load(
        &mut self,
        k: &Kernel,
        objects: &mut impl ShapeObjects,
        mark: MarkType,
        breps: &[BrepId],
        faces: &[FaceKey],
        edges: &[EdgeKey],
        vertices: &[VertexKey],
    ) -> Result<()> {
        for &id in breps {
            for vertex in p1::vertices(objects.brep(id)?, None)? {
                self.load_vertex(k, objects, mark, (id, vertex), true)?;
            }
            for edge in p1::edges(objects.brep(id)?, None)? {
                self.load_edge(k, objects, mark, (id, edge), true)?;
            }
            for face in p1::faces(objects.brep(id)?, None)? {
                self.load_face(k, objects, mark, (id, face), true)?;
            }
        }
        for &face in faces {
            if !self.load_face(k, objects, mark, face, false)? {
                continue;
            }
            for edge in face_edges(objects, face)? {
                self.load_edge(k, objects, mark, edge, true)?;
            }
            for vertex in face_vertices(objects, face)? {
                self.load_vertex(k, objects, mark, vertex, true)?;
            }
        }
        for &edge in edges {
            if !self.load_edge(k, objects, mark, edge, false)? {
                continue;
            }
            for vertex in crate::topology::edge::vertices(objects.brep(edge.0)?, edge.1, None)? {
                let vertex = required(vertex, "edge vertex list contains null")?;
                self.load_vertex(k, objects, mark, (edge.0, vertex), true)?;
            }
        }
        for &vertex in vertices {
            self.load_vertex(k, objects, mark, vertex, false)?;
        }
        Ok(())
    }
}
