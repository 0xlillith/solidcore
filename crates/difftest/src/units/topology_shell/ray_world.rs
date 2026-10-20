//! Original B-rep geometry and persistent marks for the real Shell/Shape path.
//! Main solver and `LocalCurveSector` limitations propagate without new ports.
use super::super::query_brep::query_metadata::QueryMetadata;
use solidcore::curve::curve::Curve;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::extent3d::Extent3d;
use solidcore::geom::object::{EdgeKey, FaceKey, ObjectRef, Shape, ShapeMember, ShapeObjects};
use solidcore::geom::solution_array::{ObjectKind, SolutionObject};
use solidcore::intersect::point_class::{ClassGeometry, ClassObject, ClassObjectKind};
use solidcore::intersect::topology_solver::p1::TopologyObjects;
use solidcore::kernel::Kernel;
use solidcore::math::Point3;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep_cache::CacheObject;
use solidcore::topology::topology::Topology;
use solidcore::topology::{Brep, BrepId, EdgeId, FaceId, VertexId};
use solidcore::SolidError;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum RayObject {
    Vertex(VertexId),
    Edge(EdgeId),
    Face(FaceId),
    Curve(EdgeId),
    Surface(FaceId),
}
impl SolutionObject for RayObject {
    fn kind(&self) -> ObjectKind {
        match self {
            Self::Vertex(_) => ObjectKind::Vertex,
            Self::Edge(_) => ObjectKind::Edge,
            Self::Face(_) => ObjectKind::Face,
            Self::Curve(_) => ObjectKind::Curve,
            Self::Surface(_) => ObjectKind::Surface,
        }
    }
}
impl ClassObject for RayObject {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Vertex(_) => ClassObjectKind::Vertex,
            Self::Edge(_) => ClassObjectKind::Edge,
            Self::Face(_) => ClassObjectKind::Face,
            Self::Curve(_) => ClassObjectKind::Curve,
            Self::Surface(_) => ClassObjectKind::Surface,
        }
    }
}
const RAY_WORLD_OP: &str = "difftest::topology_shell::RayWorld";

fn ray_context_error(reason: &'static str) -> SolidError {
    SolidError::invalid_input(RAY_WORLD_OP, reason)
}

/// One registered B-rep's native metadata and persistent shape/traversal marks.
/// Reuse with the same live B-rep and Kernel between calls; geometry is borrowed
/// only during a ray. `keys` must be the mapping used to apply its query metadata.
pub(crate) struct RayQueryContext {
    id: BrepId,
    records: BTreeMap<ObjectRef, (ShapeMember, Topology)>,
    registered_faces: usize,
    registered_edges: usize,
    registered_vertices: usize,
    metadata: QueryMetadata,
    keys: fn(CacheObject) -> ObjectRef,
}
impl RayQueryContext {
    pub(crate) fn new(
        brep: &Brep,
        id: BrepId,
        metadata: QueryMetadata,
        keys: fn(CacheObject) -> ObjectRef,
    ) -> Self {
        let mut context = Self {
            id,
            records: BTreeMap::new(),
            registered_faces: 0,
            registered_edges: 0,
            registered_vertices: 0,
            metadata,
            keys,
        };
        context.register_members(brep);
        context
    }

    fn register_members(&mut self, brep: &Brep) {
        let id = self.id;
        for face in brep.faces().ids().skip(self.registered_faces) {
            self.records
                .entry((self.keys)(CacheObject::Face(face)))
                .or_insert_with(|| (ShapeMember::Face((id, face)), Topology::new()));
        }
        for edge in brep.edges().ids().skip(self.registered_edges) {
            self.records
                .entry((self.keys)(CacheObject::Edge(edge)))
                .or_insert_with(|| (ShapeMember::Edge((id, edge)), Topology::new()));
        }
        for vertex in brep.vertices().ids().skip(self.registered_vertices) {
            self.records
                .entry((self.keys)(CacheObject::Vertex(vertex)))
                .or_insert_with(|| (ShapeMember::Vertex((id, vertex)), Topology::new()));
        }
        // These arenas append records and retain detached ids. Revisiting an
        // unchanged B-rep needs no owner-map lookups, including ray retries.
        self.registered_faces = brep.faces().len();
        self.registered_edges = brep.edges().len();
        self.registered_vertices = brep.vertices().len();
    }

    pub(crate) fn world<'a>(&'a mut self, brep: &'a mut Brep) -> RayWorld<'a> {
        // Edits can add members between calls. Keep existing marks and register
        // only the new topology records; do not reapply old geometry metadata.
        self.register_members(brep);
        RayWorld {
            brep,
            id: self.id,
            records: &mut self.records,
            metadata: &self.metadata,
            keys: self.keys,
        }
    }
}
pub(crate) struct RayWorld<'a> {
    brep: &'a mut Brep,
    id: BrepId,
    records: &'a mut BTreeMap<ObjectRef, (ShapeMember, Topology)>,
    metadata: &'a QueryMetadata,
    keys: fn(CacheObject) -> ObjectRef,
}
impl RayWorld<'_> {
    fn check_brep(&self, id: BrepId) -> solidcore::Result<()> {
        if id == self.id {
            Ok(())
        } else {
            Err(ray_context_error("unregistered B-rep"))
        }
    }

    fn owner_member(&self, owner: ObjectRef) -> solidcore::Result<ShapeMember> {
        self.records
            .get(&owner)
            .map(|record| record.0)
            .ok_or_else(|| ray_context_error("unregistered geometry owner"))
    }

    fn member_key(&self, member: ShapeMember) -> solidcore::Result<ObjectRef> {
        let object = self.shape_member_object(member)?;
        Ok((self.keys)(match object {
            RayObject::Vertex(id) => CacheObject::Vertex(id),
            RayObject::Edge(id) => CacheObject::Edge(id),
            RayObject::Face(id) => CacheObject::Face(id),
            _ => return Err(ray_context_error("geometry is not a topology member")),
        }))
    }
}

impl ShapeObjects for RayWorld<'_> {
    fn brep(&self, id: BrepId) -> solidcore::Result<&Brep> {
        self.check_brep(id)?;
        Ok(self.brep)
    }
    fn edge_owner(&self, owner: ObjectRef) -> solidcore::Result<EdgeKey> {
        match self.owner_member(owner)? {
            ShapeMember::Edge(edge) => Ok(edge),
            _ => Err(ray_context_error("geometry owner is not an edge")),
        }
    }
    fn face_owner(&self, owner: ObjectRef) -> solidcore::Result<FaceKey> {
        match self.owner_member(owner)? {
            ShapeMember::Face(face) => Ok(face),
            _ => Err(ray_context_error("geometry owner is not a face")),
        }
    }
    fn topology(&self, member: ShapeMember) -> solidcore::Result<&Topology> {
        self.records
            .get(&self.member_key(member)?)
            .map(|record| &record.1)
            .ok_or_else(|| ray_context_error("unregistered topology mark record"))
    }
    fn topology_mut(&mut self, member: ShapeMember) -> solidcore::Result<&mut Topology> {
        let key = self.member_key(member)?;
        self.records
            .get_mut(&key)
            .map(|record| &mut record.1)
            .ok_or_else(|| ray_context_error("unregistered topology mark record"))
    }
}

impl solidcore::topology::shell::ShellRayWorld<RayObject> for RayWorld<'_> {
    fn shell_ray_brep_mut(&mut self, id: BrepId) -> solidcore::Result<&mut Brep> {
        self.check_brep(id)?;
        Ok(self.brep)
    }
    fn shell_ray_bounds(&self, kernel: &Kernel, id: BrepId) -> solidcore::Result<Extent3d> {
        self.metadata.require_geometry()?;
        solidcore::query::brep::calculate_loose_bounding_box(kernel, self.brep(id)?, false, true)
    }
    fn shell_ray_shape(
        &mut self,
        kernel: &mut Kernel,
        faces: &[FaceKey],
    ) -> solidcore::Result<Shape> {
        if !faces.is_empty() {
            self.metadata.require_geometry()?;
        }
        Shape::new_with_objects(kernel, self, &[], faces, &[], &[])
    }
}
impl TopologyObjects<RayObject> for RayWorld<'_> {
    fn curve(&self, object: RayObject) -> solidcore::Result<&Curve> {
        let RayObject::Curve(id) = object else {
            return Err(ray_context_error("expected curve handle"));
        };
        self.brep
            .edge(id)?
            .curve
            .as_ref()
            .ok_or(SolidError::null_pointer(RAY_WORLD_OP))
    }
    fn analytic_curve_handle(&self, object: RayObject) -> Option<RayObject> {
        matches!(object, RayObject::Curve(_)).then_some(object)
    }
    fn vertex_point(&self, object: RayObject) -> solidcore::Result<Point3> {
        let RayObject::Vertex(id) = object else {
            return Err(ray_context_error("expected vertex handle"));
        };
        Ok(self.brep.vertex(id)?.point)
    }
    fn vertex_tolerance(&self, object: RayObject) -> solidcore::Result<f64> {
        let RayObject::Vertex(id) = object else {
            return Err(ray_context_error("expected vertex handle"));
        };
        Ok(self.brep.vertex(id)?.tolerance)
    }
    fn curve_edge(&self, object: RayObject) -> solidcore::Result<Option<RayObject>> {
        // Direct native owners were installed by apply_query_metadata. An
        // unrepresentable indirect chain is gated before Shape consumes it.
        TopologyObjects::curve(self, object)?
            .edge(&|_| None)?
            .map(|owner| self.edge_owner(owner).map(|(_, id)| RayObject::Edge(id)))
            .transpose()
    }
    fn edge_curve(&self, object: RayObject) -> solidcore::Result<Option<RayObject>> {
        let RayObject::Edge(id) = object else {
            return Err(ray_context_error("expected edge handle"));
        };
        Ok(self
            .brep
            .edge(id)?
            .curve
            .as_ref()
            .map(|_| RayObject::Curve(id)))
    }
    fn edge_interval(&self, object: RayObject) -> solidcore::Result<Extent1d> {
        let RayObject::Edge(id) = object else {
            return Err(ray_context_error("expected edge handle"));
        };
        Ok(self.brep.edge(id)?.interval)
    }
    fn surface_face(&self, object: RayObject) -> solidcore::Result<Option<RayObject>> {
        TopologyObjects::surface(self, object)?
            .face()?
            .map(|owner| self.face_owner(owner).map(|(_, id)| RayObject::Face(id)))
            .transpose()
    }
    fn face_surface(&self, object: RayObject) -> solidcore::Result<RayObject> {
        let RayObject::Face(id) = object else {
            return Err(ray_context_error("expected face handle"));
        };
        self.brep
            .face(id)?
            .surface
            .as_ref()
            .ok_or(SolidError::null_pointer(RAY_WORLD_OP))?;
        Ok(RayObject::Surface(id))
    }
    fn surface(&self, object: RayObject) -> solidcore::Result<&Surface> {
        let RayObject::Surface(id) = object else {
            return Err(ray_context_error("expected surface handle"));
        };
        self.brep
            .face(id)?
            .surface
            .as_ref()
            .ok_or(SolidError::null_pointer(RAY_WORLD_OP))
    }
    fn face_uv_domain(&self, object: RayObject) -> solidcore::Result<Extent2d> {
        let RayObject::Face(id) = object else {
            return Err(ray_context_error("expected face handle"));
        };
        Ok(self.brep.face(id)?.uv_domain)
    }
    fn face_tolerance(&self, object: RayObject) -> solidcore::Result<f64> {
        let RayObject::Face(id) = object else {
            return Err(ray_context_error("expected face handle"));
        };
        Ok(self.brep.face(id)?.tolerance)
    }
    fn surface_owner(&self, object: RayObject) -> solidcore::Result<Option<RayObject>> {
        TopologyObjects::surface(self, object)?
            .owner()
            .map(|owner| self.face_owner(owner).map(|(_, id)| RayObject::Face(id)))
            .transpose()
    }
    fn shape_member(&self, object: RayObject) -> ShapeMember {
        match object {
            RayObject::Vertex(id) => ShapeMember::Vertex((self.id, id)),
            RayObject::Edge(id) => ShapeMember::Edge((self.id, id)),
            RayObject::Face(id) => ShapeMember::Face((self.id, id)),
            RayObject::Curve(_) | RayObject::Surface(_) => ShapeMember::Other,
        }
    }
    fn shape_member_object(&self, member: ShapeMember) -> solidcore::Result<RayObject> {
        Ok(match member {
            ShapeMember::Vertex((brep, id)) => {
                self.check_brep(brep)?;
                self.brep.vertex(id)?;
                RayObject::Vertex(id)
            }
            ShapeMember::Edge((brep, id)) => {
                self.check_brep(brep)?;
                self.brep.edge(id)?;
                RayObject::Edge(id)
            }
            ShapeMember::Face((brep, id)) => {
                self.check_brep(brep)?;
                self.brep.face(id)?;
                RayObject::Face(id)
            }
            ShapeMember::Other => return Err(ray_context_error("not a shape member")),
        })
    }
}
impl ClassGeometry<RayObject> for RayWorld<'_> {
    fn curve(&self, object: RayObject) -> Option<&Curve> {
        TopologyObjects::curve(self, object).ok()
    }
    fn surface(&self, object: RayObject) -> Option<&Surface> {
        TopologyObjects::surface(self, object).ok()
    }
    fn volume(&self, _object: RayObject) -> Option<&solidcore::surface::volume::Volume> {
        None
    }
    fn point(&self, _object: RayObject) -> Option<Point3> {
        None
    }
    fn classification_vertex_point(&self, object: RayObject) -> Option<Point3> {
        self.vertex_point(object).ok()
    }
    fn topology_zone_tol_3d(&self, object: RayObject) -> Option<f64> {
        match object {
            RayObject::Vertex(id) => self.brep.vertex(id).ok().map(|v| v.tolerance),
            RayObject::Edge(id) => self.brep.edge(id).ok().map(|e| e.tolerance),
            RayObject::Face(id) => self.brep.face(id).ok().map(|f| f.tolerance),
            RayObject::Curve(_) | RayObject::Surface(_) => None,
        }
    }
}
