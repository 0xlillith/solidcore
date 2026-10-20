//! Classification identities and borrowed geometry for the original B-rep.
//!
//! Upstream: `source/SMLib/src/SmBrep.cpp` (`ProjectAndTrim`, `MergeCurvesOnSurface`).
//! No geometry is copied. UV trim lookups retain the real selected mate's slot.

use crate::curve::curve::Curve;
use crate::error::EntityRef;
use crate::intersect::point_class::{ClassGeometry, ClassObject, ClassObjectKind};
use crate::math::Point3;
use crate::query::face::p2::CurveOnClassifyObjects;
use crate::surface::surface::Surface;
use crate::topology::brep::p2::{MergeCurvesObjects, MergeCurvesWorld};
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, LoopId, VertexId, VertexuseId};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum CacheHandle {
    Entity(EntityRef),
    Input3d(usize),
    InputUv(usize),
    EdgeCurve(EdgeId),
    TrimCurve(EdgeuseId),
    FaceSurface(FaceId),
    Temporary(usize, ClassObjectKind),
}

impl ClassObject for CacheHandle {
    fn class_kind(&self) -> ClassObjectKind {
        match *self {
            Self::Entity(entity) => match entity {
                EntityRef::Face(_) => ClassObjectKind::Face,
                EntityRef::Edge(_) => ClassObjectKind::Edge,
                EntityRef::Vertex(_) => ClassObjectKind::Vertex,
                EntityRef::Loop(_) => ClassObjectKind::Loop,
                EntityRef::Region(_) => ClassObjectKind::Region,
                EntityRef::Faceuse(_) => ClassObjectKind::Faceuse,
                EntityRef::Edgeuse(_) => ClassObjectKind::Edgeuse,
                EntityRef::Vertexuse(_) => ClassObjectKind::Vertexuse,
                EntityRef::Loopuse(_) => ClassObjectKind::Loopuse,
                EntityRef::Brep(_) | EntityRef::Shell(_) => ClassObjectKind::Other,
            },
            Self::Input3d(_) | Self::InputUv(_) | Self::EdgeCurve(_) | Self::TrimCurve(_) => {
                ClassObjectKind::Curve
            }
            Self::FaceSurface(_) => ClassObjectKind::Surface,
            Self::Temporary(_, kind) => kind,
        }
    }
    fn temporary(index: usize, kind: ClassObjectKind) -> Option<Self> {
        Some(Self::Temporary(index, kind))
    }
}

pub(super) struct CacheMergeWorld<'a> {
    brep: &'a Brep,
    curves_3d: &'a [Curve],
    uv_curves: &'a [Curve],
}

impl ClassGeometry<CacheHandle> for CacheMergeWorld<'_> {
    fn curve(&self, object: CacheHandle) -> Option<&Curve> {
        match object {
            CacheHandle::Input3d(index) => self.curves_3d.get(index),
            CacheHandle::InputUv(index) => self.uv_curves.get(index),
            CacheHandle::EdgeCurve(edge) => self.brep.edge(edge).ok()?.curve.as_ref(),
            CacheHandle::TrimCurve(edgeuse) => {
                crate::topology::edgeuse::p1::uv_trim_curve(self.brep, edgeuse)
                    .ok()?
                    .0
            }
            _ => None,
        }
    }
    fn surface(&self, object: CacheHandle) -> Option<&Surface> {
        match object {
            CacheHandle::FaceSurface(face) => self.brep.face(face).ok()?.surface.as_ref(),
            _ => None,
        }
    }
    fn volume(&self, _object: CacheHandle) -> Option<&crate::surface::volume::Volume> {
        None
    }
    fn point(&self, _object: CacheHandle) -> Option<Point3> {
        None
    }
    fn classification_vertex_point(&self, object: CacheHandle) -> Option<Point3> {
        let CacheHandle::Entity(EntityRef::Vertex(vertex)) = object else {
            return None;
        };
        Some(self.brep.vertex(vertex).ok()?.point)
    }
    fn topology_zone_tol_3d(&self, object: CacheHandle) -> Option<f64> {
        match object {
            CacheHandle::Entity(EntityRef::Vertex(vertex)) => {
                Some(self.brep.vertex(vertex).ok()?.tolerance)
            }
            CacheHandle::Entity(EntityRef::Edge(edge)) => {
                Some(self.brep.edge(edge).ok()?.tolerance)
            }
            CacheHandle::Entity(EntityRef::Face(face)) => {
                Some(self.brep.face(face).ok()?.tolerance)
            }
            CacheHandle::Entity(EntityRef::Edgeuse(edgeuse)) => {
                let edge = self.brep.edgeuse(edgeuse).ok()?.edge?;
                Some(self.brep.edge(edge).ok()?.tolerance)
            }
            _ => None,
        }
    }
    fn topology_entity(&self, object: CacheHandle) -> Option<(&Brep, EntityRef)> {
        let CacheHandle::Entity(entity) = object else {
            return None;
        };
        Some((self.brep, entity))
    }
    fn topology_handle(&self, entity: EntityRef) -> Option<CacheHandle> {
        Some(CacheHandle::Entity(entity))
    }
    fn vertex_face_gap(
        &self,
        object: CacheHandle,
        surface: &Surface,
        face: Option<CacheHandle>,
    ) -> Option<f64> {
        let CacheHandle::Entity(EntityRef::Vertex(vertex)) = object else {
            return None;
        };
        let face = match face {
            Some(CacheHandle::Entity(EntityRef::Face(face))) => face,
            Some(_) => return None,
            None => crate::topology::brep::p1::faces(self.brep, None)
                .ok()?
                .into_iter()
                .find(|&face| {
                    self.brep
                        .face(face)
                        .ok()
                        .and_then(|f| f.surface.as_ref())
                        .is_some_and(|candidate| std::ptr::eq(candidate, surface))
                })?,
        };
        let uses = crate::topology::vertex::vertexuses_of_face(self.brep, vertex, Some(face), None)
            .ok()?;
        let edgeuse = *uses.first()?;
        if !crate::topology::vertexuse::has_vertex_face_gap(self.brep, edgeuse).ok()? {
            return Some(0.0);
        }
        let gap = self.brep.vertexuse(edgeuse).ok()?.vertex_face_gap;
        // A not-yet-computed gap cannot be replaced with zero. The existing
        // ClassGeometry contract cannot mutate its B-rep to fill that cache.
        gap.is_init().then_some(gap.gap_3d)
    }
}

impl CurveOnClassifyObjects<CacheHandle> for CacheMergeWorld<'_> {
    fn face(&self, face: FaceId) -> CacheHandle {
        CacheHandle::Entity(EntityRef::Face(face))
    }
    fn vertex(&self, vertex: VertexId) -> CacheHandle {
        CacheHandle::Entity(EntityRef::Vertex(vertex))
    }
    fn edge(&self, edge: EdgeId) -> CacheHandle {
        CacheHandle::Entity(EntityRef::Edge(edge))
    }
    fn vertexuse(&self, vertexuse: VertexuseId) -> CacheHandle {
        CacheHandle::Entity(EntityRef::Vertexuse(vertexuse))
    }
    fn edgeuse(&self, edgeuse: EdgeuseId) -> CacheHandle {
        CacheHandle::Entity(EntityRef::Edgeuse(edgeuse))
    }
    fn edge_curve(&self, edge: EdgeId) -> CacheHandle {
        CacheHandle::EdgeCurve(edge)
    }
    fn trim_curve(&self, edgeuse: EdgeuseId) -> CacheHandle {
        let selected = crate::topology::edgeuse::p1::uv_trim_curve(self.brep, edgeuse)
            .ok()
            .and_then(|(_, owner)| owner)
            .unwrap_or(edgeuse);
        CacheHandle::TrimCurve(selected)
    }
    fn face_surface(&self, face: FaceId) -> Option<CacheHandle> {
        self.brep.face(face).ok()?.surface.as_ref()?;
        Some(CacheHandle::FaceSurface(face))
    }
    fn loop_id(&self, object: CacheHandle) -> Option<LoopId> {
        match object {
            CacheHandle::Entity(EntityRef::Loop(loop_)) => Some(loop_),
            _ => None,
        }
    }
}
impl MergeCurvesObjects<CacheHandle> for CacheMergeWorld<'_> {
    fn edge_of(&self, object: CacheHandle) -> Option<EdgeId> {
        match object {
            CacheHandle::Entity(EntityRef::Edge(edge)) => Some(edge),
            _ => None,
        }
    }
    fn face_of(&self, object: CacheHandle) -> Option<FaceId> {
        match object {
            CacheHandle::Entity(EntityRef::Face(face)) => Some(face),
            _ => None,
        }
    }
}
impl MergeCurvesWorld<CacheHandle> for CacheMergeWorld<'_> {
    fn merge_brep(&self) -> &Brep {
        self.brep
    }
}

pub(super) fn merge_curves_world<'a>(
    brep: &'a Brep,
    curves_3d: &'a [Curve],
    uv_curves: &'a [Curve],
) -> (CacheMergeWorld<'a>, Vec<CacheHandle>, Vec<CacheHandle>) {
    (
        CacheMergeWorld {
            brep,
            curves_3d,
            uv_curves,
        },
        (0..curves_3d.len()).map(CacheHandle::Input3d).collect(),
        (0..uv_curves.len()).map(CacheHandle::InputUv).collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::line::Line;
    use crate::math::Vec3;

    #[test]
    fn world_resolves_original_geometry_and_keeps_handle_namespaces_distinct() {
        let mut brep = Brep::new();
        let vertex = brep.add_vertex(Vec3::new(2.0, 3.0, 5.0)).unwrap();
        let curves = [Curve::from(Line::from_points(Vec3::ZERO, Vec3::X, 3))];
        let uv = [Curve::from(Line::from_points(Vec3::ZERO, Vec3::Y, 2))];
        let (world, input, trim) = merge_curves_world(&brep, &curves, &uv);
        assert!(std::ptr::eq(world.merge_brep(), &brep));
        assert_ne!(input[0], trim[0]);
        assert!(std::ptr::eq(world.curve(input[0]).unwrap(), &curves[0]));
        assert!(std::ptr::eq(world.curve(trim[0]).unwrap(), &uv[0]));
        let handle = world.vertex(vertex);
        assert_eq!(
            world.classification_vertex_point(handle),
            Some(Vec3::new(2.0, 3.0, 5.0))
        );
        let (owner, entity) = world.topology_entity(handle).unwrap();
        assert!(std::ptr::eq(owner, &brep));
        assert_eq!(entity, EntityRef::Vertex(vertex));
        assert_ne!(
            CacheHandle::temporary(0, ClassObjectKind::Curve),
            Some(input[0])
        );
    }

    #[test]
    fn bound_entry_preserves_editing_failure_and_boolean_scope() {
        let mut brep = Brep::new();
        brep.set_editing_enabled(false);
        let (world, curves, uv) = merge_curves_world(&brep, &[], &[]);
        let mut kernel = crate::kernel::Kernel::default();
        kernel.context.set_doing_boolean(true);
        let result = crate::topology::brep::p2::merge_curves_in_world(
            &mut kernel,
            &world,
            FaceId::from_index(0),
            1e-6,
            &curves,
            Some(&uv),
            crate::topology::brep::p2::MergeCurvesOptions::default(),
        );
        assert!(matches!(result, Err(crate::SolidError::Upstream { .. })));
        assert!(kernel.context.doing_boolean());
    }
}
