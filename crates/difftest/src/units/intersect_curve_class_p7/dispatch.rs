//! Real p7 dispatch fixtures; original native PRE graph and complete PC inputs.
use super::super::boolean_merge_p1::{graph, Canon};
use super::super::topology_vertexuse::{apply_query_metadata, import};
use crate::{Case, HarnessError, Output};
use solidcore::curve::{curve::Curve, line::Line};
use solidcore::error::{EntityRef, Result as SolidResult, SolidError};
use solidcore::geom::{extent1d::Extent1d, object::ObjectRef};
use solidcore::intersect::curve_class::{p7::FaceQueryOwner, CurveClassification};
use solidcore::intersect::point_class::{
    ClassGeometry, ClassObject, ClassObjectKind, PointClassType, PointClassification,
};
use solidcore::kernel::{Kernel, KernelOptions, Tolerances};
use solidcore::math::{Point3, Vec3};
use solidcore::surface::{surface::Surface, volume::Volume};
use solidcore::topology::brep_cache::{arena_identity, CacheObject};
use solidcore::topology::{ArenaId, Brep, EdgeId, FaceId};
use solidcore_ffi::units::intersect_curve_class_p7_dispatch as ffi;

fn error(e: impl ToString) -> HarnessError {
    HarnessError::new(e.to_string())
}
fn bad(reason: &'static str) -> SolidError {
    SolidError::invalid_input("difftest::p7_dispatch", reason)
}
#[derive(Clone, Copy, Debug, PartialEq)]
enum Handle {
    Topology(EntityRef),
    Curve,
    Uv,
}
impl ClassObject for Handle {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Curve | Self::Uv => ClassObjectKind::Curve,
            Self::Topology(EntityRef::Face(_)) => ClassObjectKind::Face,
            Self::Topology(EntityRef::Edge(_)) => ClassObjectKind::Edge,
            Self::Topology(EntityRef::Vertex(_)) => ClassObjectKind::Vertex,
            Self::Topology(EntityRef::Edgeuse(_)) => ClassObjectKind::Edgeuse,
            _ => ClassObjectKind::Other,
        }
    }
}
struct Provider<'a> {
    brep: &'a mut Brep,
    face: FaceId,
    curve: Curve,
    uv: Curve,
}
impl ClassGeometry<Handle> for Provider<'_> {
    fn curve(&self, h: Handle) -> Option<&Curve> {
        match h {
            Handle::Curve => Some(&self.curve),
            Handle::Uv => Some(&self.uv),
            Handle::Topology(EntityRef::Edge(id)) => self.brep.edge(id).ok()?.curve.as_ref(),
            _ => None,
        }
    }
    fn surface(&self, h: Handle) -> Option<&Surface> {
        match h {
            Handle::Topology(EntityRef::Face(id)) => self.brep.face(id).ok()?.surface.as_ref(),
            _ => None,
        }
    }
    fn volume(&self, _: Handle) -> Option<&Volume> {
        None
    }
    fn point(&self, h: Handle) -> Option<Point3> {
        match h {
            Handle::Topology(EntityRef::Vertex(id)) => self.brep.vertex(id).ok().map(|v| v.point),
            _ => None,
        }
    }
}
impl FaceQueryOwner<Handle> for Provider<'_> {
    fn classification_brep(&self) -> SolidResult<&Brep> {
        Ok(self.brep)
    }
    fn classification_brep_mut(&mut self) -> SolidResult<&mut Brep> {
        Ok(self.brep)
    }
    fn classification_entity(&self, face: Handle, h: Handle) -> SolidResult<EntityRef> {
        if face != Handle::Topology(EntityRef::Face(self.face)) {
            return Err(bad("foreign probe face"));
        }
        match h {
            Handle::Topology(e) => Ok(e),
            _ => Err(bad("probe numerical curve is not topology")),
        }
    }
    fn classification_handle(&self, face: Handle, e: EntityRef) -> SolidResult<Handle> {
        if face != Handle::Topology(EntityRef::Face(self.face)) {
            return Err(bad("foreign probe face"));
        }
        Ok(Handle::Topology(e))
    }
    fn query_on_face<T>(
        &mut self,
        k: &Kernel,
        face: Handle,
        visit: impl FnOnce(
            &Kernel,
            &mut Brep,
            FaceId,
            &dyn Fn(CacheObject) -> ObjectRef,
            &dyn Fn(ObjectRef) -> Option<ObjectRef>,
            Option<&mut dyn solidcore::query::face::p2::FacePointCache>,
        ) -> SolidResult<T>,
    ) -> SolidResult<T> {
        if face != Handle::Topology(EntityRef::Face(self.face)) {
            return Err(bad("foreign probe face"));
        }
        // Native preparation explicitly verifies no retained trim cache exists.
        visit(k, self.brep, self.face, &arena_identity, &|_| None, None)
    }
}
fn branch(case: &Case) -> Result<u8, HarnessError> {
    #[derive(serde::Deserialize)]
    struct Input {
        dispatch_probe: u8,
    }
    let i: Input = case.parse()?;
    if i.dispatch_probe > 2 {
        return Err(error("invalid dispatch probe"));
    }
    Ok(i.dispatch_probe)
}
struct Reader<'a> {
    raw: &'a [f64],
    at: usize,
}
impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [f64], HarnessError> {
        let end = self
            .at
            .checked_add(n)
            .ok_or_else(|| error("probe envelope overflow"))?;
        let data = self
            .raw
            .get(self.at..end)
            .ok_or_else(|| error("short probe envelope"))?;
        self.at = end;
        Ok(data)
    }
    fn number(&mut self) -> Result<f64, HarnessError> {
        self.take(1)?
            .first()
            .copied()
            .ok_or_else(|| error("missing number"))
    }
    fn count(&mut self) -> Result<usize, HarnessError> {
        let v = self.number()?;
        if !v.is_finite() || v.fract() != 0.0 || v < 0.0 || v > self.raw.len() as f64 {
            return Err(error("invalid probe length"));
        }
        Ok(v as usize)
    }
    fn index(&mut self) -> Result<u32, HarnessError> {
        let v = self.number()?;
        if !v.is_finite() || v.fract() != 0.0 || v < 0.0 || v > f64::from(u32::MAX) {
            return Err(error("invalid probe index"));
        }
        Ok(v as u32)
    }
    fn point(&mut self) -> Result<Vec3, HarnessError> {
        let [x, y, z] = self.take(3)? else {
            return Err(error("probe point"));
        };
        Ok(Vec3::new(*x, *y, *z))
    }
}
fn object(out: &mut Vec<f64>, h: Option<Handle>) -> Result<(), HarnessError> {
    let pair = match h {
        None => (-1.0, -1.0),
        Some(Handle::Topology(e)) => match e {
            EntityRef::Face(id) => (3.0, f64::from(id.index())),
            EntityRef::Edge(id) => (7.0, f64::from(id.index())),
            EntityRef::Vertex(id) => (9.0, f64::from(id.index())),
            EntityRef::Edgeuse(id) => (6.0, f64::from(id.index())),
            _ => return Err(error("unexpected PC topology kind")),
        },
        _ => return Err(error("PC unexpectedly names a numerical view")),
    };
    out.extend([pair.0, pair.1]);
    Ok(())
}
fn pc(out: &mut Vec<f64>, pc: &PointClassification<Handle>) -> Result<(), HarnessError> {
    out.push(f64::from(pc.point_class().code()));
    object(out, pc.object())?;
    let face = match pc.face() {
        None => -1.0,
        Some(Handle::Topology(EntityRef::Face(id))) => f64::from(id.index()),
        _ => return Err(error("invalid PC face")),
    };
    out.extend([face, pc.deviation(), pc.t_param()]);
    let uv = pc.uv_param();
    let uvw = pc.uvw_param();
    let sss = pc.sss_point();
    out.extend([
        uv.x,
        uv.y,
        uvw.x,
        uvw.y,
        uvw.z,
        f64::from(u8::from(pc.are_parameters_set())),
        pc.pre_snap_param(),
        f64::from(u8::from(pc.sss_point_flag())),
        sss.x,
        sss.y,
        sss.z,
        pc.src_zone_tol_3d(),
    ]);
    Ok(())
}
fn pcs(
    out: &mut Vec<f64>,
    ivl: &solidcore::intersect::curve_class::CurveInterval<Handle>,
) -> Result<(), HarnessError> {
    pc(out, &ivl.start)?;
    pc(out, &ivl.mid)?;
    pc(out, &ivl.end)
}
fn output(b: &Brep, mut record: Vec<f64>) -> Result<Output, HarnessError> {
    if record.len() != 57 {
        return Err(error("probe output record size"));
    }
    let canon = Canon::of(b).map_err(error)?;
    let rows = record.get_mut(3..).ok_or_else(|| error("probe PCs"))?;
    for row in rows.as_chunks_mut::<18>().0 {
        let [_, kind, index, face, ..] = row;
        *index = match *kind as i32 {
            -1 => -1.0,
            3 => canon.faces.of_slot(*index),
            7 => canon.edges.of_slot(*index),
            9 => canon.vertices.of_slot(*index),
            6 => canon.edgeuses.of_slot(*index),
            _ => return Err(error("native PC object kind")),
        };
        *face = canon.faces.of_slot(*face);
    }
    graph(b, record)
}
pub(super) fn run_cpp(case: &Case) -> Result<Output, HarnessError> {
    let raw = ffi::probe(branch(case)?, false).map_err(error)?;
    let record = raw
        .get(..57)
        .ok_or_else(|| error("short native p7 record"))?
        .to_vec();
    let (b, _) = import(raw.get(57..).ok_or_else(|| error("missing native graph"))?)?;
    output(&b, record)
}
pub(super) fn run_rust(case: &Case) -> Result<Output, HarnessError> {
    let branch = branch(case)?;
    let raw = ffi::probe(branch, true).map_err(error)?;
    let mut rd = Reader { raw: &raw, at: 0 };
    let n = rd.count()?;
    let (mut brep, _) = import(rd.take(n)?)?;
    let n = rd.count()?;
    let metadata = apply_query_metadata(&mut brep, rd.take(n)?, &arena_identity)?;
    if let Some(out) = metadata.deferred_output() {
        return Ok(out);
    }
    let face = FaceId::from_index(rd.index()?);
    let edge_number = rd.number()?;
    let edge = if branch == 0 {
        if !edge_number.is_finite()
            || edge_number.fract() != 0.0
            || edge_number < 0.0
            || edge_number > f64::from(u32::MAX)
        {
            return Err(error("invalid real edge"));
        }
        Some(EdgeId::from_index(edge_number as u32))
    } else {
        if edge_number != -1.0 {
            return Err(error("unexpected endpoint edge"));
        }
        None
    };
    brep.face(face).map_err(error)?;
    let edge_t = rd.number()?;
    let zone = rd.number()?;
    if !zone.is_finite() || zone <= 0.0 {
        return Err(error("invalid native context tolerance"));
    }
    let tolerances = Tolerances {
        zone_tol_3d: zone,
        xsect_tol_3d: rd.number()?,
        approx_tol_3d: rd.number()?,
        ang_tol_deg: rd.number()?,
        tol_gap_gain: rd.number()?,
        model_size_estimate: rd.number()?,
        large_small_size_ratio: rd.number()?,
    };
    let p = rd.point()?;
    let q = rd.point()?;
    let [u0, v0, u1, v1] = rd.take(4)? else {
        return Err(error("probe UV endpoints"));
    };
    let mut curve = Line::from_points(p, q, 3);
    let mut uv = Line::from_points(Vec3::new(*u0, *v0, 0.0), Vec3::new(*u1, *v1, 0.0), 2);
    let interval = Extent1d::from_min_max(0.0, 1.0);
    curve.edit_parameterization(&interval).map_err(error)?;
    uv.edit_parameterization(&interval).map_err(error)?;
    let kernel = Kernel::new(&KernelOptions {
        tolerances,
        ..KernelOptions::default()
    });
    let face_handle = Handle::Topology(EntityRef::Face(face));
    let mut cc = CurveClassification::new(&kernel);
    cc.set_curve(Some(Handle::Curve), false);
    cc.set_surface1(Some(face_handle));
    cc.set_uv_curve1((branch == 2).then_some(Handle::Uv));
    cc.set_src_zone_tol_3d(zone);
    cc.set_size(&kernel, 1).map_err(error)?;
    let ivl = cc.at_mut(0).map_err(error)?;
    ivl.set_interval(&interval);
    ivl.start = PointClassification::new(zone);
    ivl.mid = PointClassification::new(zone);
    ivl.end = PointClassification::new(zone);
    if let Some(edge) = edge {
        brep.edge(edge).map_err(error)?;
        ivl.start.set_class_object(
            PointClassType::Edge,
            Some(Handle::Topology(EntityRef::Edge(edge))),
        );
        ivl.start.set_t_param(edge_t, false);
    }
    let mut initial_pcs = Vec::new();
    pcs(&mut initial_pcs, ivl)?;
    let native_pcs = rd.take(54)?;
    if initial_pcs
        .iter()
        .zip(native_pcs)
        .any(|(a, b)| a != b && !(a.is_nan() && b.is_nan()))
        || rd.at != raw.len()
    {
        return Err(error(
            "complete native PRE-PC inputs differ from Rust construction",
        ));
    }
    let mut provider = Provider {
        brep: &mut brep,
        face,
        curve: Curve::from(curve),
        uv: Curve::from(uv),
    };
    let result = cc.classify_face_interval(
        &kernel,
        &mut provider,
        0,
        true,
        face_handle,
        None,
        true,
        None,
    );
    let status = match result {
        Ok(()) => 0.0,
        Err(e) => {
            if let Some(out) = Output::unported(&e) {
                return Ok(out);
            }
            super::super::intersect_curve_class_p1::status(&e)
        }
    };
    let ivl = cc.at(0).map_err(error)?;
    if status == 0.0
        && (ivl.mid.point_class() != PointClassType::Face || ivl.mid.object() != Some(face_handle))
    {
        return Err(error("p7 dispatch did not classify the real face"));
    }
    let mut record = vec![status, ivl.interval.min(), ivl.interval.max()];
    pcs(&mut record, ivl)?;
    output(&brep, record)
}
