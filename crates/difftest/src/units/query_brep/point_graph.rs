//! Nonempty native PRE graph comparisons of the owned B-rep point classifier.
//! Genuine downstream dependencies remain deferred; no synthetic ray results.
use super::{bbox_import, bbox_input, query_metadata::query_graph_identity, BboxInput};
use crate::{Case, HarnessError, Output};
use serde::Deserialize;
use solidcore::intersect::point_class::PointClassification;
use solidcore::math::{Point3, Vec3};
use solidcore::topology::brep_cache::BrepCache;
use solidcore::topology::{Arena, Brep, BrepId};
use solidcore::EntityRef;

#[derive(Deserialize)]
struct Input {
    graph_point_classify: BboxInput,
    point: [f64; 3],
    tolerance: f64,
    boundary: bool,
    minimum_regions: u8,
}
fn native(input: &Input, export: bool) -> Result<Vec<f64>, HarnessError> {
    let (path, flat) = bbox_input(&input.graph_point_classify)?;
    solidcore_ffi::units::query_brep::classify_graph(
        export,
        &path,
        &flat,
        &input.point,
        input.tolerance,
        input.boundary,
        input.minimum_regions,
    )
    .map_err(|error| HarnessError::new(error.to_string()))
}
pub(super) fn run_cpp(case: &Case) -> Result<Output, HarnessError> {
    native(&case.parse()?, false).map(super::record)
}
pub(super) fn run_rust(case: &Case) -> Result<Output, HarnessError> {
    let input: Input = case.parse()?;
    let (graph, metadata) = bbox_import(&native(&input, true)?)?;
    if let Some(deferred) = metadata.deferred_output_for(false, false) {
        return Ok(deferred);
    }
    if graph.vertices().is_empty() && graph.edges().is_empty() && graph.faces().is_empty() {
        return Err(HarnessError::new("nonempty point fixture has no geometry"));
    }
    let mut breps: Arena<BrepId, Brep> = Arena::default();
    let id = breps
        .insert(graph)
        .map_err(|e| HarnessError::new(e.to_string()))?;
    let graph = breps
        .get_mut(id)
        .ok_or_else(|| HarnessError::new("registered B-rep missing"))?;
    let kernel = solidcore::kernel::Kernel::default();
    let geometry = || metadata.require_geometry();
    let mut slot: Option<BrepCache> = None;
    let owner_of = |_| None;
    let mut objects = solidcore::query::brep::BrepPointObjects::new(
        graph,
        id,
        &mut slot,
        &query_graph_identity,
        &owner_of,
    )
    .with_geometry_check(&geometry);
    let mut pc = PointClassification::new(17.0);
    pc.set_gap_3d(31.0);
    pc.set_t_param(47.0, false);
    pc.set_uv_param(solidcore::math::Point2::new(53.0, 59.0));
    pc.set_pre_snap_param(61.0);
    pc.set_sss_point(Vec3::new(67.0, 71.0, 73.0));
    let [x, y, z] = input.point;
    let tolerance = solidcore::geom::tol::zone_tol3d(&kernel, input.tolerance);
    let result = solidcore::query::brep::point_3d_classify(
        &kernel,
        &mut objects,
        Point3::new(x, y, z),
        tolerance,
        input.boundary,
        &mut pc,
    );
    let status = match result {
        Ok(()) => 0.0,
        Err(error) => {
            if let Some(deferred) = Output::unported(&error) {
                return Ok(deferred);
            }
            super::status(&error)
        }
    };
    let (kind, index) = match pc.object() {
        None => (-1., -1.),
        Some(EntityRef::Region(id)) => (0., f64::from(id.index())),
        Some(EntityRef::Face(id)) => (1., f64::from(id.index())),
        Some(EntityRef::Edge(id)) => (2., f64::from(id.index())),
        Some(EntityRef::Vertex(id)) => (3., f64::from(id.index())),
        Some(_) => {
            return Err(HarnessError::new(
                "point classification outside native registry",
            ))
        }
    };
    let uv = pc.uv_param();
    let sss = pc.sss_point();
    Ok(super::record(vec![
        status,
        f64::from(pc.point_class().code()),
        kind,
        index,
        pc.deviation(),
        pc.t_param(),
        uv.x,
        uv.y,
        f64::from(u8::from(pc.are_parameters_set())),
        pc.pre_snap_param(),
        f64::from(u8::from(pc.sss_point_flag())),
        sss.x,
        sss.y,
        sss.z,
        pc.src_zone_tol_3d(),
    ]))
}
