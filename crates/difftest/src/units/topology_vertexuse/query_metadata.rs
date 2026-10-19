//! Native query metadata for the shared real B-rep graph importers.
//!
//! Paired with `harness/cpp/shim/sc_brep_query_metadata.h`. Malformed metadata is
//! an error. An exact native state that the importer cannot represent is retained
//! in a report instead: known lists and owners are installed without replacing
//! geometry or silently following an indirect owner chain. Callers check topology
//! limitations before topology work, and geometric limitations only when that
//! geometry is consumed (for example via `BrepQueryContext::with_geometry_check`).

use std::collections::BTreeSet;

use solidcore::geom::a_object::BREP_TYPE;
use solidcore::geom::object::ObjectRef;
use solidcore::topology::brep_cache::CacheObject;
use solidcore::topology::{Brep, EdgeId, FaceId, FaceuseId, RegionId, VertexId};
use solidcore::SolidError;

use crate::{HarnessError, Output};

type Result<T> = std::result::Result<T, HarnessError>;

const EXPORT_LIMITATION: &str = "harness:analytic-geometry-export";

/// Representation limitations found while importing a real native B-rep.
///
/// Check topology with `deferred_output_for(false, false)` before using imported
/// lists. Curve and surface limitations can be checked later, when the operation
/// actually reads those kinds of geometry.
#[derive(Clone, Debug, Default)]
pub(crate) struct QueryMetadata {
    topology: Vec<String>,
    curves: Vec<String>,
    surfaces: Vec<String>,
}

impl QueryMetadata {
    /// A query guard for operations that may read all B-rep geometry.
    pub(crate) fn require_geometry(&self) -> solidcore::Result<()> {
        if self.topology.is_empty() && self.curves.is_empty() && self.surfaces.is_empty() {
            Ok(())
        } else {
            Err(SolidError::unported(
                "native B-rep query geometry export",
                EXPORT_LIMITATION,
            ))
        }
    }

    /// Explain all limitations when a geometry guard was actually reached.
    pub(crate) fn deferred_output(&self) -> Option<Output> {
        self.deferred_output_for(true, true)
    }

    /// Explain only limitations relevant to the requested geometric reads.
    ///
    /// Topology is always included. In particular, a vertices-only or approximate
    /// bounding-box path need not reject an unrepresented surface or curve owner.
    pub(crate) fn deferred_output_for(
        &self,
        read_curves: bool,
        read_surfaces: bool,
    ) -> Option<Output> {
        let mut reasons: Vec<&str> = self.topology.iter().map(String::as_str).collect();
        if read_curves {
            reasons.extend(self.curves.iter().map(String::as_str));
        }
        if read_surfaces {
            reasons.extend(self.surfaces.iter().map(String::as_str));
        }
        (!reasons.is_empty()).then(|| Output::deferred(EXPORT_LIMITATION, reasons.join("; ")))
    }
}

fn error(message: impl Into<String>) -> HarnessError {
    HarnessError::new(format!("B-rep query metadata: {}", message.into()))
}

enum NativeValue<T> {
    Known(T),
    Unrepresented,
}

struct Reader<'a> {
    raw: &'a [f64],
    at: usize,
}

impl Reader<'_> {
    fn scalar(&mut self, name: &str) -> Result<f64> {
        let value = *self
            .raw
            .get(self.at)
            .ok_or_else(|| error(format!("missing {name}")))?;
        self.at += 1;
        if !value.is_finite() {
            return Err(error(format!("non-finite {name}")));
        }
        Ok(value)
    }

    fn integer(&mut self, name: &str) -> Result<i64> {
        let value = self.scalar(name)?;
        // Metadata integers are indices or native types/flags. Keep conversion
        // inside f64's exact integer range as well as i64's range.
        if value.fract() != 0.0 || value.abs() > (1_u64 << 53) as f64 {
            return Err(error(format!("invalid integer {name}")));
        }
        Ok(value as i64)
    }

    fn count(&mut self, name: &str) -> Result<usize> {
        let count =
            usize::try_from(self.integer(name)?).map_err(|_| error(format!("negative {name}")))?;
        // Bound work and allocations by the remaining input before iterating.
        if count > self.raw.len().saturating_sub(self.at) {
            return Err(error(format!("truncated {name}")));
        }
        Ok(count)
    }

    fn exact_count(&mut self, expected: usize, name: &str) -> Result<()> {
        if self.count(name)? != expected {
            return Err(error(format!("{name} differs from the imported graph")));
        }
        Ok(())
    }

    fn registry_id(
        &mut self,
        maximum: usize,
        name: &str,
        limitations: &mut Vec<String>,
    ) -> Result<NativeValue<Option<u32>>> {
        match self.integer(name)? {
            -1 => Ok(NativeValue::Known(None)),
            -3 => {
                limitations.push(format!("{name} is outside the exported topology registry"));
                Ok(NativeValue::Unrepresented)
            }
            id if id >= 0 && (id as u64) < maximum as u64 => {
                let id = u32::try_from(id).map_err(|_| error(format!("{name} id overflow")))?;
                Ok(NativeValue::Known(Some(id)))
            }
            _ => Err(error(format!("invalid {name} registry id"))),
        }
    }

    fn list(
        &mut self,
        maximum: usize,
        full: bool,
        name: &str,
        limitations: &mut Vec<String>,
    ) -> Result<NativeValue<Vec<u32>>> {
        let count = self.count(name)?;
        let mut represented = true;
        if full && count != maximum {
            limitations.push(format!(
                "{name} contains {count} of {maximum} imported records; the current B-rep model requires a full list"
            ));
            represented = false;
        }
        let mut list = Vec::new();
        let mut seen = BTreeSet::new();
        for _ in 0..count {
            match self.registry_id(maximum, name, limitations)? {
                NativeValue::Known(Some(id)) => {
                    if !seen.insert(id) {
                        return Err(error(format!("duplicate id in {name}")));
                    }
                    list.push(id);
                }
                NativeValue::Known(None) => return Err(error(format!("null member in {name}"))),
                NativeValue::Unrepresented => represented = false,
            }
        }
        Ok(if represented {
            NativeValue::Known(list)
        } else {
            NativeValue::Unrepresented
        })
    }

    fn geometry(
        &mut self,
        imported_type: i64,
        expected_owner: (i64, usize),
        owner_key: impl Fn(u32) -> ObjectRef,
        name: &str,
        limitations: &mut Vec<String>,
    ) -> Result<NativeValue<Option<ObjectRef>>> {
        let (owner_type, owner_count) = expected_owner;
        let native_type = self.integer("concrete geometry type")?;
        let native_owner_type = self.integer("direct owner type")?;
        let owner_id = self.integer("direct owner id")?;
        if native_type < -1 || native_owner_type < 0 || !matches!(owner_id, -3 | -1 | 0..) {
            return Err(error(format!("invalid {name} geometry record")));
        }
        if native_type != imported_type {
            limitations.push(format!(
                "{name} concrete type {native_type} is represented as {imported_type}; exact geometry is required"
            ));
        }
        if native_owner_type == 0 && owner_id == -1 {
            return Ok(NativeValue::Known(None));
        }
        if native_type == -1 || native_owner_type == 0 || owner_id == -1 {
            return Err(error(format!("inconsistent {name} owner record")));
        }
        if native_owner_type != owner_type || owner_id == -3 {
            limitations.push(format!(
                "{name} has an unrepresented direct owner (type {native_owner_type}, registry id {owner_id}); no owner-chain substitution is allowed"
            ));
            return Ok(NativeValue::Unrepresented);
        }
        let id = usize::try_from(owner_id).map_err(|_| error("owner id overflow"))?;
        if id >= owner_count {
            return Err(error(format!("{name} owner is outside the imported graph")));
        }
        let key = owner_key(u32::try_from(id).map_err(|_| error("owner id overflow"))?);
        if key.type_code != owner_type {
            return Err(error(format!("{name} owner-key mapping changes its type")));
        }
        Ok(NativeValue::Known(Some(key)))
    }
}

/// Install the representable metadata of the same native B-rep as brep.
///
/// keys supplies the caller's exact B-rep, face and edge identities. Parsing and
/// validation finish before writes begin. Concrete geometry is never replaced;
/// unknown owners and unrepresentable lists are left unchanged and reported.
/// Existing caches must be absent or invalidated before this fixture setup.
#[allow(clippy::too_many_lines)]
pub(crate) fn apply_query_metadata(
    brep: &mut Brep,
    raw: &[f64],
    keys: &dyn Fn(CacheObject) -> ObjectRef,
) -> Result<QueryMetadata> {
    let mut reader = Reader { raw, at: 0 };
    let mut report = QueryMetadata::default();
    if reader.integer("version")? != 1 {
        return Err(error("unsupported version"));
    }
    let tolerance = reader.scalar("B-rep tolerance")?;
    let infinite = reader.registry_id(
        brep.regions().len(),
        "infinite region",
        &mut report.topology,
    )?;
    reader.exact_count(brep.regions().len(), "region count")?;
    let mut void_flags = Vec::new();
    for _ in brep.regions().ids() {
        void_flags.push(
            i32::try_from(reader.integer("region void flag")?)
                .map_err(|_| error("region void flag overflow"))?,
        );
    }
    let regions = reader.list(
        brep.regions().len(),
        false,
        "region list",
        &mut report.topology,
    )?;
    let faces = reader.list(brep.faces().len(), true, "face list", &mut report.topology)?;
    let edges = reader.list(brep.edges().len(), true, "edge list", &mut report.topology)?;
    let vertices = reader.list(
        brep.vertices().len(),
        true,
        "vertex list",
        &mut report.topology,
    )?;

    reader.exact_count(brep.edges().len(), "curve owner count")?;
    let mut curve_owners = Vec::new();
    for (id, edge) in brep.edges().iter() {
        curve_owners.push(
            reader.geometry(
                edge.curve
                    .as_ref()
                    .map_or(-1, solidcore::curve::curve::Curve::type_code),
                (solidcore::topology::edge::TYPE, brep.edges().len()),
                |id| keys(CacheObject::Edge(EdgeId::from_index(id))),
                &format!("edge {} curve", id.index()),
                &mut report.curves,
            )?,
        );
    }
    reader.exact_count(brep.faces().len(), "surface owner count")?;
    let mut surface_owners = Vec::new();
    for (id, face) in brep.faces().iter() {
        surface_owners.push(
            reader.geometry(
                face.surface
                    .as_ref()
                    .map_or(-1, solidcore::surface::surface::Surface::type_code),
                (solidcore::topology::face::p2::TYPE, brep.faces().len()),
                |id| keys(CacheObject::Face(FaceId::from_index(id))),
                &format!("face {} surface", id.index()),
                &mut report.surfaces,
            )?,
        );
    }
    reader.exact_count(brep.faces().len(), "face state count")?;
    let mut face_states = Vec::new();
    for id in brep.faces().ids() {
        let faceuse = reader.registry_id(
            brep.faceuses().len(),
            &format!("face {} primary faceuse", id.index()),
            &mut report.topology,
        )?;
        let rectangular_trim = i32::try_from(reader.integer("rectangular-trim flag")?)
            .map_err(|_| error("rectangular-trim flag overflow"))?;
        face_states.push((faceuse, rectangular_trim));
    }
    if reader.at != raw.len() {
        return Err(error("trailing values"));
    }
    let brep_key = keys(CacheObject::Brep);
    if brep_key.type_code != BREP_TYPE {
        return Err(error("B-rep key mapping changes its type"));
    }

    brep.set_key(brep_key.key);
    brep.set_tolerance(tolerance);
    if let NativeValue::Known(infinite) = infinite {
        brep.set_infinite_region(infinite.map(RegionId::from_index));
    }
    let region_ids: Vec<_> = brep.regions().ids().collect();
    for (id, is_void) in region_ids.into_iter().zip(void_flags) {
        let region = brep.region_mut(id).map_err(|e| error(e.to_string()))?;
        region.is_void = is_void;
    }
    if let NativeValue::Known(list) = regions {
        brep.set_region_list(list.into_iter().map(RegionId::from_index).collect())
            .map_err(|e| error(e.to_string()))?;
    }
    if let NativeValue::Known(list) = faces {
        brep.set_face_list(list.into_iter().map(FaceId::from_index).collect())
            .map_err(|e| error(e.to_string()))?;
    }
    if let NativeValue::Known(list) = edges {
        brep.set_edge_list(list.into_iter().map(EdgeId::from_index).collect())
            .map_err(|e| error(e.to_string()))?;
    }
    if let NativeValue::Known(list) = vertices {
        brep.set_vertex_list(list.into_iter().map(VertexId::from_index).collect())
            .map_err(|e| error(e.to_string()))?;
    }
    let edge_ids: Vec<_> = brep.edges().ids().collect();
    for (id, owner) in edge_ids.into_iter().zip(curve_owners) {
        if let NativeValue::Known(owner) = owner {
            if let Some(curve) = brep
                .edge_mut(id)
                .map_err(|e| error(e.to_string()))?
                .curve
                .as_mut()
            {
                let _ = curve.set_owner(owner, owner.map(|_| brep_key));
            }
        }
    }
    let face_ids: Vec<_> = brep.faces().ids().collect();
    for ((id, owner), (faceuse, rectangular_trim)) in
        face_ids.into_iter().zip(surface_owners).zip(face_states)
    {
        let face = brep.face_mut(id).map_err(|e| error(e.to_string()))?;
        if let NativeValue::Known(owner) = owner {
            if let Some(surface) = face.surface.as_mut() {
                surface.set_owner(owner);
            }
        }
        if let NativeValue::Known(faceuse) = faceuse {
            face.faceuse = faceuse.map(FaceuseId::from_index);
        }
        face.rectangular_trim = rectangular_trim;
    }
    Ok(report)
}
