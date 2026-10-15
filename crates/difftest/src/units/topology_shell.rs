//! WIP native shell/region query comparisons. Other owned operations are pending.
//! Native fixtures and raw graph imports reuse the edge-use fixture protocol.
//! Setup failures fail the harness; topology and native integer flags compare exactly.

use super::topology_edgeuse_p1::{fixture_input, import};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::topology::shell::{shell_brep, shell_region, wire_edges};
use solidcore::topology::{FaceuseId, Region, Shell, ShellId};
use solidcore_ffi::units::{topology_edgeuse_p1 as graph, topology_shell as ffi};

/// The shell/region unit runner (coverage is still under development).
pub struct Runner;

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Inputs {
    brep: Value,
    #[serde(default)]
    void_state: Option<i32>,
    #[serde(default)]
    wires: bool,
    #[serde(default)]
    transfer: Option<TransferInputs>,
    #[serde(default)]
    type_probe: Option<TypeProbe>,
    #[serde(default)]
    degenerate: Option<DegenerateInput>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DegenerateInput {
    tolerance: Option<Value>,
}

fn degeneracy_tolerance(input: &DegenerateInput) -> Result<Option<f64>, HarnessError> {
    let Some(value) = &input.tolerance else {
        return Ok(None);
    };
    let value = if let Some(number) = value.as_f64() {
        number
    } else {
        match value.as_str() {
            Some("nan") => f64::NAN,
            Some("inf") => f64::INFINITY,
            Some("-inf") => f64::NEG_INFINITY,
            _ => return Err(failure("invalid degeneracy tolerance")),
        }
    };
    Ok(Some(value))
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TypeProbe {
    code: i64,
    shell: Option<usize>,
}

fn type_probe_index(probe: &TypeProbe, counts: [usize; 10]) -> Result<Option<usize>, HarnessError> {
    match probe.shell {
        None => Ok(None),
        Some(index) => {
            let count = counts.get(1).copied().unwrap_or(0);
            if count == 0 {
                return Err(failure("type fixture has no existing shell"));
            }
            Ok(Some(index % count))
        }
    }
}

fn rust_type_members(
    brep: &mut solidcore::topology::Brep,
    probe: &TypeProbe,
    selected: Option<usize>,
) -> Result<Output, HarnessError> {
    use solidcore::topology::shell::ShellMembers;
    let id = if let Some(index) = selected {
        ShellId::from_index(u32::try_from(index).map_err(failure)?)
    } else {
        let region = brep
            .infinite_region()
            .ok_or_else(|| failure("fixture has no infinite region"))?;
        let id = brep
            .add_shell_record(Shell {
                region: Some(region),
                ..Shell::new()
            })
            .map_err(failure)?;
        brep.region_mut(region).map_err(failure)?.shells.push(id);
        id
    };
    brep.shell_mut(id)
        .map_err(failure)?
        .set_shell_type(probe.code);
    let shell = brep.shell(id).map_err(failure)?;
    let (kind, head, count) = match &shell.members {
        ShellMembers::Empty => (0, -1, 0),
        ShellMembers::Faceuses(ids) => (
            if ids.is_empty() { 0 } else { 16022 },
            ids.first().map_or(-1, |id| i64::from(id.index())),
            ids.len(),
        ),
        ShellMembers::Vertexuse(id) => (16006, i64::from(id.index()), 1),
        ShellMembers::Edgeuse(id) => {
            let edge = brep
                .edgeuse(*id)
                .map_err(failure)?
                .edge
                .ok_or_else(|| failure("head has no edge"))?;
            (
                16005,
                i64::from(id.index()),
                brep.edge(edge).map_err(failure)?.edgeuses.len(),
            )
        }
    };
    Ok(Output::new().ints(
        "records",
        vec![
            shell.shell_type.code(),
            flag(shell.is_vertex_shell()),
            flag(shell.is_wire_shell()),
            flag(shell.is_faceuse_shell()),
            flag(shell.is_unknown_shell()),
            shell
                .first_faceuse()
                .map_err(failure)?
                .map_or(-1, |id| i64::from(id.index())),
            shell
                .wire_edge(brep)
                .map_err(failure)?
                .map_or(-1, |id| i64::from(id.index())),
            shell
                .vertex(brep)
                .map_err(failure)?
                .map_or(-1, |id| i64::from(id.index())),
            kind,
            head,
            count as i64,
        ],
    ))
}

const TYPES: [i64; 5] = [16013, 16012, 16000, 16999, -1];

fn failure(e: impl std::fmt::Display) -> HarnessError {
    HarnessError::new(e.to_string())
}
fn flag(value: bool) -> i64 {
    i64::from(value)
}
fn list(out: &mut Vec<i64>, values: impl IntoIterator<Item = u32>) {
    let values: Vec<_> = values.into_iter().collect();
    out.push(values.len() as i64);
    out.extend(values.into_iter().map(i64::from));
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TransferInputs {
    #[serde(default)]
    split_faceuse: bool,
    #[serde(default)]
    require_faceuse_pair: bool,
    #[serde(default)]
    extra_empty: usize,
    source: usize,
    recipient: usize,
    candidates: Option<Vec<i64>>,
}

struct ResolvedTransfer {
    split_faceuse: bool,
    require_faceuse_pair: bool,
    extra_empty: usize,
    source: usize,
    recipient: usize,
    candidates: Option<Vec<Option<usize>>>,
}

fn resolve_transfer(
    input: &TransferInputs,
    counts: [usize; 10],
) -> Result<ResolvedTransfer, HarnessError> {
    if input.extra_empty > 64 {
        return Err(failure("too many empty shells"));
    }
    let shells = counts
        .get(1)
        .copied()
        .and_then(|n| n.checked_add(input.extra_empty))
        .ok_or_else(|| failure("shell count overflow"))?;
    if shells == 0 {
        return Err(failure("transfer fixture has no shells"));
    }
    let edges = counts
        .get(7)
        .copied()
        .ok_or_else(|| failure("missing edge count"))?;
    let candidates = if let Some(list) = &input.candidates {
        let mut resolved = Vec::new();
        for &entry in list {
            resolved.push(if entry == -1 {
                None
            } else if entry < 0 || edges == 0 {
                return Err(failure("invalid wire candidate selector"));
            } else {
                Some(usize::try_from(entry).map_err(failure)? % edges)
            });
        }
        Some(resolved)
    } else {
        None
    };
    Ok(ResolvedTransfer {
        split_faceuse: input.split_faceuse,
        require_faceuse_pair: input.require_faceuse_pair,
        extra_empty: input.extra_empty,
        source: input.source % shells,
        recipient: input.recipient % shells,
        candidates,
    })
}

fn rust_transfer_records(
    brep: &mut solidcore::topology::Brep,
    options: &ResolvedTransfer,
) -> Result<Output, HarnessError> {
    use solidcore::topology::shell::{faceuse_sibling_links, transfer_children, ShellMembers};
    use solidcore::topology::{EdgeId, EdgeuseOwner, VertexuseOwner};
    use solidcore::SolidError;
    let region = brep
        .infinite_region()
        .ok_or_else(|| failure("fixture has no infinite region"))?;
    for _ in 0..options.extra_empty {
        let shell = brep
            .add_shell_record(Shell {
                region: Some(region),
                ..Shell::new()
            })
            .map_err(failure)?;
        brep.region_mut(region).map_err(failure)?.shells.push(shell);
    }
    if options.split_faceuse {
        if options.extra_empty == 0 {
            return Err(failure("split fixture needs an extra shell"));
        }
        let original = ShellId::from_index(0);
        let extra_index = brep
            .shells()
            .len()
            .checked_sub(options.extra_empty)
            .ok_or_else(|| failure("extra shell index underflow"))?;
        let extra = ShellId::from_index(u32::try_from(extra_index).map_err(failure)?);
        let members = match &mut brep.shell_mut(original).map_err(failure)?.members {
            ShellMembers::Faceuses(members) if members.len() >= 2 => members,
            _ => return Err(failure("split fixture needs a face-use cycle")),
        };
        let moved = members.remove(0);
        let target = brep.shell_mut(extra).map_err(failure)?;
        target.shell_type = solidcore::topology::shell::ShellType::Faceuse;
        target.members = ShellMembers::Faceuses(vec![moved]);
        brep.faceuse_mut(moved).map_err(failure)?.shell = Some(extra);
        // Construct the two baseline cycles independently of the transfer port.
        for id in [original, extra] {
            let members = match &brep.shell(id).map_err(failure)?.members {
                ShellMembers::Faceuses(members) => members.clone(),
                _ => return Err(failure("split fixture lost its cycle")),
            };
            let mut previous = members.last().copied();
            for (&member, &next) in members.iter().zip(members.iter().cycle().skip(1)) {
                brep.faceuse_mut(member).map_err(failure)?.sibling_links =
                    Some((Some(next), previous));
                previous = Some(member);
            }
        }
    }
    let candidates = if let Some(list) = &options.candidates {
        let mut ids = Vec::new();
        for &entry in list {
            ids.push(match entry {
                None => None,
                Some(index) => Some(EdgeId::from_index(u32::try_from(index).map_err(failure)?)),
            });
        }
        Some(ids)
    } else {
        None
    };
    let source = ShellId::from_index(u32::try_from(options.source).map_err(failure)?);
    let recipient = ShellId::from_index(u32::try_from(options.recipient).map_err(failure)?);
    if options.require_faceuse_pair {
        let source_record = brep.shell(source).map_err(failure)?;
        let target_record = brep.shell(recipient).map_err(failure)?;
        if source == recipient
            || source_record.shell_type.code() != 16022
            || target_record.shell_type.code() != 16022
            || !matches!(&source_record.members,ShellMembers::Faceuses(ids) if !ids.is_empty())
            || !matches!(&target_record.members,ShellMembers::Faceuses(ids) if !ids.is_empty())
        {
            return Err(failure(
                "fixture must contain two distinct nonempty face-use shells",
            ));
        }
    }
    let status = match transfer_children(brep, source, recipient, candidates.as_deref()) {
        Ok(()) => 0,
        Err(SolidError::Upstream { code, .. }) => i64::from(code),
        Err(error) => return Output::unported(&error).ok_or_else(|| failure(error)),
    };
    let mut out = vec![
        status,
        brep.shells().len() as i64,
        brep.faceuses().len() as i64,
        brep.edgeuses().len() as i64,
        brep.vertexuses().len() as i64,
        brep.regions().len() as i64,
    ];
    for (_, shell) in brep.shells().iter() {
        let (size, kind, head) = match &shell.members {
            ShellMembers::Empty => (0, 0, -1),
            ShellMembers::Faceuses(ids) => (
                ids.len(),
                if ids.is_empty() { 0 } else { 16022 },
                ids.first().map_or(-1, |id| i64::from(id.index())),
            ),
            ShellMembers::Edgeuse(id) => {
                let edge = brep
                    .edgeuse(*id)
                    .map_err(failure)?
                    .edge
                    .ok_or_else(|| failure("wire head has no edge"))?;
                (
                    brep.edge(edge).map_err(failure)?.edgeuses.len(),
                    16005,
                    i64::from(id.index()),
                )
            }
            ShellMembers::Vertexuse(id) => (1, 16006, i64::from(id.index())),
        };
        out.extend([
            shell.shell_type.code(),
            size as i64,
            shell.region.map_or(-1, |id| i64::from(id.index())),
            kind,
            head,
        ]);
    }
    for (id, faceuse) in brep.faceuses().iter() {
        let (next, last) = faceuse_sibling_links(brep, id).map_err(failure)?;
        out.extend([
            faceuse.shell.map_or(-1, |id| i64::from(id.index())),
            next.map_or(-1, |id| i64::from(id.index())),
            last.map_or(-1, |id| i64::from(id.index())),
        ]);
    }
    for (_, edgeuse) in brep.edgeuses().iter() {
        let (kind, id) = match edgeuse.owner {
            EdgeuseOwner::Shell(id) => (16013, id.map_or(-1, |id| i64::from(id.index()))),
            EdgeuseOwner::Loopuse(id) => (16004, id.map_or(-1, |id| i64::from(id.index()))),
            EdgeuseOwner::Unknown => (16999, -1),
        };
        out.extend([kind, id]);
    }
    for (_, vertexuse) in brep.vertexuses().iter() {
        let (kind, id) = match vertexuse.owner {
            VertexuseOwner::Shell(id) => (16013, id.map_or(-1, |id| i64::from(id.index()))),
            VertexuseOwner::Loopuse(id) => (16004, id.map_or(-1, |id| i64::from(id.index()))),
            VertexuseOwner::Edgeuse(id) => (16005, id.map_or(-1, |id| i64::from(id.index()))),
            VertexuseOwner::Unknown => (16999, -1),
        };
        out.extend([kind, id]);
    }
    for (_, region) in brep.regions().iter() {
        out.extend([
            region.shells.len() as i64,
            region.shells.first().map_or(-1, |id| i64::from(id.index())),
        ]);
        list(&mut out, region.shells.iter().copied().map(ShellId::index));
    }
    for shell in brep.shells().ids() {
        let wires = wire_edges(brep, shell).map_err(failure)?;
        out.push(0);
        list(&mut out, wires.into_iter().map(EdgeId::index));
    }
    Ok(Output::new().ints("records", out))
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "topology.shell"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let inputs: Inputs = case.parse()?;
        let (path, block) = fixture_input(inputs.brep)?;
        if inputs.transfer.is_some() && (inputs.wires || inputs.void_state.is_some()) {
            return Err(failure("transfer cannot be combined with query mutations"));
        }
        if let Some(degenerate) = inputs.degenerate {
            if inputs.type_probe.is_some()
                || inputs.transfer.is_some()
                || inputs.wires
                || inputs.void_state.is_some()
            {
                return Err(failure(
                    "degeneracy cannot be combined with another operation",
                ));
            }
            let values = ffi::degeneracy(&path, &block, degeneracy_tolerance(&degenerate)?)
                .map_err(failure)?;
            return Ok(Output::new().ints(
                "records",
                values
                    .into_iter()
                    .map(|value| value as i64)
                    .collect::<Vec<_>>(),
            ));
        }
        if let Some(probe) = inputs.type_probe {
            if inputs.transfer.is_some() || inputs.wires || inputs.void_state.is_some() {
                return Err(failure(
                    "type probe cannot be combined with another operation",
                ));
            }
            let export = graph::eval(graph::Op::Export, &path, &block).map_err(failure)?;
            let counts = graph::counts(&export)
                .ok_or_else(|| failure("native type fixture build failed"))?;
            let values =
                ffi::type_members(&path, &block, type_probe_index(&probe, counts)?, probe.code)
                    .map_err(failure)?;
            return Ok(Output::new().ints("records", values.to_vec()));
        }
        let native = if let Some(transfer) = inputs.transfer {
            let export = graph::eval(graph::Op::Export, &path, &block).map_err(failure)?;
            let counts = graph::counts(&export)
                .ok_or_else(|| failure("native transfer fixture build failed"))?;
            let options = resolve_transfer(&transfer, counts)?;
            ffi::transfer(
                &path,
                &block,
                options.extra_empty,
                options.source,
                options.recipient,
                options.candidates.as_deref(),
                options.split_faceuse,
            )
        } else if inputs.wires {
            ffi::wire_edges(&path, &block)
        } else {
            ffi::query(&path, &block, inputs.void_state)
        }
        .map_err(failure)?;
        if native
            .iter()
            .any(|v| !v.is_finite() || v.fract() != 0.0 || v.abs() > (1_u64 << 53) as f64)
        {
            return Err(failure("native topology query returned an invalid integer"));
        }
        Ok(Output::new().ints(
            "records",
            native.into_iter().map(|v| v as i64).collect::<Vec<_>>(),
        ))
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let inputs: Inputs = case.parse()?;
        let (path, block) = fixture_input(inputs.brep)?;
        let export = graph::eval(graph::Op::Export, &path, &block).map_err(failure)?;
        let mut brep = import(&export)?;
        if inputs.transfer.is_some() && (inputs.wires || inputs.void_state.is_some()) {
            return Err(failure("transfer cannot be combined with query mutations"));
        }
        if let Some(degenerate) = inputs.degenerate {
            if inputs.type_probe.is_some()
                || inputs.transfer.is_some()
                || inputs.wires
                || inputs.void_state.is_some()
            {
                return Err(failure(
                    "degeneracy cannot be combined with another operation",
                ));
            }
            let tolerance = degeneracy_tolerance(&degenerate)?;
            let kernel = solidcore::kernel::Kernel::default();
            let mut out = vec![brep.shells().len() as i64, brep.regions().len() as i64];
            for id in brep.shells().ids() {
                match solidcore::topology::shell::shell_is_degenerate(&kernel, &brep, id, tolerance)
                {
                    Ok(value) => out.push(flag(value)),
                    Err(error) => return Output::unported(&error).ok_or_else(|| failure(error)),
                }
            }
            for id in brep.regions().ids() {
                match solidcore::topology::shell::region_is_degenerate(
                    &kernel, &brep, id, tolerance,
                ) {
                    Ok(value) => out.push(flag(value)),
                    Err(error) => return Output::unported(&error).ok_or_else(|| failure(error)),
                }
            }
            return Ok(Output::new().ints("records", out));
        }
        if let Some(probe) = inputs.type_probe {
            if inputs.transfer.is_some() || inputs.wires || inputs.void_state.is_some() {
                return Err(failure(
                    "type probe cannot be combined with another operation",
                ));
            }
            let counts = graph::counts(&export)
                .ok_or_else(|| failure("native type fixture build failed"))?;
            return rust_type_members(&mut brep, &probe, type_probe_index(&probe, counts)?);
        }
        if let Some(transfer) = inputs.transfer {
            let counts = graph::counts(&export)
                .ok_or_else(|| failure("native transfer fixture build failed"))?;
            return rust_transfer_records(&mut brep, &resolve_transfer(&transfer, counts)?);
        }
        if inputs.wires {
            let mut out = vec![brep.shells().len() as i64];
            for shell in brep.shells().ids() {
                match wire_edges(&brep, shell) {
                    Ok(wires) => {
                        out.push(0);
                        list(
                            &mut out,
                            wires.into_iter().map(solidcore::topology::EdgeId::index),
                        );
                    }
                    Err(error) => return Output::unported(&error).ok_or_else(|| failure(error)),
                }
            }
            return Ok(Output::new().ints("records", out));
        }
        let mut out = vec![brep.shells().len() as i64, brep.regions().len() as i64];
        for (id, shell) in brep.shells().iter() {
            out.extend([
                Shell::type_code(),
                flag(shell.is_vertex_shell()),
                flag(shell.is_wire_shell()),
                flag(shell.is_faceuse_shell()),
                flag(shell.is_unknown_shell()),
                flag(shell.is_closed(&brep, id).map_err(failure)?),
                shell_region(&brep, id)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                flag(shell_brep(&brep, id).map_err(failure)?.is_some()),
                shell
                    .first_faceuse()
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                shell
                    .wire_edge(&brep)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                shell
                    .vertex(&brep)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
            ]);
            out.extend(TYPES.into_iter().map(|ty| flag(Shell::is_kind_of(ty))));
            list(
                &mut out,
                shell
                    .get_faceuses()
                    .map_err(failure)?
                    .into_iter()
                    .map(FaceuseId::index),
            );
            list(
                &mut out,
                shell
                    .mated_shells(&brep)
                    .map_err(failure)?
                    .into_iter()
                    .map(ShellId::index),
            );
        }
        let regions: Vec<_> = brep.regions().ids().collect();
        for id in regions {
            if let Some(value) = inputs.void_state {
                brep.region_mut(id).map_err(failure)?.set_is_void(value);
            }
            let region = brep.region(id).map_err(failure)?;
            out.extend([
                Region::type_code(),
                i64::from(region.is_void()),
                flag(region.is_infinite_region(&brep, id)),
                flag(region.get_brep(&brep).is_some()),
            ]);
            out.extend(TYPES.into_iter().map(|ty| flag(Region::is_kind_of(ty))));
            list(
                &mut out,
                region
                    .get_shells(None, |_, _| Ok(false))
                    .map_err(failure)?
                    .into_iter()
                    .map(ShellId::index),
            );
        }
        Ok(Output::new().ints("records", out))
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let scale = [1e-3, 1.0, 1e3].get(rng.below(3) as usize).copied()?;
        let origin = rng.point3(-2.0, 2.0).map(|v| v * scale);
        let normal = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), 1.0];
        let wires = rng.chance(0.5);
        let void_state = if wires {
            None
        } else {
            Some(rng.below(3) as i32)
        };
        Some(
            json!({"brep": {"kind":"face", "surface":{"plane":{"o":origin,"n":normal}},
            "dom":[-scale,-scale,scale,scale]}, "void_state":void_state, "wires":wires}),
        )
    }
}
