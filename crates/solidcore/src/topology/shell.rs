//! Shell and region topology, with typed arena links and native traversal order.
//!
//! Upstream: `source/SMLib/inc/SmShell.h`, `source/SMLib/src/SmShell.cpp`,
//! `source/SMLib/inc/SmRegion.h`, `source/SMLib/src/SmRegion.cpp`,
//! and the minimum default-argument branches of `source/SMLib/src/SmEdge.cpp`
//! and `source/SMLib/src/SmVertex.cpp` used by wire-edge traversal;
//! `source/SMLib/inc/SmTopology.h` supplies literal sibling links and
//! `source/SMLib/inc/SmObject.h` defines the checked `SM_CAST_PTR` semantics.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`ShellType`], [`ShellMembers`] | `m_tShellType`, typed `m_pList` contents |
//! | [`Shell::new`], [`Region::new`], implicit `Drop` | constructors, destructors |
//! | [`Shell::get_faceuses`], [`Shell::first_faceuse`] | `GetFaceuses`, `GetFirstFaceuse` |
//! | [`Shell::wire_edge`], [`Shell::vertex`] | `GetWireEdge`, `GetVertex` |
//! | [`Shell::is_closed`], [`Shell::mated_shells`] | `IsClosed`, `GetMatedShells` |
//! | [`shell_region`], [`shell_brep`] | `GetRegion`, `GetBrep` |
//! | [`wire_edges`] | `SmShell::GetWireEdges` |
//! | [`transfer_children`] | `SmShell::TransferChildrenTo` |
//! | [`faceuse_sibling_links`] | typed `SmTopology::GetNext` / `GetLast` for face uses |
//! | `edge_vertices_for_shell`, `edge_is_wire_for_shell` | `SmEdge::GetVertices` (no attribute filter), `IsWire` |
//! | `vertex_edges_for_shell` | `SmVertex::GetVertexuses`, `GetEdges` (no face filter) |
//! | [`Region::get_shells`], [`Region::is_void`], [`Region::set_is_void`] | `GetShells`, `IsVoid`, `SetIsVoid` |
//! | [`Region::is_infinite_region`] | `IsInfiniteRegion` |
//! | [`Shell::is_inner_shell`], [`order_region_shells`] | `IsInnerShell`, `OrderShells` |
//!
//! Deviations: circular lists are typed vectors in head-first order, as in the
//! shared Brep arena model. Temporary marks use a local set without changing
//! first-seen order. Invalid native pointer casts return `InvalidInput`. Native
//! notification/logging side effects are not represented by these arena records.
//! Face-use sibling links survive a cleared owner head, as in native self-transfer.
//! There is no recursion. This is WIP; remaining operations and native gates are pending.

use crate::error::{Result, SolidError};
use crate::topology::proto_topology::p1::owning_topology;
use crate::topology::{
    Brep, EdgeId, EdgeuseId, FaceuseId, Region, RegionId, Shell, ShellId, VertexId, VertexuseId,
};
use std::collections::BTreeSet;

/// The native shell member class, independently stored from its list head.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ShellType {
    /// `SmUnknown_TYPE` (16999), the constructor value.
    #[default]
    Unknown,
    /// `SmVertexuse_TYPE` (16006).
    Vertex,
    /// `SmEdgeuse_TYPE` (16005).
    Wire,
    /// `SmFaceuse_TYPE` (16022).
    Faceuse,
    /// Other values preserved by `SetShellType` without reinterpretation.
    Other(i64),
}

impl ShellType {
    /// The exact native type code.
    #[must_use]
    pub const fn code(self) -> i64 {
        match self {
            Self::Unknown => 16999,
            Self::Vertex => 16006,
            Self::Wire => 16005,
            Self::Faceuse => 16022,
            Self::Other(code) => code,
        }
    }
    /// Retain an arbitrary native caller's shell type.
    #[must_use]
    pub const fn from_code(code: i64) -> Self {
        match code {
            16999 => Self::Unknown,
            16006 => Self::Vertex,
            16005 => Self::Wire,
            16022 => Self::Faceuse,
            other => Self::Other(other),
        }
    }
}

/// A checked representation of the shell's native list head and members.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum ShellMembers {
    /// NULL list head.
    #[default]
    Empty,
    /// Face-use circular list, from its current head.
    Faceuses(Vec<FaceuseId>),
    /// First wire edge use (the edge owns its radial list).
    Edgeuse(EdgeuseId),
    /// Vertex use of a vertex shell.
    Vertexuse(VertexuseId),
}

const fn invalid_members() -> SolidError {
    SolidError::InvalidInput {
        operation: "topology::shell",
        reason: "shell type and native member class disagree",
    }
}

impl Shell {
    /// Empty, unknown shell (`SmShell()`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    /// Native class type (`SmShell_TYPE`).
    #[must_use]
    pub const fn type_code() -> i64 {
        16013
    }
    /// Native class name.
    #[must_use]
    pub const fn class_string() -> &'static str {
        "SmShell"
    }
    /// Native type name.
    #[must_use]
    pub const fn type_string() -> &'static str {
        "SmShell_TYPE"
    }
    /// Native type/base-class membership, using the shared owning-topology port.
    #[must_use]
    pub const fn is_kind_of(code: i64) -> bool {
        code == Self::type_code() || owning_topology::is_kind_of(code)
    }
    /// Whether this is a vertex-use shell.
    #[must_use]
    pub const fn is_vertex_shell(&self) -> bool {
        self.shell_type.code() == 16006
    }
    /// Whether this is a wire-edge shell.
    #[must_use]
    pub const fn is_wire_shell(&self) -> bool {
        self.shell_type.code() == 16005
    }
    /// Whether this contains face uses.
    #[must_use]
    pub const fn is_faceuse_shell(&self) -> bool {
        self.shell_type.code() == 16022
    }
    /// Whether its type is the native unknown value.
    #[must_use]
    pub const fn is_unknown_shell(&self) -> bool {
        self.shell_type.code() == 16999
    }
    /// Change only the native type field, retaining its existing list head.
    pub fn set_shell_type(&mut self, code: i64) {
        self.shell_type = ShellType::from_code(code);
    }
    /// Face uses in native list order; empty for every other shell type.
    ///
    /// # Errors
    /// A native list head of an incompatible class.
    pub fn get_faceuses(&self) -> Result<Vec<FaceuseId>> {
        if !self.is_faceuse_shell() {
            return Ok(Vec::new());
        }
        match &self.members {
            ShellMembers::Empty => Ok(Vec::new()),
            ShellMembers::Faceuses(ids) => Ok(ids.clone()),
            _ => Err(invalid_members()),
        }
    }
    /// First face use; native logging-only errors retain the NULL result.
    ///
    /// A head-class mismatch retains native NULL rather than becoming an error.
    ///
    /// # Errors
    /// This compatibility result succeeds; native logging-only errors yield `None`.
    pub fn first_faceuse(&self) -> Result<Option<FaceuseId>> {
        if !self.is_faceuse_shell() {
            return Ok(None);
        }
        // SM_CAST_PTR checks the actual head class and returns NULL on mismatch.
        // Unlike the reinterpret-cast GetFaceuses array, this is a defined path.
        Ok(match &self.members {
            ShellMembers::Faceuses(ids) => ids.first().copied(),
            _ => None,
        })
    }
    /// First wire edge, preserving NULL for empty or other-kind shells.
    ///
    /// # Errors
    /// Dangling arena references; a head-class mismatch retains native NULL.
    pub fn wire_edge(&self, brep: &Brep) -> Result<Option<EdgeId>> {
        if !self.is_wire_shell() {
            return Ok(None);
        }
        match self.members {
            ShellMembers::Edgeuse(id) => Ok(brep.edgeuse(id)?.edge),
            _ => Ok(None),
        }
    }
    /// Vertex of a vertex shell, preserving NULL for empty or other-kind shells.
    ///
    /// # Errors
    /// Dangling arena references; a head-class mismatch retains native NULL.
    pub fn vertex(&self, brep: &Brep) -> Result<Option<VertexId>> {
        if !self.is_vertex_shell() {
            return Ok(None);
        }
        match self.members {
            ShellMembers::Vertexuse(id) => Ok(brep.vertexuse(id)?.vertex),
            _ => Ok(None),
        }
    }
    /// Whether a face in this shell has a mate belonging to a different shell.
    ///
    /// # Errors
    /// Invalid list classes or dangling face-use/mate references.
    pub fn is_closed(&self, brep: &Brep, id: ShellId) -> Result<bool> {
        if !self.is_faceuse_shell() {
            return Ok(false);
        }
        for faceuse in self.get_faceuses()? {
            if let Some(mate) = brep.faceuse(faceuse)?.mate {
                if brep.faceuse(mate)?.shell != Some(id) {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
    /// Distinct mate shells in first-seen face-use order.
    ///
    /// # Errors
    /// Missing native mates/owners or dangling arena references.
    pub fn mated_shells(&self, brep: &Brep) -> Result<Vec<ShellId>> {
        let mut output = Vec::new();
        if !self.is_faceuse_shell() {
            return Ok(output);
        }
        let mut seen = BTreeSet::new();
        for faceuse in self.get_faceuses()? {
            let mate = brep.faceuse(faceuse)?.mate.ok_or(invalid_members())?;
            let shell = brep.faceuse(mate)?.shell.ok_or(invalid_members())?;
            brep.shell(shell)?;
            if seen.insert(shell) {
                output.push(shell);
            }
        }
        Ok(output)
    }
    /// Classify whether this is an inner shell through the topology traverser.
    ///
    /// # Errors
    /// Invalid member classes or the unported topology traversal callee.
    pub fn is_inner_shell(&self) -> Result<bool> {
        let _faceuses = self.get_faceuses()?;
        // UNPORTED(topology.topology_traverser): SmTopologyTraverser::IsInnerShell
        Err(SolidError::unported(
            "SmTopologyTraverser::IsInnerShell",
            "topology.topology_traverser",
        ))
    }
}

/// Owning region (`SmShell::GetRegion`).
///
/// # Errors
/// Dangling shell id.
pub fn shell_region(brep: &Brep, id: ShellId) -> Result<Option<RegionId>> {
    Ok(brep.shell(id)?.region)
}

/// B-rep of a shell's owning region (`SmShell::GetBrep`).
///
/// # Errors
/// Dangling shell or region id.
pub fn shell_brep(brep: &Brep, id: ShellId) -> Result<Option<&Brep>> {
    match shell_region(brep, id)? {
        None => Ok(None),
        Some(region) => Ok(brep.region(region)?.in_brep.then_some(brep)),
    }
}

impl Region {
    /// Empty non-void region (`SmRegion()`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    /// Native region type.
    #[must_use]
    pub const fn type_code() -> i64 {
        16012
    }
    /// Native class name.
    #[must_use]
    pub const fn class_string() -> &'static str {
        "SmRegion"
    }
    /// Native type name.
    #[must_use]
    pub const fn type_string() -> &'static str {
        "SmRegion_TYPE"
    }
    /// Native type/base-class membership.
    #[must_use]
    pub const fn is_kind_of(code: i64) -> bool {
        code == Self::type_code() || owning_topology::is_kind_of(code)
    }
    /// Native void state, preserving `FALSE` (0), `TRUE` (1), and `UNSURE` (2).
    #[must_use]
    pub const fn is_void(&self) -> i32 {
        self.is_void
    }
    /// Set the native void flag.
    pub fn set_is_void(&mut self, value: i32) {
        self.is_void = value;
    }
    /// Owning B-rep, or NULL when the region is detached (`GetBrep`).
    #[must_use]
    pub const fn get_brep<'a>(&self, brep: &'a Brep) -> Option<&'a Brep> {
        if self.in_brep {
            Some(brep)
        } else {
            None
        }
    }

    /// Whether this is the owning B-rep's designated infinite region.
    #[must_use]
    pub fn is_infinite_region(&self, brep: &Brep, id: RegionId) -> bool {
        self.in_brep && brep.infinite_region() == Some(id)
    }
    /// Region shell list, retaining optional attribute filtering via the caller's
    /// actual attribute lookup rather than dropping or defaulting that parameter.
    ///
    /// # Errors
    /// Those of the attribute lookup.
    pub fn get_shells<F>(
        &self,
        attribute_id: Option<u64>,
        mut has_attribute: F,
    ) -> Result<Vec<ShellId>>
    where
        F: FnMut(ShellId, u64) -> Result<bool>,
    {
        let mut output = Vec::new();
        for &shell in &self.shells {
            if let Some(attribute) = attribute_id {
                if !has_attribute(shell, attribute)? {
                    continue;
                }
            }
            output.push(shell);
        }
        Ok(output)
    }
}

/// Select the first shell according to `SmRegion::OrderShells`, rotating the
/// circular head while preserving every other link's native relative order.
///
/// # Errors
/// Invalid arena references or an unported inner-shell classification.
pub fn order_region_shells(brep: &mut Brep, id: RegionId) -> Result<()> {
    let shells = brep.region(id)?.shells.clone();
    if shells.len() < 2 {
        return Ok(());
    }
    let first = shells.first().copied().ok_or(invalid_members())?;
    if brep.shell(first)?.is_faceuse_shell() {
        return Ok(());
    }
    let mut faceuse_shells = Vec::new();
    for &shell in &shells {
        if brep.shell(shell)?.is_faceuse_shell() {
            faceuse_shells.push(shell);
        }
    }
    let mut chosen = faceuse_shells.first().copied();
    if faceuse_shells.len() > 1 {
        for &shell in &faceuse_shells {
            if brep.shell(shell)?.is_closed(brep, shell)? {
                match brep.shell(shell)?.is_inner_shell() {
                    Ok(true) => {
                        chosen = Some(shell);
                        break;
                    }
                    Err(error) if error.is_unported() => return Err(error),
                    Ok(false) | Err(_) => {}
                }
            }
        }
    }
    if let Some(chosen) = chosen {
        let position = shells
            .iter()
            .position(|&shell| shell == chosen)
            .ok_or(invalid_members())?;
        brep.region_mut(id)?.shells.rotate_left(position);
    }
    Ok(())
}

/// Wire edges of a shell, preserving native stack and first-seen order.
///
/// # Errors
/// Dangling topology links, native undefined pointer accesses, or the missing
/// vertex edge-list traversal. There is no recursion.
pub fn wire_edges(brep: &Brep, id: ShellId) -> Result<Vec<EdgeId>> {
    let shell = brep.shell(id)?;
    let mut output = Vec::new();
    if shell.is_vertex_shell() {
        return Ok(output);
    }
    if shell.is_wire_shell() {
        let first = shell.wire_edge(brep)?.ok_or(invalid_members())?;
        let mut stack = vec![first];
        let mut seen_vertices = BTreeSet::new();
        let mut seen_edges = BTreeSet::new();
        while let Some(edge) = stack.pop() {
            for vertex in edge_vertices_for_shell(brep, edge)? {
                if !seen_vertices.insert(vertex) {
                    continue;
                }
                for adjacent in vertex_edges_for_shell(brep, vertex)? {
                    if seen_edges.insert(adjacent) {
                        output.push(adjacent);
                        stack.push(adjacent);
                    }
                }
            }
        }
        return Ok(output);
    }
    // The shared arena model represents SmBrep's edge list in insertion order.
    // GetBrep is required even when the edge list happens to be empty.
    let owner = shell_brep(brep, id)?.ok_or(invalid_members())?;
    for (edge, _) in owner.edges().iter() {
        if edge_is_wire_for_shell(owner, edge)? {
            let primary = crate::topology::edgeuse::p1::primary_edgeuse(owner, edge)?
                .ok_or(invalid_members())?;
            if owner.edgeuse(primary)?.owner == crate::topology::EdgeuseOwner::Shell(Some(id)) {
                output.push(edge);
            }
        }
    }
    Ok(output)
}

// Minimal no-attribute GetVertices branch needed by GetWireEdges. Endpoint
// order follows the primary edge use's orientation, not a guessed vertex pair.
fn edge_vertices_for_shell(brep: &Brep, edge: EdgeId) -> Result<Vec<VertexId>> {
    use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
    use crate::topology::edgeuse::p1;
    let Some(primary) = p1::primary_edgeuse(brep, edge)? else {
        return Ok(Vec::new());
    };
    let mate = p1::mate(brep, primary)?.ok_or(invalid_members())?;
    let first_use = p1::vertexuse(brep, primary)?.ok_or(invalid_members())?;
    let second_use = p1::vertexuse(brep, mate)?.ok_or(invalid_members())?;
    let mut first = brep.vertexuse(first_use)?.vertex;
    let mut second = brep.vertexuse(second_use)?.vertex;
    if brep.edgeuse(primary)?.orientation != OrientType::Same {
        std::mem::swap(&mut first, &mut second);
    }
    // Source adds the first pointer unconditionally. GetWireEdges immediately
    // dereferences it, so NULL is an undefined native input on this path.
    let first = first.ok_or(invalid_members())?;
    let mut output = vec![first];
    if let Some(second) = second {
        if second != first {
            output.push(second);
        }
    }
    Ok(output)
}

fn edge_is_wire_for_shell(brep: &Brep, edge: EdgeId) -> Result<bool> {
    use crate::topology::edgeuse::p1;
    match p1::primary_edgeuse(brep, edge)? {
        Some(primary) => p1::is_shell_edgeuse(brep, primary),
        None => Ok(false),
    }
}

fn vertex_edges_for_shell(brep: &Brep, vertex: VertexId) -> Result<Vec<EdgeId>> {
    use crate::topology::edgeuse::p1;
    use crate::topology::VertexuseOwner;
    let mut output = Vec::new();
    for &vertexuse in &brep.vertex(vertex)?.vertexuses {
        if let VertexuseOwner::Edgeuse(edgeuse) = brep.vertexuse(vertexuse)?.owner {
            let edgeuse = edgeuse.ok_or(invalid_members())?;
            // Upstream asks for the face even when its optional filter is NULL.
            let _face = p1::face(brep, edgeuse)?;
            let edge = brep.edgeuse(edgeuse)?.edge.ok_or(invalid_members())?;
            if !output.contains(&edge) {
                output.push(edge);
            }
        }
    }
    Ok(output)
}

/// Move a shell's vertex, wire edges and face uses to another shell.
///
/// `wire_candidates` is the native optional list of B-rep wire edges. When it
/// is present, ownership filters it in caller order; otherwise `GetWireEdges`
/// gathers the source's wires. It is never replaced by an empty default list.
/// Nullable list entries are only rejected if the native branch reads them.
///
/// # Errors
/// Native `SM_ERR` for an illegal vertex-shell move, dangling graph links, or
/// incompatible native list classes. The caller must supply connected topology.
pub fn transfer_children(
    brep: &mut Brep,
    source: ShellId,
    recipient: ShellId,
    wire_candidates: Option<&[Option<EdgeId>]>,
) -> Result<()> {
    use crate::topology::edgeuse::p1;
    use crate::topology::{EdgeuseOwner, VertexuseOwner};
    const OP: &str = "topology::shell::transfer_children";
    brep.shell(source)?;
    if shell_list_size(brep, recipient)? == 0 {
        brep.shell_mut(recipient)?.shell_type = ShellType::Unknown;
        if brep.shell(source)?.is_vertex_shell() {
            let ShellMembers::Vertexuse(vertexuse) = brep.shell(source)?.members else {
                return Err(invalid_members());
            };
            brep.vertexuse(vertexuse)?;
            let target = brep.shell_mut(recipient)?;
            target.shell_type = ShellType::Vertex;
            target.members = ShellMembers::Vertexuse(vertexuse);
            brep.vertexuse_mut(vertexuse)?.owner = VertexuseOwner::Shell(Some(recipient));
            clear_shell_members(brep, source)?;
            return Ok(());
        }
    } else {
        if brep.shell(recipient)?.is_vertex_shell() {
            return if shell_list_size(brep, source)? == 0 {
                Ok(())
            } else {
                Err(SolidError::sm_err(OP))
            };
        }
        if brep.shell(source)?.is_vertex_shell() {
            return Err(SolidError::sm_err(OP));
        }
    }

    let wires = if let Some(candidates) = wire_candidates {
        let mut wires = Vec::new();
        for &edge in candidates {
            let edge = edge.ok_or(invalid_members())?;
            let primary = p1::primary_edgeuse(brep, edge)?.ok_or(invalid_members())?;
            if brep.edgeuse(primary)?.owner == EdgeuseOwner::Shell(Some(source)) {
                wires.push(edge);
            }
        }
        wires
    } else {
        wire_edges(brep, source)?
    };
    if let Some(&first) = wires.first() {
        if brep.shell(recipient)?.is_unknown_shell() {
            let primary = p1::primary_edgeuse(brep, first)?.ok_or(invalid_members())?;
            let target = brep.shell_mut(recipient)?;
            target.members = ShellMembers::Edgeuse(primary);
            target.shell_type = ShellType::Wire;
            // UpdateListSize counts this edge's radial list; the shared arena
            // model retains that list on Edge::edgeuses rather than separately.
        }
        for edge in wires {
            // SmEdge::GetEdgeuses with neither optional filter. The shared edge
            // record already retains its native circular radial list order.
            let edgeuses = brep.edge(edge)?.edgeuses.clone();
            for edgeuse in edgeuses {
                brep.edgeuse_mut(edgeuse)?.owner = EdgeuseOwner::Shell(Some(recipient));
            }
        }
    }

    let faceuses = brep.shell(source)?.get_faceuses()?;
    if !faceuses.is_empty() {
        if !brep.shell(recipient)?.is_faceuse_shell() {
            let target = brep.shell_mut(recipient)?;
            target.shell_type = ShellType::Faceuse;
            target.members = ShellMembers::Empty;
        }
        for faceuse in faceuses {
            if brep.faceuse(faceuse)?.shell != Some(source) {
                return Err(invalid_members());
            }
            // Remove, then PostInsert: the member moves to the recipient's tail
            // without changing either list's remaining relative order.
            match &mut brep.shell_mut(source)?.members {
                ShellMembers::Faceuses(members) => {
                    let index = members
                        .iter()
                        .position(|id| *id == faceuse)
                        .ok_or(invalid_members())?;
                    members.remove(index);
                }
                _ => return Err(invalid_members()),
            }
            let target = &mut brep.shell_mut(recipient)?.members;
            match target {
                ShellMembers::Empty => *target = ShellMembers::Faceuses(vec![faceuse]),
                ShellMembers::Faceuses(members) => members.push(faceuse),
                _ => return Err(invalid_members()),
            }
            brep.faceuse_mut(faceuse)?.shell = Some(recipient);
        }
        // Self-transfer clears the shell head but leaves this face-use cycle
        // connected. Preserve those native links independently of the head.
        refresh_faceuse_sibling_links(brep, recipient)?;
    }
    clear_shell_members(brep, source)
}

fn clear_shell_members(brep: &mut Brep, id: ShellId) -> Result<()> {
    let shell = brep.shell_mut(id)?;
    shell.shell_type = ShellType::Unknown;
    shell.members = ShellMembers::Empty;
    Ok(())
}

fn shell_list_size(brep: &Brep, id: ShellId) -> Result<usize> {
    match &brep.shell(id)?.members {
        ShellMembers::Empty => Ok(0),
        ShellMembers::Faceuses(members) => Ok(members.len()),
        ShellMembers::Vertexuse(_) => Ok(1),
        ShellMembers::Edgeuse(id) => {
            let edge = brep.edgeuse(*id)?.edge.ok_or(invalid_members())?;
            Ok(brep.edge(edge)?.edgeuses.len())
        }
    }
}

/// Native next/previous face-use links, independently of a shell's list head.
pub type FaceuseSiblingLinks = (Option<FaceuseId>, Option<FaceuseId>);

/// Typed adapter for a face use's native `GetNext` and `GetLast` links.
///
/// # Errors
/// A dangling id or a constructed owner vector that does not contain its member.
pub fn faceuse_sibling_links(brep: &Brep, id: FaceuseId) -> Result<FaceuseSiblingLinks> {
    let faceuse = brep.faceuse(id)?;
    if let Some(links) = faceuse.sibling_links {
        return Ok(links);
    }
    let Some(owner) = faceuse.shell else {
        return Ok((None, None));
    };
    let ShellMembers::Faceuses(members) = &brep.shell(owner)?.members else {
        return Err(invalid_members());
    };
    let mut previous = members.last().copied();
    for (&member, &next) in members.iter().zip(members.iter().cycle().skip(1)) {
        if member == id {
            return Ok((Some(next), previous));
        }
        previous = Some(member);
    }
    Err(invalid_members())
}

fn refresh_faceuse_sibling_links(brep: &mut Brep, owner: ShellId) -> Result<()> {
    let members = brep.shell(owner)?.get_faceuses()?;
    let mut previous = members.last().copied();
    for (&member, &next) in members.iter().zip(members.iter().cycle().skip(1)) {
        brep.faceuse_mut(member)?.sibling_links = Some((Some(next), previous));
        previous = Some(member);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::Faceuse;

    #[test]
    fn region_order_rotates_the_circular_head_without_moving_other_links() {
        let mut brep = Brep::new();
        let ids = [
            ShellType::Wire,
            ShellType::Unknown,
            ShellType::Faceuse,
            ShellType::Vertex,
        ]
        .into_iter()
        .map(|shell_type| {
            brep.add_shell_record(Shell {
                shell_type,
                ..Shell::default()
            })
            .unwrap()
        })
        .collect::<Vec<_>>();
        let region = brep
            .add_region_record(Region {
                shells: ids.clone(),
                ..Region::default()
            })
            .unwrap();
        order_region_shells(&mut brep, region).unwrap();
        assert_eq!(
            brep.region(region).unwrap().shells,
            vec![ids[2], ids[3], ids[0], ids[1]]
        );
    }

    #[test]
    fn mated_shells_keep_first_seen_order_when_a_shell_reappears() {
        let mut brep = Brep::new();
        let source = brep.add_shell_record(Shell::default()).unwrap();
        let a = brep.add_shell_record(Shell::default()).unwrap();
        let b = brep.add_shell_record(Shell::default()).unwrap();
        let mut sources = Vec::new();
        for target in [b, a, b] {
            let mate = brep
                .add_faceuse_record(Faceuse {
                    shell: Some(target),
                    ..Faceuse::default()
                })
                .unwrap();
            sources.push(
                brep.add_faceuse_record(Faceuse {
                    shell: Some(source),
                    mate: Some(mate),
                    ..Faceuse::default()
                })
                .unwrap(),
            );
        }
        let shell = brep.shell_mut(source).unwrap();
        shell.shell_type = ShellType::Faceuse;
        shell.members = ShellMembers::Faceuses(sources);
        assert_eq!(
            brep.shell(source).unwrap().mated_shells(&brep).unwrap(),
            vec![b, a]
        );
        assert!(brep
            .shell(source)
            .unwrap()
            .is_closed(&brep, source)
            .unwrap());
    }

    #[test]
    fn changing_shell_type_preserves_the_list_head_and_members() {
        let members = vec![FaceuseId::from_index(2), FaceuseId::from_index(0)];
        let mut shell = Shell {
            shell_type: ShellType::Faceuse,
            members: ShellMembers::Faceuses(members.clone()),
            ..Shell::default()
        };
        shell.set_shell_type(16999);
        assert!(shell.get_faceuses().unwrap().is_empty());
        shell.set_shell_type(16022);
        assert_eq!(shell.get_faceuses().unwrap(), members);
    }
}
