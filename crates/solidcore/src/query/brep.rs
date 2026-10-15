//! B-rep queries and checked surface evaluation entry points.
//!
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`,
//! `source/SMLib/inc/SmPolyDecimate.h`, `source/SMLib/src/SmPolyDecimate.cpp`,
//! `source/SMLib/inc/SmPoly.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`evaluate_surface_point`] | `SmApiEvaluateSurfacePoint` |
//! | [`evaluate_surface_normal`] | `SmApiEvaluateSurfaceNormal` |
//! | [`evaluate_surface_derivatives`] | `SmApiEvaluateSurfaceDerivatives` |
//! | [`PolyVertexClass`] | `SmPolyVertexClass` |
//! | [`decimate_vertex_attribute`], [`DecimateVertexData`] | `SmDecimateVertexAttr` constructor/destructor |
//! | [`quadric_vertex_attribute`], [`QuadricVertexData`] | `SmQuadricVertexAttr` constructor/destructor |
//! | [`quadric_edge_attribute`], [`QuadricEdgeData`] | `SmQuadricEdgeAttr` constructor/destructor |
//! | [`crate::geom::attribute::Attribute::make_copy`] | the three `MakeCopy` overrides |
//!
//! Deviations: a nullable surface is `Option<&Surface>`; output parameters become
//! return values, so failure cannot modify caller-owned outputs. The remaining
//! B-rep functions assigned to query.brep are not implemented yet. Destructors
//! use ownership rather than poisoning fields immediately before freeing them.
//! The quadric edge constructor leaves its sort index uninitialised in C++; the
//! port represents that state by `None` instead of exposing undefined memory.

use crate::geom::vector2d::Point2;
use crate::geom::vector3d::{Point3, Vec3};
use crate::nurbs::knots::Side;
use crate::surface::surface::Surface;
use crate::{Result, SolidError};

fn checked_surface<'a>(
    surface: Option<&'a Surface>,
    uv: Point2,
    op: &'static str,
) -> Result<&'a Surface> {
    let surface = surface.ok_or_else(|| SolidError::invalid_input(op, "null surface"))?;
    let domain = surface.natural_uv_domain();
    let lo = domain.min();
    let hi = domain.max();
    // Keep the negated conjunction: a NaN parameter must fail the API range test.
    if !(uv.x >= lo.x && uv.x <= hi.x && uv.y >= lo.y && uv.y <= hi.y) {
        return Err(SolidError::invalid_input(
            op,
            "parameter outside natural domain",
        ));
    }
    Ok(surface)
}

/// Evaluate a surface point after the API's inclusive natural-domain check.
///
/// # Errors
/// Returns invalid input for a null surface or an out-of-domain (including NaN)
/// parameter, otherwise propagates the surface evaluation error unchanged.
pub fn evaluate_surface_point(surface: Option<&Surface>, uv: Point2) -> Result<Point3> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_point")?.evaluate_point(uv)
}

/// Evaluate the unit normal, with both derivative sides taken from the left.
///
/// # Errors
/// Those of [`evaluate_surface_point`] and the surface's normal evaluation.
pub fn evaluate_surface_normal(surface: Option<&Surface>, uv: Point2) -> Result<Vec3> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_normal")?.evaluate_normal(
        uv,
        Side::Left,
        Side::Left,
    )
}

/// Evaluate `(point, du, dv)`, with both derivative sides taken from the left.
///
/// # Errors
/// Those of [`evaluate_surface_point`] and the surface's first derivatives.
pub fn evaluate_surface_derivatives(
    surface: Option<&Surface>,
    uv: Point2,
) -> Result<(Point3, Vec3, Vec3)> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_derivatives")?
        .evaluate_1st_derivatives(uv, true, true)
}

/// The polygon vertex classification (`SmPolyVertexClass`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum PolyVertexClass {
    /// `SM_PVC_NOT_COMPUTED`.
    NotComputed = 0,
    /// `SM_PVC_UNKNOWN`.
    Unknown = 1,
    /// `SM_PVC_SIMPLE`.
    Simple = 2,
    /// `SM_PVC_BOUNDARY`.
    Boundary = 3,
    /// `SM_PVC_INTERIOR`.
    Interior = 4,
    /// `SM_PVC_CORNER`.
    Corner = 5,
}

/// Decimation bookkeeping (`SmDecimateVertexAttr`).
#[derive(Debug, Clone, PartialEq)]
pub struct DecimateVertexData {
    /// `m_lSortIndex`; reset to zero when copied.
    pub sort_index: u64,
    /// `m_eVertClass`.
    pub vertex_class: PolyVertexClass,
    /// `m_dVertexDistance`.
    pub vertex_distance: f64,
    /// `m_dMaxFaceError`.
    pub max_face_error: f64,
    /// `m_vAveragePlaneNormal`.
    pub average_plane_normal: Vec3,
}

/// The quadratic error form (`SmQuadricVertexAttr`).
#[derive(Debug, Clone, PartialEq)]
pub struct QuadricVertexData {
    /// Three rows of `m_Q`.
    pub matrix: [Vec3; 3],
    /// `m_b`.
    pub linear: Vec3,
    /// `m_c`.
    pub constant: f64,
}

/// Edge contraction bookkeeping (`SmQuadricEdgeAttr`).
#[derive(Debug, Clone, PartialEq)]
pub struct QuadricEdgeData {
    /// `m_vNewPoint`.
    pub new_point: Point3,
    /// `m_dError`.
    pub error: f64,
    /// `m_lSortIndex`: unset until assigned; `MakeCopy` sets it to zero.
    pub sort_index: Option<u64>,
}

/// Construct `SmDecimateVertexAttr` (upstream default behaviour: `Copy`).
#[must_use]
pub fn decimate_vertex_attribute(
    id: u64,
    behavior: crate::geom::attribute::Behavior,
) -> crate::geom::attribute::Attribute {
    use crate::geom::attribute::{Attribute, AttributeValue};
    Attribute::with_value(
        id,
        behavior,
        AttributeValue::DecimateVertex(DecimateVertexData {
            sort_index: 0,
            vertex_class: PolyVertexClass::NotComputed,
            vertex_distance: -crate::math::BIG_DOUBLE,
            max_face_error: 0.0,
            average_plane_normal: Vec3::new(0.0, 0.0, 0.0),
        }),
    )
}

/// Construct `SmQuadricVertexAttr` (upstream default behaviour: `Copy`).
#[must_use]
pub fn quadric_vertex_attribute(
    id: u64,
    behavior: crate::geom::attribute::Behavior,
) -> crate::geom::attribute::Attribute {
    use crate::geom::attribute::{Attribute, AttributeValue};
    let zero = Vec3::new(0.0, 0.0, 0.0);
    Attribute::with_value(
        id,
        behavior,
        AttributeValue::QuadricVertex(QuadricVertexData {
            matrix: [zero; 3],
            linear: zero,
            constant: 0.0,
        }),
    )
}

/// Construct `SmQuadricEdgeAttr` (upstream default behaviour: `Copy`).
#[must_use]
pub fn quadric_edge_attribute(
    id: u64,
    behavior: crate::geom::attribute::Behavior,
) -> crate::geom::attribute::Attribute {
    use crate::geom::attribute::{Attribute, AttributeValue};
    Attribute::with_value(
        id,
        behavior,
        AttributeValue::QuadricEdge(QuadricEdgeData {
            new_point: Vec3::new(0.0, 0.0, 0.0),
            error: 0.0,
            sort_index: None,
        }),
    )
}
