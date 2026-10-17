// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `blend.fillet_solver_p1`.
// Wraps source/SMLib/src/SmFilletSolver.cpp: the SmFilletSolver base (constructor, getters,
// offset-surface bookkeeping, SetupOffsetExtension, FinishConstruction, ComputeSurfaceValues,
// ComputeUVVectors, SplitFilletGeom), SmFilletSphereSolveENFO::Evaluate and
// SmRailUVCurveIntersectENFO::Evaluate.
//
// One entry point, sc_blend_fillet_solver_p1_eval(op, in, n_in, out, cap, out_len). Layouts are
// documented in crates/solidcore-ffi/src/units/blend_fillet_solver_p1.rs, which validates every
// input. SmFilletSolver is abstract: the shim derives a solver that only names itself, so every
// virtual runs its base version. Members upstream leaves uninitialised (edge uses, executive,
// orientations, saved offsets, extension flags) are set by the derived constructor.
#include "sc_shim.h"

#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmContext.h>
#include <SmCurveClass.h>
#include <SmFilletGeom.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmFilletSolver.h>
#include <SmMatrix.h>
#include <SmOffsetSurface.h>
#include <SmPlane.h>
#include <SmSurfaceIntersector.h>
#include <SmTArray.h>
#include <nurbs.h>

#include <vector>

#include "sc_edgeuse_graph.h"

namespace
{
// Ops (keep in sync with crates/solidcore-ffi/src/units/blend_fillet_solver_p1.rs).
enum : int32_t
{
    OP_STATE = 0,
    OP_EXTENSION = 1,
    OP_SPHERE = 2,
    OP_RAIL_UV = 3,
    OP_SURF_VALUES = 4,
    OP_UV_VECTORS = 5,
    OP_SPLIT = 6,
    OP_FINISH = 7,
};

struct ShimFS : SmFilletSolver
{
    ShimFS(const SmContext& ctx, double tol, double ang, double tang)
        : SmFilletSolver(ctx, tol, ang, tang)
    {
        m_pEdgeuses[0] = m_pEdgeuses[1] = NULL;
        m_pExecutive = NULL;
        m_dOrientations[0] = m_dOrientations[1] = 1.0;
        m_dSavedOffsets[0] = m_dSavedOffsets[1] = 0.0;
        m_bExtendBefore = m_bExtendAfter = FALSE;
    }
    virtual SmFilletSolverType GetSolverType() { return SM_FS_UNKNOWN; }
    virtual void GetName(TCHAR* pc, size_t n)
    {
        if (n > 0)
            pc[0] = 0;
    }
    double saved(int i) const { return m_dSavedOffsets[i]; }
    void set_orientation(int i, double o) { m_dOrientations[i] = o; }
    SmStatus finish() { return FinishConstruction(); }
};

void put_pt(std::vector<double>& o, const SmVector3d& v)
{
    o.push_back(v.x);
    o.push_back(v.y);
    o.push_back(v.z);
}

// `nx x[nx] nf f[nf] has_j rows cols` -> arrays and an optional zero-filled matrix.
struct Arrays
{
    SmTArray<double> x;
    SmTArray<double> f;
    SmMatrix* j = NULL;
    ~Arrays() { delete j; }
    void read(Rd& rd)
    {
        int64_t nx = rd.idx();
        for (int64_t i = 0; i < nx; i++)
            x.Add(rd.num());
        int64_t nf = rd.idx();
        for (int64_t i = 0; i < nf; i++)
            f.Add(rd.num());
        bool has_j = rd.flag();
        ULONG rows = (ULONG)rd.idx();
        ULONG cols = (ULONG)rd.idx();
        if (has_j)
        {
            j = new SmMatrix(rows, cols);
            for (ULONG r = 0; r < rows; r++)
                for (ULONG c = 0; c < cols; c++)
                    (*j)[r][c] = 0.0;
        }
    }
    void put(std::vector<double>& o)
    {
        for (ULONG i = 0; i < f.GetSize(); i++)
            o.push_back(f[i]);
        if (j)
            for (ULONG r = 0; r < j->GetNumRows(); r++)
                for (ULONG c = 0; c < j->GetNumColumns(); c++)
                    o.push_back((*j)[r][c]);
    }
};

// Offset surfaces made by the state script, keyed by creation order; `alive` false once
// SetOffsetsData deleted one.
struct Pool
{
    std::vector<SmOffsetSurface*> made;
    std::vector<bool> alive;
    double key(const SmOffsetSurface* s) const
    {
        for (size_t i = 0; i < made.size(); i++)
            if (made[i] == s)
                return (double)i;
        return -1.0;
    }
    SmOffsetSurface* add(SmOffsetSurface* s)
    {
        made.push_back(s);
        alive.push_back(true);
        return s;
    }
};

void put_state(std::vector<double>& o, ShimFS& fs, const Pool& pool)
{
    for (ULONG i = 0; i < 2; i++)
    {
        SmOffsetSurface* s = fs.GetSurface(i);
        o.push_back(flg(s != NULL));
        o.push_back(s ? pool.key(s) : -1.0);
        o.push_back(s ? s->GetOffsetDistance() : 0.0);
        o.push_back(fs.GetOrientation(i));
        o.push_back(fs.saved((int)i));
    }
    o.push_back(fs.GetExtensionFactor());
    o.push_back(flg(fs.GetReverseTrim()));
    o.push_back((double)fs.GetStatus());
}

int32_t run_state(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    double tol = rd.num();
    double ang = rd.num();
    double tang = rd.num();
    SmSurface* base[2] = {NULL, NULL};
    base[0] = read_surface(rd, ctx);
    base[1] = read_surface(rd, ctx);
    ShimFS* fs = new (ctx) ShimFS(ctx, tol, ang, tang);
    Pool pool;
    for (ULONG i = 0; i < 2; i++)
    {
        bool has = rd.flag();
        double off = rd.num();
        if (has)
        {
            SmOffsetSurface* s = pool.add(new (ctx) SmOffsetSurface(off, *base[i]));
            fs->SetSurface(i, s);
        }
    }
    fs->set_orientation(0, rd.num());
    fs->set_orientation(1, rd.num());

    o.push_back(fs->GetThisApproxTol3d());
    o.push_back(fs->GetThisAngTolRad());
    o.push_back(fs->GetTangencyTolerance());
    o.push_back(fs->GetConversionTol());
    o.push_back(fs->GetExtensionFactor());
    o.push_back(flg(fs->GetReverseTrim()));
    o.push_back((double)fs->GetStatus());
    o.push_back((double)fs->GetJacobianSize());
    o.push_back(flg(fs->CanMakeAnalytic()));
    o.push_back(flg(fs->OffsetRadiiCanChange()));
    SmTArray<SmFilletGeom*> geoms;
    fs->GetFilletGeoms(geoms);
    o.push_back((double)geoms.GetSize());
    o.push_back(flg(fs->GetFirstFilletGeom() == fs->GetLastFilletGeom()));
    o.push_back((double)fs->GetType());
    o.push_back(flg(fs->IsKindOf(SmFilletSolver_TYPE)));
    o.push_back(flg(fs->IsKindOf(SmObject_TYPE)));
    o.push_back(flg(fs->IsKindOf(SmFace_TYPE)));
    o.push_back(flg(fs->GetEdgeuse(0) == NULL));
    o.push_back(flg(fs->GetVertex(0) == NULL));
    o.push_back(flg(fs->GetFilletSurfaceGenerator() == NULL));
    o.push_back(flg(fs->GetFilletGeomForSurfaces() == NULL));
    o.push_back(flg(fs->GetFilletExecutive() == NULL));

    SmFilletSolverOffsetsData saved;
    int64_t k = rd.idx();
    for (int64_t c = 0; c < k; c++)
    {
        int64_t cmd = rd.idx();
        double a = rd.num();
        double b = rd.num();
        switch (cmd)
        {
        case 0: fs->SetSurfacesToZeroOffset(); break;
        case 1: fs->ReloadSurfaceOffsets(); break;
        case 2: fs->GetOffsetsData(saved); break;
        case 3:
        {
            // SetOffsetsData deletes the current surface when it differs from the saved one
            for (ULONG i = 0; i < 2; i++)
            {
                SmOffsetSurface* cur = fs->GetSurface(i);
                if (cur != saved.m_pSurfaces[i] && cur != NULL)
                {
                    double key = pool.key(cur);
                    if (key >= 0.0)
                        pool.alive[(size_t)key] = false;
                }
            }
            fs->SetOffsetsData(saved);
            break;
        }
        case 4:
        {
            ULONG rail = (ULONG)a;
            fs->SetSurface(rail, pool.add(new (ctx) SmOffsetSurface(b, *base[rail])));
            break;
        }
        case 5: fs->SetExtensionFactor(a); break;
        case 6: fs->SetReverseTrim(a != 0.0 ? TRUE : FALSE); break;
        case 7: fs->SetStatus((SmFilStatus)(int)a); break;
        default: break;
        }
        put_state(o, *fs, pool);
    }

    // the solver deletes its two surfaces; the shim deletes every other live one
    for (size_t i = 0; i < pool.made.size(); i++)
        if (pool.alive[i] && pool.made[i] != fs->GetSurface(0) && pool.made[i] != fs->GetSurface(1))
            delete pool.made[i];
    delete fs;
    delete base[0];
    delete base[1];
    return 0;
}

int32_t run_extension(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    double ext = rd.num();
    SmSurface* s = read_surface(rd, ctx);
    double u0 = rd.num();
    double v0 = rd.num();
    double u1 = rd.num();
    double v1 = rd.num();
    SmExtent2d dom(SmPoint2d(u0, v0), SmPoint2d(u1, v1));
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    fs.SetExtensionFactor(ext);
    fs.SetupOffsetExtension(s, dom);
    o.push_back(0.0);
    o.push_back(dom.GetMin().x);
    o.push_back(dom.GetMin().y);
    o.push_back(dom.GetMax().x);
    o.push_back(dom.GetMax().y);
    delete s;
    return 0;
}

int32_t run_sphere(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    SmPoint3d c = rd.vec();
    double r = rd.num();
    bool adjusts = rd.flag();
    Arrays a;
    a.read(rd);
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    SmFilletSphereSolveENFO e(c, r, fs);
    e.SetPlaneAdjusts(adjusts ? TRUE : FALSE);
    SmBoolean found = FALSE;
    SmStatus s = e.Evaluate(a.x, a.f, a.j, found);
    o.push_back(st(s));
    o.push_back(flg(found));
    a.put(o);
    return 0;
}

int32_t run_rail_uv(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    ULONG rail = (ULONG)rd.idx();
    SmBSplineCurve* crv = read_curve(rd, ctx);
    double t0 = rd.num();
    double t1 = rd.num();
    Arrays a;
    a.read(rd);
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    SmPlane plane(SmPoint3d(0.0, 0.0, 0.0), SmVector3d(0.0, 0.0, 1.0), &ctx);
    {
        SmRailUVCurveIntersectENFO e(*crv, SmExtent1d(t0, t1), plane, rail, fs);
        SmBoolean found = FALSE;
        SmStatus s = e.Evaluate(a.x, a.f, a.j, found);
        o.push_back(st(s));
        o.push_back(flg(found));
        a.put(o);
    }
    delete crv;
    return 0;
}

int32_t run_surf_values(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    SmPoint2d uv[2];
    uv[0] = rd.uv();
    uv[1] = rd.uv();
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    SmTsectPnt t;
    o.push_back(st(fs.ComputeSurfaceValues(uv, t)));
    return 0;
}

int32_t run_uv_vectors(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    SmTsectPnt t;
    for (ULONG i = 0; i < 2; i++)
    {
        t.SrfDu(i) = rd.vec();
        t.SrfDv(i) = rd.vec();
    }
    t.CrvDeriv() = rd.vec();
    for (ULONG i = 0; i < 2; i++)
        t.UVDeriv(i) = rd.uv();
    o.push_back(st(fs.ComputeUVVectors(t)));
    for (ULONG i = 0; i < 2; i++)
    {
        o.push_back(t.UVDeriv(i).x);
        o.push_back(t.UVDeriv(i).y);
    }
    return 0;
}

int32_t run_split(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    ULONG idx = (ULONG)rd.idx();
    ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
    SmTArray<SmExtent1d> i1, i2;
    SmCurveClassification cc1, cc2;
    o.push_back(st(fs.SplitFilletGeom(idx, (SmBoundaryTrimmingType)0, i1, i2, cc1, cc2)));
    SmTArray<SmFilletGeom*> geoms;
    fs.GetFilletGeoms(geoms);
    o.push_back((double)geoms.GetSize());
    return 0;
}

int32_t run_finish(Rd& rd, const SmContext& ctx, std::vector<double>& o)
{
    SmSurface* base = read_surface(rd, ctx);
    bool has0 = rd.flag();
    bool has1 = rd.flag();
    double off = rd.num();
    {
        ShimFS fs(ctx, 1.0e-6, 0.1, 0.01);
        if (has0)
            fs.SetSurface(0, new (ctx) SmOffsetSurface(off, *base));
        if (has1)
            fs.SetSurface(1, new (ctx) SmOffsetSurface(-off, *base));
        o.push_back(st(fs.finish()));
    } // the solver deletes its fillet geometries and surfaces before the base goes
    delete base;
    return 0;
}
} // namespace

SC_EXPORT int32_t sc_blend_fillet_solver_p1_eval(int32_t op, const double* in, int64_t n_in, double* out,
                                                 int64_t cap, int64_t* out_len)
{
    if (in == NULL || out_len == NULL)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = the_context();
    Rd rd{in, 0};
    std::vector<double> o;
    int32_t rc = SC_ERR_ARG;
    switch (op)
    {
    case OP_STATE: rc = run_state(rd, ctx, o); break;
    case OP_EXTENSION: rc = run_extension(rd, ctx, o); break;
    case OP_SPHERE: rc = run_sphere(rd, ctx, o); break;
    case OP_RAIL_UV: rc = run_rail_uv(rd, ctx, o); break;
    case OP_SURF_VALUES: rc = run_surf_values(rd, ctx, o); break;
    case OP_UV_VECTORS: rc = run_uv_vectors(rd, ctx, o); break;
    case OP_SPLIT: rc = run_split(rd, ctx, o); break;
    case OP_FINISH: rc = run_finish(rd, ctx, o); break;
    default: break;
    }
    (void)n_in;
    if (rc != 0)
        return rc;
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
