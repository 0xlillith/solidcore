// SPDX-License-Identifier: Apache-2.0
// query.face_p1: calls the actual SmApiEvaluateCurve entry point.
// Upstream: source/SM_API/src/SmApiCurves.cpp:1043-1080.
#include "sc_shim.h"
#include "StdAfx.h"
#include <SmApiCurves.h>
#include <SmContext.h>
#include <SmCurve.h>
#include <SmLine.h>
#include <SmHermiteCurve.h>
#include <memory>

SC_EXPORT int32_t sc_query_face_p1_curve(int32_t kind, const double* data,
                                        double parameter, int32_t mask,
                                        double* out, int32_t* operation_status)
{
    if (!data || !out || !operation_status || kind < 0 || kind > 2 || mask < 0 || mask > 7)
        return SC_ERR_ARG;
    SC_TRY
    SmContext context;
    auto vec = [&](int i) { return SmVector3d(data[i], data[i+1], data[i+2]); };
    std::unique_ptr<SmCurve> curve;
    if (kind == 1) curve.reset(new (context) SmLine(vec(0), vec(6), 3, &context));
    if (kind == 2) curve.reset(new (context) SmHermiteCurve(vec(0), vec(3), vec(6), vec(9), 3));
    SmVector3d point(11, 13, 17), d1(19, 23, 29), d2(31, 37, 41);
    SmVector3d* p = (mask & 1) ? &point : nullptr;
    SmVector3d* q = (mask & 2) ? &d1 : nullptr;
    SmVector3d* r = (mask & 4) ? &d2 : nullptr;
    *operation_status = SmApiEvaluateCurve(curve.get(), parameter, p, q, r);
    const SmVector3d values[3] = {point, d1, d2};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out[i*3+j] = values[i][j];
    return 0;
    SC_CATCH
}

// Fixture setup/export uses the shared native graph protocol. The operation
// under test remains the actual upstream method, including default arguments.
#include "sc_edgeuse_graph.h"
SC_EXPORT int32_t sc_query_face_p1_box(int32_t target, int64_t index,
    int32_t tight, int32_t use_default, const char* path, const double* input,
    int64_t n_in, double* out, int64_t cap, int64_t* out_len)
{
    if (target < 0 || target > 2 || index < 0 || !input || n_in < 1)
        return SC_ERR_ARG;
    SC_TRY
    Rd rd{input, 0}; Built built; SmStatus status=SM_SUCCESS;
    built.brep=build(rd,path,status);
    if (!built.brep || status!=SM_SUCCESS || rd.at!=n_in) return SC_ERR_ARG;
    Graph graph; walk(built.brep,graph);
    SmExtent3d box;box.SetMinMax(SmPoint3d(-17,-19,-23),SmPoint3d(29,31,37));
    if (target == 0) {
        if ((uint64_t)index >= graph.faces.size()) return SC_ERR_ARG;
        SmFace* face=graph.faces[(size_t)index];
        if (!face->GetSurface()) return SC_ERR_ARG;
        status=use_default ? face->CalculateBoundingBox(box) : face->CalculateBoundingBox(box,tight);
    } else if (target == 1) {
        if ((uint64_t)index >= graph.shells.size()) return SC_ERR_ARG;
        SmShell* shell=graph.shells[(size_t)index];
        if (shell->IsVertexShell() && !shell->GetVertex()) return SC_ERR_ARG;
        status=use_default ? shell->CalculateBoundingBox(box) : shell->CalculateBoundingBox(box,tight);
    } else {
        if ((uint64_t)index >= graph.regions.size()) return SC_ERR_ARG;
        SmRegion* region=graph.regions[(size_t)index];
        status=use_default ? region->CalculateBoundingBox(box) : region->CalculateBoundingBox(box,tight);
    }
    std::vector<double> result{st(status)};
    put_vec(result,box.GetMin());put_vec(result,box.GetMax());
    return sc_copy_out(result.data(),(int64_t)result.size(),out,cap,out_len);
    SC_CATCH
}
