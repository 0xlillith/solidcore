// SPDX-License-Identifier: Apache-2.0
// Real ClassifyFaceInterval probes. The legacy stand-in script guard is unchanged.
#include "sc_topology_fixture.h"
#include "sc_brep_query_metadata.h"
#include <SmCacheMgr.h>
#include <SmCurveClass.h>
#include <SmPointClass.h>
#include <SmTol.h>
#include <stdexcept>

namespace {
void probe_object(std::vector<double>& out,const Graph& g,const void* object) {
    if (!object) { out.insert(out.end(),{-1,-1}); return; }
    for (int kind=0;kind<10;kind++) {
        auto found=g.index[kind].find(object);
        if (found!=g.index[kind].end()) { out.push_back(kind); out.push_back(found->second); return; }
    }
    throw std::runtime_error("probe PC object not in real graph");
}
void probe_pc(std::vector<double>& out,const Graph& g,const SmPointClassification& pc) {
    out.push_back(pc.GetPointClass()); probe_object(out,g,pc.GetObject());
    out.push_back(g.find(3,pc.GetFace()));
    out.push_back(pc.GetDeviation()); out.push_back(pc.GetTParam());
    const auto uv=pc.GetUVParam(); out.insert(out.end(),{uv.x,uv.y});
    const auto uvw=pc.GetUVWParam(); out.insert(out.end(),{uvw.x,uvw.y,uvw.z});
    out.push_back(pc.AreParametersSet()?1:0); out.push_back(pc.GetPreSnapParam());
    out.push_back(pc.GetSSSPointFlag()?1:0);
    const auto sss=pc.GetSSSPoint(); out.insert(out.end(),{sss.x,sss.y,sss.z});
    out.push_back(SmTol::GetSrcZoneTol3d(&pc));
}
void probe_pcs(std::vector<double>& out,const Graph& g,const SmCurveInterval& ivl) {
    probe_pc(out,g,ivl.m_vStart); probe_pc(out,g,ivl.m_vMid); probe_pc(out,g,ivl.m_vEnd);
}
}

SC_EXPORT int32_t sc_intersect_curve_class_p7_dispatch(int32_t branch,int32_t export_only,
    double* out,int64_t cap,int64_t* len) {
    if (!len || branch<0 || branch>2 || export_only<0 || export_only>1) return SC_ERR_ARG;
    SC_TRY
    // Fresh real finite XY plane; shared constructor determines its actual UV axes.
    const double block[]={2,1,0,0,0,0,0,1,-2,-1,3,4,0,0};
    Rd rd{block,0}; Built built; SmStatus status;
    built.brep=build(rd,"",status);
    if (!built.brep || status!=SM_SUCCESS || rd.at!=static_cast<int64_t>(sizeof(block)/sizeof(block[0]))) return SC_ERR_ARG;
    SmBrep* b=built.brep; Graph g; walk(b,g,true);
    if (g.faces.size()!=1 || g.edges.empty()) return SC_ERR_ARG;
    SmFace* face=g.faces.front(); SmSurface* surface=face->GetSurface();
    if (!surface || !surface->IsKindOf(SmPlane_TYPE) || surface->GetFace()!=face) return SC_ERR_ARG;
    // Check actual cache absence, not an assumption based on a cold fixture name.
    if (SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF,surface)) return SC_ERR_ARG;
    const SmContext* ctx=face->GetContext(); if (!ctx) return SC_ERR_ARG;
    SmContext* mutable_ctx=const_cast<SmContext*>(ctx);
    if (mutable_ctx->GetUserNotifyCallback() || mutable_ctx->GetSysNotifyCallback()) return SC_ERR_ARG;
    const double zone=SmTol::GetZoneTol3d(ctx);
    const SmExtent2d domain=face->GetUVDomain();
    SmPoint2d uv0=domain.Evaluate(0.3,0.3),uv1=domain.Evaluate(0.7,0.7);
    SmPoint3d p,q; SmEdge* edge=nullptr; double edge_t=0;
    if (branch==0) {
        // A genuine boundary edge supplies a start classified to Edge. The end
        // remains Unknown; native step 1 calls EdgeCurveOnClassify before step 4.
        edge=g.edges.front(); const auto interval=edge->GetInterval(); edge_t=interval.Evaluate(0.5);
        if (!edge->GetCurve() || edge->GetCurve()->EvaluatePoint(edge_t,p)!=SM_SUCCESS) return SC_ERR_ARG;
        uv1=domain.Evaluate(0.5,0.5);
        if (surface->EvaluatePoint(uv1,q)!=SM_SUCCESS) return SC_ERR_ARG;
    } else {
        if (surface->EvaluatePoint(uv0,p)!=SM_SUCCESS || surface->EvaluatePoint(uv1,q)!=SM_SUCCESS) return SC_ERR_ARG;
    }
    if (p.DistanceBetween(q)<=zone*100) return SC_ERR_ARG;
    SmLine curve(p,q,3,ctx),uvcurve(SmPoint3d(uv0.x,uv0.y,0),SmPoint3d(uv1.x,uv1.y,0),2,ctx);
    const SmExtent1d interval(0,1);
    if (curve.EditParameterization(interval,FALSE)!=SM_SUCCESS || uvcurve.EditParameterization(interval,FALSE)!=SM_SUCCESS) return SC_ERR_ARG;
    SmCurveClassification cc; cc.SetCurve(&curve,FALSE); cc.SetSurface1(surface);
    cc.SetUVCurve1(branch==2?&uvcurve:nullptr); SmTol::SetZoneTol3d(&cc,SmZoneTol3d(zone));
    cc.SetSize(1);
    SmCurveInterval& ivl=cc.GetAt(0); ivl.SetInterval(interval);
    for (SmPointClassification* pc:{&ivl.m_vStart,&ivl.m_vMid,&ivl.m_vEnd})
        *pc=SmPointClassification(SmZoneTol3d(zone),ctx);
    if (branch==0) { ivl.m_vStart.SetClassObject(SM_PC_EDGE,edge); ivl.m_vStart.SetTParam(edge_t); }
    if (ivl.m_vMid.GetPointClass()!=SM_PC_UNKNOWN || ivl.m_vEnd.GetPointClass()!=SM_PC_UNKNOWN ||
        (branch!=0 && ivl.m_vStart.GetPointClass()!=SM_PC_UNKNOWN)) return SC_ERR_ARG;
    std::vector<double> pre; probe_pcs(pre,g,ivl);
    if (export_only) {
        std::vector<double> graph={0}; export_graph(b,g,graph);
        std::vector<double> metadata; append_query_metadata(b,g,metadata);
        std::vector<double> result={(double)graph.size()}; result.insert(result.end(),graph.begin(),graph.end());
        result.push_back((double)metadata.size()); result.insert(result.end(),metadata.begin(),metadata.end());
        result.insert(result.end(),{(double)g.find(3,face),(double)g.find(7,edge),edge_t,
            zone,SmTol::GetXSectTol3d(ctx),SmTol::GetApproxTol3d(ctx),SM_ANG_TOL_DEG,
            SmTol::GetGapGain(),SmTol::GetModelSizeEstimate(),SmTol::GetLargeSmallSizeRatio(),
            p.x,p.y,p.z,q.x,q.y,q.z,uv0.x,uv0.y,uv1.x,uv1.y});
        result.insert(result.end(),pre.begin(),pre.end());
        return sc_copy_out(result.data(),result.size(),out,cap,len);
    }
    // All three cases call the same p7 operation, not its callees in isolation.
    status=cc.ClassifyFaceInterval(0,TRUE,face,nullptr,TRUE);
    if (status==SM_SUCCESS && (ivl.m_vMid.GetPointClass()!=SM_PC_FACE || ivl.m_vMid.GetObject()!=face)) return SC_ERR_ARG;
    if (SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF,surface)) return SC_ERR_ARG;
    Graph after; walk(b,after,true);
    std::vector<double> result={(double)status,ivl.GetMin(),ivl.GetMax()};
    probe_pcs(result,after,ivl); result.push_back(0); export_graph(b,after,result);
    return sc_copy_out(result.data(),result.size(),out,cap,len);
    SC_CATCH
}
