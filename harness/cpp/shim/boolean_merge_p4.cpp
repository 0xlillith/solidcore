// SPDX-License-Identifier: Apache-2.0
// Real ManifoldBoolean entry point. Input graphs are exported before mutation;
// every post-operation snapshot walks surviving topology afresh.
#include "sc_topology_fixture.h"
#include "sc_brep_query_metadata.h"
#include <SmMerge.h>
#include <SmCacheMgr.h>
#include <SmSurfaceCache.h>
#include <memory>

namespace {
void notification_word(std::vector<double>& out, uint64_t value) {
    out.push_back((double)(value & UINT64_C(0xffffffff)));
    out.push_back((double)(value >> 32));
}
// The same default SmContext, with read-only access to its protected raw tolerance
// fields, as in geom_context.cpp. No context state is reconstructed or normalized.
struct InputContext : SmContext {
    void snapshot(std::vector<double>& out) const {
        out.push_back(1.0);
        for (auto mark : {SM_MT_MARK,SM_MT_MARK2,SM_MT_MARK3,SM_MT_MARKIO,SM_MT_MARKASSERT})
            notification_word(out,(uint64_t)GetCurrentMark(mark));
        for (auto mark : {SM_MT_MARK,SM_MT_MARK2,SM_MT_MARK3,SM_MT_MARKIO,SM_MT_MARKASSERT})
            out.push_back(IsMarkLocked(mark) ? 1.0 : 0.0);
        out.insert(out.end(),{GetDoingBoolean() ? 1.0 : 0.0,
            m_dThisModelSizeEstimate,m_dThisLargeSmallSizeRatio,(double)m_sThisZoneTol3d});
    }
};
template<class T> void notification_record(std::vector<double>& out, uint32_t kind, size_t index, T* object) {
    out.insert(out.end(),{(double)kind,(double)index});
    notification_word(out,(uint64_t)object->GetAllFlags());
    for (auto mark : {SM_MT_MARK,SM_MT_MARK2,SM_MT_MARK3,SM_MT_MARKIO,SM_MT_MARKASSERT})
        notification_word(out,(uint64_t)object->GetMark(mark));
    notification_word(out,(uint64_t)object->GetUserIndex1());
    notification_word(out,(uint64_t)object->GetUserIndex2());
    notification_word(out,(uint64_t)(uintptr_t)object->GetUserPtr1());
}
void notification_inputs(SmBrep* brep, const Graph& g, std::vector<double>& out) {
    bool attributed = brep->HasAttributes();
    std::vector<double> records;
    notification_record(records,0,0,brep);
    auto append = [&](uint32_t kind, const auto& objects) {
        for (size_t i=0;i<objects.size();++i) {
            notification_record(records,kind,i,objects[i]);
            if (objects[i]->HasAttributes()) attributed = true;
        }
    };
    append(1,g.regions); append(2,g.shells); append(3,g.faceuses); append(4,g.faces);
    append(5,g.loopuses); append(6,g.loops); append(7,g.edgeuses); append(8,g.edges);
    append(9,g.vertexuses); append(10,g.vertices);
    for (auto* edge : g.edges) if (edge->GetCurve() && edge->GetCurve()->HasAttributes()) attributed=true;
    for (auto* face : g.faces) if (face->GetSurface() && face->GetSurface()->HasAttributes()) attributed=true;
    for (auto* use : g.edgeuses) if (use->GetUVTrimCurvePointer() && use->GetUVTrimCurvePointer()->HasAttributes()) attributed=true;
    SmContext* context=const_cast<SmContext*>(brep->GetContext());
    out.insert(out.end(),{1.0,context && (context->GetUserNotifyCallback() || context->GetSysNotifyCallback()) ? 1.0:0.0,
        attributed ? 1.0:0.0,(double)(records.size()/20)});
    out.insert(out.end(),records.begin(),records.end());
}
SmBrep* make_operand(const SmContext& ctx, const double* input, SmStatus& status) {
    int kind = (int)input[0];
    status = SM_SUCCESS;
    if (kind == -1) return nullptr;
    auto* b = new (ctx) SmBrep(SM_USE_DEFAULT, SM_USE_DEFAULT, (SmContext*)&ctx);
    SmPoint3d p(input[1],input[2],input[3]);
    if (kind == 1) {
        SmShell* s = nullptr; SmVertex* v = nullptr;
        status = b->MakeShellVertex(b->GetInfiniteRegion(),p,s,v);
    } else if (kind == 2) {
        SmBSplineCurve* c = nullptr;
        status = SmBSplineCurve::CreateLineSegment(ctx,3,p,p+SmVector3d(input[4],input[5],input[6]),c);
        if (status == SM_SUCCESS) { SmEdge* e = nullptr; status = b->CreateWireEdgeFromCurve(c,c->GetNaturalInterval(),e); }
    } else if (kind == 3) {
        auto delete_unowned = [](SmPlane* s) { if (s && !s->GetOwner()) delete s; };
        std::unique_ptr<SmPlane,decltype(delete_unowned)> s(
            new (ctx) SmPlane(p,SmVector3d(0,0,1),&ctx), delete_unowned);
        SmFace* f = nullptr;
        status = b->CreateFaceFromSurface(s.get(),SmExtent2d(SmPoint2d(0,0),SmPoint2d(input[4],input[5])),f);
        // A UV-domain rejection happens before the surface gets a face owner.
        // Delete only unconsumed geometry; attached surfaces belong to the B-rep.
    }
    if (status != SM_SUCCESS) { delete b; return nullptr; }
    return b;
}

struct Access : SmMerge {
    Access(const SmContext& c, SmBrep* a, SmBrep* b, const int32_t* flags, bool manifold) : SmMerge(c) {
        m_vTI.SetPrimaryBrep(a); m_vTI.SetOtherBrep(b);
        m_vTI.SetThisApproxTol3d(1e-9); m_vTI.SetThisAngTolRad(0.3);
        m_bManifoldBoolean = manifold ? TRUE : FALSE;
        m_bIntersectionLoopClosed = TRUE;
        SetCookieCutter(flags[0] != 0); SetImprinting(flags[1] != 0);
        SetBooleanPostProcess(flags[2] != 0); SetClassifyFaces(flags[3] != 0); SetKeepOtherBrep(flags[4] != 0);
    }
    ~Access() {
        // Read current operand pointers: ManifoldBoolean may have deleted either.
        delete m_vTI.GetPrimaryBrep(); delete m_vTI.GetOtherBrep();
        for (ULONG i=0;i<m_vOneBrepPer_OrigBrepConnectedFaceSet.GetSize();++i) delete m_vOneBrepPer_OrigBrepConnectedFaceSet[i];
        for (ULONG i=0;i<m_vOneBrepPer_OtherBrepConnectedFaceSet.GetSize();++i) delete m_vOneBrepPer_OtherBrepConnectedFaceSet[i];
        m_vTI.SetPrimaryBrep(nullptr); m_vTI.SetOtherBrep(nullptr);
    }
    void snapshot(std::vector<double>& out, SmBrep* result, const SmContext& a, const SmContext& b) {
        out.insert(out.end(),{result != nullptr ? 1.0:0.0, flg(a.GetDoingBoolean()),flg(b.GetDoingBoolean()),flg(m_bIntersectionLoopClosed),
            (double)m_sIntersectionEdges.GetSize(),(double)m_sIntersectionVertices.GetSize(),
            (double)m_vOneBrepPer_OrigBrepConnectedFaceSet.GetSize(),(double)m_vOneBrepPer_OtherBrepConnectedFaceSet.GetSize()});
        for (auto* brep : {m_vTI.GetPrimaryBrep(),m_vTI.GetOtherBrep()}) {
            out.push_back(brep ? 1.0:0.0);
            if (!brep) continue;
            Graph g; walk(brep,g,true);
            std::vector<double> record{flg(brep->m_bEditingEnabled)};
            export_graph(brep,g,record);
            out.push_back((double)record.size());
            out.insert(out.end(),record.begin(),record.end());
        }
        // Counts alone cannot certify extraction after the operands are deleted.
        // Prefix each complete fresh graph with its record length for decoding.
        for (auto* components : {&m_vOneBrepPer_OrigBrepConnectedFaceSet,
                                 &m_vOneBrepPer_OtherBrepConnectedFaceSet}) {
            for (ULONG i=0; i<components->GetSize(); ++i) {
                SmBrep* brep = (*components)[i];
                std::vector<double> component;
                if (brep) {
                    Graph g; walk(brep,g,true);
                    component.push_back(flg(brep->m_bEditingEnabled));
                    export_graph(brep,g,component);
                }
                out.push_back((double)component.size());
                out.insert(out.end(),component.begin(),component.end());
            }
        }
    }
};
}

SC_EXPORT int32_t sc_boolean_merge_p4_eval(int32_t mode, const double* a, const double* b,
    int32_t operation, const int32_t* flags, int32_t manifold, int32_t distinct_contexts,
    int32_t caller_lock, double* out, int64_t cap, int64_t* len) {
    if (!a || !b || !flags || !len || (mode != -2 && mode != -1 && mode != 0)) return SC_ERR_ARG;
    SC_TRY
    InputContext ca, cb;
    const InputContext& other_context = distinct_contexts ? cb : ca;
    SmStatus status;
    std::unique_ptr<SmBrep> aa(make_operand(ca,a,status));
    if (a[0] != -1 && (status != SM_SUCCESS || !aa)) {
        double rejected[] = {-1.0}; return sc_copy_out(rejected,1,out,cap,len);
    }
    if (mode == -1) {
        std::vector<double> output{0.0};
        if (aa) { Graph g; walk(aa.get(),g,true); export_graph(aa.get(),g,output); }
        return sc_copy_out(output.data(),(int64_t)output.size(),out,cap,len);
    }
    std::unique_ptr<SmBrep> bb(make_operand(other_context,b,status));
    if (b[0] != -1 && (status != SM_SUCCESS || !bb)) {
        double rejected[] = {-1.0}; return sc_copy_out(rejected,1,out,cap,len);
    }
    if (mode == -2) {
        // Construct both original operands in the same contexts and order as eval.
        // Graph and metadata share one traversal; no Boolean output enters Rust.
        std::vector<double> output{2.0};
        // Capture both live contexts before any exporter or Access construction.
        // Null operands still have contexts; a shared pair exports the same state.
        for (const InputContext* context : {static_cast<const InputContext*>(&ca), &other_context}) {
            std::vector<double> state;
            context->snapshot(state);
            output.push_back((double)state.size());
            output.insert(output.end(),state.begin(),state.end());
        }
        for (auto* brep : {aa.get(), bb.get()}) {
            output.push_back(brep ? 1.0 : 0.0);
            if (!brep) continue;
            Graph g; walk(brep,g,true);
            std::vector<double> graph{0.0}, metadata, notifications, caches{1.0,(double)g.faces.size()};
            for (auto* face : g.faces) {
                auto* surface=face->GetSurface();
                auto* cache=surface ? static_cast<SmSurfaceCache*>(SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF,surface)) : nullptr;
                caches.insert(caches.end(),{cache ? 1.0 : 0.0,cache && cache->m_bFaceWasModified ? 1.0 : 0.0});
            }
            notification_inputs(brep,g,notifications);
            export_graph(brep,g,graph); append_query_metadata(brep,g,metadata);
            std::vector<double> record{(double)graph.size()};
            record.insert(record.end(),graph.begin(),graph.end());
            record.push_back((double)metadata.size());
            record.insert(record.end(),metadata.begin(),metadata.end());
            record.push_back((double)notifications.size());
            record.insert(record.end(),notifications.begin(),notifications.end());
            record.push_back((double)caches.size());
            record.insert(record.end(),caches.begin(),caches.end());
            output.push_back((double)record.size());
            output.insert(output.end(),record.begin(),record.end());
        }
        return sc_copy_out(output.data(),(int64_t)output.size(),out,cap,len);
    }
    // Mark lock is destroyed before contexts, after the merge; each invocation
    // owns all its state, including native contexts (no process-global fixtures).
    SmNewMarkAndLock lock;
    Access merge(ca,aa.release(),bb.release(),flags,manifold != 0);
    SmBrep* result = nullptr;
    status = merge.ManifoldBoolean((SmBooleanOperationType)operation,result,caller_lock ? &lock : nullptr);
    std::vector<double> output{st(status)};
    merge.snapshot(output,result,ca,other_context);
    return sc_copy_out(output.data(),(int64_t)output.size(),out,cap,len);
    SC_CATCH
}
