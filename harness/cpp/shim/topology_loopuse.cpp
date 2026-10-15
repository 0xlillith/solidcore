// SPDX-License-Identifier: Apache-2.0
// Real SmLoop/SmLoopuse queries on the shared native B-rep fixture.
#include "sc_edgeuse_graph.h"
#include <algorithm>

namespace {
template<class T> double loop_identity(const std::vector<T*>& objects,const T* object) {
    if(!object) return -1.0;
    auto found=std::find(objects.begin(),objects.end(),object);
    if(found==objects.end()) throw 1;
    return (double)(found-objects.begin());
}
template<class T> void loop_list(std::vector<double>& out,const SmTArray<T*>& list,const std::vector<T*>& ids) {
    out.push_back((double)list.GetSize());
    for(ULONG i=0;i<list.GetSize();++i) out.push_back(loop_identity(ids,list[i]));
}
void loop_text(std::vector<double>& out,const TCHAR* value) {
    std::vector<double> chars;
    for(const TCHAR* c=value;c && *c;++c) chars.push_back((double)*c);
    out.push_back((double)chars.size());out.insert(out.end(),chars.begin(),chars.end());
}
}

SC_EXPORT int32_t sc_topology_loopuse_query(const char* path,const double* input,int64_t count,
    int32_t orientation,double* output,int64_t capacity,int64_t* output_len) {
    if(!path || !input || count<1 || !output_len || capacity<0 || (capacity && !output)
        || orientation < -1 || orientation > 2) return SC_ERR_ARG;
    SC_TRY
    Rd reader{input,0};SmStatus status=SM_SUCCESS;
    Built built;built.brep=build(reader,path,status);
    if(!built.brep || status!=SM_SUCCESS) return SC_ERR_ARG;
    Graph graph;walk(built.brep,graph);
    std::vector<double> out{(double)graph.loopuses.size(),(double)graph.loops.size()};
    for(SmLoopuse* lu:graph.loopuses) {
        if(orientation>=0) lu->SetOrientation((SmOrientType)orientation);
        out.insert(out.end(),{(double)lu->GetType(),(double)lu->GetLoopuseType(),
            (double)lu->GetOrientation(),(double)lu->IsEdgeLoopuse(),(double)lu->IsVertexLoopuse(),(double)lu->IsClosedPtrs(),
            loop_identity(graph.faceuses,lu->GetFaceuse()),loop_identity(graph.loops,lu->GetLoop()),
            loop_identity(graph.loopuses,lu->GetOtherLoopuse()),loop_identity(graph.shells,lu->GetShell()),
            loop_identity(graph.vertexuses,lu->GetVertexuse())});
        for(SM_TYPE type:{16003L,16004L,16001L,16010L,12060L,11001L,16999L,-1L}) out.push_back((double)lu->IsKindOf(type));
        loop_text(out,lu->GetClassString());loop_text(out,lu->GetTypeString());
        SmTArray<SmEdgeuse*> eus;lu->GetEdgeuses(eus);loop_list(out,eus,graph.edgeuses);
        SmTArray<SmEdge*> edges;lu->GetEdges(edges);loop_list(out,edges,graph.edges);
        SmTArray<SmVertexuse*> vus;lu->GetVertexuses(vus);loop_list(out,vus,graph.vertexuses);
        for(ULONG i=0;i<vus.GetSize();++i) if(!vus[i]) return SC_ERR_ARG;
        SmTArray<SmVertex*> vertices;lu->GetVertices(vertices);loop_list(out,vertices,graph.vertices);
    }
    for(SmLoop* loop:graph.loops) {
        SmLoopuse* first=loop->GetLoopuse();
        if(!first || !first->GetFaceuse()) return SC_ERR_ARG;
        SmTArray<SmLoopuse*> members;first->GetFaceuse()->GetLoopuses(members);
        if(members.GetSize()==0) return SC_ERR_ARG;
        SmLoopuse *primary=nullptr,*other=nullptr;loop->GetLoopuses(primary,other);
        out.insert(out.end(),{(double)loop->GetType(),loop_identity(graph.loopuses,first),
            loop_identity(graph.loopuses,primary),loop_identity(graph.loopuses,other),
            loop_identity(graph.faces,loop->GetFace()),(double)loop->IsOuterLoop(),(double)loop->IsClosedPtrs()});
        for(SM_TYPE type:{16003L,16004L,16001L,16010L,12060L,11001L,16999L,-1L}) out.push_back((double)loop->IsKindOf(type));
        loop_text(out,loop->GetClassString());loop_text(out,loop->GetTypeString());
        SmTArray<SmEdgeuse*> eus;loop->GetEdgeuses(eus);loop_list(out,eus,graph.edgeuses);
        SmTArray<SmEdge*> edges;loop->GetEdges(edges);loop_list(out,edges,graph.edges);
        SmTArray<SmVertex*> vertices;loop->GetVertices(vertices);loop_list(out,vertices,graph.vertices);
    }
    return sc_copy_out(out.data(),(int64_t)out.size(),output,capacity,output_len);
    SC_CATCH
}
