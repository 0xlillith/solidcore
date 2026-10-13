// SPDX-License-Identifier: Apache-2.0
// Unit surface.surface_p9. Calls the pinned SmSurface methods and the file-local
// SmCSIGlobalSolver class (declaration reproduced from SmSurface.cpp:13781).
#include "sc_shim.h"
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmContext.h>
#include <SmCurveCache.h>
#include <SmSurfaceCache.h>
#include <SmGlobalSolver.h>
#include <SmSolutionArray.h>
#include <SmTree.h>
#include <nurbs.h>
#include <vector>

class SmCSIGlobalSolver : public SmGlobalSolver {
public:
    const SmSurface* m_pSurf;
    SmBoolean m_bFoundCoincidence;
    const SmExtent1d& m_crInterval;
    const SmExtent2d& m_crUVDomain;
    SmCSIGlobalSolver(SmTree*, SmTree*, SmBoolean, const SmExtent1d&, const SmExtent2d&);
    virtual SmStatus LocalSolve(SmTreeNode* [SM_GS_MAX_TREES], SmBoolean&);
    virtual SmSolutionPairType IdenticalSolutions(const SmSolution&, const SmSolution&, const SmSolution*&) const;
};

SC_EXPORT int32_t sc_surface_surface_p9_eval(int32_t op, const double* xyz, const double* bounds,
    int32_t direction, double chord, double angle, double max_side, int32_t diagonals,
    const double* line, double* output, int64_t cap, int64_t* out_len) {
    SC_TRY
    SmContext context;
    gw_CPOINT points[3][3]; gw_CPOINT* rows[3] = {points[0],points[1],points[2]};
    for (int i=0;i<3;i++) for(int j=0;j<3;j++) {
        auto& p=points[i][j]; const int at=3*(3*i+j);
        p.x=xyz[at];p.y=xyz[at+1];p.z=xyz[at+2];p.w=NL_NOW;
    }
    double knots[6]={0,0,0,1,1,1};
    SmBSplineSurface* surface=new(context) SmBSplineSurface(2,2,2,2,5,5,rows,knots,knots);
    SmExtent2d domain(bounds[0],bounds[1],bounds[2],bounds[3]);
    std::vector<double> out;
    auto status=[](SmStatus r){return r==SM_SUCCESS ? 0.0 : (double)r;};
    if(op==0) {
        SmTArray<double> splits;
        SmStatus r=surface->FindTessellationSplits((SmSurfParamType)direction,domain,chord,angle,max_side,splits,diagonals?TRUE:FALSE);
        out.push_back(status(r));
        if(r==SM_SUCCESS) for(ULONG i=0;i<splits.GetSize();i++) out.push_back(splits[i]);
    } else if(op==1) {
        SmBSplineCurve* a=NULL;SmBSplineCurve* b=NULL;
        SmStatus r=surface->CreateIsoBoundaries(context,(SmSurfParamType)direction,SmApproxTol3d(chord),a,b);
        out.push_back(status(r));
        if(r==SM_SUCCESS) {
            out.push_back((a?1.0:0.0)+(b?1.0:0.0));
            for(auto curve:{a,b}) if(curve) for(int i=0;i<=4;i++) {
                SmPoint3d p;
                SmStatus er=curve->EvaluatePoint(curve->GetNaturalInterval().Evaluate((double)i/4.0),p);
                if(er!=SM_SUCCESS) { delete a;delete b;delete surface;return SC_ERR_ARG; }
                out.push_back(p.x);out.push_back(p.y);out.push_back(p.z);
            }
        }
        delete a;delete b;
    } else if(op==2) {
        gw_CPOINT cp[2]={{line[0],line[1],line[2],NL_NOW},{line[3],line[4],line[5],NL_NOW}};
        double ck[4]={0,0,1,1};
        SmBSplineCurve* curve=new(context) SmBSplineCurve(3,1,3,1,cp,ck);
        {
            SmExtent1d interval(0.0,1.0);
            SmCurveCache cc(*curve,interval,TRUE);
            SmSurfaceCache sc(*surface,surface->GetNaturalUVDomain());
            SmStatus r=cc.Tessellate();
            if(r==SM_SUCCESS) r=sc.BuildTree(NULL);
            if(r==SM_SUCCESS) {
                SmCSIGlobalSolver solver(cc.GetTree(),sc.GetTree(),FALSE,interval,domain);
                SmSolutionArray solutions;solver.SetSolutions(&solutions);
                solver.SetSolverOperation(SM_SO_INTERSECT);solver.Set3DTolerance(1e-6);solver.SetSolutionRequested(SM_SR_ALL);
                SmTreeNode* branch[SM_GS_MAX_TREES]={cc.GetTree()->GetFirstLeafNode(),sc.GetTree()->GetFirstLeafNode(),NULL};
                SmBoolean more=FALSE;r=solver.LocalSolve(branch,more);
                out.push_back(status(r));
                if(r==SM_SUCCESS) {
                    out.push_back(more==TRUE?1.0:0.0);out.push_back((double)solutions.GetSize());
                    for(ULONG i=0;i<solutions.GetSize();i++) {
                        const SmSolution& s=solutions[i];
                        out.push_back((double)s.m_eSolutionType);out.push_back((double)s.m_lNumVariables);out.push_back((double)s.m_lNumObjects);
                        out.push_back(s.m_vStart.m_dSolutionValue);out.push_back(s.m_vEnd.m_dSolutionValue);
                        for(int j=0;j<3;j++) out.push_back(s.m_vStart[j]);
                        for(int j=0;j<3;j++) out.push_back(s.m_vEnd[j]);
                        for(int j=0;j<2;j++) out.push_back(s.m_apNodes[j]==branch[j]?0.0:-1.0);
                        for(int j=0;j<2;j++) out.push_back(s.m_apObjects[j]==curve?0.0:(s.m_apObjects[j]==surface?1.0:-1.0));
                    }
                }
            } else out.push_back(status(r));
        }
        delete curve;
    } else {delete surface;return SC_ERR_ARG;}
    delete surface;
    return sc_copy_out(out.data(),(int64_t)out.size(),output,cap,out_len);
    SC_CATCH
}

