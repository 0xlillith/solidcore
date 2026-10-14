// SPDX-License-Identifier: Apache-2.0
// surface.cone_p2: actual sm_GetEllipseTrimIntervals (SmCone.cpp:5282).
#include "sc_shim.h"
#include "StdAfx.h"
#include <SmCone.h>
#include <SmEllipse.h>
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmGeomUtility.h>
#include <vector>
#include <memory>
#include <cmath>

static void sm_GetEllipseTrimIntervals
  (const SmContext     &crContext,            // in : for temp object construction
   const SmPoint3d     &rEllipseOrigin,       // in : origin            of ellipse
   const SmVector3d    &rEllipseXAxis,        // in : xAxis             of ellipse
   const SmVector3d    &rEllipseYAxis,        // in : yAxis             of ellipse
   double               dEllipseXRadius,      // in : x radius          of ellipse
   double               dEllipseYRadius,      // in : y radius          of ellipse
   double               dEllipseStartAngDeg,  // in : start angle (deg) of ellipse
   double               dEllipseEndAngDeg,    // in : end angle (deg)   of ellipse
   const SmCone        &crCone,               // in : Cone with trim boundaries
   const SmExtent2d    &crUVDomain,           // in : nurb trim boundaries for cone
   double               dTol,                 // in : approximation tolerance
   SmTArray<SmExtent1d> &rIvls) ;             // out: array of ellipse intervals inside trimmed cone

// Verbatim internal helper from pinned SmCone.cpp:5282-5619.
void sm_GetEllipseTrimIntervals
  (const SmContext     &crContext,            // in : for temp object construction
   const SmPoint3d     &rEllipseOrigin,       // in : origin            of ellipse
   const SmVector3d    &rEllipseXAxis,        // in : xAxis             of ellipse
   const SmVector3d    &rEllipseYAxis,        // in : yAxis             of ellipse
   double               dEllipseXRadius,      // in : x radius          of ellipse
   double               dEllipseYRadius,      // in : y radius          of ellipse
   double               dEllipseStartAngDeg,  // in : start angle (deg) of ellipse
   double               dEllipseEndAngDeg,    // in : end angle (deg)   of ellipse
   const SmCone        &crCone,               // in : Cone with trim boundaries
   const SmExtent2d    &crUVDomain,           // in : nurb trim boundaries for cone
   double               dTol,                 // in : approximation tolerance
   SmTArray<SmExtent1d> &rIvls)               // out: array of ellipse intervals inside trimmed cone
{
  // gwc:NOT DONE YET - this function is not debugged and
  //     it does not yet trim to the input dEllipseStartAngDeg and
  //     dEllipseEndAngDeg.

  // init output
  rIvls.ReSet() ;

  // find ellipse/coneTrimBoundary intersections
  SmTArray<double> dXSectAngle(6,NULL,6) ;  // coneBoundary xSect Points
  SmBoolean        bInside[6] ;             // TRUE = ith Ivl is inside trim boundaries

  // ellipse locals
  ULONG ii, lCnt = 0 ;               
  SmVector3d       sEllipseZ       = rEllipseXAxis * rEllipseYAxis ;
  SmExtent1d       sEllipseIvl (dEllipseStartAngDeg, dEllipseEndAngDeg) ;
  SmEllipse        sFullEllipse(rEllipseOrigin, rEllipseXAxis, rEllipseYAxis,
                                SmExtent1d(0.0, 360.0),
                                dEllipseXRadius, dEllipseYRadius,
                                3, &crContext) ;
  // cone locals
  double           dHeight       = crCone.GetHeight();
  SmVector3d       sConeAxis     = crCone.GetPosition().GetZAxis();
//cbi change:  const SmPoint3d &rCenterBot    = crCone.GetPosition().GetOriginRef();
  SmPoint3d rCenterBot = crCone.GetPosition().GetOriginRef()
                         + crCone.GetSTEPUVDomain().GetVMin() * sConeAxis;
  SmPoint3d        sCenterTop    = rCenterBot + dHeight * sConeAxis;

  // map given NurbDomain into StepDomain - m_bInsideOut swaps Nurb GenCurve Min/Max Points
  SmExtent2d       sAnalUVDomain ;
  crCone.ConvertDomainFromNURBSToSTEP(crUVDomain, sAnalUVDomain) ;
  SmExtent1d       sSweepIvl     = sAnalUVDomain.GetUInterval() ;
  //SmBoolean        bIsClosed     = sSweepIvl.IsClosed(360.0) ;  //unused
  //double           dStartAng     = sSweepIvl.GetMin();        //unused
  //double           dEndAng       = sSweepIvl.GetMax();        //unused

  // get cone corners [StartLine start, StartLine end, EndLine start EndLine end]
  SmPoint3d sCorners[4] ;
  crCone.EvaluateSTEPPoint(sAnalUVDomain.Evaluate(0.0, 0.0), sCorners[0]) ;
  crCone.EvaluateSTEPPoint(sAnalUVDomain.Evaluate(0.0, 1.0), sCorners[1]) ;
  crCone.EvaluateSTEPPoint(sAnalUVDomain.Evaluate(1.0, 0.0), sCorners[2]) ;
  crCone.EvaluateSTEPPoint(sAnalUVDomain.Evaluate(1.0, 1.0), sCorners[3]) ;
  SmVector3d sStartLineVec = sCorners[1] - sCorners[0] ;
  SmVector3d sEndLineVec   = sCorners[3] - sCorners[2] ;

#ifdef SM_DEBUG_CODE
SmBoolean bDebugMe = FALSE;
static ULONG lCount      = 1 ; lCount++ ;
static ULONG lDebugCount = 0 ;
  if (bDebugMe || lCount == lDebugCount) 
    {
      SmBSplineCurve *pCirc1 = NULL ; 
      if(SM_IS_ZERO(crCone.GetTopRadius()))
           { SmBSplineCurve::CreateDegenerateCurve(crContext, 3, sCenterTop, pCirc1) ; }
      else { pCirc1 = new (crContext) SmCircle(sCenterTop, 
                                               crCone.GetPosition().GetXAxisRef(), 
                                               crCone.GetPosition().GetYAxisRef(),
                                               sSweepIvl, crCone.GetTopRadius(), 3, &crContext) ;
           }
      SmObjDelete sClean1(pCirc1) ;

      SmBSplineCurve *pCirc2 = NULL ; 
      if(SM_IS_ZERO(crCone.GetBotRadius()))
           { SmBSplineCurve::CreateDegenerateCurve(crContext, 3, rCenterBot, pCirc2) ; }
      else { pCirc2 = new (crContext) SmCircle(rCenterBot, 
                                               crCone.GetPosition().GetXAxisRef(), 
                                               crCone.GetPosition().GetYAxisRef(),
                                               sSweepIvl, crCone.GetBotRadius(), 3, &crContext) ;
           }
      SmObjDelete sClean2(pCirc2) ;
      SmFace *pFace = (SmFace *)crCone.GetFace() ;
      SmBrep *pBrep = pFace ? pFace->GetBrep() : NULL ;
      SmTArray<SmExtent1d> sEllipseNurbIvl ;
      sFullEllipse.ConvertIvlFromSTEPToNURBS(sEllipseIvl, sEllipseNurbIvl) ;

      smgfx_Erase();
      smgfx_SetLook(1,2, 0,0,1) ; if(pBrep) { pBrep->Draw(TRUE) ; } sm_GraphicsLoop() ;
      smgfx_SetLook(1,2, 0,1,1) ; crCone.DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 1,0,0) ; sFullEllipse.Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(6,7, 0,0,0) ; for(ii=0;ii<sEllipseNurbIvl.GetSize();ii++)
                                    { sFullEllipse.Draw(&sEllipseNurbIvl[ii]) ; sm_GraphicsLoop() ; }
      smgfx_SetLook(3,4, 0,1,0) ; (dEllipseXRadius * rEllipseXAxis).Draw(&rEllipseOrigin) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 0,0,1) ; (dEllipseYRadius * rEllipseYAxis).Draw(&rEllipseOrigin) ; sm_GraphicsLoop() ;
      smgfx_SetLook(5,6, 1,1,0) ; sStartLineVec.Draw(&sCorners[0]) ; sm_GraphicsLoop() ;
      smgfx_SetLook(5,6, 1,1,0) ; sEndLineVec.Draw(&sCorners[2]) ; sm_GraphicsLoop() ;
      smgfx_SetLook(5,6, 1,1,0) ; pCirc1->Draw() ; sm_GraphicsLoop() ;
      smgfx_SetLook(5,6, 1,1,0) ; pCirc2->Draw() ; sm_GraphicsLoop() ;
      sm_GraphicsLoop();
    }
#endif

  // get ellipse/TopCircle(max 2), ellipse/BotCircle(max 2) 
  //     ellipse/StartSeamLine(max 1), ellipse/EndSeamLine(max 1) intersections
  ULONG lNumTsect1, lNumTsect2 ;
  SmPoint3d sTsectPoints[6] ;
  double dParam[2] ;
  smgu_CirclePlaneIntersect(dTol, 
                            crCone.GetTopRadius(), sCenterTop, sConeAxis,
                            rEllipseOrigin, sEllipseZ, 
                            lNumTsect1, sTsectPoints) ;
  smgu_CirclePlaneIntersect(dTol, 
                            crCone.GetBotRadius(), rCenterBot, sConeAxis,
                            rEllipseOrigin, sEllipseZ, 
                            lNumTsect2, sTsectPoints+2) ;
  smgu_LinePlaneIntersect  (sCorners[0], sStartLineVec,
                            rEllipseOrigin, sEllipseZ,
                            dParam[0]) ;
  smgu_LinePlaneIntersect  (sCorners[2], sEndLineVec,
                            rEllipseOrigin, sEllipseZ,
                            dParam[1]) ;

  // Check for degenerate solutions, within tol. [B10]
  // (smgu_CirclePlaneIntersect() doesn't do that.)
  if ( lNumTsect1 == 2 )
  {
      if ( sTsectPoints[0].DistanceBetween( sTsectPoints[1] ) < dTol )
        { lNumTsect1 = 1; }
  }
  if ( lNumTsect2 == 2 )
  {
      if ( sTsectPoints[2].DistanceBetween( sTsectPoints[3] ) < dTol )
        { lNumTsect2 = 1; }
  }

  // circle/Plane xSect lNumTSect == 0 is a miss - nothing to do
  //                    lNumTSect == 1 is a tangent point - not put on in/out point list
  //                    lNumTSect == 2 is a trimmed section in/out boundary
  //                    lNumTSect == 3 is a coincident circle/plane 
  SM_ASSERT(   lNumTsect1 != 3 || lNumTsect2 != 3) ;
          
  // get ellipse angle in degrees for each XSectPoint and place in the in/out dXSectAngle array
  // exclude intersection points outside the cone sweep trim boundaries

  // Consider all 4 possible ConeCircle/Ellipse intersections for addition to dXSectAngle array
  SmPoint2d sConeSTEPUV ;
  SmLocationType eLoc ;
  for(ii=0;ii<4;ii++)
    {
      // only consider circle/ellipse intersections when there are two of them
      //      i.e. don't add point tangent intersections - they don't affect in/out behavior.
      //           don't add coincident intersections - they get trimmed to cone end lines later.
      // only add intersection points in the cone trim boundary
      if(   (ii < 2 && lNumTsect1 == 2)
         || (ii >=2 && lNumTsect2 == 2))
        {
          // map XSectPoint to ConeUV point
          crCone.STEPInversion(sAnalUVDomain, sTsectPoints[ii], dTol, sConeSTEPUV, eLoc) ;

          // when ConeUV is cone trim boundaries
          //      skip circle instersections on seams because the line intersections will catch seam intersections
          if(   eLoc == SM_LT_INTERIOR
             && eLoc != SM_LT_U_SEAM)
            {
              // add ellipse angleDeg for intersection to in/out angle dXSectAngle array
              sFullEllipse.STEPInversion(sTsectPoints[ii], dXSectAngle[lCnt]) ;
              lCnt++ ;
            }
        } // end skip tangent intersection cases
    } // end iter every possible Circle/Ellipse xSect

  // for StartSeamLine/ellipse xSect
  //  GWC: tolerances to think about:
  //    ContainsPeriodicValue()  increased tolerance to account for common degree domains ;
  //        to double dScaledZero = SM_EFF_ZERO * 100.0 * (1.0 + dPeriod) ;
  //        and here we use a different value for lines
  //    This computation needs both functions to have the same behavior near corners.
  double dScaledZero = SM_EFF_ZERO * (1.0 + smos_Max(sCorners[0].GetMaxDimension(),
                                                     sCorners[3].GetMaxDimension())) ;
  if(   (   lNumTsect1 == 3 || lNumTsect2 == 3)
     || (   dParam[0] >= 0.0 - dScaledZero
         && dParam[0] <= 1.0 + dScaledZero))
    {
      sTsectPoints[4] = sCorners[0] + dParam[0] * sStartLineVec;
      sFullEllipse.STEPInversion( sTsectPoints[4], dXSectAngle[lCnt] );
      lCnt += 1 ;
    }

  // for EndSeamLine/ellipse xsect
  if(   (   lNumTsect1 == 3 || lNumTsect2 == 3)
     || (   dParam[1] >= 0.0 - dScaledZero
         && dParam[1] <= 1.0 + dScaledZero))
    {
      sTsectPoints[5] = sCorners[2] + dParam[1] * sEndLineVec;
      sFullEllipse.STEPInversion( sTsectPoints[5], dXSectAngle[lCnt] );
      lCnt += 1 ;
    }

  // there should be no more than 6 intersections
  SM_ASSERT(lCnt <= 6) ;

  // arrive here when dXSectAngle array contains an ellipse angle for
  // every ConeTrimCurve/FullEllipse intersection where the FullEllipse goes in/out
  // of the ConeTrimDomain.  There should be an even number of these transitions
  // The ellipse trim range has not yet been applied.  When there are no intersections
  // the whole curve must lie outside the cone trim domain.
  if(lCnt < 2) 
    { return ; }

  // when there are intersections there should be an even number of them
  // unless ellipse intersects a corner exactly.
  //      SM_ASSERT(lCnt % 2 == 0) ;  

  // order the XSectPoints by angle
  smgu_ShellSort_double_array(dXSectAngle.GetDataArray(), lCnt) ;
  SM_ASSERT(dXSectAngle[lCnt-1] - dXSectAngle[0] <= 360.0) ;

  // now classify the segment midPoints as in/out
  // add one interval to output for every inside segment
  double dScaledAngZero = SM_EFF_ZERO * 100.0 * 360.0 ;
  ULONG iThis, iNext, iPrev;
  for( iThis=0; iThis<lCnt; iThis++ )
    {
      iNext = (iThis+1) % lCnt;

      // get interval start/stop angles
      double dStartAngDeg = dXSectAngle[iThis] ;
      double dEndAngDeg   = dXSectAngle[iNext] ;

      // the last interval includes the zero angle
      // modify the startPoint so that the interval is
      // within [-360, 360] and Start < End
      if ( iNext == 0 ) { dStartAngDeg -= 360.0; }

      // mark degenerate intervals as unknown for now - they will be classified
      // based on their neighbor classifications once they are available
      if((dEndAngDeg - dStartAngDeg < dScaledAngZero))
        { bInside[iThis] = UNSURE ; 
          continue ;
        }

      // classify nonDegenerate interval midPoint
      SmPoint3d sSegmentMidPoint ;
      double sSegmentMidAngDeg = (dEndAngDeg + dStartAngDeg) / 2.0 ;
      sFullEllipse.EvaluateSTEPPoint(sSegmentMidAngDeg, sSegmentMidPoint) ;
      crCone.STEPInversion(sAnalUVDomain, sSegmentMidPoint, dTol, sConeSTEPUV, eLoc) ;

      // Remember if segment MidPointUV is within cone trim boundaries
      bInside[iThis] = (eLoc == SM_LT_EXTERIOR) ? FALSE : TRUE ;

    } // end iter every interval to classify it

  // add one interval to output for every inside segment
  SmExtent1d sUnTrimmedIvl ;
  SmTArray<SmExtent1d> sTrimIvls ;

  for ( iThis=0; iThis<lCnt; iThis++ )
    {
      iNext = ( iThis+1 )      % lCnt;
      iPrev = ( iThis+lCnt-1 ) % lCnt;

      // get interval start/stop angles
      //   the last interval includes the zero angle
      //   modify the startPoint so that the interval is
      //   within [-360, 360] and Start < End.
      double dStartAngDeg = (iNext == 0) ? dXSectAngle[iThis] - 360.0
                                         : dXSectAngle[iThis] ;
      double dEndAngDeg   = dXSectAngle[iNext] ;

      // when the segment is inside
      if(   bInside[iThis] == TRUE   // nonDegnerate Ivl classification
         || (   bInside[iThis] == UNSURE  // Degenerate IVl is inside
             && bInside[iPrev] == FALSE   // when both neighbor intervals 
             && bInside[iNext] == FALSE)) //      are outside
        {
          // make and trim interval to ellipse interval
          // (Don't use SetMinMax, start ang can be > end ang, by machine noise.
          //   SetMinMax checks '>' without tol, and causes trouble if bad.)
          sUnTrimmedIvl.Init();
          sUnTrimmedIvl.AddValue( dStartAngDeg );
          sUnTrimmedIvl.AddValue( dEndAngDeg );

          sEllipseIvl.IntersectPeriodic(sUnTrimmedIvl, 360.0, sTrimIvls, FALSE) ;

          // add every trim invl to output
          rIvls.Append(sTrimIvls) ;                

        } // end segment inside check

    } // end iter every ellipse segment
  
  // all done

#ifdef SM_DEBUG_CODE
  if(bDebugMe || lDebugCount == lCount)
    {
      dXSectAngle.Dump() ;
      rIvls.Dump() ;

      // Build intersection ellipse - get face and Brep pointers
      SmPoint3d sPoint ;
      SmFace *pThisFace  = ((SmFace *)crCone.GetFace()) ;
      SmBrep *pThisBrep  = pThisFace  ? pThisFace->GetBrep()  : NULL ;

      // draw Breps(Blue,Green), Surfaces(cyan,yellow), faces(black), 
      // ellipse(red), TrimPoints(blue)
      smgfx_Erase() ; 
      smgfx_SetLook(1,2, 0,0,1) ; if(pThisBrep)  { pThisBrep->Draw(TRUE) ; }  sm_GraphicsLoop() ; 
      smgfx_SetLook(1,2, 0,1,1) ; crCone.DrawUV() ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 0,1,0) ; (dEllipseXRadius * rEllipseXAxis).Draw(&rEllipseOrigin) ; sm_GraphicsLoop() ;
      smgfx_SetLook(3,4, 0,0,1) ; (dEllipseYRadius * rEllipseYAxis).Draw(&rEllipseOrigin) ; sm_GraphicsLoop() ;
      smgfx_SetLook(6,7, 0,0,1) ; 

      for(ii=0;ii<lCnt;ii++) { sFullEllipse.EvaluateSTEPPoint(dXSectAngle[ii], sPoint) ;
                               smgfx_SetLook(8,9, 0,0,1) ; sPoint.Draw() ; sm_GraphicsLoop() ;
                             }
      for(ii=0;ii<4;ii++) { if(   (ii < 2 && lNumTsect1 == 2)
                               || (ii >=2 && lNumTsect2 == 2)
                               || (ii >=4))
                              { smgfx_SetLook(10,11, 1,0,0) ; sTsectPoints[ii].Draw() ; sm_GraphicsLoop() ; }
                          }
      smgfx_SetLook(10,11, 1,0,0) ; (sCorners[0] + dParam[0] * sStartLineVec).Draw() ; sm_GraphicsLoop() ; 
      smgfx_SetLook(10,11, 1,0,0) ; (sCorners[2] + dParam[1] * sEndLineVec).Draw() ; sm_GraphicsLoop() ;  
      for(ii=0;ii<rIvls.GetSize();ii++) { 
                                          ULONG jj;
                                          SmExtent1d sStepIvl = rIvls[ii] ;
                                          SmTArray<SmExtent1d> sNurbIvl ;
                                          sFullEllipse.ConvertIvlFromSTEPToNURBS(sStepIvl, sNurbIvl) ;
                                          for(jj=0;jj<sNurbIvl.GetSize();jj++)
                                            { smgfx_SetLook(8,9, 0,1,0) ; sFullEllipse.Draw(&sNurbIvl[jj]) ; sm_GraphicsLoop() ; }
                                        }
      sm_GraphicsLoop() ;
    }
#endif // SM_DEBUG_CODE

} // end sm_GetEllipseTrimIntervals

namespace {
struct TrimCone : SmCone {
    TrimCone(double b,double t,double s,double e,double h,SmBoolean swap,SmBoolean inside)
      : TrimCone(SmPoint3d(0,0,0),SmVector3d(1,0,0),SmVector3d(0,1,0),b,t,s,e,h,swap,inside) {}
    TrimCone(const SmPoint3d& o,const SmVector3d& x,const SmVector3d& y,
             double b,double t,double s,double e,double h,SmBoolean swap,SmBoolean inside)
      : SmCone(o,x,y,b,t,s,e,h,swap,inside,FALSE) {}
    TrimCone(const SmPoint3d& o,const SmVector3d& x,const SmVector3d& y,
             double radius,double semi,SmBoolean swap,SmBoolean inside)
      : SmCone(o,x,y,radius,semi,swap,inside) {}
};
}

// 18 doubles: bot,top,height,start,end,swap,inside,tilt_rad,ellipse_z,
// ellipse_rx,ellipse_ry,arc_start,arc_end,u0,v0,u1,v1,tolerance.
SC_EXPORT int32_t sc_surface_cone_p2_trim(const double* in,int64_t n_in,
    double* out,int64_t cap,int64_t* out_len)
{
    SC_TRY
    if(n_in!=18) return SC_ERR_ARG;
    SmContext context;
    std::unique_ptr<TrimCone> owned(new(context) TrimCone(in[0],in[1],in[3],in[4],in[2],in[5]?TRUE:FALSE,in[6]?TRUE:FALSE));
    TrimCone& cone=*owned;
    if(cone.MakeNurb()!=SM_SUCCESS) return SC_ERR_ARG;
    const SmPoint3d origin(0,0,in[8]);
    const SmVector3d x(std::cos(in[7]),0,std::sin(in[7])),y(0,1,0),normal=x*y;
    SmExtent2d domain(in[13],in[14],in[15],in[16]),analytic;
    if(cone.ConvertDomainFromNURBSToSTEP(domain,analytic)!=SM_SUCCESS) return SC_ERR_ARG;
    SmPoint3d corners[4];
    const double uv[4][2]={{0,0},{0,1},{1,0},{1,1}};
    for(int i=0;i<4;++i)
        if(cone.EvaluateSTEPPoint(analytic.Evaluate(uv[i][0],uv[i][1]),corners[i])!=SM_SUCCESS)
            return SC_ERR_ARG;
    // The helper ignores these statuses then reads their output. Require its
    // actual primitives to succeed before calling it, instead of excluding a
    // broad geometry family or inventing a reference result.
    for(int i=0;i<2;++i) {
        double parameter;
        if(smgu_LinePlaneIntersect(corners[2*i],corners[2*i+1]-corners[2*i],origin,normal,parameter)!=SM_SUCCESS)
            return SC_ERR_ARG;
    }
    // CirclePlaneIntersect always initializes its count, including on errors;
    // the helper can safely consume it. Only an unset line parameter is rejected.
    SmTArray<SmExtent1d> intervals;
    sm_GetEllipseTrimIntervals(context,origin,x,y,in[9],in[10],in[11],in[12],cone,domain,in[17],intervals);
    std::vector<double> values;
    values.push_back((double)intervals.GetSize());
    for(ULONG i=0;i<intervals.GetSize();++i) {
        values.push_back(intervals[i].GetMin());
        values.push_back(intervals[i].GetMax());
    }
    return sc_copy_out(values.data(),(int64_t)values.size(),out,cap,out_len);
    SC_CATCH
}

// Independent native construction and full geometric state for the remaining
// three owned functions. The input cone is never replaced by a generic surface.
#include <SmCylinder.h>
#include <SmBSplineSurface.h>
#include <SmCurve.h>
#include <nurbs.h>

namespace {
struct UnboundedCylinder : SmCylinder {
    UnboundedCylinder(const SmPoint3d& o,const SmVector3d& x,const SmVector3d& y,
                      double r,SmBoolean swap,SmBoolean inside)
      : SmCylinder(o,x,y,r,swap,inside) {}
};
struct ConeStateAccess : SmCone {
    static SmPolarConversion SmSurfOfRevolution::*polar() {return &ConeStateAccess::m_vPolarConverter;}
};
void cone_vec(std::vector<double>& out,const SmVector3d& p) {
    out.push_back(p.x);out.push_back(p.y);out.push_back(p.z);
}
bool cone_state(std::vector<double>& out,const SmCone& cone) {
    SmAxis2Placement position;double radius=0,semi=0;SmBoolean swap=FALSE,inside=FALSE;
    if(cone.GetCanonical(position,radius,semi,&swap,&inside)!=SM_SUCCESS) return false;
    out.insert(out.end(),{(double)cone.GetType(),radius,semi,(double)swap,(double)inside});
    cone_vec(out,position.GetOriginRef());cone_vec(out,position.GetXAxisRef());cone_vec(out,position.GetYAxisRef());
    for(const auto& domain:{cone.GetSTEPUVDomain(),cone.GetNaturalUVDomain()})
        out.insert(out.end(),{domain.GetMin().x,domain.GetMin().y,domain.GetMax().x,domain.GetMax().y});
    out.insert(out.end(),{cone.GetBotRadius(),cone.GetTopRadius(),cone.GetHeight()});
    const SmBSplineCurve* generator=cone.GetGenCurve();
    if(!generator) return false;
    const SmExtent1d step=generator->GetSTEPInterval(),nurb=generator->GetNaturalInterval();
    out.insert(out.end(),{(double)generator->GetType(),step.GetMin(),step.GetMax(),nurb.GetMin(),nurb.GetMax()});
    for(double fraction:{0.0,0.4,1.0}) {
        SmPoint3d point;
        if(generator->EvaluatePoint(nurb.Evaluate(fraction),point)!=SM_SUCCESS) return false;
        cone_vec(out,point);
    }
    const SmPolarConversion& polar=cone.*ConeStateAccess::polar();
    const ULONG count=polar.GetAngleCount();
    out.insert(out.end(),{cone.HavePolarConversion()?1.0:0.0,(double)count});
    if(count) out.insert(out.end(),{polar.GetFirstAngle(),polar.GetLastAngle()});
    const gw_SURFACE* surface=cone.GetGwNurbPointer();
    out.push_back(surface?1.0:0.0);
    if(!surface) return true;
    out.insert(out.end(),{(double)surface->p,(double)surface->q,(double)(surface->net->n+1),(double)(surface->net->m+1)});
    for(gw_INDEX i=0;i<=surface->net->n;++i) for(gw_INDEX j=0;j<=surface->net->m;++j) {
        const auto& p=surface->net->Pw[i][j];out.insert(out.end(),{p.x,p.y,p.z,p.w});
    }
    out.push_back((double)(surface->knu->m+1));
    for(gw_INDEX i=0;i<=surface->knu->m;++i) out.push_back(surface->knu->U[i]);
    out.push_back((double)(surface->knv->m+1));
    for(gw_INDEX i=0;i<=surface->knv->m;++i) out.push_back(surface->knv->U[i]);
    return true;
}
void cone_mutate(SmCone& cone,int mutation) {
    if(!mutation) return;
    gw_SURFACE* surface=cone.GetGwNurbPointer();
    const bool swap=cone.GetSwapUV()!=FALSE;
    if(mutation==3) {
        auto* knots=swap?surface->knv:surface->knu;
        for(gw_INDEX i=1;i<knots->m;++i) {
            if(knots->U[i]>knots->U[0] && knots->U[i]<knots->U[knots->m]) {
                knots->U[i]=(knots->U[0]+knots->U[i])/2.0;break;
            }
        }
        return;
    }
    if(mutation==4 || mutation==5) {
        for(gw_INDEX i=0;i<=surface->net->n;++i) for(gw_INDEX j=0;j<=surface->net->m;++j) {
            auto& p=surface->net->Pw[i][j];
            if(mutation==4) p.z=0.0;
            else {gw_CPOINT euclidean;TO_EUCLID(p,euclidean);p.x=euclidean.x;p.y=euclidean.y;p.z=euclidean.z;p.w=NL_NOW;}
        }
        return;
    }
    const double angle=0.05,cosine=std::cos(angle),sine=std::sin(angle);
    for(gw_INDEX i=0;i<=surface->net->n;++i) for(gw_INDEX j=0;j<=surface->net->m;++j) {
        if(swap ? i!=surface->net->n : j!=surface->net->m) continue;
        auto& p=surface->net->Pw[i][j];
        if(mutation==1) {double x=p.x,y=p.y;p.x=cosine*x-sine*y;p.y=sine*x+cosine*y;}
        else {p.x*=2.0;p.y*=2.0;p.z*=2.0;p.w*=2.0;}
    }
}
}


namespace {
// The native offset path ignores CreateCanonical's status and dereferences its
// output. Reject only a NULL result from that actual factory, before executing
// the unsafe offset call. The probe uses the source's exact truncation formulas.
bool offset_factory_defined(const SmContext& context,const SmCone& cone,double distance,double tolerance) {
    SmAxis2Placement placement;double base=0,semi=0;SmBoolean swap=FALSE,inside=FALSE;
    if(cone.GetCanonical(placement,base,semi,&swap,&inside)!=SM_SUCCESS) return false;
    const double sign=1.0*(swap?-1.0:1.0)*(inside?-1.0:1.0);
    const double sine=std::sin(semi*SM_PI/180.0),cosine=std::cos(semi*SM_PI/180.0);
    const double height=cone.GetHeight(),offset_z=-sine*sign*distance;
    const double radial=cosine*sign*distance;
    const double old_bottom=cone.GetBotRadius(),old_top=cone.GetTopRadius();
    double bottom=old_bottom+radial,top=old_top+radial;
    if(bottom<tolerance && top<tolerance) return true;
    if(bottom>=tolerance && top>=tolerance) return true;
    double new_height;SmPoint3d origin;
    const SmVector3d axis=placement.GetZAxis();
    if(bottom<tolerance) {
        const double tangent=height/old_top;new_height=top*tangent;
        origin=placement.GetOriginRef()+(height+offset_z-new_height)*axis;bottom=0.0;
    } else {
        const double tangent=height/old_bottom;new_height=bottom*tangent;
        origin=placement.GetOriginRef()+offset_z*axis;top=0.0;
    }
    SmAxis2Placement target;target.SetCanonical(origin,placement.GetXAxisRef(),placement.GetYAxisRef());
    SmCone* raw=NULL;
    SmCone::CreateCanonical(context,target,bottom,top,new_height,raw,swap,inside);
    std::unique_ptr<SmCone> probe(raw);
    return raw!=NULL;
}
}

// op 1 recognition, 2 rebuild, 3 offset. Kind bits: cylinder=1, uncached=2, unbounded=4.
// Optional final nine scalars specify origin, X axis and Y axis; real native types are used.
SC_EXPORT int32_t sc_surface_cone_p2_analytic(int32_t op,int32_t kind,int32_t mutation,
    double scale,double distance,double tolerance,const double* in,int64_t n_in,
    double* out,int64_t cap,int64_t* out_len)
{
    SC_TRY
    if(n_in!=18 && n_in!=27) return SC_ERR_ARG;
    SmContext context;
    const SmPoint3d origin=n_in==27?SmPoint3d(in[18],in[19],in[20]):SmPoint3d(0,0,0);
    const SmVector3d x=n_in==27?SmVector3d(in[21],in[22],in[23]):SmVector3d(1,0,0);
    const SmVector3d y=n_in==27?SmVector3d(in[24],in[25],in[26]):SmVector3d(0,1,0);
    const SmBoolean swap=in[5]?TRUE:FALSE,inside=in[6]?TRUE:FALSE;
    std::unique_ptr<SmCone> input;
    if(kind&4) {
        if(kind&1) input.reset(new(context) UnboundedCylinder(origin,x,y,in[0],swap,inside));
        else input.reset(new(context) TrimCone(origin,x,y,in[0],in[7]*180.0/SM_PI,swap,inside));
    } else {
        if(kind&1) input.reset(new(context) SmCylinder(origin,x,y,in[0],in[3],in[4],in[2],swap,inside,FALSE,&context));
        else input.reset(new(context) TrimCone(origin,x,y,in[0],in[1],in[3],in[4],in[2],swap,inside));
    }
    if(!(kind&2) && input->MakeNurb()!=SM_SUCCESS) return SC_ERR_ARG;
    cone_mutate(*input,mutation);
    std::unique_ptr<SmCone> result;
    const SmCone* selected=NULL;double status=0;
    if(op==1) {
        std::unique_ptr<SmBSplineSurface> pure(new(context) SmBSplineSurface(*input));
        SmCone* recognized=NULL;
        status=(double)SmCone::IsNurbSurfaceCone(context,pure.get(),recognized,scale);
        result.reset(recognized);selected=result.get();
    } else if(op==2) {
        status=(double)input->RebuildSTEPFromNURBParameters();selected=input.get();
    } else if(op==3) {
        if(!offset_factory_defined(context,*input,distance,SmApproxTol3d(tolerance))) return SC_ERR_ARG;
        SmSurface* offset=NULL;
        status=(double)input->CreateOffsetSurface(context,distance,tolerance,offset);
        if(offset && !offset->IsKindOf(SmCone_TYPE)) {delete offset;return SC_ERR_ARG;}
        result.reset(static_cast<SmCone*>(offset));selected=result.get();
    } else return SC_ERR_ARG;
    std::vector<double> values{status,selected?1.0:0.0};
    if(selected && !cone_state(values,*selected)) return SC_ERR_ARG;
    return sc_copy_out(values.data(),(int64_t)values.size(),out,cap,out_len);
    SC_CATCH
}
