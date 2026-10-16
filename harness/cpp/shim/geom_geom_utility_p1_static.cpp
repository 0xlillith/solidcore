// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `geom.geom_utility_p1`: the file-static helpers sm_AddPoint and
// sm_AddClosestPoint of upstream source/SMLib/src/SmGeomUtility.cpp.
//
// File-static helpers cannot be linked against. Copy only the two helpers,
// verbatim from pinned c9c5979 SmGeomUtility.cpp:151-248. Compiling the whole
// source in a namespace makes its class forward declarations ill-formed in Clang.
#include "sc_shim.h"
#include "StdAfx.h"
#include <SmVector3d.h>

namespace sc_gu_p1_tu {
// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0
static SmStatus sm_AddPoint
  (      ULONG     & rlNum,      // i/o: number of points in array
         SmPoint3d   aPoints[2], // i/o: array of distinct points
         double      dTol,       // in : min distance between distinct points
   const SmPoint3d & crPnt)      // in : target point to check
{
  // for every point
  double dTol2 = dTol*dTol;
  for (ULONG i=0; i<rlNum; i++)
    {
      // quit if a point within tolerance is already on the list
      if( crPnt.DistanceBetweenSquared(aPoints[i]) < dTol2)
          return SM_SUCCESS;
    }

  // error - trying to add more than 2 points to the list
  if (rlNum > 1)
    { SER(SM_ERR); }

  // add the point
  aPoints[rlNum] = crPnt;
  rlNum++;

  // all done
  return SM_SUCCESS;

} // end sm_AddPoint

/*******************************************************************//**
PURPOSE: Internal helper function for smgu_SegmentSegmentClosestPoint.

NOTES: Add points to the closest point list.
***********************************************************************/

static SmStatus sm_AddClosestPoint
  (double dDist,
   double dParam1,
   double dParam2,
   double dTolerance,
   ULONG & rlNumClosestPoints,
   double adParamSeg1[2],
   double adParamSeg2[2],
   double adDistances[2])

{
    if (rlNumClosestPoints == 0) {
        adParamSeg1[0] = dParam1;
        adParamSeg2[0] = dParam2;
        adDistances[0] = dDist;
        rlNumClosestPoints = 1;
        return SM_SUCCESS;
    }
    if (dDist > adDistances[rlNumClosestPoints-1]+dTolerance) {
        return SM_SUCCESS;  // Too big of a distance skip it
    }

    if (dDist >= adDistances[rlNumClosestPoints-1]) {
        // Just off the end either skip it or add to second slot.
        if (rlNumClosestPoints == 1) {
            adParamSeg1[1] = dParam1;
            adParamSeg2[1] = dParam2;
            adDistances[1] = dDist;
            rlNumClosestPoints ++;
        }
        return SM_SUCCESS;
    }

    // Must be TRUE if (dDist < adDistances[rlNumClosestPoints-1])

    if (dDist < adDistances[0]-dTolerance) {
        // Totally new one
        adParamSeg1[0] = dParam1;
        adParamSeg2[0] = dParam2;
        adDistances[0] = dDist;
        rlNumClosestPoints = 1;
        return SM_SUCCESS;
    }

    if (dDist > adDistances[0]) {
        // Must fall between 0 and 1
        adParamSeg1[1] = dParam1;
        adParamSeg2[1] = dParam2;
        adDistances[1] = dDist;
        rlNumClosestPoints = 2;
        return SM_SUCCESS;
    }

    // Final case is where we make 2 from 1
    rlNumClosestPoints = 2;
    adParamSeg1[1] = adParamSeg1[0];
    adParamSeg2[1] = adParamSeg2[0];
    adDistances[1] = adDistances[0];
    adParamSeg1[0] = dParam1;
    adParamSeg2[0] = dParam2;
    adDistances[0] = dDist;

    return SM_SUCCESS;

} // end sm_AddClosestPoint
} // namespace sc_gu_p1_tu

// sm_AddPoint(rlNum, aPoints, dTol, crPnt). io = [n, p0(3), p1(3)] in and out
// (n in 0..=2). Returns 0 and writes *status = the SmStatus.
SC_EXPORT int32_t sc_geom_geom_utility_p1_add_point(double* io, double tol, const double* pnt,
                                                    int32_t* status)
{
    if (io == nullptr || pnt == nullptr || status == nullptr)
        return SC_ERR_ARG;
    if (!(io[0] == 0.0 || io[0] == 1.0 || io[0] == 2.0))
        return SC_ERR_ARG;
    SC_TRY
    ULONG n = (ULONG)io[0];
    SmPoint3d pts[2] = {SmPoint3d(io[1], io[2], io[3]), SmPoint3d(io[4], io[5], io[6])};
    SmStatus st = sc_gu_p1_tu::sm_AddPoint(n, pts, tol, SmPoint3d(pnt[0], pnt[1], pnt[2]));
    *status = (int32_t)st;
    io[0] = (double)n;
    for (int i = 0; i < 2; ++i)
    {
        io[1 + 3 * i] = pts[i].x;
        io[2 + 3 * i] = pts[i].y;
        io[3 + 3 * i] = pts[i].z;
    }
    return SC_OK;
    SC_CATCH
}

// sm_AddClosestPoint(dDist, dParam1, dParam2, dTolerance, rlNum, ad1, ad2, adD).
// io = [n, p1[2], p2[2], d[2]] in and out (n in 0..=2).
SC_EXPORT int32_t sc_geom_geom_utility_p1_add_closest_point(double* io, double dist, double param1,
                                                            double param2, double tol, int32_t* status)
{
    if (io == nullptr || status == nullptr)
        return SC_ERR_ARG;
    if (!(io[0] == 0.0 || io[0] == 1.0 || io[0] == 2.0))
        return SC_ERR_ARG;
    SC_TRY
    ULONG n = (ULONG)io[0];
    double a1[2] = {io[1], io[2]}, a2[2] = {io[3], io[4]}, ad[2] = {io[5], io[6]};
    SmStatus st = sc_gu_p1_tu::sm_AddClosestPoint(dist, param1, param2, tol, n, a1, a2, ad);
    *status = (int32_t)st;
    io[0] = (double)n;
    io[1] = a1[0];
    io[2] = a1[1];
    io[3] = a2[0];
    io[4] = a2[1];
    io[5] = ad[0];
    io[6] = ad[1];
    return SC_OK;
    SC_CATCH
}
