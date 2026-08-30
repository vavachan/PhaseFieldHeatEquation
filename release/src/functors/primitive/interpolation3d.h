/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023-24 Julius Jessberger
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef INTERPOLATION_3D_H
#define INTERPOLATION_3D_H

#include "core/vector.h"

namespace olb {

/// @brief trilinear interpolation of function f
/// @tparam DataDim Dimension of function return value
/// @param distance Distance from cell
/// @param result Result is written here
/// @param f Function which is to interpolate
// Caution: Developer has to guarantee that the data are accessible.
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d(CELL& cell, RESULT& result, const Vector<V,3> distance, F f) any_platform
{
  const Vector<int,3> floorV (util::floor(distance[0]),
    util::floor(distance[1]), util::floor(distance[2]));
  Vector<Vector<int,3>,8> surroundingPoints (floorV);
  surroundingPoints[1][0] += 1;
  surroundingPoints[2][1] += 1;
  surroundingPoints[3][2] += 1;
  surroundingPoints[4][0] += 1; surroundingPoints[4][1] += 1;
  surroundingPoints[5][0] += 1; surroundingPoints[5][2] += 1;
  surroundingPoints[6][1] += 1; surroundingPoints[6][2] += 1;
  surroundingPoints[7][0] += 1; surroundingPoints[7][1] += 1; surroundingPoints[7][2] += 1;

  V resNeighbor[DataDim];
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  for (auto point : surroundingPoints) {
    const Vector<V,3> dist = distance - point;
    const V weight = (V(1) - util::abs(dist[0])) * (V(1) - util::abs(dist[1]))
      * (V(1) - util::abs(dist[2]));
    f(cell.neighbor(point), resNeighbor);
    for (unsigned i = 0; i < DataDim; ++i) {
      result[i] += weight * resNeighbor[i];
    }
  }
}

/// Interpolation between four points (linear independent case)
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d_help_li(CELL& cell, RESULT& result, const Vector<V,3> distance, F f,
  Vector<int,3> pointA, Vector<int,3> pointB, Vector<int,3> pointC, Vector<int,3> pointD) any_platform
{
  const Vector<Vector<V,3>,3> v {pointB - pointA,
    pointC - pointA,
    pointD - pointA};
  const V det = util::determinant(v[0], v[1], v[2]);
  OLB_PRECONDITION (det != 0);

  // linear interpolation
  const auto lc = util::solveLinearSystem_help(v[0], v[1], v[2],
    distance - pointA, det);
  // evaluate linear function
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  V resNeighbor[DataDim];
  f(cell.neighbor(pointA), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1 - lc[0] - lc[1] - lc[2]) * resNeighbor[i];
  }
  f(cell.neighbor(pointB), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += lc[0] * resNeighbor[i];
  }
  f(cell.neighbor(pointC), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += lc[1] * resNeighbor[i];
  }
  f(cell.neighbor(pointD), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += lc[2] * resNeighbor[i];
  }
}

/// Interpolation between four points (coplanar case)
// it is expected that point C lies diagonally opposite to pointA
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d_help_plane(CELL& cell, RESULT& result, const Vector<V,3> distance, F f,
  Vector<int,3> pointA, Vector<int,3> pointB, Vector<int,3> pointC, Vector<int,3> pointD) any_platform
{
  const Vector<Vector<V,3>,3> v {pointB - pointA,
    pointC - pointA,
    pointD - pointA};
  // bilinear interpolation in plane (constant in other direction)
  const auto normal = crossProduct3D(v[0], v[1]);
  const auto nLenSq = norm_squared(normal);
  const auto projBase = distance - pointA
    - (((distance - pointA)*normal)/nLenSq)*normal;

  const V w0 = projBase* v[0] / norm_squared(v[0]);
  const V w1 = projBase* v[2] / norm_squared(v[2]);

  // evaluate bilinear function
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  V resNeighbor[DataDim];
  f(cell.neighbor(pointA), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1-w0) * (1-w1) * resNeighbor[i];
  }
  f(cell.neighbor(pointC), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += w0 * w1 * resNeighbor[i];
  }
  f(cell.neighbor(pointB), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += w0 * (1-w1) * resNeighbor[i];
  }
  f(cell.neighbor(pointD), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1-w0) * w1 * resNeighbor[i];
  }
}

/// Linear interpolation between three points (constant in the orthogonal direction)
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d_help(CELL& cell, RESULT& result, const Vector<V,3> distance, F f,
  Vector<int,3> pointA, Vector<int,3> pointB, Vector<int,3> pointC) any_platform
{
  const Vector<V,3> v1 = pointB - pointA;
  const Vector<V,3> v2 = pointC - pointA;
  const Vector<V,3> v  = distance - pointA;
  const auto normal = crossProduct3D(v1, v2);
  const auto nLenSq = norm_squared(normal);
  const auto lc = util::solveLinearSystem(v1, v2, normal,
    v - ((v*normal)/nLenSq)*normal);  // v projected to plane
  // evaluate linear function
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  V resNeighbor[DataDim];
  f(cell.neighbor(pointA), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1 - lc[0] - lc[1]) * resNeighbor[i];
  }
  f(cell.neighbor(pointB), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += lc[0] * resNeighbor[i];
  }
  f(cell.neighbor(pointC), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += lc[1] * resNeighbor[i];
  }
}

/// Linear interpolation between two points (constant in the orthogonal directions)
// Cf. https://en.wikipedia.org/wiki/Vector_projection
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d_help(CELL& cell, RESULT& result, const Vector<V,3> distance, F f,
  Vector<int,3> pointA, Vector<int,3> pointB) any_platform
{
  const V refDistanceSq = norm_squared(Vector<V,3>(pointA - pointB));
  const V alpha = ((distance - pointB) * (pointA - pointB)) / refDistanceSq;
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  V resNeighbor[DataDim];
  f(cell.neighbor(pointA), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += alpha * resNeighbor[i];
  }
  f(cell.neighbor(pointB), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1 - alpha) * resNeighbor[i];
  }
}

/// @brief interpolation of function f
/// @tparam NEIGHBORS encodes, which neighbors are accessible: e.g. 00001100 refers
/// to the fifth and the sixth neighbor of the point of interest
/// Their order is as follows (w.r.t. the floor lattice point): 10000000=(0,0,0),
/// 1000000=(1,0,0), 100000=(0,1,0), 10000=(1,1,0), 1000=(0,0,1), 100=(1,0,1),
/// 10=(0,1,1), 1=(1,1,1)
/// @tparam DataDim Dimension of function return value
/// @param distance Distance from cell
/// @param result Result is written here
/// @param f Function which is to interpolate
/// @param missingDirections Neighbor points where function is not evaluable.
// Caution: Developer has to guarantee that function is evaluable for all other
// neighbor points on the lattice.
template <unsigned NEIGHBORS, unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate3d(CELL& cell, RESULT& result, const Vector<V,3> distance, F f) any_platform
{
  // eight neighbors: trilinear interpolation
  if constexpr (NEIGHBORS == 11111111) {
    interpolate3d<DataDim,CELL,RESULT,F,V>(cell, result, distance, f);
  }

  // four neighbors (coplanar)
  else if constexpr (NEIGHBORS == 1111) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11110000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 1010101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10101010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 110011) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11001100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10011001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 1100110) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10100101) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1011010) || (NEIGHBORS == 11011010)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11000011) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 111100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_plane<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }

  // three neighbors
  else if constexpr (NEIGHBORS == 111) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 100011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 100101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 101001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 110001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1000011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1000101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1001001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1010001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1100001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10000011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10000101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10001001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10010001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10010001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10100001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11000001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 100110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 101010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 110010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1000110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1001010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1010010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1100010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10000110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10001010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10010010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10100010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11000010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]+1), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 101100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 110100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1001100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1010100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1100100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10001100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10010100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10100100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11000100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 111000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1011000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1101000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10011000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10101000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11001000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1110000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 10110000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11010000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 11100000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  }

  // two neighbors
  else if constexpr (NEIGHBORS == 11) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 100001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1000001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10000001) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 100010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1000010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10000010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 100100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1000100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10000100) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 11000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]+1), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 101000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1001000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]+1), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10001000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 110000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]+1), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1010000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10010000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1100000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 10100000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 11000000) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  }
  // one neighbor: constant extrapolation
  else if constexpr (NEIGHBORS == 1) {
    const Vector<int,3> point (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 10) {
    const Vector<int,3> point (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 100) {
    const Vector<int,3> point (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 1000) {
    const Vector<int,3> point (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 10000) {
    const Vector<int,3> point (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 100000) {
    const Vector<int,3> point (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 1000000) {
    const Vector<int,3> point (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 10000000) {
    const Vector<int,3> point (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    f(cell.neighbor(point), result);
  }

  // four or more linearly independent neighbors: linear interpolation (use at most four points)
  else if constexpr ((NEIGHBORS % 100000 == 10111) || (NEIGHBORS % 100000 == 11111)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS % 1000000 == 100111) || (NEIGHBORS % 1000000 == 101111)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1000111) || (NEIGHBORS == 11000111) || (NEIGHBORS == 11001111) || (NEIGHBORS == 1001111)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 10000111) || (NEIGHBORS == 10001111)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS % 100000 == 11011) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 101011) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1001011) || (NEIGHBORS == 11001011)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10001011) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS % 100000 == 11101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 101101) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1001101) || (NEIGHBORS == 11001101)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10001101) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS % 100000 == 11110) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 101110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1001110) || (NEIGHBORS == 11001110)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10001110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 11101000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11100100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11100010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11100001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 11011000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11010100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 11010010) || (NEIGHBORS == 11110010)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11010001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 10111000) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10110100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10110010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10110001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 1111000) || (NEIGHBORS == 11111000)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1110100) || (NEIGHBORS == 11110100)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 1110010) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1110001) || (NEIGHBORS == 11110001)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 10101100) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10101001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 10100110) || (NEIGHBORS == 11100110)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10100011) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 11001010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 10011010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 1101010) || (NEIGHBORS == 11101010)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 111010) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 1011100) || (NEIGHBORS == 1111100) || (NEIGHBORS == 11111100) || (NEIGHBORS == 11011100)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1010110) || (NEIGHBORS == 11010110)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 1011001) || (NEIGHBORS == 11011001)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1010011) || (NEIGHBORS == 1110011) || (NEIGHBORS == 11110011) || (NEIGHBORS == 11010011)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 11000101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1100101) || (NEIGHBORS == 11100101)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr ((NEIGHBORS == 10010101) || (NEIGHBORS == 11010101)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 110101) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 11001001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS == 11000110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 10011100) || (NEIGHBORS == 10111100)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1101100) || (NEIGHBORS == 11101100)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 110110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr (NEIGHBORS % 1000000 == 111001) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1100011) || (NEIGHBORS == 11100011)) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 10010011) || (NEIGHBORS == 10110011)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }
  else if constexpr (NEIGHBORS == 10010110) {
    const Vector<int,3> pointA (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  } else if constexpr ((NEIGHBORS == 1101001) || (NEIGHBORS == 11101001)) {
    const Vector<int,3> pointA (util::floor(distance[0])+1, util::floor(distance[1]), util::floor(distance[2]));
    const Vector<int,3> pointB (util::floor(distance[0]), util::floor(distance[1])+1, util::floor(distance[2]));
    const Vector<int,3> pointC (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2])+1);
    const Vector<int,3> pointD (util::floor(distance[0])+1, util::floor(distance[1])+1, util::floor(distance[2])+1);
    interpolate3d_help_li<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC, pointD);
  }

  else {
    std::string s = "Invalid NEIGHBORS" + std::to_string(NEIGHBORS);
    throw std::invalid_argument(s);
  }
}

/// Trilinear interpolation of FIELD
/// see above for explanation of the parameters
template <unsigned NEIGHBORS, typename FIELD, typename CELL, typename RESULT, typename V = typename CELL::value_t>
void interpolate3d(CELL& cell, RESULT& result, const Vector<V,3> distance) any_platform
{
  FieldD<V, typename CELL::descriptor_t, FIELD> v;
  constexpr unsigned DataDim = CELL::descriptor_t::template size<FIELD>();
  auto f = [DataDim](const auto& c, RESULT& res) -> void {
    for (unsigned i = 0; i < DataDim; ++i) {
      res[i] = c.template getFieldComponent<FIELD>(i);
    }
  };
  interpolate3d<NEIGHBORS,DataDim,CELL,RESULT,decltype(f),V>(cell, result, distance, f);
}

/// Trilinear interpolation of FIELD
template <typename FIELD, typename CELL, typename RESULT, typename V = typename CELL::value_t>
void interpolate3d(CELL& cell, RESULT& result, const Vector<V,3> distance) any_platform
{
  interpolate3d<11111111,FIELD,CELL,RESULT,V>(cell, result, distance);
}

}

#endif