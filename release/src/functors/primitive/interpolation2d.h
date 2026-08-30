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

#ifndef INTERPOLATION_2D_H
#define INTERPOLATION_2D_H

#include "core/vector.h"

namespace olb {

/// Bilinear interpolation between four neighboring mesh points
/// @tparam DataDim Dimension of function return value
/// @param result Result is written here
/// @param distance Distance from cell (in lattice units)
/// @param f Function which is to interpolate
// Caution: Developer has to guarantee that the data are accessible.
// It is expected that the point of interest (=cell+distance) lies nearby cell
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate2d(CELL& cell, RESULT& result, const Vector<V,2> distance, F f) any_platform
{
  const Vector<int,2> floorV (util::floor(distance[0]), util::floor(distance[1]));
  Vector<Vector<int,2>,4> surroundingPoints (floorV);
  surroundingPoints[1][0] += 1;
  surroundingPoints[2][1] += 1;
  surroundingPoints[3][0] += 1; surroundingPoints[3][1] += 1;

  V resNeighbor[DataDim];
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  for (auto point : surroundingPoints) {
    const Vector<V,2> dist = distance - point;
    const V weight = (V(1) - util::abs(dist[0])) * (V(1) - util::abs(dist[1]));
    f(cell.neighbor(point), resNeighbor);
    for (unsigned i = 0; i < DataDim; ++i) {
      result[i] += weight * resNeighbor[i];
    }
  }
}

/// Linear interpolation between three points
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate2d_help(CELL& cell, RESULT& result, const Vector<V,2> distance, F f,
  Vector<int,2> pointA, Vector<int,2> pointB, Vector<int,2> pointC) any_platform
{
  // express point as a linear combination of the two sides
  const auto lc = util::solveLinearSystem(Vector<V,2>(pointB - pointA),
    Vector<V,2>(pointC - pointA),
    distance - pointA);
  // evaluate linear function
  V resNeighbor[DataDim];
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
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

/// Linear interpolation between two points (constant in the orthogonal direction)
template <unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate2d_help(CELL& cell, RESULT& result, const Vector<V,2> distance, F f,
  Vector<int,2> pointA, Vector<int,2> pointB) any_platform
{
  const V refDistanceSq = norm_squared(Vector<V,2>(pointA - pointB));
  const V alpha = ((distance - pointB) * (pointA - pointB)) / refDistanceSq;
  V resNeighbor[DataDim];
  for (unsigned i=0; i<DataDim; ++i) {
    result[i] = V();
  }
  f(cell.neighbor(pointA), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += alpha * resNeighbor[i];
  }
  f(cell.neighbor(pointB), resNeighbor);
  for (unsigned i = 0; i < DataDim; ++i) {
    result[i] += (1 - alpha) * resNeighbor[i];
  }
}

/// @brief linear interpolation/ extrapolation of function f
/// @tparam DataDim Dimension of function return value
/// @tparam NEIGHBORS encodes, which neighbors are accessible: e.g. 1100 refers
/// to the first two neighbors of the point of interest (in lexicographic order).
/// Hence, their order is as follows (w.r.t. the floor lattice point): 1000=(0,0),
/// 100=(0,1), 10=(1,0), 1=(1,1).
/// @param result Result is written here
/// @param distance Distance from cell (in lattice units)
/// @param f Function which is to interpolate
// Caution: Developer has to guarantee that the data are accessible.
// It is expected that the point of interest (=cell+distance) lies nearby cell
template <unsigned NEIGHBORS, unsigned DataDim, typename CELL, typename RESULT, typename F, typename V = typename CELL::value_t>
void interpolate2d(CELL& cell, RESULT& result, const Vector<V,2> distance, F f) any_platform
{
  // four neighbors: bilinear interpolation
  if constexpr (NEIGHBORS == 1111) {
    interpolate2d<DataDim,CELL,RESULT,F,V>(cell, result, distance, f);
  }
  // three neighbors: linear interpolation
  else if constexpr (NEIGHBORS == 1110) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1]));
    const Vector<int,2> pointC (util::floor(distance[0]), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1101) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1]));
    const Vector<int,2> pointC (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 1011) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]), util::floor(distance[1])+1);
    const Vector<int,2> pointC (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  } else if constexpr (NEIGHBORS == 111) {
    const Vector<int,2> pointA (util::floor(distance[0]+1), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]), util::floor(distance[1])+1);
    const Vector<int,2> pointC (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB, pointC);
  }
  else if constexpr (NEIGHBORS == 1100) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1]));
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 11) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1])+1);
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1010) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 101) {
    const Vector<int,2> pointA (util::floor(distance[0]+1), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 1001) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1]));
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1])+1);
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  } else if constexpr (NEIGHBORS == 110) {
    const Vector<int,2> pointA (util::floor(distance[0]), util::floor(distance[1])+1);
    const Vector<int,2> pointB (util::floor(distance[0]+1), util::floor(distance[1]));
    interpolate2d_help<DataDim,CELL,RESULT,F,V>(cell, result, distance, f, pointA, pointB);
  }
  // one neighbor: constant extrapolation
  else if constexpr (NEIGHBORS == 1) {
    const Vector<int,2> point (util::floor(distance[0])+1, util::floor(distance[1])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 10) {
    const Vector<int,2> point (util::floor(distance[0]), util::floor(distance[1])+1);
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 100) {
    const Vector<int,2> point (util::floor(distance[0])+1, util::floor(distance[1]));
    f(cell.neighbor(point), result);
  } else if constexpr (NEIGHBORS == 1000) {
    const Vector<int,2> point (util::floor(distance[0]), util::floor(distance[1]));
    f(cell.neighbor(point), result);
  }
  else {
    throw std::invalid_argument("Invalid NEIGHBORS");
  }
}

/// @brief linear interpolation/ extrapolation of FIELD
/// @tparam NEIGHBORS encodes, which neighbors are accessible: e.g. 1100 refers
/// to the first two neighbors of the point of interest (in lexicographic order)
/// Hence, their order is as follows (w.r.t. the floor lattice point): 1000=(0,0),
/// 100=(0,1), 10=(1,0), 1=(1,1).
/// @param distance Distance between point of interest and cell (in lattice units)
/// @param result Result is written here
// Caution: Developer has to guarantee that the data are accessible.
// It is expected that the point of interest lies nearby cell
template <unsigned NEIGHBORS, typename FIELD, typename CELL, typename RESULT, typename V = typename CELL::value_t>
void interpolate2d(CELL& cell, RESULT& result, const Vector<V,2> distance) any_platform
{
  FieldD<V, typename CELL::descriptor_t, FIELD> v;
  constexpr unsigned DataDim = CELL::descriptor_t::template size<FIELD>();
  auto f = [DataDim](const auto& c, RESULT& res) -> void {
    for (unsigned i = 0; i < DataDim; ++i) {
      res[i] = c.template getFieldComponent<FIELD>(i);
    }
  };
  interpolate2d<NEIGHBORS,DataDim,CELL,RESULT,decltype(f),V>(cell, result, distance, f);
}

/// @brief bilinear interpolation/ extrapolation of FIELD
/// @param distance Distance between point of interest and cell (in lattice units)
/// @param result Result is written here
// Caution: Developer has to guarantee that the data are accessible.
// It is expected that the point of interest lies nearby cell
template <typename FIELD, typename CELL, typename RESULT, typename V = typename CELL::value_t>
void interpolate2d(CELL& cell, RESULT& result, const Vector<V,2> distance) any_platform
{
  interpolate2d<1111,FIELD,CELL,RESULT,V>(cell, result, distance);
}


}

#endif