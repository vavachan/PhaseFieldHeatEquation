/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Dennis Teutscher
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

#ifndef OLB_LINE_3D_HH
#define OLB_LINE_3D_HH

#include "line3D.h"
#include "core/olbDebug.h"
#include "utilities/vectorHelpers.h"

namespace olb {


template <typename T>
Line3D<T>& Line3D<T>::originAt(const Vector<T,3>& o)
{
  origin[0] = o[0] - 2*std::numeric_limits<T>::epsilon()*util::fabs(o[0]);
  origin[1] = o[1] - 2*std::numeric_limits<T>::epsilon()*util::fabs(o[1]);
  origin[2] = o[2] - 2*std::numeric_limits<T>::epsilon()*util::fabs(o[2]);

  return *this;
}

template <typename T>
Line3D<T>& Line3D<T>::centeredIn(const Cuboid3D<T>& cuboid)
{
  const Vector<T,3>& cuboidOrigin = cuboid.getOrigin();
  const Vector<int,3>& extend     = cuboid.getExtent();
  const T deltaR = cuboid.getDeltaR();

  origin[0] = (cuboidOrigin[0] + 0.5 * deltaR * extend[0]);
  origin[1] = (cuboidOrigin[1] + 0.5 * deltaR * extend[1]);
  origin[2] = (cuboidOrigin[2] + 0.5 * deltaR * extend[2]);
  origin[0] -= 2*std::numeric_limits<T>::epsilon()*util::fabs(origin[0]);
  origin[1] -= 2*std::numeric_limits<T>::epsilon()*util::fabs(origin[1]);
  origin[2] -= 2*std::numeric_limits<T>::epsilon()*util::fabs(origin[2]);

  return *this;
}

template <typename T>
Line3D<T>& Line3D<T>::parallelTo(const Vector<T,3>& direction)
{
  u = direction;
  // help vektor which is not parallel to u
  Vector<T,3> helper = (u[0] == 0 && u[1] == 0) ? Vector<T,3>({1, 0, 0}) : Vector<T,3>({0, 0, 1});

  normal = normalize(crossProduct(u, helper));

  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u, normal)));
  return *this;
}

template <typename T>
Line3D<T>& Line3D<T>::normalTo(const Vector<T,3>& n)
{
  normal = normalize(n);
  // Check if normal is aligned with any of the axes
  if (util::nearZero(normal[0]) && util::nearZero(normal[1])) {
    // Normal is along the z-axis
    u = {T(1), T(0), T(0)};
  } else if (util::nearZero(normal[0]) && util::nearZero(normal[2])){
    // Normal is along the y-axis
    u = {T(1), T(0), T(0)};
  } else if (util::nearZero(normal[1]) && util::nearZero(normal[2])) {
    // Normal is along the x-axis
    u = {T(0), T(1), T(0)};
  } else {
    // Use cross product with a vector not parallel to the normal to find u
    Vector<T,3> arbitraryVec = (util::nearZero(normal[0]) && util::nearZero(normal[1]))
                                ? Vector<T,3>{T(1), T(0), T(0)}
                                : Vector<T,3>{T(0), T(1), T(0)};
    u = crossProduct(normal, arbitraryVec);// Calculate orthogonal vector
  }
  u = normalize(u);
  OstreamManager clout ("Test");
  clout <<"Vector"<<u <<std::endl;
  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u, normal)));
  return *this;
}

template <typename T>
bool Line3D<T>::isParallelToX() const
{
  return util::nearZero(util::dotProduct3D(normal, {1,0,0}));
}

template <typename T>
bool Line3D<T>::isParallelToY() const
{
  return util::nearZero(util::dotProduct3D(normal, {0,1,0}));
}
template <typename T>
bool Line3D<T>::isParallelToZ() const
{
  return util::nearZero(util::dotProduct3D(normal, {0,0,1}));
}


}

#endif
