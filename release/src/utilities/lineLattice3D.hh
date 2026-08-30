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

#ifndef LINE_LATTICE_3D_HH
#define LINE_LATTICE_3D_HH

#include "lineLattice3D.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T>
int LineLattice3D<T>::computeMaxLatticeDistance() const
{
  OstreamManager clout ("ComputeMaxLatticeDistance");
  const Cuboid3D<T>&  cuboid = _geometry.getMotherCuboid();
  const Vector<T,3>   origin = cuboid.getOrigin();
  const Vector<int,3> extend = cuboid.getExtent();
  const T             deltaR = cuboid.getDeltaR();

  T maxPhysDistance = T();
  T tmp;
  Vector<T,3> tmpO;
  Vector<T,3> tmpE;

  for (int iDim=0; iDim<3; ++iDim) {
    tmpO[iDim] = origin[iDim] - _origin[iDim];
    tmpE[iDim] = origin[iDim] + extend[iDim]*deltaR - _origin[iDim];
  }

  tmp = util::sqrt(tmpE[0]*tmpE[0] + tmpE[1]*tmpE[1]+tmpE[2]*tmpE[2]);
  if (maxPhysDistance < tmp) {
    maxPhysDistance = tmp;
  }
  return int(maxPhysDistance/_h) + 1;
}

template <typename T>
void LineLattice3D<T>::constructCuboid(int maxLatticeDistance)
{
  int min = -maxLatticeDistance;
  int max = maxLatticeDistance;

  for ( int i = -maxLatticeDistance; i < maxLatticeDistance; ++i ) {
    if (auto iC = _geometry.getC(getPhysR(i))) {
      min = i;
      break;
    }
  }
  for ( int i = maxLatticeDistance; i > -maxLatticeDistance; --i ) {
    if (auto iC = _geometry.getC(getPhysR(i))) {
      max = i;
      break;
    }
  }

  _n = max - min + 1;
  _origin = _origin + T(min)*_u;
}

template <typename T>
void LineLattice3D<T>::setToResolution(int resolution)
{
  T newH = _n*_h/(T) resolution;
  _n = resolution;
  _h = newH;
  _u = normalize(_u, _h);
}

template<typename T>
LineLattice3D<T>::LineLattice3D(
  CuboidDecomposition<T,3>& geometry, Line3D<T> line3D)
  : _geometry(geometry),
    _line3D(line3D),
    _origin(line3D.origin),
    _u(line3D.u),
    _h(geometry.getDeltaR())
{
  _u = normalize(_u, _h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _line3D.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);
}

template<typename T>
LineLattice3D<T>::LineLattice3D(
  CuboidDecomposition<T,3>& geometry, Line3D<T> line3D, int resolution)
  : _geometry(geometry),
    _line3D(line3D),
    _origin(line3D.origin),
    _u(line3D.u),
    _h(geometry.getDeltaR())
{
  _u = normalize(_u, _h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _line3D.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);

  if ( resolution > 0 ) {
    setToResolution(resolution);
  }
}

template<typename T>
LineLattice3D<T>::LineLattice3D(
  CuboidDecomposition<T,3>& geometry, Line3D<T> line3D, T h)
  : _geometry(geometry),
    _line3D(line3D),
    _origin(line3D.origin),
    _u(line3D.u),
    _h(h)
{
  if ( util::nearZero(_h) ) {
    _h = _geometry.getDeltaR();
  }

  _u = normalize(_u, _h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _line3D.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);
}

template <typename T>
const Line3D<T>& LineLattice3D<T>::getLine3D() const
{
  return _line3D;
}

template <typename T>
Vector<T,3> LineLattice3D<T>::getPhysR(const int& n) const
{
  return Vector<T,3> {
    _origin[0] + T(n)*_u[0],
    _origin[1] + T(n)*_u[1],
    _origin[2] + T(n)*_u[2]
  };
}

template <typename T>
int LineLattice3D<T>::getN() const
{
  return _n;
}

template <typename T>
T LineLattice3D<T>::getPhysSpacing() const
{
  return _h;
}

template <typename T>
Vector<T,3> LineLattice3D<T>::getPhysOrigin() const
{
  return _origin;
}

template <typename T>
Vector<T,3> LineLattice3D<T>::getVectorU() const
{
  return _u;
}


}

#endif
