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

#ifndef OLB_LINE_LATTICE_3D_H
#define OLB_LINE_LATTICE_3D_H

#include "core/vector.h"
#include "geometry/cuboidDecomposition.h"
#include "line3D.h"

namespace olb {

/// Parametrization of a line3D lattice (i.e. a line lattice).
/**
 * This class provides a common interface for describing how to discretize the intersection
 * of a line3D given by Line3D<T> and the mother cuboid of CuboidDecomposition2D<T>.
 **/
template <typename T>
class LineLattice3D {
private:
  CuboidDecomposition<T,3>&  _geometry;

  /// \return max possible distance
  int computeMaxLatticeDistance() const;
  /// Compute _line3D.origin, _n so that the cuboid is right inside the geometry
  void constructCuboid(int maxLatticeDistance);
  /// Update _h, _n _line3D.u so that the length matches the given resolution
  void setToResolution(int resolution);

protected:
  const Line3D<T> _line3D;

  /// Origin vector of the lattice
  /**
   * Note that this origin is set to a outermost point of the intersection between
   * cuboid geometry and line3D. Thus it is different from the Line3D<T>
   * origin vector in the general case.
   **/
  Vector<T,3> _origin;
  /// Direction vector of the lattice, normalized to grid width _h
  Vector<T,3> _u;

  /// Distance between discrete lattice points
  T _h;
  /// Number of lattice points in the direction of _u
  int _n;

public:
  /// Constructor for automatic discretization.
  /**
   * i.e. the grid width is set to CuboidDecomposition2D<T>::getDeltaR.
   **/
  LineLattice3D(CuboidDecomposition<T,3>& geometry,
                Line3D<T>      line3D);
  /// Constructor for discretization of a given resolution.
  LineLattice3D(CuboidDecomposition<T,3>& geometry,
                Line3D<T>      line3D,
                int                  resolution);
  /// Constructor for discretization of a given grid width.
  LineLattice3D(CuboidDecomposition<T,3>& geometry,
                Line3D<T>      line3D,
                T                    h);

  LineLattice3D(const LineLattice3D&) = default;

  const Line3D<T>& getLine3D() const;

  /// Transform 1d lattice coordinates to their physical 3d location
  Vector<T,3> getPhysR(const int& n) const;

  /// \return _n
  int getN() const;
  /// \return _h
  T getPhysSpacing() const;

  /// \return _origin
  Vector<T,3> getPhysOrigin() const;
  /// \return _u
  Vector<T,3> getVectorU() const;

};

}

#endif
