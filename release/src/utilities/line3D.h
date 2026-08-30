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

#ifndef OLB_LINE_3D_H
#define OLB_LINE_3D_H

#include "core/vector.h"
#include "geometry/cuboid.h"

namespace olb {

/// Definition of a analytical line embedded in 3D space
/**
 * Line3D defines a line using its origin and a direction vector.
 **/
template <typename T>
struct Line3D {
  Vector<T,3> origin;
  Vector<T,3> u;
  Vector<T,3> normal;

  Line3D() = default;

  /// Center the line at the given origin vector
  /// \return Line3D reference for further construction
  Line3D& originAt(const Vector<T,3>& origin);
  /// Center the line relative to the given cuboid
  /// \return Line3D reference for further construction
  Line3D& centeredIn(const Cuboid3D<T>& cuboid);
  /// Set the direction of the line parallel to a vector
  /// \return Line3D reference for further construction
  Line3D& parallelTo(const Vector<T,3>& direction);
  /// Calculate the direction vector of the line to be orthogonal to the given normal
  /// \return Line3D reference for further construction
  Line3D& normalTo(const Vector<T,3>& normal);

  /// \return true iff normal is orthogonal to X axis
  bool isParallelToX() const;
  /// \return true iff normal is orthogonal to Y axis
  bool isParallelToY() const;

  bool isParallelToZ() const;
};

}

#endif
