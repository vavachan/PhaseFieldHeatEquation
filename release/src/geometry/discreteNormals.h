/**  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Mathias Krause, Simon Zimny, Adrian Kummerlaender
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
**/

#ifndef GEOMETRY_DISCRETE_NORMAL_H
#define GEOMETRY_DISCRETE_NORMAL_H

namespace olb {

template <typename T> class BlockIndicatorF2D;
template <typename T> class BlockIndicatorF3D;

/// Type associated with a discrete normal vector
enum class DiscreteNormalType : int {
  Flat           = 0, /// Normal detected as flat plane
  ExternalCorner = 1, /// Normal detected as external corner
  InternalCorner = 2, /// Normal detected as internal corner
  ExternalEdge   = 3, /// Normal detected as external edge (only 3D)
  InternalEdge   = 4  /// Normal detected as internal edge (only 3D)
};

template <concepts::Descriptor DESCRIPTOR>
using DiscreteNormal = Vector<int,DESCRIPTOR::d>;

/// Returns type (e.g. edge / corner) and discrete normal in 2D
template<concepts::BaseType T>
std::pair<DiscreteNormalType,Vector<int,2>> computeBoundaryTypeAndNormal(
  BlockIndicatorF2D<T>& fluidI,
  BlockIndicatorF2D<T>& outsideI,
  Vector<int,2> latticeR);

/// Returns type (e.g. edge / corner) and discrete normal in 3D
template<concepts::BaseType T>
std::pair<DiscreteNormalType,Vector<int,3>> computeBoundaryTypeAndNormal(
  BlockIndicatorF3D<T>& fluidI,
  BlockIndicatorF3D<T>& outsideI,
  Vector<int,3> latticeR);

}

#endif
