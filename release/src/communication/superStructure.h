/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod, Albert Mink, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_STRUCTURE_H
#define SUPER_STRUCTURE_H

#include "utilities/aliases.h"

#include "geometry/cuboidDecomposition.h"

#include "communication/loadBalancer.h"

namespace olb {


template<typename T, unsigned D>
class SuperStructure {
protected:
  /// The grid structure is referenced here
  CuboidDecomposition<T,D>& _cuboidDecomposition;
  /// Distribution of the cuboids of the cuboid structure
  LoadBalancer<T>& _loadBalancer;
  /// Size of ghost cell layer (must be greater than 1 and
  /// greater_overlapBC, default =1)
  int _overlap;
  /// class specific output stream
  mutable OstreamManager clout;

public:
  using value_t = T;

  /// Virtual Destructor for inheritance
  virtual ~SuperStructure() {};

  /// Construction of a super structure
  SuperStructure(CuboidDecomposition<T,D>& cuboidDecomposition,
                 LoadBalancer<T>& loadBalancer, int overlap = 2);

  /// Read and write access to cuboid geometry
  CuboidDecomposition<T,D>& getCuboidDecomposition();
  /// Read only access to cuboid geometry
  const CuboidDecomposition<T,D>& getCuboidDecomposition() const;

  /// Read and write access to the overlap
  int getOverlap();
  /// Read only access to the overlap
  int getOverlap() const;

  /// Read and write access to the load balancer
  LoadBalancer<T>& getLoadBalancer();
  /// Read only access to the load balancer
  LoadBalancer<T> const& getLoadBalancer() const;

  virtual void communicate() { };

  /// Iterate over discrete physical locations
  template <typename F>
  void forCorePhysLocations(F f) const;

  /// Iterate over discrete physical locations between min and max
  template <typename F>
  void forCorePhysLocations(PhysR<T,D> min, PhysR<T,D> max, F f) const;

  /// Iterate over spatial locations
  /// NOTE: Based on physical locations (as opposed to its blockStructure version)
  template <typename F>
  void forCoreSpatialLocations(F f) const;

  /// Iterate over spatial locations between min and max
  /// NOTE: Based on physical locations (as opposed to its blockStructure version)
  template <typename F>
  void forCoreSpatialLocations(PhysR<T,D> min, PhysR<T,D> max, F f) const;

};


}

#endif
