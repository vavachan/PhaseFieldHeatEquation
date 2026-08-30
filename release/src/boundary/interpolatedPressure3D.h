/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

//This file contains the Interpolated Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef INTERPOLATED_PRESSURE_3D_H
#define INTERPOLATED_PRESSURE_3D_H

#include "setBoundary.h"
#include "interpolatedPressure.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct InterpolatedPressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 2;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    if (n[0] != 0 && n[0] == -1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<0,-1>>>{};
    }
    else if (n[0] != 0 && n[0] == 1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<0,1>>>{};
    }
    else if (n[1] != 0 && n[1] == -1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<1,-1>>>{};
    }
    else if (n[1] != 0 && n[1] == 1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<1,1>>>{};
    }
    else if (n[2] != 0 && n[2] == -1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<2,-1>>>{};
    }
    else if (n[2] != 0 && n[2] == 1) {
      return meta::id<typename MixinDynamics::template exchange_momenta<momenta::BasicDirichletPressureBoundaryTuple<2,1>>>{};
    }
    else {
      return std::nullopt;
    }

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForDirectionOrientation<T,DESCRIPTOR,PlaneFdBoundaryProcessor3D
          >(n);

  default:
    return std::nullopt;
  }
}

};

}

}

#endif
