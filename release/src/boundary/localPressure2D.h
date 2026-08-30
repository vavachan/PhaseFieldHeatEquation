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

//This file contains the Local Pressure Boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_LOCAL_PRESSURE_BOUNDARY_2D_H
#define SET_LOCAL_PRESSURE_BOUNDARY_2D_H

#include "setBoundary2D.h"
#include "localPressure.h"

namespace olb {
namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct LocalPressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 0;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type){
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
     >::construct(n);
  default:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return std::nullopt;
}

};

}

}


#endif
