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

//This file contains the ZouHe Pressure Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef OLB_ZOUHE_PRESSURE_2D_H
#define OLB_ZOUHE_PRESSURE_2D_H

#include "zouHePressure.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct ZouHePressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}
std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
      return boundaryhelper::DirectionOrientationMixinDynamicsForDirectionOrientationMomenta
        <T,DESCRIPTOR,ZouHeDynamics,MixinDynamics,momenta::BasicDirichletPressureBoundaryTuple>
        ::construct(n);
  default:
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
