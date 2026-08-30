/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Luiz Eduardo Czelusniak, Tim Bingert
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

//This file contains the Regularized Temperature Boundary 2D
//This is a new version of the Boundary, which only contains free floating functions

#ifndef OLB_REGULARIZED_TEMPERATURE_2D_H
#define OLB_REGULARIZED_TEMPERATURE_2D_H

#include "setBoundary.h"
#include "regularizedTemperature.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct RegularizedTemperature<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
    CombinedAdvectionDiffusionRLBdynamics,MixinDynamics,momenta::RegularizedTemperatureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::NormalMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
    AdvectionDiffusionCornerDynamics2D,MixinDynamics,momenta::EquilibriumBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::InternalCorner:
    return std::nullopt;

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
