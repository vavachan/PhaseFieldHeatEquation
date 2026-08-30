/*  This file is part of the OpenLB library
 *
 *   Copyright (C) 2024 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef OLB_ADVECTION_DIFFUSION_DIRICHLET_3D_H
#define OLB_ADVECTION_DIFFUSION_DIRICHLET_3D_H

#include "setBoundary.h"
#include "advectionDiffusionDirichlet.h"
#include "dynamics/advectionDiffusionBoundaries.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct AdvectionDiffusionDirichlet<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
            AdvectionDiffusionBoundariesDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple
    >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::NormalMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          AdvectionDiffusionCornerDynamics3D,MixinDynamics,momenta::EquilibriumBoundaryTuple
    >::construct(n);

  case DiscreteNormalType::ExternalEdge:
    return boundaryhelper::NormalSpecialMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          AdvectionDiffusionEdgesDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple
    >::construct(n);

  default:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
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
