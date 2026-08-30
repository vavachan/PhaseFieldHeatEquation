/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef OLB_BOUNDARY_ZERO_DISTRIBUTION_H
#define OLB_BOUNDARY_ZERO_DISTRIBUTION_H

#include "setBoundary.h"

namespace olb {

template <int... NORMAL>
struct ZeroDistributionBoundaryO {

static constexpr OperatorScope scope = OperatorScope::PerCell;

int getPriority() const {
  return 0;
}

template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
  using V = typename CELL::value_t;
  using DESCRIPTOR = typename CELL::descriptor_t;
  FieldD<V,DESCRIPTOR,descriptors::POPULATION> pop{};
  for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    const bool reset = descriptors::c<DESCRIPTOR>(iPop)*Vector{NORMAL...} > 0;
    pop[iPop] = !reset *  cell[iPop]
              +  reset * -descriptors::t<V,DESCRIPTOR>(iPop);
  }
  cell.template setField<descriptors::POPULATION>(pop);
}

};

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
struct ZeroDistribution {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 0;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  return std::nullopt;
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return boundaryhelper::promiseForNormal<PostProcessorPromise<T,DESCRIPTOR>,
                                          ZeroDistributionBoundaryO>(-1*n);
}

};

}

}

#endif
