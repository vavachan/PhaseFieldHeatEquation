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

#ifndef OLB_BOUNDARY_EXTERNAL_FIELD_H
#define OLB_BOUNDARY_EXTERNAL_FIELD_H

#include "setBoundary.h"

namespace olb {

template <typename FIELD_A, typename FIELD_B>
struct ExternalFieldO {

template <int... NORMAL>
struct normal {

static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

using parameters = typename meta::list<descriptors::LATTICE_TIME>;

int getPriority() const {
  return 0;
}

template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& params) any_platform {
  const std::size_t iT = params.template get<descriptors::LATTICE_TIME>();
  auto neighbor = cell.neighbor({NORMAL...});
  if (iT % 2 == 0) {
    cell.template setField<FIELD_A>(neighbor.template getField<FIELD_B>());
  } else {
    cell.template setField<FIELD_B>(neighbor.template getField<FIELD_A>());
  }
}

};

};

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR,
          typename FIELD_A, typename FIELD_B>
struct ExternalField {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  return std::nullopt;
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return boundaryhelper::promiseForNormal<PostProcessorPromise<T,DESCRIPTOR>,
                                          ExternalFieldO<FIELD_A,FIELD_B>::template normal>(-1*n);
}

};

}

}

#endif
