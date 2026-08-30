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

//This file contains the Slip Boundary
//This is an onLattice boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SLIP_3D_H
#define SLIP_3D_H

#include "slip.h"
#include "setBoundary.h"

namespace olb {

namespace boundary {

template <int NX, int NY, int NZ> struct FullSlipO3D;

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 3)
struct FullSlip<T,DESCRIPTOR> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  if (n != 0) {
    return std::nullopt;
  } else {
    return meta::id<olb::BounceBack<T,DESCRIPTOR>>();
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  if (n != 0) {
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FullSlipO3D>(n);
  }
  return std::nullopt;
}

};

}

}

#endif
