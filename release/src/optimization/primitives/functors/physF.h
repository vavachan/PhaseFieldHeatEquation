/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef PHYS_F_H
#define PHYS_F_H

#include "../concept.h"

namespace olb {

namespace functors {

/// Computes the physical velocity
struct VelocityF {
  using parameters = meta::list<descriptors::CONVERSION>;

  using result_t = descriptors::VELOCITY;
  using fields_t = meta::list<descriptors::VELOCITY>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    V conversion = parameters.template get<descriptors::CONVERSION>();
    V u[DESCRIPTOR::d]{0};
    cell.computeU(u);

    FieldD<V,DESCRIPTOR,descriptors::VELOCITY> physU;
    for (auto iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      physU[iDim] = u[iDim] * conversion;
    }
    return physU;
  }
};

}

}

#endif
