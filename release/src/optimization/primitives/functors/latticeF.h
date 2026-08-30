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

#ifndef LATTICE_F_H
#define LATTICE_F_H

#include "../concept.h"

namespace olb {

namespace functors {

/// Computes the shifted populations
struct PopulationF {
  using parameters = meta::list<>;

  using result_t = descriptors::POPULATION;
  using fields_t = meta::list<descriptors::POPULATION>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    auto pop = cell.template getField<descriptors::POPULATION>();
    for (auto iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      pop[iPop] += descriptors::t<V,DESCRIPTOR>(iPop);
    }
    return pop;
  }
};

}

}

#endif
