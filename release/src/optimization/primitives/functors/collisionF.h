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

#ifndef COLLISION_F_H
#define COLLISION_F_H

namespace olb {

namespace functors {

/// Execute a collision for given dynamics for a cell
template <typename DYNAMICS>
struct CollisionF {
  using parameters = typename DYNAMICS::parameters;
  using result_t = typename descriptors::POPULATION;
  using fields_t = meta::list<>; // Is this correct?

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    // Instantiate copy of passed cell
    CellD<V,DESCRIPTOR> copyCell;
    DESCRIPTOR::fields_t::for_each([&](auto id) {
      using FIELD = typename decltype(id)::type;
      copyCell.template setField<FIELD>(FieldD<V,DESCRIPTOR,FIELD>(cell.template getField<FIELD>()));
    });

    // Apply dynamics collision on copy cell
    typename DYNAMICS::CollisionO{}.apply(copyCell, parameters);
    return copyCell.template getField<descriptors::POPULATION>();
  }
};

}

}

#endif
