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

#ifndef COPY_FIELDS_O_H
#define COPY_FIELDS_O_H

#include <core/operator.h>
#include <core/superLatticeCoupling.h>

namespace olb {

namespace couplers {

/// Simple coupling operator for copying field data between two SuperLattice instances.
template <typename FROM, typename TO>
struct CopyFieldsO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  template <typename CELLS>
  void apply(CELLS& cells) any_platform
  {
    // Using data type and descriptor from destination lattice
    using V = typename CELLS::template value_t<names::Lattice2>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Lattice2>::descriptor_t;

    auto& cell_from = cells.template get<names::Lattice>();
    auto& cell_to = cells.template get<names::Lattice2>();

    FieldD<V,DESCRIPTOR,TO> content = cell_from.template getField<FROM>();
    cell_to.template setField<TO>(content);
  }
};

}

/// Copies field data stored in FROM in lattice_from to TO in lattice_to
template <typename FROM, typename TO,
          typename T_FROM, typename DESCRIPTOR_FROM,
          typename T_TO, typename DESCRIPTOR_TO
>
void copyFields(SuperLattice<T_FROM,DESCRIPTOR_FROM>& lattice_from,
                SuperLattice<T_TO,DESCRIPTOR_TO>& lattice_to) {
  SuperLatticeCoupling coupling(couplers::CopyFieldsO<FROM,TO>{},
                       names::Lattice{}, lattice_from,
                       names::Lattice2{}, lattice_to);
  coupling.execute();
}

}
#endif
