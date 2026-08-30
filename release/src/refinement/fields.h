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

#ifndef REFINEMENT_FIELDS_H
#define REFINEMENT_FIELDS_H

#include "descriptor/fields.h"

namespace olb {

namespace fields::refinement {

struct CELL_ID_COARSE : public descriptors::TYPED_FIELD_BASE<CellID,1> { };
struct CELL_ID_FINE : public descriptors::TYPED_FIELD_BASE<CellID,1> { };

struct NORMAL : public descriptors::FIELD_BASE<0,1> { };

struct PREV_RHO  : public descriptors::FIELD_BASE<1> { };
struct PREV_U    : public descriptors::FIELD_BASE<0,1> { };
struct PREV_FNEQ : public descriptors::FIELD_BASE<0,0,1> { };

/// Field for storing indices to context data neighbors
struct CONTEXT_NEIGHBORS : public descriptors::NEIGHBOR_FIELD {
  template <typename T>
  using value_type = CellID;

  template <typename T>
  using column_type = AbstractColumn<int>;

  template <typename DESCRIPTOR>
  requires (DESCRIPTOR::d == 2)
  static constexpr auto count() any_platform {
    return descriptors::q<descriptors::D2Q9<>>()-1;
  }

  template <typename DESCRIPTOR>
  requires (DESCRIPTOR::d == 2)
  static constexpr auto c(unsigned i) any_platform {
    return descriptors::c<descriptors::D2Q9<>>(i+1);
  }

  template <typename DESCRIPTOR>
  requires (DESCRIPTOR::d == 3)
  static constexpr auto count() any_platform {
    return descriptors::q<descriptors::D3Q27<>>()-1;
  }

  template <typename DESCRIPTOR>
  requires (DESCRIPTOR::d == 3)
  static constexpr auto c(unsigned i) any_platform {
    return descriptors::c<descriptors::D3Q27<>>(i+1);
  }

};

}

namespace refinement {

template <typename DESCRIPTOR>
using DATA_DESCRIPTOR = descriptors::LATTICE_DESCRIPTOR<DESCRIPTOR::d, DESCRIPTOR::q,
  fields::refinement::CELL_ID_COARSE,
  fields::refinement::CELL_ID_FINE
>;

}

}

#endif
