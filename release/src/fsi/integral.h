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

#ifndef FSI_INTEGRAL_H
#define FSI_INTEGRAL_H

#include "fields.h"

namespace olb {

/// Operator for integrating per-element momentum exchange forces in a HLBM-FSI context
template <typename... FIELDS>
struct IntegratePorousElementFieldsO {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  using parameters = meta::list<
    fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>,
    fields::array_of<FIELDS>...,
    fields::fsi::REDUCED_ELEMENTS_COUNT
  >;

  int getPriority() const {
    return 2;
  }

  template <typename BLOCK>
  struct type {
    void setup(BLOCK& blockLattice) { }
    void apply(BLOCK& blockLattice) {
      throw std::runtime_error("IntegratePorousElementFieldsO not implemented");
    }
  };

  template <typename BLOCK>
  void setup(BLOCK& blockLattice) {
    type<BLOCK>{}.setup(blockLattice);
  }

  template <typename BLOCK>
  void apply(BLOCK& blockLattice) {
    type<BLOCK>{}.apply(blockLattice);
  }

};

}

#endif
