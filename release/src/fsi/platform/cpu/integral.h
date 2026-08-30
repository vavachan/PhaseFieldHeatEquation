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

#ifndef FSI_PLATFORM_CPU_INTEGRAL_H
#define FSI_PLATFORM_CPU_INTEGRAL_H

#include "fsi/integral.h"

namespace olb {

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
requires (isPlatformCPU(PLATFORM))
struct IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> {

void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>();
}

void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& blockLattice) {
  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>().parameters;

  const auto& tag = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();

  std::map<int, ParametersD<T,DESCRIPTOR,FIELDS...>> data;
  // Sum up per-element fields
  blockLattice.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
    const std::size_t iCell = blockLattice.getCellId(latticeR);
    if (tag[0][iCell] > 0) {
      auto& fields = data[tag[0][iCell]];
      ((fields.template set<FIELDS>(fields.template get<FIELDS>() + blockLattice.template getField<FIELDS>().getRow(iCell))), ...);
    }
  });

  auto reducedElementTag = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();

  std::vector<int> tags;
  tags.reserve(data.size());
  for (int tag : std::views::keys(data)) {
    tags.emplace_back(tag);
  }
  std::sort(std::begin(tags), std::end(tags));

  for (std::size_t iReduced=0; iReduced < tags.size(); ++iReduced) {
    reducedElementTag[iReduced] = tags[iReduced];
    const auto& fields = data[tags[iReduced]];
    auto set = [&]<typename FIELD>(meta::id<FIELD>, const FieldD<T,DESCRIPTOR,FIELD>& value) {
      auto reduced = parameters.template get<fields::array_of<FIELD>>();
      if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
        reduced[iReduced] = value[0];
      } else {
        for (unsigned iD=0; iD < value.d; ++iD) {
          reduced[iD][iReduced] = value[iD];
        }
      }
    };
    ((set(meta::id<FIELDS>{}, fields.template get<FIELDS>())), ...);
  }

  parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(tags.size());
}

};

}

#endif
