/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef FIELD_PARAMETERS_D_HH
#define FIELD_PARAMETERS_D_HH

#include "fieldParametersD.h"

namespace olb {

template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, Platform PLATFORM, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,PLATFORM,PARAMETERS>::getNblock() const
{
  return PARAMETERS::size;
}

template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, Platform PLATFORM, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,PLATFORM,PARAMETERS>::getSerializableSize() const
{
  std::size_t size = 0;
  decltype(parameters)::fields_t::for_each([&size](auto field) {
    using field_t = typename decltype(field)::type;
    size += FieldD<T,DESCRIPTOR,field_t>{}.getSerializableSize();
  });
  return size;
}

template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, Platform PLATFORM, typename PARAMETERS>
bool* ConcreteParametersD<T,DESCRIPTOR,PLATFORM,PARAMETERS>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;
  decltype(parameters)::fields_t::for_each([&](auto field) {
    using field_t = typename decltype(field)::type;
    if constexpr (DESCRIPTOR::template size<field_t>() == 1) {
      registerVar(iBlock, sizeBlock, currentBlock, dataPtr,
                  parameters.template get<field_t>(), loadingMode);
    } else {
      registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr,
                                      parameters.template get<field_t>(), loadingMode);
    }
  });
  return dataPtr;
}


}

#endif
