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

#ifndef CORE_FIELDS_H
#define CORE_FIELDS_H

#include "descriptor/fields.h"

#include <source_location>

namespace olb {

template <typename FIELD>
constexpr std::string_view getFieldName() {
#ifndef USING_LEGACY_CODEGEN
  const std::string_view raw = std::source_location::current().function_name();
  return std::string_view(raw.cbegin() + raw.find_first_of('=')+2,
                          raw.cbegin() + std::min(raw.find_first_of(']'),
                                                  raw.find_first_of(';')));
#else
  return std::string_view();
#endif
}

namespace fields {

/// Returns name of FIELD for human consumption
template <typename FIELD>
std::string name() {
  auto raw = getFieldName<FIELD>();
  if (raw.starts_with("olb::")) {
    raw = std::string_view(raw.cbegin() + 5,
                           raw.cend());
  }
  // Call root namespace functions in order to display full sub-namespace
  // alongside FIELD name
  return std::string(raw);
}
template <typename FIELD>
using array_of = descriptors::POINTER_FIELD_BASE<FIELD>;

template <typename FIELD>
using identity = FIELD;

/// Copy of FIELD without inheritance (used for checkpointing)
template <typename FIELD>
struct copy_of {
  template <typename T>
  using value_type =  typename FIELD::template value_type<T>;

  template <typename T>
  using column_type = typename FIELD::template column_type<T>;

  template <unsigned D, unsigned Q>
  static constexpr unsigned size() {
    return FIELD::template size<D,Q>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return FIELD::template getInitialValue<T,DESCRIPTOR>();
  }

  static constexpr bool isSerializable() {
    return FIELD::isSerializable();
  }
};

}

}

#endif
