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

#ifndef UTILITIES_TYPE_MAP_H
#define UTILITIES_TYPE_MAP_H

#include "core/meta.h"

namespace olb {

namespace meta {

template <typename...>
struct unzip_flattened_keys;

template <typename KEY, typename VALUE, typename... TAIL>
struct unzip_flattened_keys<KEY,VALUE,TAIL...> {
  static_assert(sizeof...(TAIL) % 2 == 0, "TAIL must be valid map size");
  using type = typename unzip_flattened_keys<TAIL...>::type::template push<KEY>;
};

template <>
struct unzip_flattened_keys<> {
  using type = list<>;
};

template <typename...>
struct unzip_flattened_values;

template <typename KEY, typename VALUE, typename... TAIL>
struct unzip_flattened_values<KEY,VALUE,TAIL...> {
  static_assert(sizeof...(TAIL) % 2 == 0, "TAIL must be valid map size");
  using type = typename unzip_flattened_values<TAIL...>::type::template push<VALUE>;
};

template <>
struct unzip_flattened_values<> {
  using type = list<>;
};

template <typename KEYS, typename VALUES>
struct plain_map {
  using keys_t   = KEYS;
  using values_t = VALUES;

  static_assert(keys_t::size == values_t::size,
                "Count of keys and values must match");

  static constexpr unsigned size = keys_t::size;

  template <typename KEY>
  using value = typename values_t::template get<(keys_t::template index<KEY>())>;

  template <template<typename> typename F>
  using map_values = plain_map<keys_t, typename values_t::template map<F>>;

  template <typename F>
  static constexpr void for_each(F f) {
    KEYS::for_each([&](auto key) {
      f(key, meta::id<value<typename decltype(key)::type>>{});
    });
  }

};

/// Map of types
template <typename... KVs>
using map = plain_map<
  typename unzip_flattened_keys<KVs...>::type,
  typename unzip_flattened_values<KVs...>::type
>;

}

}

#endif
