/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-22 Adrian Kummerlaender, Julius Jessberger
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

#ifndef TYPE_INDEXED_CONTAINERS_H
#define TYPE_INDEXED_CONTAINERS_H

#include <vector>
#include <optional>
#include <tuple>

#include "core/meta.h"
#include "utilities/typeMap.h"

namespace olb {

namespace utilities {

template<typename KEYS, typename VALUE>
class FixedTypeIndexedMap {
private:
  std::array<VALUE,KEYS::size> _index;

public:
  template<typename TYPE>
  static constexpr bool provides() {
    return KEYS::template contains<TYPE>();
  }

  template<typename TYPE>
  const VALUE& get() const {
    return _index[KEYS::template index<TYPE>()];
  }
  template<typename TYPE>
  VALUE& get() {
    return _index[KEYS::template index<TYPE>()];
  }

  template<typename TYPE>
  void set(VALUE value) {
    _index[KEYS::template index<TYPE>()] = value;
  }

};

/// (Time) efficient mapping between TYPEs and VALUEs
/**
 * CONTEXT should be set to the containing type / some one-of identifier
 * to restrict the index space to the types that are actually used in the
 * specific context.
 **/
template<typename VALUE, typename CONTEXT=void>
class TypeIndexedMap {
private:
  std::vector<std::optional<VALUE>> _index;

  static std::size_t next_index();
  template <typename TYPE>
  static std::size_t get_index();

public:
  TypeIndexedMap() {
    _index.resize(1);
  }

  template<typename TYPE>
  bool provides() const;

  template<typename TYPE>
  std::size_t index() const;

  template<typename TYPE>
  const VALUE& get() const;
  template<typename TYPE>
  VALUE& get();

  template<typename TYPE>
  void set(VALUE value);

};

template <typename VALUE, typename CONTEXT>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::next_index()
{
  std::size_t curr;
  #ifdef PARALLEL_MODE_OMP
  #pragma omp critical
  #endif
  {
    static std::size_t index = 0;
    curr = index++;
  }
  return curr;
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::get_index()
{
  static std::size_t index = next_index();
  return index;
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
bool TypeIndexedMap<VALUE,CONTEXT>::provides() const
{
  const std::size_t index = get_index<TYPE>();
  if (index >= _index.size()) {
    return false;
  } else {
    return _index[index].has_value();
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::index() const
{
  return get_index<TYPE>();
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
const VALUE& TypeIndexedMap<VALUE,CONTEXT>::get() const
{
  const std::size_t index = get_index<TYPE>();
  if (auto& value = _index[index]) {
    return value.value();
  } else {
    throw std::out_of_range("VALUE for TYPE doesn't exist as it was not set.");
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
VALUE& TypeIndexedMap<VALUE,CONTEXT>::get()
{
  const std::size_t index = get_index<TYPE>();
  if (auto& value = _index[index]) {
    return value.value();
  } else {
    throw std::out_of_range("VALUE for TYPE doesn't exist as it was not set.");
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
void TypeIndexedMap<VALUE,CONTEXT>::set(VALUE value)
{
  const std::size_t index = get_index<TYPE>();
  if (index >= _index.size()) {
    _index.resize(2*index);
  }
  _index[index] = value;
}


/// Mapping between KEYs and instances of type VALUEs
/**
  * \arg MAP Type map structured as `meta::plain_map`
  **/
template <typename MAP>
struct TypeIndexedTuple {
  /// Storage of values in MAP-order
  using tuple_t = typename MAP::values_t::template decompose_into<std::tuple>;
  using map_t = MAP;
  tuple_t tuple;

  TypeIndexedTuple() any_platform = default;
  TypeIndexedTuple(TypeIndexedTuple&) any_platform = default;
  TypeIndexedTuple(TypeIndexedTuple&& rhs) any_platform = default;

  template <typename... KEYS>
  struct generator {
    template <typename F>
    auto operator()(F f) any_platform {
      // Don't use std::make_tuple here due to weird incompatability with CUDA 11.7+
      return std::tuple<decltype(f(meta::id<KEYS>{}))...>(f(meta::id<KEYS>{})...);
    }
  };

  template <typename F>
  TypeIndexedTuple(F f) any_platform:
    tuple{typename MAP::keys_t::template decompose_into<generator>{}(f)} { }

  template <typename KEY>
  using value_t = typename MAP::template value<KEY>;

  /// Return true iff MAP contains KEY
  template <typename KEY>
  constexpr bool contains(meta::id<KEY> = meta::id<KEY>{}) any_platform {
    return MAP::keys_t::template contains<KEY>();
  }

  /// Access Ith element
  template <unsigned I>
  constexpr auto& get() any_platform {
    return std::get<I>(tuple);
  }
  /// Access Ith element
  template <unsigned I>
  constexpr const auto& get() const any_platform {
    return std::get<I>(tuple);
  }

  /// Access by KEY of MAP
  template <typename KEY>
  constexpr auto& get(meta::id<KEY> = meta::id<KEY>{}) any_platform {
    return get<(MAP::keys_t::template index<KEY>())>();
  }
  /// Access by KEY of MAP
  template <typename KEY>
  constexpr const auto& get(meta::id<KEY> = meta::id<KEY>{}) const any_platform {
    return get<(MAP::keys_t::template index<KEY>())>();
  }

  /// Set value assigned to KEY of MAP
  template <typename KEY>
  constexpr void set(typename MAP::template value<KEY> value) any_platform {
    get<KEY>() = value;
  }

  template <typename KEY>
  constexpr void set(meta::id<KEY>, typename MAP::template value<KEY> value) any_platform {
    get<KEY>() = value;
  }

  /// Calls f(name, value) for all pairs
  template <typename F>
  constexpr void for_each(F f) any_platform {
    MAP::keys_t::for_each([&](auto key) {
      f(key, get(key));
    });
  }

  /// Returns new TypeIndexedTuple with same keys but values generated by f(key)
  template <typename F>
  constexpr auto exchange_values(F f) any_platform {
    return TypeIndexedTuple<meta::plain_map<
      typename MAP::keys_t,
      typename MAP::keys_t::template map_to_callable_result<F>
    >>(f);
  }

};

template <typename MAP>
using TypeIndexedSharedPtrTuple = TypeIndexedTuple<
  typename MAP::template map_values<std::shared_ptr>
>;

}

}

#endif
