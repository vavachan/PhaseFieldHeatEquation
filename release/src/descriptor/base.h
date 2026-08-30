/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender, Mathias J. Krause
 *                2021 Adrian Kummerlaender, Nicolas Hafen, Mathias J. Krause
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

#ifndef DESCRIPTOR_BASE_H
#define DESCRIPTOR_BASE_H

#include "tag.h"
#include "fields.h"

#include "core/meta.h"

namespace olb {

namespace descriptors {

template <typename FIELD>
using is_data_field = std::integral_constant<bool, !is_tag_field<FIELD>::value>;

/// \defgroup descriptor
//@{

/// Tuple of abstract field declarations
template <typename... FIELDS>
struct FIELD_TUPLE : public meta::list<FIELDS...> {
  /// Deleted constructor to enforce pure usage as type and prevent implicit narrowing conversions
  FIELD_TUPLE() = delete;

  using fields_t = meta::list<FIELDS...>;

  template <template<typename...> class COLLECTION>
  using decompose_into = COLLECTION<FIELDS...>;

  /// Number of fields
  static constexpr std::size_t field_count = sizeof...(FIELDS);

  /// Returns FIELD_TUPLE with subset of FIELDS meeting COND
  template <template<typename> class COND>
  using filter = typename meta::filter_t<COND, FIELDS...>
    ::template decompose_into<FIELD_TUPLE>;

  /// Returns whether WANTED_FIELD is contained in FIELDS
  template <typename WANTED_FIELD>
  static constexpr bool provides() {
    return meta::contains<WANTED_FIELD,FIELDS...>();
  }

  /// Returns whether field is contained in FIELDS
  static bool provides(std::type_index field) {
    return meta::contains<FIELDS...>(field);
  }

};

/// Base descriptor of a d-dimensional system
template <unsigned D, typename... FIELDS>
struct SPATIAL_DESCRIPTOR : public FIELD_TUPLE<FIELDS...> {
  /// Number of dimensions
  static constexpr int d = D;
  /// Number of velocities
  static constexpr int q = 0;

  /// Returns size of given FIELD
  template <typename FIELD>
  static constexpr std::size_t size() {
    return FIELD::template size<SPATIAL_DESCRIPTOR>();
  }
};

struct LATTICE_DESCRIPTOR_BASE { };

/// Base descriptor of a D-dimensional lattice with Q directions and a list of additional fields
template <unsigned D, unsigned Q, typename... FIELDS>
struct LATTICE_DESCRIPTOR : public LATTICE_DESCRIPTOR_BASE
                          , public FIELD_TUPLE<FIELDS...> {
  /// Number of dimensions
  static constexpr int d = D;
  /// Number of velocities
  static constexpr int q = Q;

  using all_t = meta::list<FIELDS...>;
  using fields_t = meta::filter_t<is_data_field, FIELDS...>;
  using tags_t   = meta::filter_t<is_tag_field, FIELDS...>;

  /// Returns size of given FIELD
  template <typename FIELD>
  static constexpr std::size_t size() {
    return FIELD::template size<LATTICE_DESCRIPTOR>();
  }

  /// Tag that describes the category of the descriptor
  /**
   * e.g. tag::RTLBM for RTLBM descriptors or tag::DEFAULT if no category tag-field is found.
   *
   * This is needed to enable tag-dispatching of descriptor function overloads.
   **/
  using category_tag = typename tags_t::template first_with_base_or_fallback<tag::CATEGORY, tag::DEFAULT>;

};

/// Pair of base value and descriptor type
template <typename T, typename DESCRIPTOR>
struct VALUED_DESCRIPTOR {
  using value_t = T;
  using descriptor_t = DESCRIPTOR;
};

/// Return pair of base value and descriptor type declared by TYPE
/**
 * Used to extract descriptor and associated value type from cell and lattice types.
 **/
template <typename TYPE>
using extract_valued_descriptor_t = VALUED_DESCRIPTOR<
  typename std::remove_reference_t<TYPE>::value_t,
  typename std::remove_reference_t<TYPE>::descriptor_t
>;


/// Base descriptor of a particle system
template <unsigned D, typename... FIELDS>
struct PARTICLE_DESCRIPTOR : public FIELD_TUPLE<FIELDS...> {
  /// Number of dimensions
  static constexpr int d = D;
  /// Number of velocities
  static constexpr int q = 0;

  /// Returns size of given FIELD
  template <typename FIELD>
  static constexpr std::size_t size() {
    return FIELD::template size<PARTICLE_DESCRIPTOR>();
  }

  /// Returns FIELD by provided BASE
  /**
   * returns void, if FIELD_TUPE does not provide BASE field
   **/
  template<typename BASE>
  using derivedField = std::conditional_t<
    meta::contains<
      typename meta::list_item_with_base_default_base<BASE,FIELDS...>::type,
      FIELDS...
    >(),
    typename meta::list_item_with_base_default_base<BASE,FIELDS...>::type,
    void
  >;

  /// Returns whether WANTED_FIELD (last in NESTED_FIELDS) is provided by CONCRETE_FIELD_TUPLE
  /**
   * Offers the same funcionality as it overloaded pendant in in FIELD_TUPLE, however allows beeing called
   * with multiple BASE types (as currently used in the particle framework)
   **/
  template <typename ...NESTED_FIELDS>
  static constexpr bool providesNested() {
    return meta::derived_type_in_nested<PARTICLE_DESCRIPTOR, NESTED_FIELDS...>::contains();
  }

};

//@}

}

}

#endif
