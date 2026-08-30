/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef CORE_CONCEPTS_H
#define CORE_CONCEPTS_H

#if (defined __cpp_concepts && (__cplusplus >= __cpp_concepts))

#include <concepts>
#include <type_traits>

#include "baseType.h"
#include "vector.h"
#include "fieldArrayD.h"
#include "operatorScope.h"

#include "descriptor/descriptor.h"

#include "utilities/aDiff.h"

namespace olb {

namespace concepts {

/// Foundational arithmetic type for all computations
template <typename T>
concept BaseType = requires(T x, T y) {
  // Basic arithmetic
  x + y; x - y; x * y; x / y;
  // in-place arithmetic
  x += y; x -= y; x *= y; x /= y;
  // Establish guaranteed core set of additional functions
  util::sqrt(x);
  //util::abs(x);
  //util::sin(x);
  //util::cos(x);
  //util::pow(x,y);
};

template <typename DESCRIPTOR>
concept Descriptor = requires () {
  { DESCRIPTOR::d } -> std::convertible_to<std::size_t>;
  requires std::is_base_of_v<meta::list_base, DESCRIPTOR>;
};

namespace placeholder {

using ValueType = float;
using Descriptor = descriptors::D2Q9<>;

using Descriptor2D = descriptors::D2Q9<>;
using Descriptor3D = descriptors::D3Q19<>;

}

/// Descriptor of DdQq lattice and associated constants
template <typename DESCRIPTOR>
concept LatticeDescriptor = requires () {
  requires Descriptor<DESCRIPTOR>;
  /// Count of discrete velocities is defined
  { DESCRIPTOR::q } -> std::convertible_to<std::size_t>;
  /// Population field is provided
  requires (DESCRIPTOR::template provides<descriptors::POPULATION>());
  /// Discrete velocities are defined as D-dimensional vectors
  { descriptors::c<DESCRIPTOR>(unsigned{}) } -> std::same_as<Vector<int,DESCRIPTOR::d>>;
  /// Discrete velocities are associated with weight constants
  { descriptors::t<placeholder::ValueType,DESCRIPTOR>(unsigned{}) } -> std::same_as<placeholder::ValueType>;
  /// Indices of opposing discrete velocities are defined
  { descriptors::opposite<DESCRIPTOR>(unsigned{}) } -> std::same_as<int>;
  /// Inverse speed of sound squared is defined
  { descriptors::invCs2<placeholder::ValueType,DESCRIPTOR>() } -> std::same_as<placeholder::ValueType>;
};

template <typename FIELD>
concept Field = requires() {
  // Column type is declared
  //FIELD::template column_type<placeholder::ValueType>;
  // Size is exposed
  { FIELD::template size<placeholder::Descriptor>() } -> std::convertible_to<std::size_t>;
  // Initial value is defined
  { FIELD::template getInitialValue<placeholder::ValueType,placeholder::Descriptor>() } -> std::same_as<
    Vector<typename FIELD::template value_type<placeholder::ValueType>,
           placeholder::Descriptor::template size<FIELD>()>
  >;
  // (De)Serializability is defined
  { FIELD::isSerializable() } -> std::same_as<bool>;
};

namespace placeholder {

struct Field {
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractColumn<T>;

  static constexpr bool isSerializable() {
    return true;
  }

  template <typename DESCRIPTOR>
  static constexpr std::size_t size() {
    return 3;
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,size<DESCRIPTOR>()>{};
  }
};

}

/// Basic cell exposing value-typed population references
template <typename CELL>
concept MinimalCell = requires(CELL cell) {
  // Population references are exposed via operator[]
  { cell[unsigned{}] } -> std::convertible_to<const typename CELL::value_t&>;
};

/// Cell exposing populations and associated fields
template <typename CELL>
concept Cell = requires(CELL cell) {
  // Population access
  requires MinimalCell<CELL>;
  // CELL::value_t is valid BaseType
  requires BaseType<typename CELL::value_t>;
  // CELL::descriptor_t is a DdQq lattice descriptor
  requires LatticeDescriptor<typename CELL::descriptor_t>;
  // Access to associated fields (using guaranteed POPULATION field as placeholder)
  // get by-value
  { cell.template getField<descriptors::POPULATION>() } -> std::same_as<
    FieldD<typename CELL::value_t,typename CELL::descriptor_t,descriptors::POPULATION>
  >;
  // set by-value
  //cell.template setField<descriptors::POPULATION>(
  //  FieldD<typename CELL::value_t,typename CELL::descriptor_t,descriptors::POPULATION>{});
  // get component by-value
  { cell.template getFieldComponent<descriptors::POPULATION>(unsigned{}) } -> std::convertible_to<typename CELL::value_t>;
  // get by-pointer
  { cell.template getFieldPointer<descriptors::POPULATION>()[unsigned{}] } -> std::convertible_to<typename CELL::value_t>;
};

/// Parameters tuple for use in collision, non-local operators
template <typename PARAMETERS>
concept Parameters = requires (typename PARAMETERS::template include_fields<placeholder::Field> parameters) {
  // PARAMETERS::value_t is valid BaseType
  requires BaseType<typename PARAMETERS::value_t>;
  // PARAMETERS::descriptor_t is a descriptor
  //requires LatticeDescriptor<typename PARAMETERS::descriptor_t>;

  // Access to parameter fields
  { parameters.template provides<placeholder::Field>() } -> std::same_as<bool>;
  { parameters.template get<placeholder::Field>() } -> std::convertible_to<
    FieldD<typename PARAMETERS::value_t,typename PARAMETERS::descriptor_t,placeholder::Field>
  >;
  // set by-value
  parameters.template set<placeholder::Field>(
    FieldD<typename PARAMETERS::value_t,typename PARAMETERS::descriptor_t,placeholder::Field>{});
};

namespace placeholder {

struct MinimalCell {
  using value_t = ValueType;
  using descriptor_t = Descriptor;

  ValueType f[Descriptor::q];
  ValueType& operator[](unsigned iPop) { return f[iPop]; }
};

struct Cell : public MinimalCell {
  template <typename FIELD>
  auto getField() { return FieldD<ValueType,Descriptor,FIELD>(this->f); }
  template <typename FIELD>
  auto getFieldComponent(unsigned iPop) { return this->f[iPop]; }
  template <typename FIELD>
  auto& getFieldPointer() { return this->f; }
  template <typename FIELD>
  void setField(FieldD<ValueType,Descriptor,FIELD>&&) { }
};

struct Parameters2D  {
  using value_t = ValueType;
  using descriptor_t = Descriptor2D;

  template <typename...>
  using include_fields = Parameters2D;

  template <typename FIELD>
  bool provides() const { return false; }
  template <typename FIELD>
  auto get() {
    if constexpr (descriptor_t::template size<FIELD>() == 1) {
      return typename FieldD<value_t,descriptor_t,FIELD>::value_t{};
    } else {
      return FieldD<value_t,descriptor_t,FIELD>{};
    }
  }
  template <typename FIELD>
  void set(FieldD<ValueType,Descriptor2D,FIELD>&&) { }
};

struct Parameters3D  {
  using value_t = ValueType;
  using descriptor_t = Descriptor3D;

  template <typename...>
  using include_fields = Parameters3D;

  template <typename FIELD>
  bool provides() const { return false; }
  template <typename FIELD>
  auto get() {
    if constexpr (descriptor_t::template size<FIELD>() == 1) {
      return typename FieldD<value_t,descriptor_t,FIELD>::value_t{};
    } else {
      return FieldD<value_t,descriptor_t,FIELD>{};
    }
  }
  template <typename FIELD>
  void set(FieldD<ValueType,Descriptor3D,FIELD>&&) { }
};

using Parameters = Parameters2D;

}

/// Full cell exposing populations, moments, dynamics and associated fields
template <typename CELL>
concept DynamicCell = requires(CELL cell, typename CELL::value_t value) {
  // Field access
  requires Cell<CELL>;

  // Computation of moments
  { cell.computeRho() } -> std::convertible_to<typename CELL::value_t>;
  cell.computeU(&value);
  cell.computeJ(&value);
  cell.computeRhoU(value, &value);
  cell.computeStress(&value);
  cell.computeAllMomenta(value,&value,&value);

  // Definition of moments
  cell.defineRho(value);
  cell.defineU(&value);
  cell.defineRhoU(value,&value);
  cell.defineAllMomenta(value,&value,&value);
  cell.inverseShiftRhoU(value,&value);

};

/// Cell-wise operator with scope, priority and apply template
template <typename OPERATOR>
concept CellOperator = requires(OPERATOR op) {
  // Operator scope is defined
  { OPERATOR::scope } -> std::convertible_to<OperatorScope>;
  // Operator scope is cell-wise
  requires OPERATOR::scope == OperatorScope::PerCell
        || OPERATOR::scope == OperatorScope::PerCellWithParameters;
  // Operator priority is defined
  { op.getPriority() } -> std::same_as<int>;
  // Todo: Operator provides apply method template
};

/// Block-wise operator with scope, priority and apply template
template <typename OPERATOR>
concept BlockOperator = requires(OPERATOR op) {
  // Operator scope is defined
  { OPERATOR::scope } -> std::convertible_to<OperatorScope>;
  // Operator scope is block-wise
  requires OPERATOR::scope == OperatorScope::PerBlock;
  // Operator priority is defined
  { op.getPriority() } -> std::same_as<int>;
  // Todo: Operator provides apply method template
};

}

}

#else // Fallback if concepts are not available (standard < C++20)

#error C++20 concept support is required

#endif

#endif
