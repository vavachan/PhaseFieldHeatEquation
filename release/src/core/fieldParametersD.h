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

#ifndef FIELD_PARAMETERS_D_H
#define FIELD_PARAMETERS_D_H

#include "fieldArrayD.h"
#include "utilities/typeMap.h"

#include "core/concepts.h"

#include <stdexcept>

namespace olb {

/// Storage of a single FIELD-valued parameter
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, concepts::Field FIELD>
class ParameterD {
private:
  std::conditional_t<
    DESCRIPTOR::template size<FIELD>() == 1,
    typename FIELD::template value_type<T>,
    FieldD<T,DESCRIPTOR,FIELD>
  > _data;

public:
  using data_t = decltype(_data);

  ParameterD() = default;

  template <typename V>
  ParameterD(const ParameterD<V,DESCRIPTOR,FIELD>& rhs) any_platform:
    _data(rhs.read()) { }

  const auto& read() const any_platform {
    return _data;
  }

  auto& read() any_platform {
    return _data;
  }

  void write(FieldD<T,DESCRIPTOR,FIELD>&& value) any_platform {
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      _data = value[0];
    } else {
      _data = value;
    }
  }

  void write(const FieldD<T,DESCRIPTOR,FIELD>& value) any_platform {
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      _data = value[0];
    } else {
      _data = value;
    }
  }
};

/// Dynamic access interface for FIELD-valued parameters
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR>
struct AbstractParameters {
  virtual ~AbstractParameters() = default;

  template <concepts::Field FIELD>
  bool provides() const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return true;
    } else {
      return false;
    }
  };

  template <concepts::Field FIELD>
  auto get() const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->read();
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

  template <concepts::Field FIELD>
  auto getOrFallback(typename ParameterD<T,DESCRIPTOR,FIELD>::data_t&& fallback) const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->read();
    } else {
      return fallback;
    }
  }

  template <concepts::Field FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD>&& value) {
    if (auto concrete = dynamic_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->write(std::forward<decltype(value)>(value));
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

  template <concepts::Field FIELD>
  void set(const FieldD<T,DESCRIPTOR,FIELD>& value) {
    if (auto concrete = dynamic_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->write(value);
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

};

/// Set of FIELD-valued parameters
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, concepts::Field... FIELDS>
struct ParametersD final : public AbstractParameters<T,DESCRIPTOR>
                         , public ParameterD<T,DESCRIPTOR,FIELDS>...
{
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  using fields_t = meta::list<FIELDS...>;

  template <concepts::Field... Fs>
  using include_fields = ParametersD<T,DESCRIPTOR,FIELDS...,Fs...>;
  /// Return ParametersD containing FIELDS in addition to all entries of FIELD_LIST
  template <typename FIELD_LIST>
  using include = typename FIELD_LIST::template decompose_into<include_fields>;

  ParametersD() = default;

  template <typename V>
  ParametersD(const ParametersD<V,DESCRIPTOR,FIELDS...>& rhs) any_platform :
    ParameterD<T,DESCRIPTOR,FIELDS>(static_cast<const ParameterD<T,DESCRIPTOR,FIELDS>&>(rhs))...
  { }

  template <concepts::Field FIELD>
  bool provides() const any_platform {
    return fields_t::template contains<FIELD>();
  }

  template <concepts::Field FIELD>
  const auto& get() const any_platform {
    static_assert(fields_t::template contains<FIELD>(),
                  "FIELD must be available");
    return static_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)->read();
  }

  template <typename V>
  ParametersD<V,DESCRIPTOR,FIELDS...> copyAs() const any_platform {
    return ParametersD<V,DESCRIPTOR,FIELDS...>(*this);
  }

  template <concepts::Field FIELD>
  auto& get() any_platform {
    static_assert(fields_t::template contains<FIELD>(),
                  "FIELD must be available");
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->read();
  }

  template <concepts::Field FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD>&& value) any_platform {
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->write(
      std::forward<decltype(value)>(value));
  };

  template <concepts::Field FIELD>
  void set(const FieldD<T,DESCRIPTOR,FIELD>& value) any_platform {
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->write(value);
  };

};

/// Builds a ParametersD from field-value pairs
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, typename... MAP>
auto makeParametersD(MAP&&... args) {
  auto tuple = std::tie(args...);
  using result_t = typename ParametersD<T,DESCRIPTOR>::template include<typename meta::map<MAP...>::keys_t>;
  result_t params;
  result_t::fields_t::for_each([&](auto id) {
    using field_t = typename decltype(id)::type;
    constexpr unsigned idx = result_t::fields_t::template index<field_t>();
    // If parameter is pointer-valued…
    if constexpr (std::is_pointer_v<typename field_t::template value_type<T>>) {
      auto& fieldArray = std::get<2*idx+1>(tuple);
      // …and AbstractFieldArrayD is provided
      if constexpr (std::is_base_of_v<AbstractFieldArrayBase,
                                      typename std::remove_reference_t<decltype(fieldArray)>>) {
        using arg_field_t = typename std::remove_reference_t<decltype(fieldArray)>::field_t;
        using arg_descriptor_t = typename std::remove_reference_t<decltype(fieldArray)>::descriptor_t;
        static_assert(arg_descriptor_t::template size<arg_field_t>() == DESCRIPTOR::template size<field_t>(),
                      "Pointer field and field array must match dimension");
        callUsingConcretePlatform<ConcretizableFieldArrayD<T,arg_descriptor_t,arg_field_t>>(
          fieldArray.getPlatform(),
          &fieldArray.asColumnVectorBase(),
          [&](auto* concreteFieldArray) {
            if constexpr (std::remove_reference_t<decltype(*concreteFieldArray)>::platform == Platform::GPU_CUDA) {
              FieldD<T,DESCRIPTOR,field_t> fieldPtr{};
              for (unsigned iD=0; iD < concreteFieldArray->d; ++iD) {
                fieldPtr[iD] = concreteFieldArray->operator[](iD).deviceData();
              }
              params.template set<field_t>(fieldPtr);
            } else {
              FieldD<T,DESCRIPTOR,field_t> fieldPtr{};
              for (unsigned iD=0; iD < concreteFieldArray->d; ++iD) {
                fieldPtr[iD] = concreteFieldArray->operator[](iD).data();
              }
              params.template set<field_t>(fieldPtr);
            }
          }
        );
      } else {
        params.template set<field_t>(fieldArray);
      }
    // Default case, just try to set what is given
    } else {
      params.template set<field_t>(std::get<2*idx+1>(tuple));
    }
  });
  return params;
}


/// Deduce ParametersD of OPERATOR w.r.t. T and DESCRIPTOR
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, typename OPERATOR>
using ParametersOfOperatorD = typename ParametersD<T,DESCRIPTOR>::template include<
  typename OPERATOR::parameters
>;

/// Deduce ParametersD of DYNAMICS w.r.t. its value type and descriptor
template <typename DYNAMICS>
using ParametersOfDynamicsD = typename ParametersD<
  typename DYNAMICS::value_t,
  typename DYNAMICS::descriptor_t
>::template include<
  typename DYNAMICS::parameters
>;

/// Abstract base of ConcreteParametersD
/**
 * Used for platform-agnostic access to concrete parameter storage.
 **/
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR>
struct AbstractedConcreteParameters {
  virtual AbstractParameters<T,DESCRIPTOR>& asAbstract() = 0;

  virtual void setProcessingContext(ProcessingContext context) = 0;
};

/// Concrete storage of ParametersD in olb::Data
/**
 * AbstractParameters resp. ParametersD are not directly used
 * in order to preserve a minimal cross-device implementation
 * of this critical data structure.
 **/
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR, Platform PLATFORM, typename PARAMETERS>
struct ConcreteParametersD final : public AbstractedConcreteParameters<T,DESCRIPTOR>
                                 , public Serializable {
  typename ParametersD<T,DESCRIPTOR>::template include<PARAMETERS> parameters;

  ConcreteParametersD(std::size_t): // TODO: Implement more generic non-cellwise field allocation in Data
    parameters{}
  { }

  /// Return abstract interface to host-side parameters
  AbstractParameters<T,DESCRIPTOR>& asAbstract() override {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) override { }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


}

#endif
