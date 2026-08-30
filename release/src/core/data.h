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

#ifndef CORE_DATA_H
#define CORE_DATA_H

#include "platform/platform.h"

#include "fieldArrayD.h"
#include "fieldParametersD.h"
#include "utilities/typeIndexedContainers.h"

#include <optional>

namespace olb {

template <typename T, Platform PLATFORM> class ConcreteBlockMask;
template <typename T, typename DESCRIPTOR> class BlockLattice;
template <typename T, typename DESCRIPTOR> class BlockD;
template <typename T, typename DESCRIPTOR> class Cell;

/// TODO Make nice
struct field_type_of_field_array_tag { };

/// Describe FieldArray of a FIELD in Data
/**
 * i.e. Data::getField<Array<FIELD>> will expose FieldArrayD<FIELD>
 **/
template <typename FIELD>
struct Array final : public field_type_of_field_array_tag {
  using field_t = FIELD;

  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>;
};

/// Describe paramaters of OPERATOR in Data
/**
 * i.e. Data::getField<OperatorParameters<DYNAMICS>> will expose ConcreteParametersD<DYNAMICS::parameters>
 **/
template <typename OPERATOR>
struct OperatorParameters {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteParametersD<T,DESCRIPTOR,PLATFORM,typename OPERATOR::parameters>;
};

/// Describe mask of DYNAMICS in Data
/**
 * i.e. Data::getField<DynamicsMask<DYNAMICS>> will expose ConcreteBlockMask
 **/
template <typename DYNAMICS>
struct DynamicsMask {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteBlockMask<T,PLATFORM>;
};

namespace fields {

/// Helper to determine if FIELD contains FieldArray
template <typename FIELD>
constexpr bool isArray() {
  return std::is_base_of_v<field_type_of_field_array_tag, FIELD>;
}

}

/// Helper for conveying the ability for operations on FIELD_TYPE
/**
 * e.g. enable field allocation and communication setup without exposing the concrete
 * field type (which would necessitate templatization)
 *
 * c.f. DynamicsPromise and PostProcessorPromise
 **/
template<typename T, typename DESCRIPTOR>
class FieldTypePromise {
private:
  /// Name of stored FIELD_TYPE
  const std::string _name;

  /// Dimension of FIELD if it is a FieldArray
  std::optional<unsigned> _dim = std::nullopt;

  /// Function for allocating FIELD_TYPE without concrete knowledge of FIELD_TYPE
  std::function<void(BlockLattice<T,DESCRIPTOR>&)> _allocateInLattice;
  /// Function for allocating FIELD_TYPE without concrete knowledge of FIELD_TYPE
  std::function<void(BlockD<T,DESCRIPTOR>&)> _allocateInBlockD;

  std::function<void(Cell<Expr,DESCRIPTOR>&, std::vector<Expr>&)> _setPlaceholderExpr;
  std::function<std::vector<Expr>(Cell<Expr,DESCRIPTOR>&)> _getPlaceholderExpr;

public:
  template <typename FIELD_TYPE>
  FieldTypePromise(meta::id<FIELD_TYPE>):
    _name{fields::name<FIELD_TYPE>()},
    _allocateInLattice([](BlockLattice<T,DESCRIPTOR>& block) {
      if constexpr (concepts::LatticeDescriptor<DESCRIPTOR>) {
        block.template getData<FIELD_TYPE>();
      }
    }),
    _allocateInBlockD([](BlockD<T,DESCRIPTOR>& block) {
       block.template getData<FIELD_TYPE>();
    }),
    _setPlaceholderExpr([](Cell<Expr,DESCRIPTOR>& cell, std::vector<Expr>& expr) {
      if constexpr (concepts::LatticeDescriptor<DESCRIPTOR>) {
        if constexpr (fields::isArray<FIELD_TYPE>()) {
          if constexpr (std::is_same_v<Expr, typename FieldD<Expr,DESCRIPTOR,typename FIELD_TYPE::field_t>::value_t>) {
            FieldD<Expr,DESCRIPTOR,typename FIELD_TYPE::field_t> placeholder{};
            for (unsigned iD=0; iD < placeholder.d; ++iD) {
              placeholder[iD] = expr[iD];
            }
            cell.template setField<typename FIELD_TYPE::field_t>(placeholder);
          } else {
            throw std::runtime_error("Can not set placeholder for non-T types");
          }
        } else {
          throw std::runtime_error("Can not set placeholder for non-FieldArray types");
        }
      }
    }),
    _getPlaceholderExpr([](Cell<Expr,DESCRIPTOR>& cell) {
      std::vector<Expr> expressions;
      if constexpr (concepts::LatticeDescriptor<DESCRIPTOR>) {
        if constexpr (fields::isArray<FIELD_TYPE>()) {
          if constexpr (std::is_same_v<Expr, typename FieldD<Expr,DESCRIPTOR,typename FIELD_TYPE::field_t>::value_t>) {
            FieldD<Expr,DESCRIPTOR,typename FIELD_TYPE::field_t> placeholder{};
            placeholder = cell.template getFieldPointer<typename FIELD_TYPE::field_t>();
            for (unsigned iD=0; iD < placeholder.d; ++iD) {
              expressions.push_back(placeholder[iD]);
            }
          } else {
            throw std::runtime_error("Can not get placeholder for non-T types");
          }
        } else {
          throw std::runtime_error("Can not get placeholder for non-FieldArray types");
        }
      }
      return expressions;
    })
  {
    if constexpr (fields::isArray<FIELD_TYPE>()) {
      _dim = std::optional<unsigned>{DESCRIPTOR::template size<typename FIELD_TYPE::field_t>()};
    }
  }

  std::string name() const {
    return _name;
  }

  std::optional<unsigned> dimension() const {
    return _dim;
  }

  void ensureAvailabilityIn(BlockLattice<T,DESCRIPTOR>& block) {
    _allocateInLattice(block);
  }

  void ensureAvailabilityIn(BlockD<T,DESCRIPTOR>& block) {
    _allocateInBlockD(block);
  }

  void setPlaceholderExpression(Cell<Expr,DESCRIPTOR>& cell, std::vector<Expr>& expr) {
    _setPlaceholderExpr(cell, expr);
  }

  std::vector<Expr> getPlaceholderExpression(Cell<Expr,DESCRIPTOR>& cell) {
    return _getPlaceholderExpr(cell);
  }

  bool operator<(const FieldTypePromise<T,DESCRIPTOR>& rhs) const {
    return name() < rhs.name();
  }

};

/// Container for arbitrary data instances
/**
 * Enables runtime-allocated field storage in Data
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class AnyFieldTypeD {
private:
  const std::string _name;
  const FieldTypePromise<T,DESCRIPTOR> _promise;

  std::unique_ptr<Serializable> _data;
  std::function<void(ProcessingContext)> _setProcessingContext;

public:
  template<typename FIELD_TYPE, typename... ARGS>
  AnyFieldTypeD(meta::id<FIELD_TYPE> id, ARGS&&... args):
    _name(fields::name<FIELD_TYPE>()),
    _promise(id),
    _data(new typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>(
      std::forward<decltype(args)>(args)...)),
    _setProcessingContext([&](ProcessingContext context) {
      static_cast<typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>*>(_data.get())->setProcessingContext(context);
    })
  { }

  std::string getName() const {
    return _name;
  }

  FieldTypePromise<T,DESCRIPTOR> getPromise() const {
    return _promise;
  }

  template<typename FIELD_TYPE>
  const auto* as() const {
    return static_cast<const typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>*>(_data.get());
  }
  template<typename FIELD_TYPE>
  auto* as() {
    return static_cast<typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>*>(_data.get());
  }

  Serializable* asSerializable() {
    return _data.get();
  }

  template<typename TYPE>
  std::optional<TYPE*> tryAs() {
    if (TYPE* ptr = dynamic_cast<TYPE*>(_data.get())) {
      return std::optional<TYPE*>(ptr);
    } else {
      return std::nullopt;
    }
  }

  void setProcessingContext(ProcessingContext context) {
    _setProcessingContext(context);
  }

};

/// Efficient indexing of dynamically allocated data fields
/**
 * May be specialized to manage platform specific access structures
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class FieldTypeRegistry {
private:
  utilities::TypeIndexedMap<AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>*,FieldTypeRegistry> _index;

public:
  /// Returns true iff FIELD_TYPE is registered
  /**
   * FIELD_TYPE is registered if track<FIELD_TYPE> was called previously
   **/
  template <typename FIELD_TYPE>
  bool provides() const {
    return _index.template provides<FIELD_TYPE>();
  }

  /// Return read-only pointer to FIELD_TYPE data
  /**
   * FIELD_TYPE must be registered
   **/
  template <typename FIELD_TYPE>
  const AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>* get() const {
    return _index.template get<FIELD_TYPE>();
  }
  /// Return pointer to FIELD_TYPE data
  /**
   * FIELD_TYPE must be registered
   **/
  template <typename FIELD_TYPE>
  AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  /// Track newly allocated FIELD_TYPE
  template <typename FIELD_TYPE>
  void track(AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>* fieldType) {
    _index.template set<FIELD_TYPE>(fieldType);
  }

};


template<typename T, typename DESCRIPTOR> class Data;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteData;

/// Curried ConcreteData template for use in callUsingConcretePlatform
template<typename T, typename DESCRIPTOR>
struct ConcretizableData {

using value_t = T;

using base_t = Data<T,DESCRIPTOR>;

template <Platform PLATFORM>
using type = ConcreteData<T,DESCRIPTOR,PLATFORM>;

};

/// Platform-indepentent interface to ConcreteData
template<typename T, typename DESCRIPTOR>
class Data : public Serializable {
protected:
  const Platform _platform;

  /// Calls f on concrete implementation of the Data interface
  template<typename F>
  auto dispatch(F&& f) const {
    return callUsingConcretePlatform<ConcretizableData<T,DESCRIPTOR>>(
      _platform,
      this,
      f);
  }
  template<typename F>
  auto dispatch(F&& f) {
    return callUsingConcretePlatform<ConcretizableData<T,DESCRIPTOR>>(
      _platform,
      this,
      f);
  }

public:
  Data(Platform platform):
    _platform{platform}
  { }
  virtual ~Data() = default;

  /// Set processing context
  /**
   * This is currently used to trigger data transfers between host
   * and GPU data for Platform::GPU_CUDA.
   **/
  virtual void setProcessingContext(ProcessingContext) = 0;

  /// Return platform of concrete implementation
  Platform getPlatform() const {
    return _platform;
  }

  /// Return whether FIELD_TYPE is available / has been allocated
  template<typename FIELD_TYPE>
  bool provides() {
    return dispatch([&](const auto* data) -> bool {
      return data->template provides<FIELD_TYPE>();
    });
  }

  /// Allocate and return FIELD_TYPE data
  /**
   * Must only be called once for each FIELD_TYPE
   **/
  template <typename FIELD_TYPE, typename... ARGS>
  auto& allocate(ARGS&&... args) {
    return *dispatch([&](auto* data) {
      return &(data->template allocate<FIELD_TYPE>(std::forward<ARGS&&>(args)...).asAbstract());
    });
  }

  /// Return reference to data of FIELD_TYPE
  template <typename FIELD_TYPE>
  const auto& get() const {
    return *dispatch([&](const auto* data) {
      return &(data->template get<FIELD_TYPE>().asAbstract());
    });
  }
  /// Return reference to data of FIELD_TYPE
  template <typename FIELD_TYPE>
  auto& get() {
    return *dispatch([&](auto* data) {
      return &(data->template get<FIELD_TYPE>().asAbstract());
    });
  }

  virtual void resize(std::size_t newSize) = 0;

  template <typename OPERATOR>
  void apply() {
    dispatch([&](auto* data) {
      data->template apply<OPERATOR>();
    });
  }

  /// Sets FIELD_TYPE serialization state to active
  /**
   * By default Data doesn't (de)serialize any fields but
   * ConcreteBlockLattice sets all descriptor-declared
   * field arrays to active.
   **/
  template <typename FIELD_TYPE>
  void setSerialization(bool active) {
    dispatch([&](auto* data) {
      data->template setSerialization<FIELD_TYPE>(active);
    });
  }

};

/// Storage of any FIELD_TYPE data on PLATFORM
/**
 * Manages access and serialization of dynamically allocated FIELD_TYPE data.
 * This is the core storage class used for ConcreteBlockLattice on any platform
 * supported by OpenLB.
 *
 * See Data<T,DESCRIPTOR> for the platform-agnostic interface.
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteData final : public Data<T,DESCRIPTOR> {
private:
  /// Map for dynamic-dispatch of field accesses
  std::map<std::type_index, std::unique_ptr<AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>>> _map;
  /// Efficient static-dispatch access to fields
  FieldTypeRegistry<T,DESCRIPTOR,PLATFORM> _registry;
  /// Map for ordering serialization of fields
  std::map<std::string, Serializable*> _serializationNames;

public:
  /// Construct empty
  ConcreteData():
    Data<T,DESCRIPTOR>(PLATFORM)
  { }

  /// Construct Array<FIELD> of size for each FIELD in DESCRIPTOR and mark for serialization
  ConcreteData(std::size_t size):
    Data<T,DESCRIPTOR>(PLATFORM)
  {
    DESCRIPTOR::fields_t::for_each([&](auto id) {
      using field = typename decltype(id)::type;
      using field_type = Array<field>;
      allocate<field_type>(size);
      setSerialization<field_type>(true);
    });
  }

  /// Returns true iff FIELD_TYPE is allocated and can be accessed
  template <typename FIELD_TYPE>
  bool provides() const {
    return _registry.template provides<FIELD_TYPE>();
  }

  /// Allocate and return FIELD_TYPE data
  template <typename FIELD_TYPE, typename... ARGS>
  auto& allocate(ARGS&&... args) {
    _map[typeid(FIELD_TYPE)] = std::make_unique<AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>>(
      meta::id<FIELD_TYPE>(), std::forward<decltype(args)>(args)...);

    AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>* newField = _map[typeid(FIELD_TYPE)].get();
    _registry.template track<FIELD_TYPE>(newField);
    return *newField->template as<FIELD_TYPE>();
  }

  /// Return reference to data of FIELD_TYPE
  template <typename FIELD_TYPE>
  const auto& get() const {
    return *_registry.template get<FIELD_TYPE>()->template as<FIELD_TYPE>();
  }
  /// Return reference to data of FIELD_TYPE
  template <typename FIELD_TYPE>
  auto& get() {
    return *_registry.template get<FIELD_TYPE>()->template as<FIELD_TYPE>();
  }

  /// Call f for each managed AnyFieldType
  template <typename F>
  void forEach(F f) {
    for (auto& [_, field] : _map) {
      f(*field);
    }
  }

  /// Call f for each managed AnyFieldType of TYPE
  /**
   * Used to e.g. iterate over all field arrays
   **/
  template <typename TYPE, typename F>
  void forEachCastable(F f) {
    forEach([&f](AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>& any) {
      if (auto casted = any.template tryAs<TYPE>()) {
        f(*casted);
      }
    });
  }

  void resize(std::size_t newSize) override {
    forEachCastable<ColumnVectorBase>([&](auto columnVector) {
      columnVector->resize(newSize);
    });
  }

  /// Set processing context of all managed fields
  void setProcessingContext(ProcessingContext context) override {
    forEach([context](AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>& any) {
      any.setProcessingContext(context);
    });
  }

  /// Expose FieldTypeRegistry for device-support
  auto& getRegistry() {
    return _registry;
  }

  /// Sets FIELD_TYPE serialization state to active
  template <typename FIELD_TYPE>
  void setSerialization(bool active) {
    if (active) {
      _serializationNames[meta::name<FIELD_TYPE>()] = _registry.template get<FIELD_TYPE>()->asSerializable();
    } else {
      _serializationNames.erase(meta::name<FIELD_TYPE>());
    }
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

/// Constructs Data accessible on all locally used platforms
template <typename T, typename DESCRIPTOR, typename... ARGS>
std::unique_ptr<Data<T,DESCRIPTOR>> makeSharedData(
  LoadBalancer<T>& loadBalancer, ARGS&&... args)
{
  #if defined(PLATFORM_GPU_CUDA)
  if (loadBalancer.isLocal(Platform::GPU_CUDA)) {
    return std::unique_ptr<Data<T,DESCRIPTOR>>(
      new ConcreteData<T,DESCRIPTOR,Platform::GPU_CUDA>(std::forward<ARGS&&>(args)...));
  } else {
    return std::unique_ptr<Data<T,DESCRIPTOR>>(
      new ConcreteData<T,DESCRIPTOR,Platform::CPU_SISD>(std::forward<ARGS&&>(args)...));
  }
  #else
  return std::unique_ptr<Data<T,DESCRIPTOR>>(
    new ConcreteData<T,DESCRIPTOR,Platform::CPU_SISD>(std::forward<ARGS&&>(args)...));
  #endif
}

}

#endif
