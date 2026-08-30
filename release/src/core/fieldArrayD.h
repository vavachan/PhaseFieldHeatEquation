/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef FIELD_ARRAY_D_H
#define FIELD_ARRAY_D_H

#include "meta.h"
#include "vector.h"
#include "columnVector.h"
#include "serializer.h"
#include "descriptor/base.h"

#include "platform/platform.h"

#include <memory>
#include <tuple>

namespace olb {


template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename FIELD>
class FieldArrayD;

/// Base to identify whether a type is some kind of AbstractFieldArrayD without caring about the specifics
struct AbstractFieldArrayBase { };

/// Platform-agnostic interface to concrete host-side field arrays
template<typename T, typename DESCRIPTOR, typename FIELD>
class AbstractFieldArrayD : public AbstractFieldArrayBase {
private:
  virtual const typename FIELD::template column_type<T>& getAbstractColumn(unsigned iDim) const = 0;
  virtual       typename FIELD::template column_type<T>& getAbstractColumn(unsigned iDim)       = 0;

public:
  class const_ptr;
  class ptr;

  using descriptor_t = DESCRIPTOR;
  using field_t = FIELD;

  /// Returns Platform of concrete FieldArrayD
  virtual Platform getPlatform() const = 0;
  /// Cast to concrete FieldArrayD<PLATFORM>
  template <Platform PLATFORM>
  FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& asConcrete() {
    if (auto* ptr = dynamic_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>*>(this);
        ptr != nullptr) {
      return *ptr;
    } else {
      throw std::invalid_argument("FieldArrayD is not of PLATFORM");
    }
  }

  ColumnVectorBase& asColumnVectorBase() {
    if (auto* ptr = dynamic_cast<ColumnVectorBase*>(this);
        ptr != nullptr) {
      return *ptr;
    } else {
      throw std::invalid_argument("FieldArrayD is not based on ColumnVectorBase. Something is wrong.");
    }
  }

  /// Sets processing context (Platform dependent)
  virtual void setProcessingContext(ProcessingContext context) = 0;
  /// Resize to size newCount
  virtual void resize(std::size_t newCount) = 0;

  const auto& operator[](unsigned iDim) const
  {
    return getAbstractColumn(iDim);
  }

  auto& operator[](unsigned iDim)
  {
    return getAbstractColumn(iDim);
  }

  auto get(std::size_t i) const
  {
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      return operator[](0)[i];
    } else {
      return Vector<typename FIELD::template value_type<T>,
                    DESCRIPTOR::template size<FIELD>()>([this,i](unsigned iD) {
        return operator[](iD)[i];
      });
    }
  }

  void set(std::size_t i, const FieldD<T,DESCRIPTOR,FIELD>& data)
  {
    for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
      operator[](iD)[i] = data[iD];
    }
  }

  const_ptr getPointer(std::size_t i) const
  {
    return const_ptr(*this, i);
  }

  ptr getPointer(std::size_t i)
  {
    return ptr(*this, i);
  }

};

/// Read-only proxy for accessing a column vector entry
template<typename T, typename DESCRIPTOR, typename FIELD>
class AbstractFieldArrayD<T,DESCRIPTOR,FIELD>::const_ptr
  : public ScalarVector<const typename FIELD::template value_type<T>,
                        DESCRIPTOR::template size<FIELD>(),
                        const_ptr> {
private:
  const AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& _data;
  std::size_t _index;

  friend typename ScalarVector<const typename FIELD::template value_type<T>,
                               DESCRIPTOR::template size<FIELD>(),
                               const_ptr>::type;

protected:

public:
  const_ptr(const AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& data, std::size_t index):
    _data(data),
    _index(index) { }

  const_ptr(const_ptr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  const typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) const
  {
    return &_data[iDim][_index];
  }

  std::size_t getIndex() const
  {
    return _index;
  }

  void setIndex(std::size_t index)
  {
    _index = index;
  }

};

/// Proxy for accessing a column vector entry
template<typename T, typename DESCRIPTOR, typename FIELD>
class AbstractFieldArrayD<T,DESCRIPTOR,FIELD>::ptr
  : public ScalarVector<typename FIELD::template value_type<T>,
                        DESCRIPTOR::template size<FIELD>(),
                        ptr> {
private:
  AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& _data;
  std::size_t _index;

  friend typename ScalarVector<typename FIELD::template value_type<T>,
                               DESCRIPTOR::template size<FIELD>(),
                               ptr>::type;

protected:

public:
  ptr(AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& data, std::size_t index):
    _data(data),
    _index(index) { }

  ptr(ptr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  const typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) const
  {
    return &_data[iDim][_index];
  }

  typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim)
  {
    return &_data[iDim][_index];
  }

  template <typename U, typename IMPL>
  ptr& operator=(const GenericVector<U,DESCRIPTOR::template size<FIELD>(),IMPL>& rhs)
  {
    for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
      this->operator[](iD) = rhs[iD];
    }
    return *this;
  }

  std::size_t getIndex() const
  {
    return _index;
  }

  void setIndex(std::size_t index)
  {
    _index = index;
  }

};


/// SoA storage for instances of a single FIELD
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename FIELD>
class FieldArrayD final : public ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                                              DESCRIPTOR::template size<FIELD>()>
                        , public AbstractFieldArrayD<T,DESCRIPTOR,FIELD>
{
private:
  const typename FIELD::template column_type<T>& getAbstractColumn(unsigned iDim) const override
  {
    return this->operator[](iDim);
  }

  typename FIELD::template column_type<T>& getAbstractColumn(unsigned iDim) override
  {
    return this->operator[](iDim);
  }

public:
  using field_t = FIELD;
  using value_type = typename FIELD::template value_type<T>;
  using column_type = typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type;

  constexpr static Platform platform = PLATFORM;

  using ColumnVector<column_type,DESCRIPTOR::template size<FIELD>()>::operator[];

  FieldArrayD(std::size_t count):
    ColumnVector<column_type,
                 DESCRIPTOR::template size<FIELD>()>(count)
  {
    const auto initial = FIELD::template getInitialValue<T,DESCRIPTOR>();
    for (std::size_t i=0; i < count; ++i) {
      this->getRowPointer(i) = initial;
    }
  }

  const AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& asAbstract() const
  {
    return static_cast<const AbstractFieldArrayD<T,DESCRIPTOR,FIELD>&>(*this);
  }

  AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& asAbstract()
  {
    return static_cast<AbstractFieldArrayD<T,DESCRIPTOR,FIELD>&>(*this);
  }

  Platform getPlatform() const override
  {
    return PLATFORM;
  }

  void setProcessingContext(ProcessingContext context) override
  {
    for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
      this->operator[](iDim).setProcessingContext(context);
    }
  }

  void resize(std::size_t newCount) override {
    const std::size_t oldCount = this->_count;
    for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
      operator[](iDim).resize(newCount);
    }
    this->_count = newCount;
    if (oldCount < newCount) {
      const auto initial = FIELD::template getInitialValue<T,DESCRIPTOR>();
      for (std::size_t i=oldCount; i < newCount; ++i) {
        this->getRowPointer(i) = initial;
      }
    }
  }

  /// Return copy of FIELD data for cell iCell
  auto getField(std::size_t iCell) const
  {
    return this->getRow(iCell);
  }

  /// Set FIELD data at cell iCell
  void setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v)
  {
    this->setRow(iCell, v);
  }

  auto getFieldPointer(std::size_t iCell) const
  {
    return this->getRowPointer(iCell);
  }
  auto getFieldPointer(std::size_t iCell)
  {
    return this->getRowPointer(iCell);
  }

};

/// Curried FieldArrayD template for use in callUsingConcretePlatform
template <typename T, typename DESCRIPTOR, typename FIELD>
struct ConcretizableFieldArrayD {
  using value_t = T;

  using base_t = ColumnVectorBase;

  template <Platform PLATFORM>
  using type = FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>;

};

template<typename T> class LoadBalancer;

/// Constructs FieldArrayD accessible on all locally used platforms
/**
 * For data that is shared between blocks, e.g. reference porosity fields
 **/
template <typename T, typename DESCRIPTOR, typename FIELD>
std::unique_ptr<AbstractFieldArrayD<T,DESCRIPTOR,FIELD>> makeSharedAbstractFieldArrayD(
  LoadBalancer<T>& loadBalancer, std::size_t count=1)
{
  #if defined(PLATFORM_GPU_CUDA)
  if (loadBalancer.isLocal(Platform::GPU_CUDA)) {
    return std::unique_ptr<AbstractFieldArrayD<T,DESCRIPTOR,FIELD>>(
      new FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>(count));
  } else {
    return std::unique_ptr<AbstractFieldArrayD<T,DESCRIPTOR,FIELD>>(
      new FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,FIELD>(count));
  }
  #else
  return std::unique_ptr<AbstractFieldArrayD<T,DESCRIPTOR,FIELD>>(
    new FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,FIELD>(count));
  #endif
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename FIELD>
class ConcreteCommunicatable<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>> final : public Communicatable {
private:
  FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& _communicatee;

public:
  ConcreteCommunicatable(FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& communicatee):
    _communicatee{communicatee} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    return (ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                 DESCRIPTOR::template size<FIELD>()>
    >(_communicatee).size(indices));
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::uint8_t* curr = buffer;
    curr += ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<FIELD>()>
    >(_communicatee).serialize(indices, curr);
    return curr - buffer;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    const std::uint8_t* curr = buffer;
    curr += ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<FIELD>()>
    >(_communicatee).deserialize(indices, curr);
    return curr - buffer;
  }
};

/// Storage for a fixed set of FIELDS
/**
 * Actual field data is stored by individual FieldArrayD instances
 *
 * This class is not tied to block-structured resp. list-like data concepts
 * and may be used for both lattice and (resizable) particle data.
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
class MultiFieldArrayD : public Serializable {
private:
  /// Current row count of each field array
  std::size_t _count;
  /// FieldArrayD instances for FIELDS
  std::tuple<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS>...> _static;

public:
  using fields_t = meta::list<FIELDS...>;

  MultiFieldArrayD(std::size_t count=1):
    _count(count),
    // Trickery to construct each member of _static with `count`.
    // Uses the comma operator in conjunction with type dropping.
    _static((std::void_t<FIELDS>(), count)...)
  { }

  template <typename FIELD>
  const FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& get(meta::id<FIELD> field = meta::id<FIELD>()) const
  {
    static_assert(meta::contains<FIELD,FIELDS...>(), "FIELD not contained in FIELDS");
    return std::get<(fields_t::template index<FIELD>())>(_static);
  }

  template <typename FIELD>
  FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& get(meta::id<FIELD> field = meta::id<FIELD>())
  {
    static_assert(meta::contains<FIELD,FIELDS...>(), "FIELD not contained in FIELDS");
    return std::get<(fields_t::template index<FIELD>())>(_static);
  }

  /// Return copy of FIELD data for cell iCell
  template <typename FIELD>
  auto getField(std::size_t iCell) const
  {
    return get<FIELD>().getRow(iCell);
  }

  /// Set FIELD data at cell iCell
  template <typename FIELD>
  void setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v)
  {
    get<FIELD>().setRow(iCell, v);
  }

  template <typename FIELD>
  const typename FIELD::template value_type<T>& getFieldComponent(std::size_t iCell, unsigned iDim) const;
  template <typename FIELD>
  typename FIELD::template value_type<T>& getFieldComponent(std::size_t iCell, unsigned iDim);

  template <typename FIELD>
  auto getFieldPointer(std::size_t iCell) const
  {
    return get<FIELD>().getRowPointer(iCell);
  }
  template <typename FIELD>
  auto getFieldPointer(std::size_t iCell)
  {
    return get<FIELD>().getRowPointer(iCell);
  }

  /// Apply generic expression to each FIELD array
  template <typename F>
  void forFields(F f) const;
  template <typename F>
  void forFields(F f);
  /// Apply generic lambda expression to each FIELD of a cell
  template <typename F>
  void forFieldsAt(std::size_t idx, F f);

  /// Change number of rows
  /**
   * Drops the last (newCount-_count) rows when shrinking
   **/
  void resize(std::size_t newCount)
  {
    forFields([newCount](auto& fieldArray) {
      fieldArray.resize(newCount);
    });
    _count = newCount;
  }

  /// Swap contents of rows i and j
  void swap(std::size_t i, std::size_t j)
  {
    forFields([i,j](auto& fieldArray) {
      fieldArray.swap(i,j);
    });
  }

  void setProcessingContext(ProcessingContext context) {
    forFields([context](auto& fieldArray) {
      fieldArray.setProcessingContext(context);
    });
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

//template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename FIELD>
//ConcreteCommunicatable(FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&) -> ConcreteCommunicatable<
//  ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
//               DESCRIPTOR::template size<FIELD>()>
//>;


template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
class ConcreteCommunicatable<MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>> final : public Communicatable {
private:
  MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>& _communicatee;

public:
  ConcreteCommunicatable(MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>& communicatee):
    _communicatee{communicatee} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    return (ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename FIELDS::template column_type<T>,PLATFORM>::type,
                 DESCRIPTOR::template size<FIELDS>()>
    >(_communicatee.template get<FIELDS>()).size(indices) + ... + 0);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::uint8_t* curr = buffer;
    meta::list<FIELDS...>::for_each([&](auto field) {
      using FIELD = typename decltype(field)::type;
      curr += ConcreteCommunicatable<
        ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                     DESCRIPTOR::template size<FIELD>()>
      >(_communicatee.get(field)).serialize(indices, curr);
    });
    return curr - buffer;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    const std::uint8_t* curr = buffer;
    meta::list<FIELDS...>::for_each([&](auto field) {
      using FIELD = typename decltype(field)::type;
      curr += ConcreteCommunicatable<
        ColumnVector<typename ImplementationOf<typename FIELD::template column_type<T>,PLATFORM>::type,
                     DESCRIPTOR::template size<FIELD>()>
      >(_communicatee.get(field)).deserialize(indices, curr);
    });
    return curr - buffer;
  }
};


/// Storage for dynamic field groups (Prototype for ParticleSystem)
template <typename T, typename GROUPS>
class DynamicFieldGroupsD;
template <typename T, typename... GROUPS>
class DynamicFieldGroupsD<T, meta::list<GROUPS...>> {
private:
  template <typename GROUP>
  struct GroupedFieldArrayD {
    template <typename... FIELDS>
    using curried = MultiFieldArrayD<T,GROUP,Platform::CPU_SISD,FIELDS...>;
    using type = typename GROUP::template decompose_into<curried>;
  };

  std::tuple<typename GroupedFieldArrayD<GROUPS>::type...> _data;

  std::size_t _count;

public:
  DynamicFieldGroupsD(std::size_t count):
    _data((std::void_t<GROUPS>(), count)...), _count(count) { }

  template <typename GROUP>
  auto& get() {
    if constexpr (meta::contains<GROUP,GROUPS...>()) {
      return std::get<descriptors::getIndexInFieldList<GROUP,GROUPS...>()>(_data);
    } else {
      throw std::invalid_argument("This DynamicFieldGroupsD does not provide GROUP.");
    }
  }

  std::size_t count(){ return _count; }

  void resize(std::size_t newCount) {
    (get<GROUPS>().resize(newCount), ...);
    _count = newCount;
  }

  void swap(std::size_t i, std::size_t j)
  {
    (get<GROUPS>().swap(i,j), ...);
  }

  std::size_t constexpr getSerializableSize()
  {
    return (get<GROUPS>().getSerializableSize() + ... + 0);
  }

};




//Communicatable for GroupedFieldArrays
template<typename DATA, typename... GROUPS>
class GroupedDataCommunicatable{
private:
  DATA& _data;

public:
  GroupedDataCommunicatable( DATA& data )
    : _data(data) {}

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    std::size_t size = 0;
    meta::list<GROUPS...>::for_each([&](auto group) {
      using GROUP = typename decltype(group)::type;
      // Workaround for Clang argument deduction bug
      auto communicatable = ConcreteCommunicatable(_data.template get<GROUP>());
      size += communicatable.size(indices);
    });
    return size;
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices, std::uint8_t* buffer) const
  {
    std::uint8_t* curr = buffer;
    meta::list<GROUPS...>::for_each([&](auto group) {
      using GROUP = typename decltype(group)::type;
      // Workaround for Clang argument deduction bug
      auto communicatable = ConcreteCommunicatable(_data.template get<GROUP>());
      curr += communicatable.serialize(indices, curr);
    });
    return curr - buffer;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices, const std::uint8_t* buffer)
  {
    const std::uint8_t* curr = buffer;
    meta::list<GROUPS...>::for_each([&](auto group) {
      using GROUP = typename decltype(group)::type;
      // Workaround for Clang argument deduction bug
      auto communicatable = ConcreteCommunicatable(_data.template get<GROUP>());
      curr += communicatable.deserialize(indices, curr);
    });
    return curr - buffer;
  }
};

/// Declare GroupedDataCommunicatable containing each GROUP in DESCRIPTOR::fields_t
template <typename DATA, typename DESCRIPTOR>
struct GroupedDataCommunicatableHelper {
  template<typename... GROUPS>
  using CurriedFieldGroupsCommunicatable = GroupedDataCommunicatable<DATA,GROUPS...>;

  using type = typename DESCRIPTOR::template decompose_into<CurriedFieldGroupsCommunicatable>;
};

/// Declare MultiFieldArrayD containing each field in DESCRIPTOR::fields_t
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct MultiFieldArrayForDescriptorHelper {
  template<typename... FIELDS>
  using CurriedMultiFieldArrayD = MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>;

  using type = typename DESCRIPTOR::fields_t::template decompose_into<CurriedMultiFieldArrayD>;
};

/// MultiFieldArrayD containing each field in DESCRIPTOR::fields_t
template <typename T, typename DESCRIPTOR, Platform PLATFORM=Platform::CPU_SISD>
using MultiFieldArrayForDescriptorD = typename MultiFieldArrayForDescriptorHelper<T,DESCRIPTOR,PLATFORM>::type;


}

#endif
