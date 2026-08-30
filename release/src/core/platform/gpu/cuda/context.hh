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

#ifndef GPU_CUDA_CONTEXT_HH
#define GPU_CUDA_CONTEXT_HH

#include "context.h"

#include "registry.hh"

namespace olb {

namespace gpu {

namespace cuda {

namespace concepts {

template <typename CONTAINER>
concept ConcreteFieldsContainer = requires (CONTAINER& container) {
  requires olb::concepts::BaseType<typename CONTAINER::value_t>;
  requires olb::concepts::Descriptor<typename CONTAINER::descriptor_t>;
  requires (CONTAINER::platform == Platform::GPU_CUDA);
};

}

/// Return pointers to on-device data for FIELD on lattice
template <concepts::ConcreteFieldsContainer CONTAINER, typename FIELD>
static auto getDeviceFieldPointer(CONTAINER& lattice) {
  using T = typename CONTAINER::value_t;
  using DESCRIPTOR = typename CONTAINER::descriptor_t;
  auto& fieldArray = lattice.template getField<FIELD>();
  return Vector<typename FIELD::template value_type<T>*,DESCRIPTOR::template size<FIELD>()>([&](unsigned iD) {
    return fieldArray[iD].deviceData();
  });
}

/// Return tuple of pointers to on-device data for multiple FIELDS on lattice
template <concepts::ConcreteFieldsContainer CONTAINER, typename... FIELDS>
static auto getDeviceFieldPointers(typename meta::list<FIELDS...>,
                                   CONTAINER& lattice) {
  return std::make_tuple(getDeviceFieldPointer<CONTAINER,FIELDS>(lattice)...);
}

/// Structure for passing pointers to on-device data into CUDA kernels
template <typename T, typename DESCRIPTOR>
class DeviceContext {
private:
  const std::size_t _nCells;

  /// Pointers to descriptor-declared field data
  decltype(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(),
                                  std::declval<ConcreteBlockD<T,DESCRIPTOR,Platform::GPU_CUDA>&>()))
  _staticFieldsD;
  /// Pointer to type-erased field data
  /**
   * Maintained by FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>
   **/
  void*** _customFieldsD;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  template <concepts::ConcreteFieldsContainer CONTAINER>
  DeviceContext(CONTAINER& lattice) __host__:
    _nCells(lattice.getNcells()),
    _staticFieldsD(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(), lattice)),
    _customFieldsD(lattice.getDataRegistry().deviceData())
  {
    static_assert(std::is_same_v<typename CONTAINER::value_t, T>,
                  "CONTAINER must share value_t");
    static_assert(std::is_same_v<typename CONTAINER::descriptor_t, DESCRIPTOR>,
                  "CONTAINER must share descriptor_t");
  }

  DeviceContext(ConcreteBlockD<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) __host__:
    _nCells(lattice.getNcells()),
    _staticFieldsD(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(), lattice)),
    _customFieldsD(lattice.getDataRegistry().deviceData())
  { }

  std::size_t getNcells() const any_platform {
    return _nCells;
  }

  template <typename FIELD>
  typename FIELD::template value_type<T>** getField() __device__ {
    if constexpr (DESCRIPTOR::template provides<FIELD>()) {
      return std::get<(DESCRIPTOR::fields_t::template index<FIELD>())>(_staticFieldsD).data();
    } else {
      return reinterpret_cast<typename FIELD::template value_type<T>**>(
        _customFieldsD[field_type_index<FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>, Array<FIELD>>]);
    }
  }

};

/// Pointer to row of a D-dimensional field
template <typename T, typename DESCRIPTOR, typename FIELD>
class FieldPtr : public ScalarVector<typename FIELD::template value_type<T>,
                                     DESCRIPTOR::template size<FIELD>(),
                                     FieldPtr<T,DESCRIPTOR,FIELD>> {
private:
  DeviceContext<T,DESCRIPTOR>& _data;
  std::size_t _index;

  friend typename ScalarVector<typename FIELD::template value_type<T>,
                               DESCRIPTOR::template size<FIELD>(),
                               FieldPtr>::type;

protected:

public:
  FieldPtr(DeviceContext<T,DESCRIPTOR>& data, std::size_t index) __device__:
    _data(data),
    _index(index) { }

  FieldPtr(FieldPtr<T,DESCRIPTOR,FIELD>&& rhs) __device__:
    _data(rhs._data),
    _index(rhs._index) { }

  const typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) const __device__
  {
    return _data.template getField<FIELD>()[iDim] + _index;
  }

  typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) __device__
  {
    return _data.template getField<FIELD>()[iDim] + _index;
  }

  template <typename U, typename IMPL>
  FieldPtr<T,DESCRIPTOR,FIELD>& operator=(const GenericVector<U,DESCRIPTOR::template size<FIELD>(),IMPL>& rhs) __device__
  {
    for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

};

/// Device-side implementation of the data-only Cell concept for collision steps
template <typename T, typename DESCRIPTOR>
class DataOnlyCell {
protected:
  DeviceContext<T,DESCRIPTOR>& _data;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  DataOnlyCell(DeviceContext<T,DESCRIPTOR>& data, CellID iCell) __device__:
    _data(data),
    _iCell(iCell)
  { }

  void setCellId(std::size_t iCell) __device__ {
    _iCell = iCell;
  }

  template <typename FIELD>
  typename FIELD::template value_type<T> getFieldComponent(unsigned iD) __device__ {
    return _data.template getField<FIELD>()[iD][_iCell];
  }

  value_t& operator[](int iPop) __device__ {
    return _data.template getField<descriptors::POPULATION>()[iPop][_iCell];
  }

  template <typename FIELD>
  auto getField() const __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    if constexpr (descriptor_t::template size<FIELD>() == 1) {
      return fieldArray[0][_iCell];
    } else {
      return FieldD<value_t,descriptor_t,FIELD>([&](unsigned iD) {
        return fieldArray[iD][_iCell];
      });
    }
  }

  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer() __device__ {
    return FieldPtr<T,DESCRIPTOR,FIELD>(_data, _iCell);
  }

  template <typename FIELD>
  auto getFieldComponent(unsigned iD) const __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    return fieldArray[iD][_iCell];
  }

  template <typename FIELD>
  void setField(FieldD<value_t,descriptor_t,FIELD>&& value) __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

  template <typename FIELD>
  void setField(const FieldD<value_t,descriptor_t,FIELD>& value) __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

};


template <typename T, typename DESCRIPTOR>
class Cell;

/// Device-side view of a block lattice
/**
 * Used for non-local operators such as post processors
 **/
template <typename T, typename DESCRIPTOR>
class DeviceBlockLattice final : public DeviceContext<T,DESCRIPTOR> {
private:
  LatticeR<DESCRIPTOR::d> _core;
  LatticeR<DESCRIPTOR::d> _projection;
  int _padding;

public:
  DeviceBlockLattice(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) __host__:
    DeviceContext<T,DESCRIPTOR>(lattice),
    _core(lattice.getExtent()),
    _padding{lattice.getPadding()}
  {
    auto size = lattice.getExtent() + 2*_padding;
    if constexpr (DESCRIPTOR::d == 3) {
      _projection = {size[1]*size[2], size[2], 1};
    } else {
      _projection = {size[1], 1};
    }
  }

  CellID getCellId(LatticeR<DESCRIPTOR::d> loc) const __device__ {
    return (loc + _padding) * _projection;
  }

  CellDistance getNeighborDistance(LatticeR<DESCRIPTOR::d> dir) const __device__ {
    return dir * _projection;
  }

  Cell<T,DESCRIPTOR> get(CellID iCell) __device__ {
    return Cell<T,DESCRIPTOR>(*this, iCell);
  }

  Cell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d> loc) __device__ {
    return get(getCellId(loc));
  }

  LatticeR<DESCRIPTOR::d> getLatticeR(CellID iCell) const __device__ {
    if constexpr (DESCRIPTOR::d == 3) {
      const std::int32_t iX = iCell / _projection[0];
      const std::int32_t remainder = iCell % _projection[0];
      const std::int32_t iY = remainder / _projection[1];
      const std::int32_t iZ = remainder % _projection[1];
      return LatticeR<DESCRIPTOR::d>{iX - _padding, iY - _padding, iZ - _padding};
    } else {
      const std::int32_t iX = iCell / _projection[0];
      const std::int32_t iY = iCell % _projection[0];
      return LatticeR<DESCRIPTOR::d>{iX - _padding, iY - _padding};
    }
  };

  bool isInside(LatticeR<DESCRIPTOR::d> loc) const __device__ {
    return loc >= -_padding && loc < _core + _padding;
  }

  bool isPadding(LatticeR<DESCRIPTOR::d> loc) const __device__ {
    return isInside(loc) && !(loc >= 0 && loc < _core);
  }

  bool isPadding(CellID iCell) const __device__ {
    return isPadding(getLatticeR(iCell));
  }

};

/// Device-side implementation of the Cell concept for post processors
/**
 * Adds neighborhood and dynamically-dispatched momenta access to DataOnlyCell
 **/
template <typename T, typename DESCRIPTOR>
class Cell final : public DataOnlyCell<T,DESCRIPTOR> {
protected:
  DeviceBlockLattice<T,DESCRIPTOR>& getBlock() __device__ {
    return static_cast<DeviceBlockLattice<T,DESCRIPTOR>&>(this->_data);
  }

  DeviceBlockLattice<T,DESCRIPTOR>& getBlock() const __device__ {
    return static_cast<DeviceBlockLattice<T,DESCRIPTOR>&>(this->_data);
  }

public:
  Cell(DeviceBlockLattice<T,DESCRIPTOR>& data, CellID iCell) __device__:
    DataOnlyCell<T,DESCRIPTOR>(data, iCell)
  { }

  Cell<T,DESCRIPTOR> neighbor(LatticeR<DESCRIPTOR::d> offset) __device__ {
    return {getBlock(), CellID(this->_iCell + getBlock().getNeighborDistance(offset))};
  }

  Dynamics<T,DESCRIPTOR>& getDynamics() __device__ {
    return *this->template getField<DYNAMICS<T,DESCRIPTOR>>();
  }

  std::size_t getCellId() const __device__ {
    return this->_iCell;
  }

  void setCellId(std::size_t iCell) __device__ {
    this->_iCell = iCell;
  }

  void setLatticeR(LatticeR<DESCRIPTOR::d> latticeR) __device__ {
    this->_iCell = getBlock().getCellId(latticeR);
  }

  bool isPadding() const __device__ {
    return getBlock().isPadding(this->_iCell);
  }

  T computeRho() __device__ {
    return getDynamics().computeRho(*this);
  }
  void computeU(T* u) __device__ {
    getDynamics().computeU(*this, u);
  }
  void computeJ(T* j) __device__ {
    getDynamics().computeJ(*this, j);
  }
  void computeRhoU(T& rho, T* u) __device__ {
    getDynamics().computeRhoU(*this, rho, u);
  }
  void computeStress(T* pi) __device__ {
    T rho, u[DESCRIPTOR::d] { };
    getDynamics().computeRhoU(*this, rho, u);
    getDynamics().computeStress(*this, rho, u, pi);
  }
  void computeAllMomenta(T& rho, T* u, T* pi) __device__ {
    getDynamics().computeAllMomenta(*this, rho, u, pi);
  }

  void defineRho(T& rho) __device__ {
    return getDynamics().defineRho(*this, rho);
  }
  void defineU(T* u) __device__ {
    return getDynamics().defineU(*this, u);
  }
  void defineRhoU(T& rho, T* u) __device__ {
    return getDynamics().defineRhoU(*this, rho, u);
  }
  void defineAllMomenta(T& rho, T* u, T* pi) __device__ {
    return getDynamics().defineAllMomenta(*this, rho, u, pi);
  }
  void definePopulations(const T* f) __device__ {
    for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      this->operator[](iPop) = f[iPop];
    }
  }

  void iniEquilibrium(T rho, T* u) __device__ {
    getDynamics().iniEquilibrium(*this, rho, u);
  }
  void iniRegularized(T rho, T* u, T* pi) __device__ {
    getDynamics().iniRegularized(*this, rho, u, pi);
  }
  void inverseShiftRhoU(T& rho, T* u) __device__ {
    getDynamics().inverseShiftRhoU(*this, rho, u);
  }

};

}

}

template <typename T, typename DESCRIPTOR, typename PARAMETERS>
ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS>::ConcreteParametersD(std::size_t):
  _deviceParameters{gpu::cuda::device::malloc<ParametersD>(1)},
  parameters{}
{
  gpu::cuda::device::copyToDevice(&parameters,
                                  _deviceParameters.get(),
                                  sizeof(ParametersD));
}

template <typename T, typename DESCRIPTOR, typename PARAMETERS>
void ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS>::setProcessingContext(
  ProcessingContext context)
{
  if (context == ProcessingContext::Simulation) {
    gpu::cuda::device::copyToDevice(&parameters,
                                    _deviceParameters.get(),
                                    sizeof(ParametersD));
  }
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS>::getNblock() const
{
  return decltype(parameters)::fields_t::size;
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS>::getSerializableSize() const
{
  std::size_t size = 0;
  decltype(parameters)::fields_t::for_each([&size](auto field) {
    using field_t = typename decltype(field)::type;
    size += FieldD<T,DESCRIPTOR,field_t>{}.getSerializableSize();
  });
  return size;
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
bool* ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS>::getBlock(
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
