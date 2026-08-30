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

#ifndef GPU_CUDA_REGISTRY_HH
#define GPU_CUDA_REGISTRY_HH

#include "registry.h"

namespace olb {

namespace gpu {

namespace cuda {

/// Mapping of TYPE in CONTEXT to runtime-fixed index
/**
 * Used for dynamic access to field arrays
 **/
template <typename CONTEXT, typename TYPE>
__constant__ std::size_t field_type_index;

/// Host-side version of gpu::cuda::AnyDeviceFieldArrayD
struct FieldArrayPointer {
  gpu::cuda::device::unique_ptr<void*> data;
  const unsigned column_count;
  const unsigned element_size;
};

}

}

template <typename T, typename DESCRIPTOR>
struct FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::Data {
  Data():
    indexOnHost(2, nullptr)
  { }

  std::map<std::type_index,gpu::cuda::FieldArrayPointer> fieldArrayPointers;

  std::vector<void**>                   indexOnHost;
  gpu::cuda::device::unique_ptr<void**> indexOnDevice;

};

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::FieldTypeRegistry():
  _data(new Data()),
  _modified{true} { }

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::~FieldTypeRegistry()
{ }

template <typename T, typename DESCRIPTOR>
template <typename FIELD_TYPE>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::track(
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* fieldType)
{
  _index.template set<FIELD_TYPE>(fieldType);

  // Copy host-side field index to device-side field index
  std::size_t index = _index.template index<FIELD_TYPE>();
  cudaMemcpyToSymbol(gpu::cuda::field_type_index<FieldTypeRegistry,FIELD_TYPE>, &index, sizeof(std::size_t));
  gpu::cuda::device::check();
  if (index >= _data->indexOnHost.size()) {
    _data->indexOnHost.resize(2*index);
  }

  // Update device-side pointers to Array<FIELD> field types
  using ConcreteFieldType = typename FIELD_TYPE::template type<T,DESCRIPTOR,Platform::GPU_CUDA>;
  if constexpr (std::is_base_of_v<ColumnVectorBase, ConcreteFieldType>) {
    using field_t = typename ConcreteFieldType::field_t;
    auto& fieldArray = *fieldType->template as<FIELD_TYPE>();

    std::array<void*,DESCRIPTOR::template size<field_t>()> componentPointers;
    for (unsigned iD=0; iD < componentPointers.size(); ++iD) {
      componentPointers[iD] = fieldArray[iD].deviceData();
    }

    auto componentPointersOnDevice = gpu::cuda::device::malloc<void*>(DESCRIPTOR::template size<field_t>());
    gpu::cuda::device::copyToDevice(componentPointers.data(),
                                    componentPointersOnDevice,
                                    componentPointers.size()*sizeof(void*));

    _data->indexOnHost[index] = componentPointersOnDevice;
    _data->fieldArrayPointers.emplace(typeid(field_t), gpu::cuda::FieldArrayPointer{
      gpu::cuda::device::unique_ptr<void*>(componentPointersOnDevice),
      componentPointers.size(),
      sizeof(typename field_t::template value_type<T>)
    });
  }

  _modified = true;
}

template<typename T, typename DESCRIPTOR>
void*** FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceData()
{
  if (_modified) {
    _data->indexOnDevice = gpu::cuda::device::malloc<void**>(_data->indexOnHost.size());
    gpu::cuda::device::copyToDevice(_data->indexOnHost.data(),
                                    _data->indexOnDevice.get(),
                                    _data->indexOnHost.size()*sizeof(void**));
    _modified = false;
  }
  return _data->indexOnDevice.get();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::refreshDeviceFieldArray(
  FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>& fieldArray)
{
  std::array<void*,DESCRIPTOR::template size<FIELD>()> componentPointers;
  for (unsigned iD=0; iD < componentPointers.size(); ++iD) {
    componentPointers[iD] = fieldArray[iD].deviceData();
  }
  auto& ptr = _data->fieldArrayPointers.at(typeid(FIELD));
  gpu::cuda::device::copyToDevice(componentPointers.data(),
                                  ptr.data.get(),
                                  componentPointers.size()*sizeof(void*));
}

template <typename T, typename DESCRIPTOR>
gpu::cuda::AnyDeviceFieldArrayD
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceFieldArray(
  std::type_index field)
{
  auto& ptr = _data->fieldArrayPointers.at(field);
  return gpu::cuda::AnyDeviceFieldArrayD{
    ptr.data.get(),
    ptr.column_count,
    ptr.element_size
  };
}

template <typename T, typename DESCRIPTOR>
std::vector<gpu::cuda::AnyDeviceFieldArrayD>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceFieldArrays(
  const std::vector<std::type_index>& fields)
{
  std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFields;
  deviceFields.reserve(fields.size());
  for (std::type_index field : fields) {
    deviceFields.emplace_back(deviceFieldArray(field));
  }
  return deviceFields;
}

}

#endif
