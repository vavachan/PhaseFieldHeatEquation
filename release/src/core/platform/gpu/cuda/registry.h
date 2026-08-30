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

#ifndef GPU_CUDA_REGISTRY_H
#define GPU_CUDA_REGISTRY_H

#include "utilities/typeIndexedContainers.h"

#include <map>
#include <memory>

namespace olb {

namespace gpu {

namespace cuda {

/// Type-erased pointer to FieldArrayD device data
/**
 * Used for dynamic communication buffer (de)serialization
 **/
struct AnyDeviceFieldArrayD {
  void**   data;
  unsigned column_count;
  unsigned element_size;

  std::uint8_t* operator[](unsigned iD) __device__ {
    return reinterpret_cast<std::uint8_t*>(data[iD]);
  }
};

}

}

/// Maintain on-device structure for dynamic field access
/**
 * Used to enable access to fields that are dynamically allocated
 * instead of being declared by the descriptor.
 **/
template<typename T, typename DESCRIPTOR>
class FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  /// Host-side field type index
  utilities::TypeIndexedMap<AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>*,FieldTypeRegistry> _index;

  struct Data;
  std::unique_ptr<Data> _data;

  bool _modified;

public:
  FieldTypeRegistry();
  ~FieldTypeRegistry();

  template <typename FIELD_TYPE>
  bool provides() const {
    return _index.template provides<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  const AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* get() const {
    return _index.template get<FIELD_TYPE>();
  }
  template <typename FIELD_TYPE>
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  void track(AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* fieldType);

  void*** deviceData();

  template <typename FIELD>
  void refreshDeviceFieldArray(FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>& fieldArray);

  gpu::cuda::AnyDeviceFieldArrayD deviceFieldArray(std::type_index field);
  std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFieldArrays(
    const std::vector<std::type_index>& fields);

};

}

#endif
