/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
 *
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

#ifndef GPU_CUDA_COLUMN_HH
#define GPU_CUDA_COLUMN_HH

#include "column.h"
#include "device.h"

#include <thrust/gather.h>
#include <thrust/scatter.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

#include <cstring>

namespace olb {

namespace gpu {

namespace cuda {

template<typename T>
struct Column<T>::Data {
  thrust::host_vector<T>   host;
  thrust::device_vector<T> device;

  Data(std::size_t count):
    host(count),
    device(count)
  { }

  void resize(std::size_t newCount)
  {
    host.resize(newCount);
    device.resize(newCount);
  }

};

template<typename T>
Column<T>::Column(std::size_t count):
  _count(count),
  _data(new Data(count))
{ }

template<typename T>
Column<T>::Column(Column<T>&& rhs):
  _count(rhs._count),
  _data(rhs._data.release())
{ }

template<typename T>
Column<T>::Column(const Column<T>& rhs):
  _count(rhs._count),
  _data(new Data(*rhs._data))
{ }

template<typename T>
Column<T>::~Column()
{ }

template<typename T>
void Column<T>::clear()
{
  _data->device.clear();
  _data->host.clear();
  _count = _data->host.size();
}

template<typename T>
void Column<T>::resize(std::size_t newCount)
{
  _data->resize(newCount);
  _count = newCount;
}

template<typename T>
void Column<T>::push_back(T value)
{
  _data->host.push_back(value);
  _count = _data->host.size();
}

template<typename T>
void Column<T>::deduplicate()
{
  std::sort(_data->host.begin(), _data->host.end());
  _data->host.erase(std::unique(_data->host.begin(), _data->host.end()), _data->host.end());
}

template<typename T>
const T& Column<T>::operator[](std::size_t i) const
{
  return _data->host[i];
}

template<typename T>
T& Column<T>::operator[](std::size_t i)
{
  return _data->host[i];
}

template<typename T>
std::size_t Column<T>::size() const
{
  return _count;
}

template<typename T>
const T* Column<T>::data() const
{
  return _data->host.data();
}

template<typename T>
T* Column<T>::data()
{
  return _data->host.data();
}

template<typename T>
const T* Column<T>::deviceData() const
{
  return _data->device.data().get();
}

template<typename T>
T* Column<T>::deviceData()
{
  return _data->device.data().get();
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context)
{
  if (_count != _data->device.size()) {
    _data->device.resize(_count);
  }
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(_data->device.data().get(), _data->host.data(), size()*sizeof(T));
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_data->host.data(), _data->device.data().get(), size()*sizeof(T));
    return;
  }
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context, device::Stream& stream)
{
  if (_count != _data->device.size()) {
    _data->device.resize(_count);
  }
  switch (context) {
  case ProcessingContext::Evaluation:
    device::asyncCopyToHost(stream, _data->device.data().get(), _data->host.data(), size()*sizeof(T));
    return;
  case ProcessingContext::Simulation:
    device::asyncCopyToDevice(stream, _data->host.data(), _data->device.data().get(), size()*sizeof(T));
    return;
  }
}

template<typename T>
std::size_t Column<T>::getNblock() const
{
  return 2;
}

template<typename T>
std::size_t Column<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::size_t);
}

template<typename T>
bool* Column<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  if (loadingMode && iBlock == 1) {
    resize(_count);
  }
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *data(), _count);

  return dataPtr;
}


template<typename T>
struct CyclicColumn<T>::Data {
  const std::size_t size;

  std::unique_ptr<T[]> host;

  CUmemGenericAllocationHandle handle;
  CUmemAllocationProp prop{};
  CUmemAccessDesc access{};
  CUdeviceptr ptr;

  Data(std::size_t count):
    host(new T[count] { }),
    size(count * sizeof(T))
  {
    const int device = device::get();

    prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    prop.location.id = device;
    device::check(cuMemAddressReserve(&ptr, 2 * size, 0, 0, 0));

    // per-population handle until cuMemMap accepts non-zero offset
    device::check(cuMemCreate(&handle, size, &prop, 0));
    device::check(cuMemMap(ptr,        size, 0, handle, 0));
    device::check(cuMemMap(ptr + size, size, 0, handle, 0));

    access.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    access.location.id = device;
    access.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
    device::check(cuMemSetAccess(ptr, 2 * size, &access, 1));
  }

  ~Data() {
    device::check(cuMemUnmap(ptr, 2 * size));
    device::check(cuMemRelease(handle));
    device::check(cuMemAddressFree(ptr, 2 * size));
  }

};

template<typename T>
CyclicColumn<T>::CyclicColumn(std::size_t count):
  _count(device::getPageAlignedCount<T>(count)),
  _size(_count * sizeof(T)),
  _data(new Data(_count)),
  _shift(0)
{
  _deviceBase = reinterpret_cast<T*>(_data->ptr);
  _devicePopulation = _deviceBase;
}

template<typename T>
CyclicColumn<T>::~CyclicColumn()
{ }

template<typename T>
const T& CyclicColumn<T>::operator[](std::size_t i) const
{
  return _data->host[i];
}

template<typename T>
T& CyclicColumn<T>::operator[](std::size_t i)
{
  return _data->host[i];
}

template<typename T>
void CyclicColumn<T>::setProcessingContext(ProcessingContext context)
{
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(_devicePopulation, _data->host.get(), _size);
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_data->host.get(), _devicePopulation, _size);
    return;
  }
}

template<typename T>
std::size_t CyclicColumn<T>::getNblock() const
{
  return 3;
}

template<typename T>
std::size_t CyclicColumn<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::ptrdiff_t) + sizeof(std::size_t);
}

template<typename T>
bool* CyclicColumn<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _shift);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_data->host.get(), _count);

  return dataPtr;
}


}

}


template <typename T>
std::size_t ConcreteCommunicatable<gpu::cuda::Column<T>>::serialize(ConstSpan<CellID> indices,
                                                                    std::uint8_t* buffer) const
{
  thrust::gather(thrust::device,
                 thrust::device_pointer_cast(indices.begin()),
                 thrust::device_pointer_cast(indices.end()),
                 thrust::device_pointer_cast(_column.deviceData()),
                 thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::cuda::Column<T>>::deserialize(ConstSpan<CellID> indices,
                                                                      const std::uint8_t* buffer)
{
  thrust::scatter(thrust::device,
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                  thrust::device_pointer_cast(indices.begin()),
                  thrust::device_pointer_cast(_column.deviceData()));
  return indices.size() * sizeof(T);
}


template <typename T>
std::size_t ConcreteCommunicatable<gpu::cuda::CyclicColumn<T>>::serialize(
  ConstSpan<CellID> indices,
  std::uint8_t* buffer) const
{
  thrust::gather(thrust::device,
                 thrust::device_pointer_cast(indices.begin()),
                 thrust::device_pointer_cast(indices.end()),
                 thrust::device_pointer_cast(_column.deviceData()),
                 thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::cuda::CyclicColumn<T>>::deserialize(
  ConstSpan<CellID> indices,
  const std::uint8_t* buffer)
{
  thrust::scatter(thrust::device,
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                  thrust::device_pointer_cast(indices.begin()),
                  thrust::device_pointer_cast(_column.deviceData()));
  return indices.size() * sizeof(T);
}


}

#endif
