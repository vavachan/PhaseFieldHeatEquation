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

#ifndef GPU_CUDA_DEVICE_HH
#define GPU_CUDA_DEVICE_HH

#include "device.h"

#include <stdexcept>

#include <cuda.h>

namespace olb {

namespace gpu {

namespace cuda {

namespace device {

int getCount() {
  int devices{};
  cudaGetDeviceCount(&devices);
  return devices;
}

void check() {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(error));
  }
}

/// Check CUDA driver errors
void check(CUresult result) {
  if (result != CUDA_SUCCESS) {
    const char* description{};
    cuGetErrorString(result, &description);
    throw std::runtime_error(std::string(description));
  }
}

void synchronize() {
  if (getCount() > 0) {
    cudaDeviceSynchronize();
    check();
  }
}

int get() {
  int device{};
  cudaGetDevice(&device);
  check();
  return device;
}

void copyToHost(void* src, void* dst, std::size_t count) {
  cudaMemcpy(dst, src, count, cudaMemcpyDeviceToHost);
  check();
}

void copyToDevice(void* src, void* dst, std::size_t count) {
  cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice);
  check();
}

template <typename T>
T* malloc(std::size_t size) {
  T* ptr{};
  cudaMalloc(&ptr, size*sizeof(T));
  check();
  return ptr;
}

std::size_t getDevicePageSize() {
  std::size_t granularity = 0;
  CUmemAllocationProp prop = {};
  prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
  prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
  prop.location.id = device::get();
  check(cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM));
  return granularity;
}


template <typename T>
unique_ptr<T>::~unique_ptr() {
  if (_ptr != nullptr) {
    cudaFree(_ptr);
    check();
  }
}

template <typename T>
unique_ptr<T>& unique_ptr<T>::operator=(unique_ptr<T>&& rhs) {
  if (_ptr != nullptr) {
    cudaFree(_ptr);
    check();
  }
  _ptr = rhs.release();
  return *this;
}


Stream::Stream(unsigned int flags) {
  cudaStreamCreateWithFlags(&_stream, flags);
  check();
}

Stream::~Stream() {
  cudaStreamDestroy(_stream);
  check();
}

void Stream::synchronize() {
  cudaStreamSynchronize(_stream);
  check();
}

void asyncCopyToHost(Stream& stream, void* src, void* dst, std::size_t count) {
  cudaMemcpyAsync(dst, src, count, cudaMemcpyDeviceToHost, stream.get());
}

void asyncCopyToDevice(Stream& stream, void* src, void* dst, std::size_t count) {
  cudaMemcpyAsync(dst, src, count, cudaMemcpyHostToDevice, stream.get());
}

}

}

}

}

#endif
