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

#ifndef GPU_CUDA_DEVICE_H
#define GPU_CUDA_DEVICE_H

#include <cuda_runtime_api.h>

namespace olb {

namespace gpu {

namespace cuda {

/// Basic wrappers of common CUDA functions
namespace device {

/// Return number of available devices
int getCount();

/// Check errors
void check();

/// Synchronize device
void synchronize();

/// Get current device
int get();

/// Copy data from device to host
void copyToHost(void* src, void* dst, std::size_t count);
/// Copy data from host to device
void copyToDevice(void* src, void* dst, std::size_t count);

/// Allocate data on device
template <typename T>
T* malloc(std::size_t size);

/// Returns device memory page size
std::size_t getDevicePageSize();
/// Returns `count` rounded up to be a multiple of `getDevicePageSize`
template <typename T>
std::size_t getPageAlignedCount(std::size_t count)
{
  const std::size_t page_size = getDevicePageSize();
  const std::size_t size = ((count * sizeof(T) - 1) / page_size + 1) * page_size;
  const std::size_t volume = size / sizeof(T);

  if (size % page_size != 0) {
    throw std::invalid_argument("Buffer size must be multiple of PAGE_SIZE");
  }

  return volume;
};

/// Managed pointer for device-side memory
template <typename T>
class unique_ptr {
private:
  T* _ptr;

public:
  unique_ptr():
    _ptr{nullptr} { }
  unique_ptr(T* ptr):
    _ptr{ptr} { }
  unique_ptr(unique_ptr&& rhs):
    _ptr{rhs.release()} { }

  ~unique_ptr();

  unique_ptr& operator=(unique_ptr&& rhs);

  void reset(T* ptr) {
    operator=(unique_ptr<T>(ptr));
  }

  unique_ptr& operator=(T* ptr) {
    reset(ptr);
    return *this;
  }

  T* release() {
    T* ptr = _ptr;
    _ptr = nullptr;
    return ptr;
  }

  const T* get() const {
    return _ptr;
  }
  T* get() {
    return _ptr;
  }

};


/// Basic wrapper for device stream
class Stream {
private:
  cudaStream_t _stream;

public:
  Stream(unsigned int flags);
  ~Stream();

  cudaStream_t& get() {
    return _stream;
  }

  void synchronize();

};

/// Copy data from device to host (async)
void asyncCopyToHost(Stream& stream, void* src, void* dst, std::size_t count);
/// Copy data from host to device (async)
void asyncCopyToDevice(Stream& stream, void* src, void* dst, std::size_t count);

}

}

}

}

#endif
