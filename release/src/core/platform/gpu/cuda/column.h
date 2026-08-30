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

#ifndef GPU_CUDA_COLUMN_H
#define GPU_CUDA_COLUMN_H

#include <memory>

#include "core/platform/platform.h"
#include "core/serializer.h"
#include "communication/communicatable.h"

namespace olb {

namespace gpu {

namespace cuda {

namespace device {

struct Stream;

}

/// Plain column for CUDA GPU targets
template<typename T>
class Column final : public AbstractColumn<T>
                   , public Serializable {
private:
  std::size_t _count;

  struct Data;
  std::unique_ptr<Data> _data;

public:
  using value_t = T;

  Column(std::size_t count);
  Column(Column<T>&& rhs);
  Column(const Column<T>& rhs);
  ~Column();

  const T& operator[](std::size_t i) const override;
  T& operator[](std::size_t i) override;

  std::size_t size() const override;

  const T* data() const;
  T* data();

  const T* deviceData() const;
  T* deviceData();

  /// Reset size to zero
  void clear();
  void resize(std::size_t newCount);
  void push_back(T value);
  /// Combined ascending sort and removal of duplicate entries
  void deduplicate();

  void setProcessingContext(ProcessingContext);
  void setProcessingContext(ProcessingContext, device::Stream&);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


/// Virtual memory based cyclic column for usage in ColumnVector
/**
 * Column type used for propagatable population storage using
 * the virtual memory PS pattern.
 **/
template<typename T>
class CyclicColumn final : public AbstractCyclicColumn<T>
                         , public Serializable {
private:
  const std::ptrdiff_t _count;
  const std::size_t    _size;

  struct Data;
  std::unique_ptr<Data> _data;

  T* _deviceBase;
  T* _devicePopulation;

  std::ptrdiff_t _shift;

public:
  using value_t = T;

  CyclicColumn(std::size_t count);
  ~CyclicColumn();

  const T& operator[](std::size_t i) const override;
  T& operator[](std::size_t i) override;

  const T* deviceData() const
  {
    return _devicePopulation;
  }

  T* deviceData()
  {
    return _devicePopulation;
  }

  std::size_t size() const override
  {
    return _count;
  }

  void refresh()
  {
    _devicePopulation = _deviceBase + _shift;
  }

  void rotate(std::ptrdiff_t offset)
  {
    _shift -= offset;
    if (_shift >= _count) {
      _shift -= _count;
    }
    else if (_shift < 0) {
      _shift += _count;
    }
    refresh();
  }

  void resize(std::size_t count)
  {
    throw std::logic_error("Cyclic column can not be resized");
  }

  void setProcessingContext(ProcessingContext);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  void postLoad() override { refresh(); }

};

}

}


/// Declare gpu::cuda::Column as the AbstractColumn implementation for GPU CUDA targets
template <typename T>
struct ImplementationOf<AbstractColumn<T>,Platform::GPU_CUDA> {
  using type = gpu::cuda::Column<T>;
};

/// Declare gpu::cuda::CyclicColumn as the AbstractCyclicColumn implementation for GPU CUDA targets
template <typename T>
struct ImplementationOf<AbstractCyclicColumn<T>,Platform::GPU_CUDA> {
  using type = gpu::cuda::CyclicColumn<T>;
};


/// Communicatable implementation for a single gpu::cuda::Column
/**
 * This is only implemented to match the interface of CPU communication.
 * The actual communication buffers for communicating block lattice data
 * are handled using custom kernel functions in the the Platform::GPU_CUDA
 * specialization of ConcreteBlockCommunicator.
 **/
template <typename T>
class ConcreteCommunicatable<gpu::cuda::Column<T>> final : public Communicatable {
private:
  gpu::cuda::Column<T>& _column;

public:
  ConcreteCommunicatable(gpu::cuda::Column<T>& column):
    _column{column} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    return indices.size() * sizeof(T);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const;

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer);

};

/// Communicatable implementation for a single gpu::cuda::CyclicColumn
template <typename T>
class ConcreteCommunicatable<gpu::cuda::CyclicColumn<T>> final : public Communicatable {
private:
  gpu::cuda::CyclicColumn<T>& _column;

public:
  ConcreteCommunicatable(gpu::cuda::CyclicColumn<T>& column):
    _column{column} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    return indices.size() * sizeof(T);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const;

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer);

};

}

#endif
