/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Yuji (Sam) Shimojima
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

#ifndef GPU_CUDA_RUSUCTION_OPERATORS_HH
#define GPU_CUDA_RUSUCTION_OPERATORS_HH

#include "fieldReduction.h"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace olb {

namespace gpu {

namespace cuda {

// ref : https://github.com/NVIDIA/cuda-samples/blob/master/Samples/2_Concepts_and_Techniques/reduction/reduction_kernel.cu
// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template <class T>
struct SharedMemory {
  __device__ inline operator T*()
  {
    extern __shared__ int __smem[];
    return (T*)__smem;
  }

  __device__ inline operator const T*() const
  {
    extern __shared__ int __smem[];
    return (T*)__smem;
  }
};

// specialize for double to avoid unaligned memory
// access compile errors
template <>
struct SharedMemory<double> {
  __device__ inline operator double*()
  {
    extern __shared__ double __smem_d[];
    return (double*)__smem_d;
  }

  __device__ inline operator const double*() const
  {
    extern __shared__ double __smem_d[];
    return (double*)__smem_d;
  }
};
} // namespace cuda
} // namespace gpu
/// ref: https://github.com/NVIDIA/cuda-samples/blob/master/Samples/2_Concepts_and_Techniques/reduction/reduction_kernel.cu
template <typename T, typename DESCRIPTOR, typename REDUCTION_OP, typename CONDITION>
void reduceKernelInBlock(thrust::device_ptr<const T> field, T* g_odata,
                         gpu::cuda::DeviceBlockLattice<T, DESCRIPTOR> lattice, CellID size) __global__
{
  namespace cg = cooperative_groups; // ref :: https://developer.nvidia.com/blog/cooperative-groups/
  CONDITION mask {};
  // Handle to thread block group
  cg::thread_block cta   = cg::this_thread_block();
  T*               sdata = gpu::cuda::SharedMemory<T>();

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  unsigned int                   tid   = threadIdx.x;
  CellID                         iCell = blockIdx.x * (blockDim.x * 2) + threadIdx.x; //for 2element per 1thread
  gpu::cuda::Cell<T, DESCRIPTOR> cell(lattice, iCell);

  T mySum = (iCell < size) && mask(cell) && (cell.template getField<field::reduction::TAG_CORE>() == (int)1)
                ? field[iCell]
                : (T)0.0;

  //for 2elements per 1thread
  if (iCell + blockDim.x < size) {
    gpu::cuda::Cell<T, DESCRIPTOR> cell2(lattice, iCell + blockDim.x);

    if (mask(cell2) && (cell2.template getField<field::reduction::TAG_CORE>() == (int)1)) {
      mySum = REDUCTION_OP {}(mySum, thrust::raw_pointer_cast(field)[iCell + blockDim.x]);
    }
  }
  sdata[tid] = mySum;
  cg::sync(cta);

  // do reduction in shared mem
  for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (tid < stride) {
      sdata[tid] = mySum = REDUCTION_OP {}(mySum, sdata[tid + stride]);
    }

    cg::sync(cta);
  }

  // write result for this block to global mem
  if (tid == 0)
    g_odata[blockIdx.x] = mySum;
  return;
}

template <typename T, typename REDUCTION_OP>
void reduceKernelInGrid(const T* partialSums, int partialCount, T* d_result) __global__
{
  namespace cg           = cooperative_groups; // ref :: https://developer.nvidia.com/blog/cooperative-groups/
  cg::thread_block cta   = cg::this_thread_block();
  T*               sdata = gpu::cuda::SharedMemory<T>();

  unsigned int tid = threadIdx.x;
  unsigned int i   = threadIdx.x + blockDim.x * blockIdx.x;

  T mySum = 0.0;

  for (int idx = i; idx < partialCount; idx += blockDim.x) {
    mySum = REDUCTION_OP {}(mySum, partialSums[idx]);
  }
  sdata[tid] = mySum;
  cg::sync(cta);

  for (unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (tid < stride) {
      sdata[tid] = mySum = REDUCTION_OP{}(mySum, sdata[tid + stride]);
    }
    cg::sync(cta);
  }

  if (tid == 0) {
    d_result[0] = sdata[0];
  }
  return;
}

template <typename T, typename DESCRIPTOR, typename REDUCTION_OP, typename CONDITION>
void reductionFunctionDevice(thrust::device_ptr<const T>                              devPtr,
                             ConcreteBlockLattice<T, DESCRIPTOR, Platform::GPU_CUDA>& blockLattice, T* paramInDevice)
{

  gpu::cuda::DeviceBlockLattice<T, DESCRIPTOR> lattice(blockLattice);
  const auto                                   nCells     = lattice.getNcells();
  const auto                                   block_size = 32;
  const auto block_count = (nCells + (block_size * 2) - 1) / (block_size * 2); //for 2elements per 1thread

  thrust::device_vector<T> d_partial(block_count, (T)0.0);

  reduceKernelInBlock<T, DESCRIPTOR, REDUCTION_OP, CONDITION><<<block_count, block_size, block_size * sizeof(T)>>>(
      devPtr, thrust::raw_pointer_cast(d_partial.data()), lattice, nCells);
  gpu::cuda::device::check();

  reduceKernelInGrid<T, REDUCTION_OP><<<1, block_size, block_size * sizeof(T)>>>(
      thrust::raw_pointer_cast(d_partial.data()), block_count, paramInDevice);
  gpu::cuda::device::check();

  return;
}

template <typename FIELD, typename REDUCTION_OP, typename CONDITION>
template <typename T, typename DESCRIPTOR>
void BlockLatticeFieldReductionO<FIELD, REDUCTION_OP, CONDITION>::type<ConcreteBlockLattice<
    T, DESCRIPTOR, Platform::GPU_CUDA>>::setup(ConcreteBlockLattice<T, DESCRIPTOR, Platform::GPU_CUDA>& blockLattice)
{
  blockLattice.template getData<OperatorParameters<BlockLatticeFieldReductionO>>();
}

template <typename FIELD, typename REDUCTION_OP, typename CONDITION>
template <typename T, typename DESCRIPTOR>
void BlockLatticeFieldReductionO<FIELD, REDUCTION_OP, CONDITION>::type<ConcreteBlockLattice<
    T, DESCRIPTOR, Platform::GPU_CUDA>>::apply(ConcreteBlockLattice<T, DESCRIPTOR, Platform::GPU_CUDA>& blockLattice)
{
  auto&       parameters = blockLattice.template getData<OperatorParameters<BlockLatticeFieldReductionO>>().parameters;
  const auto& blockField = blockLattice.template getField<FIELD>();

  FieldD<T, DESCRIPTOR, fields::array_of<FIELD>> elementField =
      parameters.template get<fields::array_of<FIELD>>(); //not auto. because of considering 1 dimensional field

  /*note
   -blockField[iD].deviceData(); is return  row pointer T* in device.
   -thrust::device_pointer_cast(elementField[iD])[0] is that thrust::device_pointer_cast is casting to thrust::device_ptr<T> from T* in device.
   -thrust::device_pointer_cast(
          blockField[iD].deviceData()); is that thrust::device_pointer_cast is casting to thrust::device_ptr<T> from T* in device.
   -thrust::raw_pointer_cast(d_partial.data());, with thrust::device_vector<T> d_partial(block_count, (T)0.0);, is that thrust::raw_pointer_cast is casting to T* from thrust::raw_pointer_cast in device.

  */

  for (unsigned iD = 0; iD < blockField.d; iD++) {
    cudaPointerAttributes attributes;
    cudaError_t           error = cudaPointerGetAttributes(&attributes, elementField[iD]);
    if (error == cudaSuccess) {
      if (attributes.devicePointer != nullptr) {
      }
      else if (attributes.hostPointer != nullptr) {
        std::cout << "elementField is on the host" << std::endl;
      }
    }
    else {
      std::cerr << "Error in cudaPointerGetAttributes: " << cudaGetErrorString(error) << std::endl;
    }
    auto devPtr = thrust::device_pointer_cast(
        blockField[iD].deviceData()); // thrust::device_pointer_cast is casting to thrust::device_ptr<T> from T*.
    if (devPtr == nullptr) {
      std::cerr << "Invalid device pointer at index " << iD << std::endl;
      return;
    }
    reductionFunctionDevice<T, DESCRIPTOR, REDUCTION_OP, CONDITION>(devPtr, blockLattice, elementField[iD]);
  }
}

} // namespace olb

#endif //GPU_CUDA_RUSUCTION_OPERATORS_HH
