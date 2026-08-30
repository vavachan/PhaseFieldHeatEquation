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

#ifndef GPU_CUDA_OPERATOR_HH
#define GPU_CUDA_OPERATOR_HH

#include "operator.h"

#include "context.hh"
#include "dynamics.hh"

#include "solver/names.h"

namespace olb {

/// Mask of non-overlap block subdomain
struct CollisionSubdomainMask;

/// Implementations of GPU specifics
namespace gpu {

/// Implementations of Nvidia CUDA specifics
namespace cuda {

/// Masked application of DYNAMICS::collide for use in kernel::call_operators
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class MaskedCollision {
private:
  /// Pointer to on-device parameters structure to be passed to DYNAMICS::collide
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS> _parameters;
  /// Pointer to on-device mask array
  bool* _mask;

  CellStatistic<T> apply(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    return DYNAMICS().collide(cell, _parameters);
  }

public:
  /// Constructor (commonly called on the host side)
  /**
   * See e.g. getFusedCollisionO
   **/
  MaskedCollision(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>& parameters, bool* mask) any_platform:
    _parameters{parameters},
    _mask{mask}
  { }

  /// Chainable call operator for use in kernel::call_operators
  /**
   * Returns true iff MaskedCollision applies to iCell, enabling easy chaining
   * by a fold expression to yield a fused collision kernel.
   **/
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (_mask[iCell]) {
      apply(lattice, iCell);
      return true;
    }
    return false;
  }

  /// Chainable call operator with statistics storage
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    if (_mask[iCell]) {
      statistic = apply(lattice, iCell);
      return true;
    }
    return false;
  }

};

/// List-based application of DYNAMICS::collide for use in kernel::call_list_operators
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ListedCollision {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS> _parameters;

public:
  ListedCollision(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>& parameters) __host__:
    _parameters{parameters}
  { }

  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    DYNAMICS().collide(cell, _parameters);
    return true;
  }

  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    statistic = DYNAMICS().collide(cell, _parameters);
    return true;
  }

};


/// Masked application of OPERATOR::apply
template <typename OPERATOR>
class MaskedPostProcessor {
private:
  /// Pointer to on-device mask array
  bool* _mask;

public:
  MaskedPostProcessor(bool* mask) any_platform:
    _mask{mask}
  { }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (_mask[iCell]) {
      Cell<T,DESCRIPTOR> cell(lattice, iCell);
      OPERATOR().apply(cell);
      return true;
    }
    return false;
  }

};

/// List-based application of OPERATOR::apply
/**
 * Most common approach to calling post processors on device data
 **/
template <typename OPERATOR>
struct ListedPostProcessor {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell);
    return true;
  }

};

/// List-based application of OPERATOR::apply with parameters
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ListedPostProcessorWithParameters {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,OPERATOR> _parameters;

public:
  ListedPostProcessorWithParameters(ParametersOfOperatorD<T,DESCRIPTOR,OPERATOR>& parameters) __host__:
    _parameters{parameters}
  { }

  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell, _parameters);
    return true;
  }

};

/// Unrestricted application of COUPLING::apply
template <typename COUPLER>
struct UnmaskedCoupling {
  template <typename CONTEXT>
  bool operator()(CONTEXT& lattices,
                  CellID iCell) __device__ {
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      return Cell{lattices.get(name), iCell};
    });
    COUPLER().apply(cells);
    return true;
  }
};

/// Unrestricted application of COUPLING::apply with parameters
template <typename COUPLER, typename COUPLEES>
class UnmaskedCouplingWithParameters {
private:
  typename COUPLER::parameters::template decompose_into<
    AbstractCouplingO<COUPLEES>::ParametersD::template include_fields
  > _parameters;

public:
  template <typename PARAMETERS>
  UnmaskedCouplingWithParameters(PARAMETERS& parameters) any_platform:
    _parameters{parameters}
  { }

  template <typename CONTEXT>
  bool operator()(CONTEXT& lattices,
                  CellID iCell) __device__ {
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      return Cell{lattices.get(name), iCell};
    });
    COUPLER().apply(cells, _parameters);
    return true;
  }

};

/// Unrestricted application of COUPLING::apply with parameters
template <typename COUPLER, typename COUPLEES>
class ParticleCouplingWithParameters {
private:
  typename AbstractCouplingO<COUPLEES>::ParametersD::template include<
    typename meta::merge<
      typename COUPLER::parameters,
      meta::list<
        fields::converter::PHYS_DELTA_X,
        fields::BLOCK_LOWER
      >
    >
  > _parameters;

public:
  template <typename PARAMETERS>
  ParticleCouplingWithParameters(PARAMETERS& parameters) any_platform:
    _parameters{parameters}
  { }

  template <typename CONTEXT>
  bool operator()(CONTEXT& lattices,
                  CellID iParticle) __device__ {
    using V = typename CONTEXT::template value_t<names::Points>::value_t;
    using DESCRIPTOR = typename CONTEXT::template value_t<names::Points>::descriptor_t;
    const auto originR = _parameters.template get<fields::BLOCK_LOWER>();
    const auto deltaX = _parameters.template get<fields::converter::PHYS_DELTA_X>();
    const Vector<V*,DESCRIPTOR::d> physRs(lattices.get(meta::id<names::Points>{})
                                                  .template getField<fields::PHYS_R>());
    const PhysR<V,DESCRIPTOR::d> physR{physRs, iParticle};
    const LatticeR<DESCRIPTOR::d> latticeR = util::floor((physR - originR) / deltaX + V{0.5});
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      if constexpr (std::is_same_v<typename decltype(name)::type, names::Points>) {
        return DataOnlyCell{lattices.get(name), iParticle};
      } else {
        const auto iCell = lattices.get(name).getCellId(latticeR);
        return Cell{lattices.get(name), iCell};
      }
    });
    COUPLER().apply(cells, _parameters);
    return true;
  }

};

/// Helper for constructing fused collision operators
/**
 * This is a convenient way for potentially improving performance by
 * injecting application knowledge. E.g. if the lattice contains
 * primarily BGK and BounceBack dynamics this can be declared using:
 *
 * \code{.cpp}
 * superLattice.forBlocksOnPlatform<Platform::GPU_CUDA>([](auto& block) {
 *   block.setCollisionO(
 *     gpu::cuda::getFusedCollisionO<T,DESCRIPTOR,
 *                                   BGKdynamics<T,DESCRIPTOR>,
 *                                   BounceBack<T,DESCRIPTOR>,
 *                                   BounceBackVelocity<T,DESCRIPTOR>>());
 * });
 * \endcode
 **/
template <typename T, typename DESCRIPTOR, typename... DYNAMICS>
std::function<void(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>&)>
getFusedCollisionO() {
  return [](ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) {
    bool* subdomain = block.template getData<CollisionSubdomainMask>().deviceData();
    DeviceContext<T,DESCRIPTOR> lattice(block);
    if (block.statisticsEnabled()) {
      call_operators_with_statistics(
        lattice,
        subdomain,
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{
          block.template getData<OperatorParameters<DYNAMICS>>().parameters,
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    } else {
      call_operators(
        lattice,
        subdomain,
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{
          block.template getData<OperatorParameters<DYNAMICS>>().parameters,
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    }
  };
}


/// CUDA kernels to execute collisions and post processors
namespace kernel {

/// CUDA kernel for applying purely local collision steps
template <typename CONTEXT, typename... OPERATORS>
void call_operators(CONTEXT lattice, bool* subdomain, OPERATORS... ops) __global__ {
  const CellID iCell = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iCell < lattice.getNcells()) || !subdomain[iCell]) {
    return;
  }
  (ops(lattice, iCell) || ... );
}

/// CUDA kernel for applying purely local collision steps while tracking statistics
/**
 * Statistics data is reduced by StatisticsPostProcessor
 **/
template <typename CONTEXT, typename... OPERATORS>
void call_operators_with_statistics(CONTEXT lattice, bool* subdomain, OPERATORS... ops) __global__ {
  const CellID iCell = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iCell < lattice.getNcells()) || !subdomain[iCell]) {
    return;
  }
  typename CONTEXT::value_t** statistic = lattice.template getField<descriptors::STATISTIC>();
  int* statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
  CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};
  if ((ops(lattice, iCell, cellStatistic) || ... )) {
    if (cellStatistic) {
      statisticGenerated[iCell] = 1;
      statistic[0][iCell] = cellStatistic.rho;
      statistic[1][iCell] = cellStatistic.uSqr;
    } else {
      statisticGenerated[iCell] = 0;
      statistic[0][iCell] = 0;
      statistic[1][iCell] = 0;
    }
  }
}

/// CUDA kernel for applying generic OPERATORS with OperatorScope::PerCell or ListedCollision
template <typename CONTEXT, typename... OPERATORS>
void call_list_operators(CONTEXT lattice,
                         const CellID* indices, std::size_t nIndices,
                         OPERATORS... ops) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  (ops(lattice, indices[iIndex]) || ... );
}

/// CUDA kernel for applying ListedCollision
/**
 * Statistics data is reduced by StatisticsPostProcessor
 **/
template <typename CONTEXT, typename... OPERATORS>
void call_list_operators_with_statistics(CONTEXT lattice,
                                         const CellID* indices, std::size_t nIndices,
                                         OPERATORS... ops) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  typename CONTEXT::value_t** statistic = lattice.template getField<descriptors::STATISTIC>();
  int* statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
  CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};
  if ((ops(lattice, indices[iIndex], cellStatistic) || ... )) {
    if (cellStatistic) {
      statisticGenerated[indices[iIndex]] = 1;
      statistic[0][indices[iIndex]] = cellStatistic.rho;
      statistic[1][indices[iIndex]] = cellStatistic.uSqr;
    } else {
      statisticGenerated[indices[iIndex]] = 0;
      statistic[0][indices[iIndex]] = 0;
      statistic[1][indices[iIndex]] = 0;
    }
  }
}

/// CUDA kernel for applying UnmaskedCoupling(WithParameters)
template <typename CONTEXTS, typename... OPERATORS>
void call_coupling_operators(CONTEXTS lattices, bool* subdomain, OPERATORS... ops) __global__ {
  const CellID iCell = blockIdx.x * blockDim.x + threadIdx.x;
  const auto nCells = lattices.template get<0>().getNcells();
  if (!(iCell < nCells) || !subdomain[iCell]) {
    return;
  }
  (ops(lattices, iCell) || ... );
}

/// CUDA kernel for applying UnmaskedCoupling(WithParameters)
template <typename CONTEXTS, typename... OPERATORS>
void call_particle_coupling_operators(CONTEXTS lattices,
                                      std::size_t nParticles,
                                      OPERATORS... ops) __global__ {
  const CellID iParticle = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iParticle < nParticles)) {
    return;
  }
  (ops(lattices, iParticle) || ... );
}


/// CUDA kernel for constructing on-device ConcreteDynamics
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename PARAMETERS=typename DYNAMICS::ParametersD>
void construct_dynamics(void* target, PARAMETERS* parameters) __global__ {
  new (target) ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(parameters);
}

}

/// Apply masked collision operators to lattice
/**
 * ARGS are instances of MaskedCollision or DynamicDispatchCollision
 **/
template <typename CONTEXT, typename... ARGS>
void call_operators(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice (async)
template <typename CONTEXT, typename... ARGS>
void async_call_operators(cudaStream_t stream, CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators<CONTEXT,ARGS...><<<block_count,block_size,0,stream>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice while tracking statistics
/**
 * ARGS are instances of MaskedCollision or DynamicDispatchCollision
 **/
template <typename CONTEXT, typename... ARGS>
void call_operators_with_statistics(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators_with_statistics<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice while tracking statistics (async)
template <typename CONTEXT, typename... ARGS>
void async_call_operators_with_statistics(cudaStream_t stream, CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators_with_statistics<CONTEXT,ARGS...><<<block_count,block_size,0,stream>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply operators to listed cell indices
/**
 * Used to call post processors in ConcreteBlockO with OperatorScope::PerCell
 **/
template <typename CONTEXT, typename... ARGS>
void call_list_operators(CONTEXT& lattice,
                         const gpu::cuda::Column<CellID>& cells,
                         ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators<CONTEXT,ARGS...><<<block_count, block_size>>>(
    lattice,
    cells.deviceData(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply operators to listed cell indices (async version)
template <typename CONTEXT, typename... ARGS>
void async_call_list_operators(cudaStream_t stream,
                               CONTEXT& lattice,
                               const gpu::cuda::Column<CellID>& cells,
                               ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators<<<block_count,block_size,0,stream>>>(
    lattice,
    cells.deviceData(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply ListedCollision with statistics (async version)
template <typename CONTEXT, typename... ARGS>
void async_call_list_operators_with_statistics(cudaStream_t stream,
                                               CONTEXT& lattice,
                                               const gpu::cuda::Column<CellID>& cells,
                                               ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators_with_statistics<<<block_count,block_size,0,stream>>>(
    lattice,
    cells.deviceData(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply coupling on subdomain
template <typename CONTEXT, typename... ARGS>
void call_coupling_operators(CONTEXT& lattices, bool* subdomain, ARGS&&... args) {
  const auto nCells = lattices.template get<0>().getNcells();
  const auto block_size = 32;
  const auto block_count = (nCells + block_size - 1) / block_size;
  kernel::call_coupling_operators<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattices, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply coupling on particles and lattices
template <typename CONTEXT, typename... ARGS>
void call_particle_coupling_operators(CONTEXT& lattices, ARGS&&... args) {
  const auto nParticles = lattices.get(meta::id<names::Points>{}).getNcells();
  const auto block_size = 32;
  const auto block_count = (nParticles + block_size - 1) / block_size;
  kernel::call_particle_coupling_operators<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattices, nParticles, std::forward<decltype(args)>(args)...);
  device::check();
}

}

}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::ConcreteBlockCollisionO():
  _dynamics(new DYNAMICS()),
  _parameters(nullptr),
  _mask(nullptr),
  _cells(0),
  _modified(true),
  _stream(cudaStreamDefault)
{ }

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::applyDominant(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
  ConcreteBlockMask<T,Platform::GPU_CUDA>&               subdomain)
{
  using namespace gpu::cuda;
  DeviceContext<T,DESCRIPTOR> lattice(block);
  if (block.statisticsEnabled()) {
    call_operators_with_statistics(
      lattice,
      subdomain.deviceData(),
      MaskedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters, _mask->deviceData()},
      DynamicDispatchCollision{});
  } else {
    call_operators(
      lattice,
      subdomain.deviceData(),
      MaskedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters, _mask->deviceData()},
      DynamicDispatchCollision{});
  }
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::applyIndividual(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
  ConcreteBlockMask<T,Platform::GPU_CUDA>&               subdomain)
{
  using namespace gpu::cuda;
  DeviceContext<T,DESCRIPTOR> lattice(block);
  // Primitive heuristic for preferring mask-based to list-based dispatch
  if (_mask->weight() > 0.5*subdomain.weight()) {
    if (block.statisticsEnabled()) {
      async_call_operators_with_statistics(
        _stream.get(),
        lattice,
        subdomain.deviceData(),
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters, _mask->deviceData()});
    } else {
      async_call_operators(
        _stream.get(),
        lattice,
        subdomain.deviceData(),
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters, _mask->deviceData()});
    }

  // Use list of cell indices
  } else {
    // Update cell list from mask
    if (_modified) {
      _cells.clear();
      for (CellID iCell=0; iCell  < block.getNcells(); ++iCell) {
        if (_mask->operator[](iCell)) {
          _cells.push_back(iCell);
        }
      }
      _cells.setProcessingContext(ProcessingContext::Simulation);
      _modified = false;
    }

    if (block.statisticsEnabled()) {
      async_call_list_operators_with_statistics(
        _stream.get(),
        lattice,
        _cells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters});
    } else {
      async_call_list_operators(
        _stream.get(),
        lattice,
        _cells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters});
    }
  }
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::setup(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block)
{
  // Fetch pointers to DYNAMICS-specific parameter and mask data
  _parameters = &block.template getData<OperatorParameters<DYNAMICS>>();
  _mask = &block.template getData<DynamicsMask<DYNAMICS>>();

  {
    // Construct on-device dynamics proxy for dynamic dispatch
    _deviceDynamics = gpu::cuda::device::malloc<gpu::cuda::ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>>(1);
    gpu::cuda::kernel::construct_dynamics<T,DESCRIPTOR,DYNAMICS><<<1,1>>>(
      _deviceDynamics.get(),
      _parameters->deviceData());
    gpu::cuda::device::check();

    // Fetch pointer to on-device dynamic-dispatch field
    _dynamicsOfCells = block.template getField<gpu::cuda::DYNAMICS<T,DESCRIPTOR>>()[0].data();
  }
}

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCell>::ConcreteBlockO():
  _cells(0),
  _modified{false},
  _stream{cudaStreamDefault}
{ }

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
void ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCell>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block)
{
  if (_cells.size() > 0) {
    if (_modified) {
      _cells.deduplicate();
      _cells.setProcessingContext(ProcessingContext::Simulation);
      _modified = false;
    }
    gpu::cuda::DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
    gpu::cuda::async_call_list_operators(_stream.get(),
                                         lattice,
                                         _cells,
                                         gpu::cuda::ListedPostProcessor<OPERATOR>{});
  }
}

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCellWithParameters>::ConcreteBlockO():
  _cells(0),
  _modified{false},
  _stream{cudaStreamDefault}
{ }

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
void ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCellWithParameters>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block)
{
  if (_cells.size() > 0) {
    if (_modified) {
      _cells.deduplicate();
      _cells.setProcessingContext(ProcessingContext::Simulation);
      _modified = false;
    }
    using namespace gpu::cuda;
    DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
    async_call_list_operators(_stream.get(),
                              lattice,
                              _cells,
                              ListedPostProcessorWithParameters<T,DESCRIPTOR,OPERATOR>{_parameters->parameters});
  }
}


template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCell>::set(
  CellID iCell, bool state)
{
  if (!_mask) {
    _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                               Platform::GPU_CUDA>>(
      _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
    );
  }
  _mask->set(iCell, state);
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCell>::execute()
{
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    return gpu::cuda::DeviceBlockLattice{*_lattices.get(name)};
  });
  if (_mask) {
    _mask->setProcessingContext(ProcessingContext::Simulation);
    gpu::cuda::call_coupling_operators(
      deviceLattice, _mask->deviceData(),
      gpu::cuda::UnmaskedCoupling<COUPLER>{});
  } else {
    auto& mask = _lattices.template get<0>()->template getData<CollisionSubdomainMask>();
    gpu::cuda::call_coupling_operators(
      deviceLattice, mask.deviceData(),
      gpu::cuda::UnmaskedCoupling<COUPLER>{});
  }
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>::set(
  CellID iCell, bool state)
{
  if (!_mask) {
    _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                               Platform::GPU_CUDA>>(
      _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
    );
  }
  _mask->set(iCell, state);
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>::execute()
{
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    return gpu::cuda::DeviceBlockLattice{*_lattices.get(name)};
  });
  if (_mask) {
    _mask->setProcessingContext(ProcessingContext::Simulation);
    gpu::cuda::call_coupling_operators(
      deviceLattice, _mask->deviceData(),
      gpu::cuda::UnmaskedCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
  } else {
    auto& mask = _lattices.template get<0>()->template getData<CollisionSubdomainMask>();
    gpu::cuda::call_coupling_operators(
      deviceLattice, mask.deviceData(),
      gpu::cuda::UnmaskedCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
  }
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockPointCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>::set(
  CellID iCell, bool state)
{
  throw std::bad_function_call();
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockPointCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>::execute()
{
  if (_lattices.get(meta::id<names::Points>{})->getNcells() == 0) {
    return;
  }

  auto* particleLattice = _lattices.get(meta::id<names::Points>{});
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    if constexpr (std::is_same_v<typename decltype(name)::type, names::Points>) {
      return gpu::cuda::DeviceContext{*_lattices.get(name)};
    } else {
      return gpu::cuda::DeviceBlockLattice{*_lattices.get(name)};
    }
  });

  gpu::cuda::call_particle_coupling_operators(
    deviceLattice,
    gpu::cuda::ParticleCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
}

}

#endif
