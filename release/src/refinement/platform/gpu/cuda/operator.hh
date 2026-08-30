/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef REFINEMENT_PLATFORM_GPU_CUDA_OPERATOR_HH
#define REFINEMENT_PLATFORM_GPU_CUDA_OPERATOR_HH

#include "operator.h"

namespace olb {

namespace gpu::cuda {

namespace kernel {

template <typename COARSE, typename FINE, typename CONTEXT, typename PARAMETERS, typename OPERATOR>
void call_refinement_coupling_operator(COARSE cLattice,
                                       FINE fLattice,
                                       CONTEXT data,
                                       PARAMETERS parameters,
                                       std::size_t n) __global__ {
  const CellID i = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(i < n)) {
    return;
  }

  auto cellIdCoarse = data.template getField<fields::refinement::CELL_ID_COARSE>();
  auto cellIdFine = data.template getField<fields::refinement::CELL_ID_FINE>();

         if constexpr (OPERATOR::scope == refinement::OperatorScope::PerCoarseCell) {
    if (cellIdCoarse[0][i] != 0) {
      Cell cCell{cLattice, cellIdCoarse[0][i]};
      Cell fCell{fLattice, cellIdFine[0][i]};
      DataOnlyCell cData{data, i};
      OPERATOR().apply(cCell, fCell, cData, parameters);
    }
  } else if constexpr (OPERATOR::scope == refinement::OperatorScope::PerFineCell) {
    if constexpr (OPERATOR::data::template contains<fields::refinement::CONTEXT_NEIGHBORS>()) {
      Cell fCell{fLattice, cellIdFine[0][i]};
      refinement::CoarseCell cCell{Cell{cLattice, cellIdCoarse[0][i]},
                                   fLattice.getLatticeR(cellIdFine[0][i])};
      refinement::ContextDataWithNeighbors cData{DataOnlyCell{data, i}};
      OPERATOR().apply(cCell, fCell, cData, parameters);
    } else {
      Cell fCell{fLattice, cellIdFine[0][i]};
      refinement::CoarseCell cCell{Cell{cLattice, cellIdCoarse[0][i]},
                                   fLattice.getLatticeR(cellIdFine[0][i])};
      refinement::ContextData cData{DataOnlyCell{data, i}};
      OPERATOR().apply(cCell, fCell, cData, parameters);
    }
  }
}

}

template <typename COARSE, typename FINE, typename CONTEXT, typename PARAMETERS, typename OPERATOR>
void call_refinement_coupling_operator(COARSE& cLattice, FINE& fLattice,
                                       CONTEXT& data,
                                       PARAMETERS& parameters,
                                       meta::id<OPERATOR>) {
  const std::size_t n = data.getNcells();
  if (n > 0) {
    const auto block_size = 32;
    const auto block_count = (n + block_size - 1) / block_size;
    kernel::call_refinement_coupling_operator<COARSE,FINE,CONTEXT,PARAMETERS,OPERATOR>
      <<<block_count,block_size>>>(cLattice, fLattice, data, parameters, n);
    device::check();
  }
}

}

template <typename T, typename DESCRIPTOR, typename OPERATOR>
void ConcreteBlockRefinementO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR>::apply(
  ConcreteBlockRefinementContextD<T,DESCRIPTOR,Platform::GPU_CUDA>& context)
{
  auto& cLattice = context.getCoarseLattice();
  auto& fLattice = context.getFineLattice();
  auto& data = context.getConcreteData();

  auto parameters = data.template getData<OperatorParameters<OPERATOR>>().parameters;

  gpu::cuda::DeviceBlockLattice deviceCoarseLattice{cLattice};
  gpu::cuda::DeviceBlockLattice deviceFineLattice{fLattice};
  gpu::cuda::DeviceContext deviceData{data};

  gpu::cuda::call_refinement_coupling_operator(deviceCoarseLattice,
                                               deviceFineLattice,
                                               deviceData,
                                               parameters,
                                               meta::id<OPERATOR>{});
}

}

#endif
