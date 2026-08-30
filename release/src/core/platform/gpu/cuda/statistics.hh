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

#ifndef GPU_CUDA_STATISTICS_HH
#define GPU_CUDA_STATISTICS_HH

#include "statistics.h"

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

namespace olb {

namespace gpu {

namespace cuda {

/// Plain pair type with single-value constructor for use in gpu::cuda::maximum_and_plus
template <typename T, typename U>
struct pair {
  T first;
  U second;

  pair() any_platform:
    first{},
    second{}
  { }

  pair(T init) any_platform:
    first{init},
    second{init}
  { }

  pair(T a, U b) any_platform:
    first{a},
    second{b}
  { }

};

/// Function object for simulateneously computing maximum and sum in a single thrust::reduce
template <typename T>
struct maximum_and_plus {
  using first_argument_type  = pair<T,T>;
  using second_argument_type = pair<T,T>;
  using result_type          = pair<T,T>;

  __device__ result_type operator()(const first_argument_type&  lhs,
                                    const second_argument_type& rhs) const {
    return {lhs.first < rhs.first ? rhs.first : lhs.first,
            lhs.second + rhs.second};
  }
};

}

}

template <typename T, typename DESCRIPTOR>
void StatisticsPostProcessor::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& blockLattice)
{
  if (!blockLattice.statisticsEnabled()) {
    return;
  }

  /// Todo: Reimplement these four separate reductions as single kernel
  auto& statisticGenerated = blockLattice.template getField<descriptors::STATISTIC_GENERATED>();
  auto& statistic = blockLattice.template getField<descriptors::STATISTIC>();
  std::size_t nCells = thrust::reduce(thrust::device,
                                      statisticGenerated[0].deviceData(),
                                      statisticGenerated[0].deviceData() + blockLattice.getNcells(),
                                      T{0},
                                      thrust::plus<T>());
  T rhoSum = thrust::reduce(thrust::device,
                            statistic[0].deviceData(),
                            statistic[0].deviceData() + blockLattice.getNcells(),
                            T{0},
                            thrust::plus<T>());
  auto [maxU, uSqrSum]  = thrust::reduce(thrust::device,
                                         statistic[1].deviceData(),
                                         statistic[1].deviceData() + blockLattice.getNcells(),
                                         gpu::cuda::pair{std::numeric_limits<T>::min(), T{0}},
                                         gpu::cuda::maximum_and_plus<T>());

  typename LatticeStatistics<T>::Aggregatable statistics{
    .nCells = nCells,
    .avRho = rhoSum,
    .avEnergy = uSqrSum,
    .maxU = maxU
  };
  blockLattice.getStatistics().incrementStats(statistics);

  blockLattice.getStatistics().reset();
}


}

#endif
