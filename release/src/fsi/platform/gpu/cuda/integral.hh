/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef FSI_PLATFORM_GPU_CUDA_INTEGRAL_HH
#define FSI_PLATFORM_GPU_CUDA_INTEGRAL_HH

#include "integral.h"

#include <thrust/partition.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/tuple.h>


namespace olb {

namespace gpu {

namespace cuda {

template <typename F, std::size_t... INDICES>
auto make_thrust_tuple_f(F&& f, std::index_sequence<INDICES...>) {
  return thrust::make_tuple(f(INDICES)...);
}

template <unsigned D, typename F>
auto make_thrust_tuple_f(F&& f) {
  return make_thrust_tuple_f<F>(std::forward<F&&>(f), std::make_index_sequence<D>{});
}

template <typename A, std::size_t... AIs>
auto thrust_tuple_cat(A& a, std::index_sequence<AIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...);
}

template <typename A, typename B, std::size_t... AIs, std::size_t... BIs>
auto thrust_tuple_cat(A& a, B& b, std::index_sequence<AIs...>, std::index_sequence<BIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...,
                            thrust::get<BIs>(b)...);
}

template <typename A, typename B, typename C, std::size_t... AIs, std::size_t... BIs, std::size_t... CIs>
auto thrust_tuple_cat(A& a, B& b, C& c, std::index_sequence<AIs...>, std::index_sequence<BIs...>, std::index_sequence<CIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...,
                            thrust::get<BIs>(b)...,
                            thrust::get<CIs>(c)...);
}

template <typename... ARGS>
auto thrust_tuple_cat(ARGS&&... args) {
  return thrust_tuple_cat(args...,
                          std::make_index_sequence<thrust::tuple_size<ARGS>::value>{}...);
}

template <typename T, unsigned D>
auto make_thrust_tuple_of(Vector<T,D>& f) {
  return gpu::cuda::make_thrust_tuple_f<f.d>([&](unsigned iDim) {
    return f[iDim];
  });
}

template <typename T>
auto make_thrust_tuple_of(T& f) {
  return thrust::make_tuple(f);
}

template <typename VECTOR>
auto make_thrust_tuple_of_device_data(VECTOR& f) {
  return gpu::cuda::make_thrust_tuple_f<f.d>([&](unsigned iDim) {
    return f[iDim].deviceData();
  });
}

struct first_tuple_element_non_zero {
  template <typename TUPLE>
  bool operator()(const TUPLE& t) __host__ __device__ {
    return thrust::get<0>(t) != 0;
  }
};

struct sum_of_thrust_tuple {
  template <typename TUPLE, std::size_t... INDICES>
  auto operator()(const TUPLE& t1, const TUPLE& t2,
                  std::index_sequence<INDICES...>) __host__ __device__ {
    return thrust::make_tuple((thrust::get<INDICES>(t1) + thrust::get<INDICES>(t2))...);
  }

  template <typename TUPLE>
  auto operator()(const TUPLE& t1, const TUPLE& t2) __host__ __device__ {
    return operator()(t1, t2, std::make_index_sequence<thrust::tuple_size<TUPLE>::value>{});
  }
};

}

}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::setup(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& blockLattice)
{
  blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>();
}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& blockLattice)
{
  auto& tag = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();

  auto begin = thrust::make_zip_iterator(gpu::cuda::thrust_tuple_cat(
                 gpu::cuda::make_thrust_tuple_of_device_data(tag),
                 gpu::cuda::make_thrust_tuple_of_device_data(
                   blockLattice.template getField<FIELDS>())...
               ));
  auto end = thrust::partition(thrust::device,
                               begin,
                               begin + blockLattice.getNcells(),
                               gpu::cuda::first_tuple_element_non_zero{});
  std::size_t nActiveCells = end - begin;

  if (nActiveCells > 0) {
    // Group per-cell field values by FSI tag
    thrust::sort_by_key(thrust::device,
                        tag[0].deviceData(),
                        tag[0].deviceData() + nActiveCells,
                        thrust::make_zip_iterator(gpu::cuda::thrust_tuple_cat(
                          gpu::cuda::make_thrust_tuple_of_device_data(
                            blockLattice.template getField<FIELDS>())...
                        )));
  }

  // Sum up grouped force values
  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>().parameters;
  if (nActiveCells > 0) {
    auto reducedElementTag = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
    const std::size_t nKeys = thrust::reduce_by_key(
      thrust::device,
      tag[0].deviceData(),
      tag[0].deviceData() + nActiveCells,
      thrust::make_zip_iterator(gpu::cuda::thrust_tuple_cat(
        gpu::cuda::make_thrust_tuple_of_device_data(
          blockLattice.template getField<FIELDS>())...
      )),
      reducedElementTag,
      thrust::make_zip_iterator(gpu::cuda::thrust_tuple_cat(
        gpu::cuda::make_thrust_tuple_of(parameters.template get<fields::array_of<FIELDS>>())...
      )),
      thrust::equal_to<T>{},
      gpu::cuda::sum_of_thrust_tuple{}
    ).first - reducedElementTag;
    parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(nKeys);
  } else {
    parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(0);
  }

  if (nActiveCells > 0) {
    thrust::fill(thrust::device,
                 tag[0].deviceData(),
                 tag[0].deviceData() + nActiveCells,
                 0);
  }
}

}

#endif
