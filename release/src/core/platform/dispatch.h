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

#ifndef PLATFORM_DISPATCH_H
#define PLATFORM_DISPATCH_H

#include <stdexcept>

#include "platform.h"
#include "core/expr.h"

/// Top level namespace for all of OpenLB
namespace olb {

/// Dispatcher for concrete platform access
/**
 * See e.g. ConcretizableBlockLattice usage in BlockLattice::getField
 **/
template <typename CONCRETIZABLE, typename F>
inline auto callUsingConcretePlatform(Platform platform, typename CONCRETIZABLE::base_t* ptr, F f)
{
  if constexpr (std::is_same_v<typename CONCRETIZABLE::value_t,Expr>) {
    if (platform != Platform::CPU_SISD) {
      throw std::domain_error("Only CPU_SISD supports Expr value type");
    }
    return f(static_cast<typename CONCRETIZABLE::template type<Platform::CPU_SISD>*>(ptr));
  } else {
    switch (platform) {
#ifdef PLATFORM_CPU_SISD
    case Platform::CPU_SISD:
      return f(static_cast<typename CONCRETIZABLE::template type<Platform::CPU_SISD>*>(ptr));
#endif
#ifdef PLATFORM_CPU_SIMD
    case Platform::CPU_SIMD:
      return f(static_cast<typename CONCRETIZABLE::template type<Platform::CPU_SIMD>*>(ptr));
#endif
#ifdef PLATFORM_GPU_CUDA
    case Platform::GPU_CUDA:
      return f(static_cast<typename CONCRETIZABLE::template type<Platform::GPU_CUDA>*>(ptr));
#endif
    default:
      throw std::invalid_argument("Invalid PLATFORM");
    }
  }
}

/// Dispatcher for read-only concrete platform access
/**
 * See e.g. ConcretizableBlockLattice usage in BlockLattice::getField
 **/
template <typename CONCRETIZABLE, typename F>
inline auto callUsingConcretePlatform(Platform platform, const typename CONCRETIZABLE::base_t* ptr, F f)
{
  return callUsingConcretePlatform<CONCRETIZABLE>(platform,
                                                  const_cast<typename CONCRETIZABLE::base_t*>(ptr),
                                                  [&](auto* concrete) {
    return f(static_cast<typename std::add_const_t<decltype(concrete)>>(concrete));
  });
}

template <typename F>
inline void callUsingConcretePlatform(Platform platform, F f)
{
  switch (platform) {
#ifdef PLATFORM_CPU_SISD
  case Platform::CPU_SISD:
    f(meta::value<Platform::CPU_SISD>{});
    break;
#endif
#ifdef PLATFORM_CPU_SIMD
  case Platform::CPU_SIMD:
    f(meta::value<Platform::CPU_SIMD>{});
    break;
#endif
#ifdef PLATFORM_GPU_CUDA
  case Platform::GPU_CUDA:
    f(meta::value<Platform::GPU_CUDA>{});
    break;
#endif
  default:
    throw std::invalid_argument("Invalid PLATFORM");
  }
}

template <typename CONCRETIZABLE, typename... ARGS>
typename CONCRETIZABLE::base_t* constructUsingConcretePlatform(Platform platform, ARGS&&... args)
{
  if constexpr (std::is_same_v<typename CONCRETIZABLE::value_t,Expr>) {
    if (platform != Platform::CPU_SISD) {
      throw std::domain_error("Only CPU_SISD supports Expr value type");
    }
    return new typename CONCRETIZABLE::template type<Platform::CPU_SISD>(std::forward<decltype(args)>(args)...);
  } else {
    switch (platform) {
#ifdef PLATFORM_CPU_SISD
    case Platform::CPU_SISD:
      return new typename CONCRETIZABLE::template type<Platform::CPU_SISD>(std::forward<decltype(args)>(args)...);
#endif
#ifdef PLATFORM_CPU_SIMD
    case Platform::CPU_SIMD:
      return new typename CONCRETIZABLE::template type<Platform::CPU_SIMD>(std::forward<decltype(args)>(args)...);
#endif
#ifdef PLATFORM_GPU_CUDA
    case Platform::GPU_CUDA:
      return new typename CONCRETIZABLE::template type<Platform::GPU_CUDA>(std::forward<decltype(args)>(args)...);
#endif
    default:
      throw std::invalid_argument("Invalid PLATFORM");
    }
  }
}

}

#endif
