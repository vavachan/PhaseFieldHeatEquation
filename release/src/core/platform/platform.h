/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef PLATFORM_H
#define PLATFORM_H

#include <cstdint>
#include <stdexcept>

/// Top level namespace for all of OpenLB
namespace olb {

/// OpenLB execution targets
enum struct Platform : std::uint8_t {
  CPU_SISD, /// Basic scalar CPU
  CPU_SIMD, /// Vector CPU (AVX2 / AVX-512 collision)
  GPU_CUDA, /// GPU using CUDA
};

/// Verifies requirements for using PLATFORM
template <Platform PLATFORM>
void checkPlatform();

/// OpenLB processing contexts
/**
 * Currently of no relevance for CPU_SISD and CPU_SIMD target platforms
 *
 * Used to control synchronization between mirrored device and host data
 * for non-host processed block lattices.
 *
 * Preliminary for first GPU release.
 **/
enum struct ProcessingContext {
  Evaluation, /// Data available on host for e.g. functor evaluation
  Simulation  /// Data available on device for evolving the simulation
};

namespace stage {

template <ProcessingContext CONTEXT>
struct PreContextSwitchTo { };

}

/// Define preprocessor macros for device-side functions, constant storage
#ifdef __CUDACC__
  #define any_platform __device__ __host__
  #ifdef __CUDA_ARCH__
    #define platform_constant constexpr __constant__
    #define platform_constant_definition constexpr __constant__
  #else
    #define platform_constant constexpr
    #define platform_constant_definition constexpr
  #endif
#else
  #define any_platform
  #define platform_constant constexpr
  #define platform_constant_definition constexpr
#endif

/// Returns true if platform is equal to Platform::CPU_*
constexpr bool isPlatformCPU(Platform platform) {
  switch (platform) {
#ifdef PLATFORM_CPU_SISD
  case Platform::CPU_SISD:
    return true;
#endif
#ifdef PLATFORM_CPU_SIMD
  case Platform::CPU_SIMD:
    return true;
#endif
#ifdef PLATFORM_GPU_CUDA
  case Platform::GPU_CUDA:
    return false;
#endif
  default:
    throw std::invalid_argument("Invalid PLATFORM");
  }
}

}

#endif
