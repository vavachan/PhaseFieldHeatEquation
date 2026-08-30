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

#ifndef HETEROGENEOUS_LOAD_BALANCER_H
#define HETEROGENEOUS_LOAD_BALANCER_H

#include "communication/loadBalancer.h"

#include <random>
#include <algorithm>
#include <queue>

namespace olb {

/// Load balancer for heterogeneous CPU-GPU systems
/**
 * Assigns largest cuboids to GPUs until given largeBlockFraction is reached.
 * Remaining (small) cuboids are assigned to CPUs of their GPU-placed neighbors.
 *
 * This balancer should only be used in conjunction with a heterogeneous
 * cuboid decomposition and appropriate, system specific, CPU_SIMD+OpenMP
 * configuration.
 **/
template<typename T>
class HeterogeneousLoadBalancer final : public LoadBalancer<T> {
private:
  std::vector<Platform> _platform;

public:
  HeterogeneousLoadBalancer(CuboidDecomposition<T,3>& cGeometry, T largeBlockFraction=0.9):
    LoadBalancer<T>(0)
  {
    OstreamManager clout(std::cout, "HeterogeneousLoadBalancer");

    std::vector<int> rankBuffer(cGeometry.size(), 0);
    std::vector<int> locBuffer(cGeometry.size(), 0);
    std::vector<std::uint8_t> platformBuffer(cGeometry.size(), static_cast<std::uint8_t>(Platform::CPU_SISD));

    const auto nRank = singleton::mpi().getSize();
    const auto iRank = singleton::mpi().getRank();

    auto& cuboids = cGeometry.cuboids();
    // Prioritize assignment of largest cuboids
    std::sort(cuboids.begin(), cuboids.end(),
              [&](const auto& lhs, const auto& rhs) {
                return lhs.getLatticeVolume() > rhs.getLatticeVolume();
              });

    // Distribute cuboids to ranks on rank 0
    if (iRank == 0) {
      // Total volume for tracking targeted GPU block fraction
      std::size_t totalVolume = std::accumulate(cuboids.begin(), cuboids.end(),
                                                std::size_t{0},
                                                [](const auto& currVolume, const auto& rhs) -> std::size_t {
                                                  return currVolume + rhs.getLatticeVolume();
                                                });

      std::map<int,int> nLoc;
      int jRank = 0;
      int iCuboid = 0;

      std::size_t globalAssignedVolume = 0;
      std::vector<std::size_t> localAssignedVolume(nRank, 0);

      // Assign largest cuboids to GPUs until desired fraction is reached
      do {
        jRank = std::distance(localAssignedVolume.begin(),
                              std::min_element(localAssignedVolume.begin(),
                                               localAssignedVolume.end()));
        rankBuffer[iCuboid] = jRank;
        platformBuffer[iCuboid] = static_cast<std::uint8_t>(Platform::GPU_CUDA);
        locBuffer[iCuboid]  = nLoc[jRank]++;
        localAssignedVolume[jRank] += cuboids[iCuboid].getLatticeVolume();
        globalAssignedVolume += cuboids[iCuboid].getLatticeVolume();
        clout << iCuboid << ", " << jRank << ", assignedVolumeFraction=" << (1.*globalAssignedVolume) / totalVolume << std::endl;
        iCuboid += 1;
      } while (globalAssignedVolume < largeBlockFraction*totalVolume);

      // Compute GPU rank affinity of remaining cuboids
      std::vector<int> preferredRank(cGeometry.size(), -1);
      for (int jCuboid=iCuboid; jCuboid < cGeometry.size(); ++jCuboid) {
        std::vector<int> neighbours;
        cGeometry.getNeighbourhood(jCuboid, neighbours, 1);
        std::vector<int> neighboursInRank(nRank, 0);
        int nGpuNeighbors = 0;
        for (int neighborC : neighbours) {
          if (platformBuffer[neighborC] == static_cast<std::uint8_t>(Platform::GPU_CUDA)) {
            neighboursInRank[rankBuffer[neighborC]] += 1;
            nGpuNeighbors += 1;
          }
        }
        if (nGpuNeighbors > 0) {
          preferredRank[jCuboid] = std::distance(neighboursInRank.begin(), std::max_element(neighboursInRank.begin(),
                                                                                            neighboursInRank.end()));
          clout << "Preferred rank of " << jCuboid << " is " << preferredRank[jCuboid] << std::endl;
        }
      }

      // Distribute remaining GPU-neighboring blocks over ranks using CPU_SIMD platform
      for (int jCuboid=iCuboid; jCuboid < cGeometry.size(); ++jCuboid) {
        if (preferredRank[jCuboid] != -1) {
          auto iRank = preferredRank[jCuboid];

          rankBuffer[jCuboid] = iRank;
          platformBuffer[jCuboid] = static_cast<std::uint8_t>(Platform::CPU_SIMD);
          locBuffer[jCuboid]  = nLoc[iRank]++;
        }
      }

      // Compute fallback CPU rank affinity of remaining cuboids
      std::vector<int> preferredRankFallback(cGeometry.size(), -1);
      for (int jCuboid=iCuboid; jCuboid < cGeometry.size(); ++jCuboid) {
        if (preferredRank[jCuboid] == -1) {
          std::vector<int> neighbours;
          cGeometry.getNeighbourhood(jCuboid, neighbours, 1);
          std::vector<int> neighboursInRank(nRank, 0);
          int nCpuNeighbors = 0;
          for (int neighborC : neighbours) {
            if (platformBuffer[neighborC] != static_cast<std::uint8_t>(Platform::GPU_CUDA)) {
              neighboursInRank[rankBuffer[neighborC]] += 1;
              nCpuNeighbors += 1;
            }
          }
          if (nCpuNeighbors > 0) {
            preferredRankFallback[jCuboid] = std::distance(neighboursInRank.begin(), std::max_element(neighboursInRank.begin(),
                                                                                                      neighboursInRank.end()));
            clout << "Preferred rank of " << jCuboid << " is " << preferredRankFallback[jCuboid] << std::endl;
          }
        }
      }

      // Distribute remaining blocks over ranks using CPU_SIMD platform
      for (int jCuboid=iCuboid; jCuboid < cGeometry.size(); ++jCuboid) {
        if (preferredRank[jCuboid] == -1) {
          auto iRank = preferredRankFallback[jCuboid];

          rankBuffer[jCuboid] = iRank;
          platformBuffer[jCuboid] = static_cast<std::uint8_t>(Platform::CPU_SIMD);
          locBuffer[jCuboid]  = nLoc[iRank]++;
        }
      }

    }

    #ifdef PARALLEL_MODE_MPI
    // Broadcast assignments to all processes
    singleton::mpi().bCast(rankBuffer.data(), rankBuffer.size());
    singleton::mpi().bCast(locBuffer.data(), locBuffer.size());
    singleton::mpi().bCast(platformBuffer.data(), platformBuffer.size());
    #endif

    // Update internal LoadBalancer structure to match given assignment
    for (int iCuboid=0; iCuboid < cGeometry.size(); ++iCuboid) {
      this->_rank[iCuboid] = rankBuffer[iCuboid];
      this->_loc[iCuboid] = locBuffer[iCuboid];

      if (rankBuffer[iCuboid] == iRank) {
        this->_glob.resize(std::max(int(this->_glob.size()), this->_loc[iCuboid]+1));
        this->_platform.resize(this->_glob.size());

        this->_glob[this->_loc[iCuboid]] = iCuboid;
        this->_platform[this->_loc[iCuboid]] = static_cast<Platform>(platformBuffer[iCuboid]);

        this->_size = this->_glob.size();
      }
    }
  }

  Platform platform(int loc) const override {
    return _platform[loc];
  }

  void setPlatform(int loc, Platform platform) {
    _platform[loc] = platform;
  }

};

/// Load balancer for heterogeneous CPU-GPU systems
/**
 * Assigns largest cuboids to GPUs until given largeBlockFraction is reached.
 * Remaining (small) cuboids are assigned to ranks without GPUs.
 *
 * This balancer should only be used in conjunction with a heterogeneous
 * cuboid decomposition and appropriate, system specific, CPU_SIMD+OpenMP
 * configuration.
 **/
template<typename T>
class OrthogonalHeterogeneousLoadBalancer final : public LoadBalancer<T> {
private:
  std::vector<Platform> _platform;

public:
  OrthogonalHeterogeneousLoadBalancer(CuboidDecomposition<T,3>& cGeometry, T largeBlockFraction=0.9):
    LoadBalancer<T>(0)
  {
    OstreamManager clout(std::cout, "OrthogonalHeterogeneousLoadBalancer");

    std::vector<int> rankBuffer(cGeometry.size(), 0);
    std::vector<int> locBuffer(cGeometry.size(), 0);
    std::vector<std::uint8_t> platformBuffer(cGeometry.size(), static_cast<std::uint8_t>(Platform::CPU_SISD));

    const auto nRank = singleton::mpi().getSize();
    const auto iRank = singleton::mpi().getRank();

    auto& cuboids = cGeometry.cuboids();
    // Prioritize assignment of largest cuboids
    std::sort(cuboids.begin(), cuboids.end(),
              [&](const auto& lhs, const auto& rhs) {
                return lhs.getLatticeVolume() > rhs.getLatticeVolume();
              });

    #if defined(PARALLEL_MODE_MPI) || defined(PLATFORM_GPU_CUDA)
    int localPreferredPlatform = static_cast<int>(Platform::CPU_SIMD);
    #endif
    #ifdef PLATFORM_GPU_CUDA
    if (gpu::cuda::device::getCount() > 0) {
      localPreferredPlatform = static_cast<int>(Platform::GPU_CUDA);
    }
    #endif
    std::vector<int> preferredPlatform(cuboids.size(), -1);
    #ifdef PARALLEL_MODE_MPI
    singleton::mpi().gather(&localPreferredPlatform,  1,
                            preferredPlatform.data(), 1);
    #endif

    // Distribute cuboids to ranks on rank 0
    if (iRank == 0) {
      std::set<int> cpuRanks;
      std::set<int> gpuRanks;
      for (int jRank=0; jRank < singleton::mpi().getSize(); ++jRank) {
        switch (static_cast<Platform>(preferredPlatform[jRank])) {
          case Platform::CPU_SIMD:
            cpuRanks.insert(jRank);
            break;
          case Platform::GPU_CUDA:
            gpuRanks.insert(jRank);
            break;
          default:
            break;
       }
      }

      // Total volume for tracking targeted GPU block fraction
      std::size_t totalVolume = std::accumulate(cuboids.begin(), cuboids.end(),
                                                std::size_t{0},
                                                [](const auto& currVolume, const auto& rhs) -> std::size_t {
                                                  return currVolume + rhs.getLatticeVolume();
                                                });

      std::map<int,int> nLoc;
      int jRank = 0;
      int iCuboid = 0;

      std::size_t globalAssignedVolume = 0;
      std::vector<std::size_t> localAssignedVolume(nRank, 0);
      // Prevent GPU assignment to CPU ranks
      for (int jRank : cpuRanks) {
        localAssignedVolume[jRank] = std::numeric_limits<std::size_t>::max();
      }

      // Assign largest cuboids to GPUs until desired fraction is reached
      do {
        jRank = std::distance(localAssignedVolume.begin(),
                              std::min_element(localAssignedVolume.begin(),
                                               localAssignedVolume.end()));
        rankBuffer[iCuboid] = jRank;
        platformBuffer[iCuboid] = static_cast<std::uint8_t>(Platform::GPU_CUDA);
        locBuffer[iCuboid]  = nLoc[jRank]++;
        localAssignedVolume[jRank] += cuboids[iCuboid].getLatticeVolume();
        globalAssignedVolume += cuboids[iCuboid].getLatticeVolume();
        clout << iCuboid << ", " << jRank << ", assignedVolumeFraction=" << (1.*globalAssignedVolume) / totalVolume << std::endl;
        iCuboid += 1;
      } while (globalAssignedVolume < largeBlockFraction*totalVolume);

      for (int iRank : gpuRanks) {
        clout << "assignedVolumeOfRank[" << iRank << "]=" << localAssignedVolume[iRank] / double(totalVolume) << std::endl;
      }

      // Prevent CPU assignment to GPU ranks
      for (int jRank : cpuRanks) {
        localAssignedVolume[jRank] = 0;
      }
      for (int jRank : gpuRanks) {
        localAssignedVolume[jRank] = std::numeric_limits<std::size_t>::max();
      }

      for (int jCuboid=iCuboid; jCuboid < cGeometry.size(); ++jCuboid) {
        jRank = std::distance(localAssignedVolume.begin(),
                              std::min_element(localAssignedVolume.begin(),
                                               localAssignedVolume.end()));
        rankBuffer[jCuboid] = jRank;
        platformBuffer[jCuboid] = static_cast<std::uint8_t>(Platform::CPU_SIMD);
        locBuffer[jCuboid]  = nLoc[jRank]++;
        localAssignedVolume[jRank] += cuboids[jCuboid].getLatticeVolume();
        clout << jCuboid << ", " << jRank << std::endl;
      }

      for (int iRank : cpuRanks) {
        clout << "assignedVolumeOfRank[" << iRank << "]=" << localAssignedVolume[iRank] / double(totalVolume) << std::endl;
      }
    }

    #ifdef PARALLEL_MODE_MPI
    // Broadcast assignments to all processes
    singleton::mpi().bCast(rankBuffer.data(), rankBuffer.size());
    singleton::mpi().bCast(locBuffer.data(), locBuffer.size());
    singleton::mpi().bCast(platformBuffer.data(), platformBuffer.size());
    #endif

    // Update internal LoadBalancer structure to match given assignment
    for (int iCuboid=0; iCuboid < cGeometry.size(); ++iCuboid) {
      this->_rank[iCuboid] = rankBuffer[iCuboid];
      this->_loc[iCuboid] = locBuffer[iCuboid];

      if (rankBuffer[iCuboid] == iRank) {
        this->_glob.resize(std::max(int(this->_glob.size()), this->_loc[iCuboid]+1));
        this->_platform.resize(this->_glob.size());

        this->_glob[this->_loc[iCuboid]] = iCuboid;
        this->_platform[this->_loc[iCuboid]] = static_cast<Platform>(platformBuffer[iCuboid]);

        this->_size = this->_glob.size();
      }
    }
  }

  Platform platform(int loc) const override {
    return _platform[loc];
  }

  void setPlatform(int loc, Platform platform) {
    _platform[loc] = platform;
  }

};

}

#endif
