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

#ifndef RANDOM_LOAD_BALANCER_H
#define RANDOM_LOAD_BALANCER_H

#include "communication/loadBalancer.h"

#include <random>
#include <algorithm>

namespace olb {

/// Basic Random Load Balancer
/**
 * Randomly assigns cuboids to ranks while approximately
 * preserving equal number of blocks (not cells) per rank.
 **/
template<typename T>
struct RandomLoadBalancer final : public LoadBalancer<T> {
  RandomLoadBalancer(int nCuboid, int nRank, int iRank):
    LoadBalancer<T>(0)
  {
    std::vector<int> rankBuffer(nCuboid, 0);
    std::vector<int> locBuffer(nCuboid, 0);

    // Randomly distribute cuboids to ranks on rank 0
    if (iRank == 0) {
      std::vector<int> cuboids(nCuboid, 0);
      std::iota(cuboids.begin(),
                cuboids.end(),
                0);
      std::random_device randomDevice;
      std::default_random_engine randomEngine(randomDevice());
      std::shuffle(cuboids.begin(), cuboids.end(), randomEngine);

      std::map<int,int> nLoc;
      int jRank = 0;
      for (int iCuboid : cuboids) {
        if (jRank >= nRank) {
          jRank -= nRank;
        }
        rankBuffer[iCuboid] = jRank;
        locBuffer[iCuboid]  = nLoc[jRank]++;
        jRank += 1;
      }
    }

    #ifdef PARALLEL_MODE_MPI
    // Broadcast assignments to all processes
    singleton::mpi().bCast(rankBuffer.data(), rankBuffer.size());
    singleton::mpi().bCast(locBuffer.data(), locBuffer.size());
    #endif

    // Update internal LoadBalancer structure to match given assignment
    for (int iCuboid=0; iCuboid < nCuboid; ++iCuboid) {
      this->_rank[iCuboid] = rankBuffer[iCuboid];
      this->_loc[iCuboid] = locBuffer[iCuboid];
      if (rankBuffer[iCuboid] == singleton::mpi().getRank()) {
        this->_glob.resize(std::max(int{this->_glob.size()}, this->_loc[iCuboid]+1));
        this->_glob[this->_loc[iCuboid]] = iCuboid;
        this->_size = this->_glob.size();
      }
    }
  }

  RandomLoadBalancer(CuboidDecomposition<T,3>& cGeometry):
    RandomLoadBalancer(cGeometry.size(),
                       singleton::mpi().getSize(),
                       singleton::mpi().getRank())
  { }

};

}

#endif
