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

#ifndef REFINEMENT_LOAD_BALANCER_H
#define REFINEMENT_LOAD_BALANCER_H

#include "communication/loadBalancer.h"
#include "geometry/cuboidDecomposition.h"

namespace olb {

/// Load balancer for refined lattice hierarchies
/**
 * Ensures that fine cuboids are assigned to the same
 * rank and platform as their corresponding coarse parents.
 *
 * Coarse parent load balancer must be aware of additional
 * weight due to refinement for a good balancing.
 **/
template<typename T, unsigned D>
class RefinedLoadBalancer final : public LoadBalancer<T> {
private:
  std::map<int,int> _coarseLoc;

public:
  RefinedLoadBalancer(CuboidDecomposition<T,D>& coarseGeometry,
                      LoadBalancer<T>& coarseBalancer,
                      CuboidDecomposition<T,D>& fineGeometry)
    : LoadBalancer<T>(0)
  {
    OstreamManager clout(std::cout, "RefinedLoadBalancer");
    clout.setMultiOutput(true);
    std::vector<int> rankBuffer(fineGeometry.size(), 0);
    std::vector<int> locBuffer(fineGeometry.size(), 0);
    std::map<int,int> nLoc;
    for (unsigned iC=0; iC < fineGeometry.size(); ++iC) {
      const auto& fineCuboid = fineGeometry.get(iC);
      if (auto latticeR = coarseGeometry.getLatticeR(fineCuboid.getOrigin())) {
        const int coarseC = (*latticeR)[0];
        if (coarseBalancer.isLocal(coarseC)) {
          rankBuffer[iC] = coarseBalancer.rank(coarseC);
          locBuffer[iC] = nLoc[rankBuffer[iC]]++;
          this->_platform[iC] = coarseBalancer.platform(coarseBalancer.loc(coarseC));
          clout << "Matching fineC=" << iC << " to coarseC=" << coarseC << " on rank=" << rankBuffer[iC] << std::endl;
          _coarseLoc[locBuffer[iC]] = coarseBalancer.loc(coarseC);
        }
      }
    }

    #ifdef PARALLEL_MODE_MPI
    std::vector<int> globalRankBuffer(fineGeometry.size(), 0);
    std::vector<int> globalLocBuffer(fineGeometry.size(), 0);
    // Broadcast assignments to all processes
    singleton::mpi().allreduce(rankBuffer.data(), globalRankBuffer.data(), rankBuffer.size(), MPI_SUM);
    singleton::mpi().allreduce(locBuffer.data(), globalLocBuffer.data(), locBuffer.size(), MPI_SUM);
    #else
    auto& globalRankBuffer = rankBuffer;
    auto& globalLocBuffer = locBuffer;
    #endif

    // Update internal LoadBalancer structure to match given assignment
    for (unsigned iC=0; iC < fineGeometry.size(); ++iC) {
      this->_rank[iC] = globalRankBuffer[iC];
      this->_loc[iC] = globalLocBuffer[iC];
      if (globalRankBuffer[iC] == singleton::mpi().getRank()) {
        this->_glob.resize(std::max(static_cast<int>(this->_glob.size()), this->_loc[iC]+1));
        this->_glob[this->_loc[iC]] = iC;
        this->_size = this->_glob.size();
      }
    }
    clout.setMultiOutput(false);
  }

  /// Returns coarse loc for given fine loc
  int cloc(int fineLoc) const {
    return _coarseLoc.at(fineLoc);
  }

};

}

#endif
