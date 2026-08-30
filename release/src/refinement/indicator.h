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

#ifndef REFINEMENT_INDICATOR_H
#define REFINEMENT_INDICATOR_H

#include "loadBalancer.h"
#include "geometry/cuboidDecomposition.h"

namespace olb {

/// Indicator to identify a cell layer at given distance from the cuboid decomposition's frontier
/**
 * Used to identify regions where coupling between refined resolutions should take place.
 **/
template <typename T, typename DESCRIPTOR>
struct SuperIndicatorDomainFrontierDistanceF : public SuperIndicatorF<T,DESCRIPTOR::d> {
  SuperIndicatorDomainFrontierDistanceF(CuboidDecomposition<T,DESCRIPTOR::d>& cDecomposition,
                                        SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                                        int distance)
  : SuperIndicatorF<T,DESCRIPTOR::d>(sGeometry)
  {
    auto& load = sGeometry.getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockIndicatorFfromCallableF<T,DESCRIPTOR::d>(sGeometry.getBlock(iC),
                                                          [&,iC,distance](LatticeR<DESCRIPTOR::d> latticeR) -> bool {
          const auto& block = sGeometry.getBlock(iC);
          const T deltaX = cDecomposition.getMotherCuboid().getDeltaR();
          if (block.getNeighborhoodRadius(latticeR) < 1) {
            return false;
          }
          auto physR = cDecomposition.get(load.glob(iC)).getPhysR(latticeR);
          const T outsideTestDistance = 1 + distance + 0.75;
          const T  insideTestDistance = 1 + distance;
          bool outsideTestOnce = false;
          bool insideTestAlways = true;
          for (int iDir=0; iDir < fields::refinement::CONTEXT_NEIGHBORS::count<DESCRIPTOR>(); ++iDir) {
            const auto c_i = fields::refinement::CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iDir);
            outsideTestOnce  |= !cDecomposition.isInside(physR + c_i * outsideTestDistance * deltaX);
            insideTestAlways &=  cDecomposition.isInside(physR + c_i *  insideTestDistance * deltaX);
          }
          return outsideTestOnce && insideTestAlways;
        })
      );
    }
  }

  using SuperIndicatorF<T,DESCRIPTOR::d>::operator();

  bool operator() (bool output[], const int input[]) override {
    auto& load = this->_superGeometry.getLoadBalancer();
    if (load.isLocal(input[0])) {
      return this->getBlockF(load.loc(input[0]))(output,&input[1]);
    } else {
      return false;
    }
  }

};

}

#endif
