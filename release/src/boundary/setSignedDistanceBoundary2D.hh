/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert
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

//This file contains the Signed Distance Function Boundary
#ifndef SET_SIGNED_DISTANCE_BOUNDARY_2D_HH
#define SET_SIGNED_DISTANCE_BOUNDARY_2D_HH

#include "setSignedDistanceBoundary2D.h"

namespace olb {
///Initialising the setSignedDistanceBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSignedDistanceBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material)
{
  setSignedDistanceBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setSignedDistanceBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSignedDistanceBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  int _overlap = 1;
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setSignedDistanceBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc),indicator->getBlockIndicatorF(iCloc));
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  auto& communicator = sLattice.getCommunicator(stage::IterativePostProcess());
  communicator.template requestField<descriptors::PSI>();

  SuperIndicatorBoundaryNeighbor<T,DESCRIPTOR::d> neighborIndicator(std::forward<decltype(indicator)>(indicator), _overlap);
  communicator.requestOverlap(_overlap, neighborIndicator);
  communicator.exchangeRequests();
}

/// Set Signed Distance Function boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setSignedDistanceBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator)
{
  using namespace boundaryhelper;
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = 1;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      T discreteNormalSum=0;
      for (int iD=0; iD<DESCRIPTOR::d; iD++) {
        discreteNormalSum += abs(discreteNormal[iD+1]);
      }
      if (discreteNormalSum == 0) {
        block.addPostProcessor(typeid(stage::IterativePostProcess),{iX,iY},meta::id<normGradPsi>{});
      } else {
        block.addPostProcessor(typeid(stage::IterativePostProcess),{iX,iY},
          promisePostProcessorForNormal<T,DESCRIPTOR,normGradPsiBoundary2D>(
          Vector<int,2>(discreteNormal.data() + 1)));
      }
    }
  });
}

}//namespace olb

#endif
