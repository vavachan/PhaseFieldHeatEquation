/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

//This file contains the LocalConvection Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_LOCAL_CONVECTION_BOUNDARY_3D_HH
#define SET_LOCAL_CONVECTION_BOUNDARY_3D_HH

#include "setLocalConvectionBoundary3D.h"

namespace olb {

///Initialising the setLocalConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material,
                                T* uAv)
{
  setLocalConvectionBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), uAv);
}

///Initialising the setLocalConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T* uAv)
{
  int _overlap = 0;  // TODO: Is intended?
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setLocalConvectionBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc),
        uAv);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  // TODO: Is communication really needed for this BC?
  auto& communicator = sLattice.getCommunicator(stage::PostStream());
  communicator.template requestField<descriptors::POPULATION>();

  SuperIndicatorBoundaryNeighbor<T,DESCRIPTOR::d> neighborIndicator(std::forward<decltype(indicator)>(indicator), _overlap);
  communicator.requestOverlap(_overlap, neighborIndicator);
  communicator.exchangeRequests();
}


/// Add convection boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setLocalConvectionBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, T* uAv)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

      if (discreteNormal[0] == 0) {//set postProcessors for indicated boundary cells
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          postProcessor = nullptr;
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          postProcessor = nullptr;
        }
        if (postProcessor) {
          _block.addPostProcessor(*postProcessor);
        }
      }
    }
  });
}

}//namespace olb
#endif

