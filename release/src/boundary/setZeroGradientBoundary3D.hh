/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Fedor Bukreev, Adrian Kummerländer
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

//This file contains the ZeroGradientBoundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_ADVECTION_DIFFUSION_ZERO_GRADIENT_BOUNDARY_3D_HH
#define SET_ADVECTION_DIFFUSION_ZERO_GRADIENT_BOUNDARY_3D_HH

#include "setZeroGradientBoundary3D.h"

namespace olb {

///Initialising the setZeroGradientBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroGradientBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material)
{
  setZeroGradientBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setZeroGradientBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroGradientBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setZeroGradientBoundary");

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setZeroGradientBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc));
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  int _overlap = 2;
  auto& communicator = sLattice.getCommunicator(stage::PostStream());
  communicator.template requestField<descriptors::POPULATION>();

  SuperIndicatorBoundaryNeighbor<T,DESCRIPTOR::d> neighborIndicator(std::forward<decltype(indicator)>(indicator), _overlap);
  communicator.requestOverlap(_overlap, neighborIndicator);
  communicator.exchangeRequests();
}


/// Set ZeroGradientBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroGradientBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = 2;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
        for(int iPop = 0; iPop < DESCRIPTOR::q; iPop++){
          auto c = descriptors::c<DESCRIPTOR>(iPop);
          T k = c[0]*(-discreteNormal[1]) + c[1]*(-discreteNormal[2]) + c[2]*(-discreteNormal[3]);
          if(k > T{0.} && blockGeometryStructure.getMaterial(iX+c[0], iY+c[1], iZ+c[2]) == 1 ){
            _block.get(iX+c[0], iY+c[1], iZ+c[2]).template getFieldPointer<descriptors::NEIGHBOR>()[0] = T{1.};
          }
          if(k > T{0.} && blockGeometryStructure.getMaterial(iX+2*c[0], iY+2*c[1], iZ+2*c[2]) == 1 ){
            _block.get(iX+2*c[0], iY+2*c[1], iZ+2*c[2]).template getFieldPointer<descriptors::NEIGHBOR>()[0] = T{1.};
          }
        }
        _block.addPostProcessor(
            typeid(stage::PostCollide), {iX,iY,iZ},
            meta::id<ZeroGradientLatticePostProcessor3D<T,DESCRIPTOR>>{}
            );
        _block.template defineDynamics<NoCollideDynamics>({iX, iY, iZ});
      } else {
        _block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
    }
  });
}


}//namespace olb
#endif
