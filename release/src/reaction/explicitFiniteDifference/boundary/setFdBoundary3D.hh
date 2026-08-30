/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, Davide Dapelo
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

//This file contains the Interpolated Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FD_BOUNDARY_DEV03_HH
#define SET_FD_BOUNDARY_DEV03_HH

namespace olb {

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setFdBoundary3D function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                SuperGeometry<T,3>& superGeometry, int material)
{
  setFdBoundary3D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setFdBoundary3D function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setFdBoundary3D");
  int _overlap = indicator->getSuperGeometry().getOverlap();
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setFdBoundary3D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>(sLattice.getBlock(iC), indicator->getBlockIndicatorF(iC),includeOuterCells);
    /// Adds needed Cells to the Communicator _commBC in SuperLattice
    //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary3D.h/hh
    addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
  }
}

////////// BlockLattice Domain  /////////////////////////////////////////

//set FdBoundary on indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary3D(BlockLattice<T,DESCRIPTOR>& block,
                                BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      std::vector<int> discreteNormal(4,0);
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY, iZ);
      block.get(iX, iY, iZ).template setField<descriptors::NORMAL_X>(-discreteNormal[1]);
      block.get(iX, iY, iZ).template setField<descriptors::NORMAL_Y>(-discreteNormal[2]);
      block.get(iX, iY, iZ).template setField<descriptors::NORMAL_Z>(-discreteNormal[3]);
    }
  });
  block.addPostProcessor ( typeid(stage::PostStream), indicator,
                           meta::id<FdBoundaryPostProcessor3D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>> {} );
}

}
#endif
