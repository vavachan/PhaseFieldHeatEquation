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
#ifndef SET_FD_BOUNDARY_2D_DEV03_HH
#define SET_FD_BOUNDARY_2D_DEV03_HH

namespace olb {

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setFdBoundary2D function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary2D(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material)
{
  setFdBoundary2D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setFdBoundary2D function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary2D(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = indicator->getSuperGeometry().getOverlap();
  OstreamManager clout(std::cout, "setOnBCInterpolatedBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFdBoundary2D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>(sLattice.getBlock(iCloc),
        indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary2D.h/hh
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);


}
////////// BlockLattice Domain  /////////////////////////////////////////

/// Set Interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdBoundary2D(BlockLattice<T,DESCRIPTOR>& block,
                                BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      block.get(iX, iY).template setField<descriptors::NORMAL_X>(-discreteNormal[1]);
      block.get(iX, iY).template setField<descriptors::NORMAL_Y>(-discreteNormal[2]);
    }
  });
  block.addPostProcessor ( typeid(stage::PostStream), indicator,
                           meta::id<FdBoundaryPostProcessor2D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMETERS,FIELD,SOURCE>>{} );
}


}//namespace olb

#endif
