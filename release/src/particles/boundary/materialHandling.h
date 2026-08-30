/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Nicolas Hafen, Mathias J. Krause
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



#ifndef PARTICLE_MATERIAL_HANDLING_H
#define PARTICLE_MATERIAL_HANDLING_H


namespace olb {

namespace particles {

namespace boundaries {



//Return whether location is in direct vicinity of specified materials
template<typename T,unsigned D>
bool materialVicinity( SuperIndicatorMaterial<T,D>& materialIndicator, LatticeR<D+1>& latticeR )
{
  static_assert((D==2 || D==3), "Dimension unknown!");
  static_assert(D!=2, "2D VERSION NOT IMPLEMENTED YET");
  if constexpr(D==2){
    //TODO: implement 2D version
  } else if constexpr(D==3){
    return (materialIndicator(latticeR[0],latticeR[1]+1,latticeR[2],  latticeR[3] ) ||
            materialIndicator(latticeR[0],latticeR[1]+1,latticeR[2],  latticeR[3]+1 ) ||
            materialIndicator(latticeR[0],latticeR[1]+1,latticeR[2]+1,latticeR[3] ) ||
            materialIndicator(latticeR[0],latticeR[1]+1,latticeR[2]+1,latticeR[3]+1 ) ||
            materialIndicator(latticeR[0],latticeR[1],latticeR[2],    latticeR[3] ) ||
            materialIndicator(latticeR[0],latticeR[1],latticeR[2],    latticeR[3]+1 ) ||
            materialIndicator(latticeR[0],latticeR[1],latticeR[2]+1,  latticeR[3] ) ||
            materialIndicator(latticeR[0],latticeR[1],latticeR[2]+1,  latticeR[3]+1 ));
  }
}


template<typename T, typename PARTICLETYPE>
bool checkMaterialVicinity(
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator,
  Particle<T,PARTICLETYPE>& particle )
{
  constexpr unsigned D = PARTICLETYPE::d;
  //Retrieve particle position
  auto position = access::getPosition( particle );
  //Retrieve super geometry
  auto& sGeometry = materialIndicator.getSuperGeometry();
  //Get lattice position
  auto latticeR = sGeometry.getCuboidDecomposition().getFloorLatticeR(position);
  bool foundPos = latticeR.operator bool();
  if (!foundPos) {
    std::cerr << "Position (" << position << ") not found!" << std::endl;
  } else {
    //Check whether responsible iC
    bool isLocal = sGeometry.getLoadBalancer().isLocal((*latticeR)[0]);
    if(isLocal){
      //Check if in vicinity of specified material numbers
      bool vicinity = materialVicinity<T,D>(materialIndicator, *latticeR);
      return vicinity;
    }
  }
  return false;
}


} //namespace boundaries

} //namespace particles

} //namespace olb


#endif
