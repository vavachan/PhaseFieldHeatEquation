/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 František Prinz, Nicolas Hafen, Mathias J. Krause
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

#ifndef LATTICE_SCHILLER_NAUMANN_DRAG_FORCE_HH
#define LATTICE_SCHILLER_NAUMANN_DRAG_FORCE_HH

namespace olb {



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeSchillerNaumannDragForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeSchillerNaumannDragForce(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  PhysR<T,DESCRIPTOR::d>cellMin, PhysR<T,DESCRIPTOR::d> cellMax,
  Vector<bool,DESCRIPTOR::d> periodic,
  std::size_t iP0 )
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice, converter,
                                   (DESCRIPTOR::d)*(particleSystem.size()-iP0)),
    _blockGeometry(blockGeometry), _blockLattice(blockLattice),
    _particleSystem(particleSystem),
    _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic), _iP0(iP0)
{
  this->getName() = "physSchillerNaumannDragForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeSchillerNaumannDragForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned D = DESCRIPTOR::d;
  const int serialSize = D;

  using namespace descriptors;
    static_assert(PARTICLETYPE::template providesNested<MOBILITY, FLUIDVEL>(), "Field MOBILITY:FLUIDVEL has to be provided");

  //Retrieve particle quantities
  Vector<T,D> position = particle.template getField<GENERAL,POSITION>();
  Vector<T,D> velocity = particle.template getField<MOBILITY,VELOCITY>();
  T radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
  T mass = particle.template getField<PHYSPROPERTIES,MASS>();
  T* positionArray = position.data();

  //TODO: check, whether creation can be avoided by using a second _blockF list
  const auto& cuboid = _blockGeometry.getCuboid();
  BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> blockInterpPhysVelF(
    _blockLattice, this->_converter, cuboid);

  //Calculate general constants
  T part_density = mass/ (4./3.*M_PI*radius*radius*radius);

  T Coeff = 18* this->_converter.getPhysViscosity()*this->_converter.getPhysDensity()/(part_density*2*radius*2*radius);

  //Calculate particle coefficiants


  T fluidVelArray[D] = {0.};

  //Check whether inside cuboid (when not parallelized)
  bool inside = true;
  if constexpr ( !particles::access::providesParallelization<PARTICLETYPE>() ){
    inside = cuboid.isInside(position);
  }
  if (inside){

    //Retrieve interpolated velocity at position
    blockInterpPhysVelF(fluidVelArray, positionArray);
 particle.template setField<MOBILITY, FLUIDVEL>(fluidVelArray);
 T relVel = std::sqrt((velocity[0]-fluidVelArray[0])*(velocity[0]-fluidVelArray[0])+(velocity[1]-fluidVelArray[1])*(velocity[1]-fluidVelArray[1])+(velocity[2]-fluidVelArray[2])*(velocity[2]-fluidVelArray[2]));

 T Re_p = relVel*2*radius/this->_converter.getPhysViscosity();

 T C_d;
 if (Re_p < 0.0000000001)
   C_d=0;
 else
 {
 if (Re_p <1000)
  C_d = 24./Re_p*(1+0.15*olb::util::pow(Re_p, 0.687));
 else
  C_d = 0.44;
 }

    //Calculate SchillerNaumann force
    T tmpForce[D] = {0.};
    for (int iDim = 0; iDim < PARTICLETYPE::d; ++iDim) {
      tmpForce[iDim] = C_d*Coeff*(fluidVelArray[iDim] - velocity[iDim])*mass/24.*Re_p;
    }

    particle.template setField<FORCING,FORCE> (tmpForce);
  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeSchillerNaumannDragForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  using namespace descriptors;
  // iterate over all particles in _indicator //TODO: add periodic treatment analogous to momentum exchange
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); iP++) {
    auto particle = _particleSystem.get(iP);
    evaluate(output, particle, iP);
  }
  return true;
}



}
#endif
