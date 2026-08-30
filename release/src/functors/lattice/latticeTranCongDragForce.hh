/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 František Prinz, Nicolas Hafen,
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

#ifndef LATTICE_TRAN_CONG_DRAG_FORCE_HH
#define LATTICE_TRAN_CONG_DRAG_FORCE_HH

namespace olb {


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeTranCongDragForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeTranCongDragForce(
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
    _cellMin(cellMin), _cellMax(cellMax),_periodic(periodic), _iP0(iP0)
{
  this->getName() = "physTranCongDragForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeTranCongDragForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned Dd = DESCRIPTOR::d;
  const int serialSize = Dd;
  using namespace descriptors;
    static_assert(PARTICLETYPE::template providesNested<MOBILITY, FLUIDVEL>(), "Field MOBILITY:FLUIDVEL has to be provided");

  //Retrieve particle quantities
  Vector<T,Dd> position = particle.template getField<GENERAL,POSITION>();
  Vector<T,Dd> velocity = particle.template getField<MOBILITY,VELOCITY>();
  T radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
  T mass = particle.template getField<PHYSPROPERTIES,MASS>();
  T* positionArray = position.data();

  //TODO: check, whether creation can be avoided by using a second _blockF list
  const auto& cuboid = _blockGeometry.getCuboid();
  BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> blockInterpPhysVelF(
    _blockLattice, this->_converter, cuboid);

  //Calculate general constants
  T part_density = particle.template getField<PHYSPROPERTIES,DENSITY>();;

  T fluidVelArray[Dd] = {0.};

  //Check whether inside cuboid (when not parallelized)
  bool inside = true;
  if constexpr ( !particles::access::providesParallelization<PARTICLETYPE>() ){
    inside = cuboid.isInside(position);
  }
  if (inside){

    //Retrieve interpolated velocity at position
    blockInterpPhysVelF(fluidVelArray, positionArray);
 particle.template setField<MOBILITY, FLUIDVEL>(fluidVelArray);
 T relVel = util::sqrt((velocity[0]-fluidVelArray[0])*(velocity[0]-fluidVelArray[0])+(velocity[1]-fluidVelArray[1])*(velocity[1]-fluidVelArray[1])+(velocity[2]-fluidVelArray[2])*(velocity[2]-fluidVelArray[2]));

  Vector<T,Dd> direction  ((velocity[0]-fluidVelArray[0])/relVel,(velocity[1]-fluidVelArray[1])/relVel,(velocity[2]-fluidVelArray[2])/relVel);
 T Re_p = relVel*2*radius/this->_converter.getPhysViscosity();//needed
 T C_d;
  if (Re_p < 0.000001)
   C_d = 0.;
  else
   C_d = 24./Re_p*particle.template getField<NUMERICPROPERTIES, DA_DN>()*(1.+0.15/util::sqrt(particle.template getField<NUMERICPROPERTIES,CIRCULARITY>())*util::pow(particle.template getField<NUMERICPROPERTIES, DA_DN>()*Re_p,0.687))+
 (0.42*particle.template getField<NUMERICPROPERTIES, DA_DN>()*particle.template getField<NUMERICPROPERTIES, DA_DN>())/(util::sqrt(particle.template getField<NUMERICPROPERTIES, CIRCULARITY>())*(1+4.25*10000*util::pow(particle.template getField<NUMERICPROPERTIES,DA_DN>()*Re_p,-1.16)));


    //Calculate TranCong force
    T tmpForce[Dd] = {0.};
    for (int iDim = 0; iDim < PARTICLETYPE::d; ++iDim) {
      tmpForce[iDim] =-C_d*particle.template getField<NUMERICPROPERTIES, PROJECTED_SURFACE>()*this->_converter.getPhysDensity()*relVel*relVel*direction[iDim];

    }
    particle.template setField<FORCING,FORCE> (tmpForce);
  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeTranCongDragForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  using namespace descriptors;
  using namespace particles;
  // iterate over all particles in _indicator //TODO: add periodic treatment analogous to momentum exchange

      forParticlesInParticleSystem<T,PARTICLETYPE,conditions::valid_particles>( _particleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){

      evaluate(output, particle, particle.getId());


    });
  return true;
}



}
#endif
