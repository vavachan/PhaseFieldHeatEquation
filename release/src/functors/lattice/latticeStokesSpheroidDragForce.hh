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

#ifndef LATTICE_STOKES_SPHEROID_DRAG_FORCE_HH
#define LATTICE_STOKES_SPHEROID_DRAG_FORCE_HH

#include "../../particles/functions/eulerRotation.h"

namespace olb {


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeStokesSpheroidDragForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeStokesSpheroidDragForce(
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
  this->getName() = "physStokesSpheroidDragForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeStokesSpheroidDragForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned D = DESCRIPTOR::d;

  using namespace descriptors;
  using namespace eler;
  static_assert(PARTICLETYPE::template providesNested<MOBILITY, FLUIDVEL>(), "Field MOBILITY:FLUIDVEL has to be provided");
  //Retrieve particle quantities
  Vector<T,D> position = particle.template getField<GENERAL,POSITION>();
  Vector<T,D> velocity = particle.template getField<MOBILITY,VELOCITY>();
  T* positionArray = position.data();

  //TODO: check, whether creation can be avoided by using a second _blockF list
  const auto& cuboid = _blockGeometry.getCuboid();
  BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> blockInterpPhysVelF(
    _blockLattice, this->_converter, cuboid);

  //Calculate general constants
  T Const = this->_converter.getPhysViscosity()*this->_converter.getPhysDensity()*M_PI*particle.template getField<PHYSPROPERTIES,RADIUS>();
  //Calculate particle coefficiants

  T fluidVelArray[D] = {0.};

  //Check whether inside cuboid (when not parallelized)
  bool inside = true;
  if constexpr ( !particles::access::providesParallelization<PARTICLETYPE>() ){
    inside = cuboid.isInside(position);
  }
  if (inside){
   //Force computation
    //Retrieve interpolated velocity at position
    blockInterpPhysVelF(fluidVelArray, positionArray);
 particle.template setField<MOBILITY, FLUIDVEL>(fluidVelArray);

 Vector<T,3> fvel = particle.template getField<MOBILITY, FLUIDVEL>();

 Vector<T,9> rotmat = particle.template getField<NUMERICPROPERTIES,ROT_MATRIX>();

 Vector<T,9> Khathat (0.);
Khathat = eler::matrixMultiply(eler::inverseTransformMatrix(rotmat),
eler::matrixMultiply(particle.template getField<NUMERICPROPERTIES,TRANSLATIONAL_DYADIC>(), rotmat));
Vector<T,3> Khathatvelo (0.);
Khathatvelo = eler::matrixVectorMultiply(Khathat, fvel - velocity);
Vector<T,3> force = Const*Khathatvelo;

    T dynvisc = this->_converter.getPhysViscosity()*this->_converter.getPhysDensity();


    //Torque computation
    Vector <T,3> torq (0.);
    Vector <T,3> posVector (0.);
     posVector    = positionArray;
     Vector<T,3> x (1.,0.,0.);
     Vector<T,3> y (0.,1.,0.);
     Vector<T,3> z (0.,0.,1.);
     Vector<T,3> gradx (0.);
     Vector<T,3> grady (0.);
     Vector<T,3> gradz (0.);


    veloGradientInterpFD<T,DESCRIPTOR>(gradx,blockInterpPhysVelF, x,posVector, this->_converter.getPhysDeltaX()/100.);
    veloGradientInterpFD<T,DESCRIPTOR>(grady,blockInterpPhysVelF, y,posVector, this->_converter.getPhysDeltaX()/100.);
    veloGradientInterpFD<T,DESCRIPTOR>(gradz,blockInterpPhysVelF, z,posVector, this->_converter.getPhysDeltaX()/100.);


    /* / 0 1 2 \  vx/dx vx/dy vx/dz
       | 3 4 5 |  vy/dx vy/dy vy/dz
       \ 6 7 8 /  vz/dx vz/dy vz/dz  */
    Vector <T,9> gradvel (0.);
    gradvel[0] = gradx[0];
    gradvel[3] = gradx[1];
    gradvel[6] = gradx[2];
    gradvel[1] = grady[0];
    gradvel[4] = grady[1];
    gradvel[7] = grady[2];
    gradvel[2] = gradz[0];
    gradvel[5] = gradz[1];
    gradvel[8] = gradz[2];

    Vector<T,3> abg0 =particle.template getField<NUMERICPROPERTIES,ABG0> ();
    T beta = particle.template getField<NUMERICPROPERTIES,BETA> ();
    Vector<T,3> ang_vel = particle.template getField<NUMERICPROPERTIES,ANG_VELOCITY> ();
    T rad = particle.template getField<PHYSPROPERTIES,RADIUS> ();

    gradvel = matrixMultiply(rotmat,matrixMultiply(gradvel, inverseTransformMatrix(rotmat)));

    torq[0] = (16*M_PI*dynvisc*util::pow(rad,3)*beta)/
    (3*(abg0[1]+util::pow(beta,2)*abg0[2]))*
    ((1- util::pow(beta,2))*(0.5*(gradvel[7]+gradvel[5]))+ (1+ util::pow(beta,2))*
    (0.5*(gradvel[7]-gradvel[5])- ang_vel[0]));

    torq[1] = (16*M_PI*dynvisc*util::pow(rad,3)*beta)/
    (3*(abg0[0]+util::pow(beta,2)*abg0[2]))*
    ((util::pow(beta,2)-1)*(0.5*(gradvel[2]+gradvel[6]))+ (1+ util::pow(beta,2))*
    (0.5*(gradvel[2]-gradvel[6])-ang_vel[1]));

    torq[2] = (32*M_PI*dynvisc*util::pow(rad,3)*beta)/
    (3*(abg0[0]+abg0[1]))*
    (0.5*(gradvel[3]-gradvel[1])-ang_vel[2]);

    particle.template setField<FORCING,TORQUE> ( particle.template getField<FORCING,TORQUE>()+torq);
    particle.template setField<FORCING,FORCE> ( particle.template getField<FORCING,FORCE>()+force);



  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeStokesSpheroidDragForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  using namespace descriptors;
  using namespace particles;
  // iterate over all particles in _indicator //TODO: add periodic treatment analogous to momentum exchange
    forParticlesInParticleSystem<T,PARTICLETYPE,conditions::valid_particles>( _particleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){

      evaluate(output, particle, particle.getId());


    });
   // forParticlesInSuperParticleSystem(

  return true;
}



}
#endif
