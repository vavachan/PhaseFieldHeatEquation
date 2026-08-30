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


#ifndef LATTICE_SPHEROID_Lift_FORCE_HH
#define LATTICE_SPHEROID_Lift_FORCE_HH

#include "../../particles/functions/eulerRotation.h"

namespace olb {


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
 BlockLatticeSpheroidLiftForce<T, DESCRIPTOR, PARTICLETYPE>:: BlockLatticeSpheroidLiftForce(
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
  this->getName() = "physStokesSpheroidLiftForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeSpheroidLiftForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
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

  //Calculate particle coefficiants


  T fluidVelArray[D] = {0.};

  //Check whether inside cuboid (when not parallelized)
  bool inside = true;
  if constexpr ( !particles::access::providesParallelization<PARTICLETYPE>() ){
    inside = cuboid.checkPoint(position);
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
Khathatvelo = eler::matrixVectorMultiply(Khathat, fvel - velocity);//in Schachar Beerman is K

 /*   / 0 1 2 \
 *   | 3 4 5 |
 *   \ 6 7 8 /
 */

Vector<T,3> F12, F13, F23, F21, F31, F32;
Vector<T,3> force = Vector<T,3> (0.,0.,0.);
Vector<T,9> B12, B13, B23, B21, B31, B32;
Vector<T,9> invB12, invB13, invB23, invB21, invB31, invB32;
B12 = Vector<T,9> (1.,0.,0.,0.,1.,0.,0.,0.,1.);
invB12 = B12;
B13 = Vector<T,9> (1.,0.,0.,0.,0.,-1.,0.,1.,0.);
invB13 = Vector<T,9> (1.,0.,0.,0.,0.,1.,0.,-1.,0.);
B21 = Vector<T,9> (0.,0.,1.,1.,0.,0.,0.,-1.,0.);
invB21 = Vector<T,9>(0.,1.,0.,0.,0.,-1.,1.,0.,0.);
B23 = Vector<T,9>(0.,0.,1.,1.,0.,0.,0.,1.,0.);
invB23 = Vector<T,9>(0.,1.,0.,0.,0.,1.,1.,0.,0.);
B31 = Vector<T,9>(0.,1.,0.,0.,0.,1.,1.,0.,0.);
invB31 = Vector<T,9>(0.,0.,1.,1.,0.,0.,0.,1.,0.);
B32 = Vector<T,9>(0.,0.,-1.,0.,1.,0.,1.,0.,0.);
invB32 = Vector<T,9>(0.,0.,1.,0.,1.,0.,-1.,0.,0.);

Vector<T,9> L (0.0501, 0.0329, 0.0, 0.0182, 0.0173, 0.00, 0.00, 0.00, 0.0373);


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

using namespace eler;
using namespace util;

T lim = std::numeric_limits<T>::min();

T Const = M_PI*M_PI*particle.template getField<PHYSPROPERTIES,RADIUS>()*particle.template getField<PHYSPROPERTIES,RADIUS>()*this->_converter.getPhysViscosity()*this->_converter.getPhysDensity()/sqrt(this->_converter.getPhysViscosity());//this->_converter.getPhysViscosity()*this->_converter.getPhysDensity()*M_PI*particle.template getField<PHYSPROPERTIES,RADIUS>();

if(std::fabs(gradvel[1]) < lim )
F12 = Vector<T,3> (0.,0.,0.);
else
F12 = Const*gradvel[1]/sqrt(std::fabs(gradvel[1]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B12,matrixMultiply(L,matrixMultiply(invB12,Khathat)))),fvel - velocity);
if(std::fabs(gradvel[2]) < lim )
F13 = Vector<T,3> (0.,0.,0.);
else
F13 = Const*gradvel[2]/sqrt(std::fabs(gradvel[2]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B13,matrixMultiply(L,matrixMultiply(invB13,Khathat)))),fvel - velocity) ;
if(std::fabs(gradvel[3]) < lim )
F21 = Vector<T,3> (0.,0.,0.);
else
F21 = Const*gradvel[3]/sqrt(std::fabs(gradvel[3]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B21,matrixMultiply(L,matrixMultiply(invB21,Khathat)))),fvel - velocity) ;
if(std::fabs(gradvel[5]) < lim )
F23 = Vector<T,3> (0.,0.,0.);
else
F23 = Const*gradvel[5]/sqrt(std::fabs(gradvel[5]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B23,matrixMultiply(L,matrixMultiply(invB23,Khathat)))),fvel - velocity) ;
if(std::fabs(gradvel[6]) < lim )
F31 = Vector<T,3> (0.,0.,0.);
else
F31 = Const*gradvel[6]/sqrt(std::fabs(gradvel[6]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B31,matrixMultiply(L,matrixMultiply(invB31,Khathat)))),fvel - velocity) ;
if(std::fabs(gradvel[7]) < lim )
F32 = Vector<T,3> (0.,0.,0.);
else
F32 = Const*gradvel[7]/sqrt(std::fabs(gradvel[7]))*matrixVectorMultiply(matrixMultiply(Khathat,matrixMultiply(B32,matrixMultiply(L,matrixMultiply(invB32,Khathat)))),fvel - velocity) ;


force = F12 + F13 + F21 + F23 + F31 + F32;

    particle.template setField<FORCING,FORCE> ( particle.template getField<FORCING,FORCE>()+force);

  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeSpheroidLiftForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
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
