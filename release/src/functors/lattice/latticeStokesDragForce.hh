/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef LATTICE_STOKES_DRAG_FORCE_HH
#define LATTICE_STOKES_DRAG_FORCE_HH

#include "latticeStokesDragForce.h"

namespace olb {


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, bool serialize>
BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE, serialize>::BlockLatticeStokesDragForce(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  PhysR<T,DESCRIPTOR::d>cellMin, PhysR<T,DESCRIPTOR::d> cellMax,
  Vector<bool,DESCRIPTOR::d> periodic,
  std::size_t iP0,
  const std::set<int>& ignoredMaterials,
  const F f)
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice, converter,
                                   (DESCRIPTOR::d)*(particleSystem.size()-iP0)),
    _blockGeometry(blockGeometry), _blockLattice(blockLattice),
    _particleSystem(particleSystem),
    _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic),
    _iP0(iP0), _ignoredMaterials(ignoredMaterials), _f(f)
{
  this->getName() = "physStokesDragForce";
  //Calculate precalculated constants
  _delTinv = 1./this->_converter.getPhysDeltaT();
  _C1 = 6. * M_PI
           * converter.getPhysViscosity()
           * converter.getPhysDensity()
           * converter.getConversionFactorTime();
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, bool serialize>
void BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE, serialize>::evaluate(
  T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned D = DESCRIPTOR::d;
  const int serialSize = D;

  using namespace descriptors;

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

  //Calculate particle coefficiants
  T c = _C1 * radius * 1./mass;
  T C2 = 1. / (1. + c);

  T fluidVelArray[D] = {0.};

  //Check whether inside cuboid (when not parallelized)
  bool inside = true;
  if constexpr ( !particles::access::providesParallelization<PARTICLETYPE>() ){
    inside = cuboid.isInside(position);
  }
  if (inside){

    //Retrieve interpolated velocity at position
    blockInterpPhysVelF(fluidVelArray, positionArray);

    //Record interpolated fluid field, if MOBILITY,FLUIDVEL provided
    if constexpr( PARTICLETYPE::template providesNested<MOBILITY,FLUIDVEL>() ){
      particle.template setField<MOBILITY, FLUIDVEL>(fluidVelArray);
    }

    //Calculate stokes force
    Vector<T,D> tmpForce(0.);
    for (int iDim = 0; iDim < PARTICLETYPE::d; ++iDim) {
      tmpForce[iDim] = mass * _delTinv
        * ((c * fluidVelArray[iDim] + velocity[iDim]) * C2 - velocity[iDim]);
    }

    _f(particle, position, tmpForce, Vector<T,utilities::dimensions::convert<D>::rotation>(T{0}));

    //Add force and torque to output or apply directly
    //INFO: serialization only possible to provide analogy to momentumExchange!
    if constexpr (serialize){
      std::size_t iPeval = iP-_iP0; //Shifted output index, if iP!=0
      for (unsigned iDim=0; iDim<D; iDim++) {
        output[iDim+iPeval*serialSize] += tmpForce[iDim];
      }
    } else {
      //Apply tmpForce as absolute force (no increment, as no empty output[] available)
      particle.template setField<FORCING,FORCE>( tmpForce );
    }
  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, bool serialize>
bool BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE, serialize>::operator()(T output[], const int input[])
{
  using namespace descriptors;
  // iterate over all particles in _indicator //TODO: add periodic treatment analogous to momentum exchange
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); iP++) {
    auto particle = _particleSystem.get(iP);
    if (particles::access::isValid(particle)){
      evaluate(output, particle, iP);
    }
  }
  return true;
}



}
#endif
