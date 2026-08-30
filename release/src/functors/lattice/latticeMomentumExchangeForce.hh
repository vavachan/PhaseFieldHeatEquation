/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Robin Trunk, Mathias J. Kraus
 *                2021 Nicolas Hafen
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

#ifndef LATTICE_MOMENTUM_EXCHANGE_FORCE_HH
#define LATTICE_MOMENTUM_EXCHANGE_FORCE_HH

#include "latticeMomentumExchangeForce.h"
#include "particles/resolved/blockLatticeInteraction.hh"

namespace olb {


template <typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename BLOCKFUNCTOR>
SuperLatticeParticleForce<T,DESCRIPTOR,PARTICLETYPE,BLOCKFUNCTOR>::SuperLatticeParticleForce(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  Vector<bool,DESCRIPTOR::d> periodic, std::size_t iP0,
  const std::set<int>& ignoredMaterials, const F f )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice,converter,
                                    (DESCRIPTOR::d+utilities::dimensions::convert<DESCRIPTOR::d>::rotation+1)*(particleSystem.size()-iP0))
{
  constexpr unsigned D = DESCRIPTOR::d;

  this->getName() = "physParticleForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  const PhysR<T,D> min = particles::communication::getCuboidMin<T,D>(
      superGeometry.getCuboidDecomposition());
  const PhysR<T,D> max = particles::communication::getCuboidMax<T,D>(
      superGeometry.getCuboidDecomposition(), min);

  for (int iC = 0; iC < maxC; ++iC) {
    this->_blockF.emplace_back( new BLOCKFUNCTOR(
                                  this->_sLattice.getBlock(iC),
                                  superGeometry.getBlockGeometry(iC),
                                  particleSystem, converter, min, max,
                                  periodic, iP0, ignoredMaterials, f));
  }
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename BLOCKFUNCTOR>
bool SuperLatticeParticleForce<T,DESCRIPTOR,PARTICLETYPE,BLOCKFUNCTOR>::operator() (T output[],
    const int input[])
{
  for (int iS=0; iS<this->getTargetDim(); ++iS) {
    output[iS] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->getBlockF(iC)(output, &input[1]);
  }

#ifdef PARALLEL_MODE_MPI
  for (int iS = 0; iS < this->getTargetDim(); ++iS) {
    singleton::mpi().reduceAndBcast(output[iS], MPI_SUM);
  }
#endif
  return true;

}


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeMomentumExchangeForce(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  PhysR<T,DESCRIPTOR::d>cellMin, PhysR<T,DESCRIPTOR::d> cellMax,
  Vector<bool,DESCRIPTOR::d> periodic, std::size_t iP0,
  const std::set<int>& ignoredMaterials,  const F f )
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice, converter,
                                   (DESCRIPTOR::d+utilities::dimensions::convert<DESCRIPTOR::d>::rotation+1)*(particleSystem.size()-iP0)),
    _blockGeometry(blockGeometry), _blockLattice(blockLattice),
    _particleSystem(particleSystem),
    _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic),
    _iP0(iP0), _ignoredMaterials(ignoredMaterials), _f(f)
{
  this->getName() = "physMomentumExchangeForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned D = DESCRIPTOR::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation; //TODO: prevent rotation calculatiion, when no angle
  constexpr int serialSize = D+Drot+1;

  using namespace descriptors;
  const PhysR<T, D> position = particles::access::getPosition(particle);
  const T circumRadius = particles::access::getRadius(particle);

  //For all cells in block particle intersection
  int numVoxels = 0;
  constexpr int padding = 0;
  const std::size_t iPeval = iP - _iP0; //Shifted output index, if iP!=0
  particles::forSpatialLocationsInBlockParticleIntersection( _blockGeometry,
      _blockLattice, padding, position, circumRadius,
  [&](const LatticeR<D>& latticeRinner) {

    if (_ignoredMaterials.find(_blockGeometry.getMaterial(latticeRinner)) == _ignoredMaterials.end()) {
      Vector<T,D>tmpForce(0.);
      PhysR<T,D> lever(0.);
      if ( particles::resolved::momentumExchangeAtSurfaceLocation(tmpForce.data(),
            lever, latticeRinner, this->_blockGeometry,
            this->_blockLattice, this->_converter, particle) ){

        // count cells of considered object
        ++numVoxels;

        //Shift centre of mass (for oblique rotation) if offset provided
        if constexpr ( PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>() ){
          lever -= Vector<T,D> ( particle.template getField<SURFACE,COR_OFFSET>() );
        }

        //Calculate torque
        const Vector<T,Drot> torque = particles::dynamics::torque_from_force<D,T>::calculate( tmpForce, lever );

        Vector<T,D> physRv;
        this->_blockGeometry.getPhysR(physRv, latticeRinner);
        _f(particle, physRv.data(), tmpForce, torque);

        //Add force and torque to output
        for (unsigned iDim=0; iDim<D; ++iDim) {
          output[iDim+iPeval*serialSize] += tmpForce[iDim];
        }
        for (unsigned iRot=0; iRot<Drot; ++iRot) {
          output[(iRot+D)+iPeval*serialSize] += torque[iRot];
        }
      }
    }
  });

  output[D+Drot+iPeval*serialSize] += numVoxels;
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  constexpr unsigned D = DESCRIPTOR::d;

  using namespace descriptors;

  //Evaluate periodicity
  bool anyPeriodic = false;
  for (unsigned iDim=0; iDim<D; ++iDim) {
    anyPeriodic = anyPeriodic || _periodic[iDim];
  }

  // iterate over all particles in _particleSystem
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); ++iP) {
    auto particle = _particleSystem.get(iP);
    if (particles::access::isValid(particle)){
      // If any direction requires periodic treatment
      if (anyPeriodic) {
        const T circumRadius = particles::access::getRadius(particle);
        const PhysR<T,D> position = particles::access::getPosition(particle);

        bool surfaceOutOfGeometry = false;
        PhysR<T,D> ghostPos;
        particles::checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos,
          _cellMin, _cellMax, position, circumRadius, _periodic);

        if (!surfaceOutOfGeometry) {
          evaluate(output, particle, iP);
        }
        else {
          // Calculate force for the ghost particle
          particle.template setField<GENERAL,POSITION>(ghostPos);
          evaluate(output, particle, iP);
          // Calculate force of actual particle
          particle.template setField<GENERAL,POSITION>(position);
          evaluate(output, particle, iP);
        }
      }
      else {
        evaluate(output, particle, iP);
      }
    }
  }
  return true;
}


///The following are functors that work in the traditional (output[], input[]) sense,
///They can therefore be used e.g. in the vtk writer as well


template<typename T,typename DESCRIPTOR,typename PARTICLETYPE, bool useTorque>
SuperLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE,useTorque>::SuperLatticeMomentumExchangeForceLocal(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  const SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice, converter,
      useTorque ? utilities::dimensions::convert<DESCRIPTOR::d>::rotation : DESCRIPTOR::d)
{
  this->getName() = "localMomentumExchange";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; ++iC) {
    this->_blockF.emplace_back(new BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE,useTorque>(
      this->_sLattice.getBlock(iC),
      superGeometry.getBlockGeometry(iC),
      particleSystem,
      this->_converter));
  }
}

template<typename T,typename DESCRIPTOR,typename PARTICLETYPE, bool useTorque>
SuperLatticeMomentumExchangeForceLocalParallel<T,DESCRIPTOR,PARTICLETYPE,useTorque>::SuperLatticeMomentumExchangeForceLocalParallel(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  const SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice, converter,
      useTorque ? utilities::dimensions::convert<DESCRIPTOR::d>::rotation : DESCRIPTOR::d)
{
  this->getName() = "localMomentumExchange";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; ++iC) {

    auto bParticleSystems = sParticleSystem.getBlockParticleSystems();
    auto& particleSystem = *bParticleSystems[iC];

    this->_blockF.emplace_back(new BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE,useTorque>(
      this->_sLattice.getBlock(iC),
      superGeometry.getBlockGeometry(iC),
      particleSystem,
      this->_converter));
  }
}


template <typename T, typename DESCRIPTOR,typename PARTICLETYPE, bool useTorque>
BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE,useTorque>::BlockLatticeMomentumExchangeForceLocal(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice,converter,
      useTorque ? utilities::dimensions::convert<DESCRIPTOR::d>::rotation : DESCRIPTOR::d),
    _blockLattice(blockLattice),
    _blockGeometry(blockGeometry),
    _particleSystem(particleSystem)
{
  this->getName() = "localMomentumExchange";
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE, bool useTorque>
bool BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE,useTorque>::operator() (T output[], const int input[])
{
  using namespace descriptors;
  const unsigned D = DESCRIPTOR::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;

  for (unsigned iDim=0; iDim<D; ++iDim) {
    output[iDim] = 0.;
  }

  // iterate over all particles in _particleSystem
  int _iP0 = 0;
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); ++iP) {
    auto particle = _particleSystem.get(iP);
    auto circumRadius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();
    auto position = particle.template getField<GENERAL,POSITION>();

    LatticeR<DESCRIPTOR::d> start, end;
    if ( particles::getBlockParticleIntersection(this->_blockGeometry, T(1.)/this->_converter.getPhysDeltaX(), start, end,
                                      position, circumRadius) ){

      Vector<T,D>tmpForce(0.);
      PhysR<T,D> lever(0.);
      if ( particles::resolved::momentumExchangeAtSurfaceLocation(tmpForce.data(),
            lever, input, this->_blockGeometry,
            this->_blockLattice, this->_converter, particle) ){

        //Shift centre of mass (for oblique rotation) if offset provided
        if constexpr ( PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>() ){
          lever -= Vector<T,D> ( particle.template getField<SURFACE,COR_OFFSET>() );
        }

        //Calculate torque
        const auto torque = particles::dynamics::torque_from_force<D,T>::calculate( tmpForce, lever );

        //Add force and torque to output
        if constexpr(!useTorque){
          for (unsigned iDim=0; iDim<D; ++iDim) {
            output[iDim] = tmpForce[iDim];
          }
        } else {
          for (unsigned iRot=0; iRot<Drot; ++iRot) {
            output[iRot] = torque[iRot];
          }
        }
      }
    }
  }
  return true;
}



}
#endif
