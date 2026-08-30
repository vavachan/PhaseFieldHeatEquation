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


/* This file contains functions used for the calculation of particle related dynamics.
 *
*/

#ifndef PARTICLE_DYNAMICS_FUNCTIONS_H
#define PARTICLE_DYNAMICS_FUNCTIONS_H


#include "particles/contact/contactContainer.h"
#include "particles/contact/wall.h"
#include <cassert>

namespace olb {

namespace particles {

template <unsigned D>
constexpr bool isPeriodic(const Vector<bool, D>& periodic)
{
  bool isPeriodic = periodic[0] || periodic[1];
  if constexpr (D == 3) {
    isPeriodic = isPeriodic || periodic[2];
  }
  return isPeriodic;
}

namespace defaults{
  template <unsigned D>
  const auto periodicity = [](){
    if constexpr (D==3) {
      return Vector<bool,3>(false,false,false);
    } else {
      return Vector<bool,2>(false,false);
    }
  };
}

namespace contact {
template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void communicateContacts(ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer);
}

namespace dynamics {



//Calculate torque from force and lever
template<unsigned D, typename T>
struct torque_from_force {
  static constexpr Vector<T,utilities::dimensions::convert<D>::rotation> calculate(
    Vector<T,D> force, PhysR<T,D> lever )
  {
    if constexpr (D==2){
      return Vector<T,utilities::dimensions::convert<D>::rotation>( crossProduct( lever, force ) );
    } else {
      return crossProduct( lever, force );
    }
  }
};


////////////// Dimension sensitive Functions ////////////


//Calculate local velocity
template <typename T, typename PARTICLETYPE>
constexpr Vector<T,PARTICLETYPE::d> calculateLocalVelocity(Particle<T,PARTICLETYPE>& particle, const PhysR<T,PARTICLETYPE::d>& input)
{
  using namespace descriptors;

  const PhysR<T,PARTICLETYPE::d> position =
    particle.template getField<GENERAL,POSITION>();
  const Vector<T,PARTICLETYPE::d> velocity =
    particle.template getField<MOBILITY,VELOCITY>();
  const Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>
    angVelocity(particle.template getField<MOBILITY,ANG_VELOCITY>());

  return util::calculateLocalVelocity(position, velocity, angVelocity, input);
}


////////////// Force and Utility Functions ////////////

/// Unserialize force field provieded by force integration functor (e.g. momentumExchange)
template<typename T, typename PARTICLETYPE>
void unserializeForceTorqueVoxels( Vector<T,PARTICLETYPE::d>& force,
                                   Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& torque,
                                   T serializedFunctorForceField[], int iP )
{
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  constexpr int serialSize = D+Drot+1;

  const int idxForce = iP*(serialSize);
  const int idxTorque = idxForce+D;
  const int idxTorqueEnd = idxTorque+Drot;

  //Get force
  //TODO: include funcionality to OLB-Vector someday to avoid std::vector cast
  force = std::vector<T>(  serializedFunctorForceField+idxForce,
                           serializedFunctorForceField+idxTorque);
  //Get torque
  torque = std::vector<T>( serializedFunctorForceField+idxTorque,
                           serializedFunctorForceField+idxTorqueEnd );
}


/// Unserialize force field provieded by force integration functor (e.g. stokesDragForce)
template<typename T, typename PARTICLETYPE>
void unserializeForce( Vector<T,PARTICLETYPE::d>& force,
                       T serializedFunctorForceField[], int iP )
{
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr int serialSize = D;

  const int idxForce = iP*(serialSize);
  const int idxForceEnd = idxForce+D;

  //Get force
  //TODO: include funcionality to OLB-Vector someday to avoid std::vector cast
  force = std::vector<T>(  serializedFunctorForceField+idxForce,
                           serializedFunctorForceField+idxForceEnd);
}

/// Apply boundary force provided by force functor to the particle center as torque and force
template<typename T, typename PARTICLETYPE, typename FORCEFUNCTOR>
void applySerializableParticleForce( FORCEFUNCTOR& forceF, ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t iP0=0 )
{
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;

  //Create serialized force field and dummy input
  int input[1]{};
  T serializedFunctorForceField[forceF.getTargetDim()]; //TODO: could be decreased by only considering valid_particles

  //Initialize serialized force field when directly using block functors
  // -as this is usually done inside the super functor
  if constexpr (std::is_base_of<BlockF<T,PARTICLETYPE::d>, FORCEFUNCTOR>::value){
    for (int iS=0; iS<forceF.getTargetDim(); ++iS) {
      serializedFunctorForceField[iS] = 0.;
    }
  }

  //Retrieve boundary force field
  forceF(serializedFunctorForceField, input);

  //Loop over particles and apply individual force and torque contribution
  //TODO: for parallized particles, this represents a redundant loop over the particle system
  for (std::size_t iP=iP0; iP<particleSystem.size(); iP++) {
    auto particle = particleSystem.get(iP);

    Vector<T,D> force;
    std::size_t iPeval = iP-iP0; //Shift, if iP!=0
    if constexpr ( PARTICLETYPE::template providesNested<descriptors::FORCING,descriptors::TORQUE>() ) {
      Vector<T,Drot> torque;
      unserializeForceTorqueVoxels<T,PARTICLETYPE>( force, torque, serializedFunctorForceField, iPeval );
      particle.template setField<descriptors::FORCING,descriptors::TORQUE>(
        utilities::dimensions::convert<D>::serialize_rotation(torque) );
    }
    else {
      unserializeForce<T,PARTICLETYPE>( force, serializedFunctorForceField, iPeval );
    }

    //DEBUG OUTPUT
//    int rank = singleton::mpi().getRank();
//    std::cout << "  force(pSys=" << &particleSystem << ",rank=" << rank << ")=" << force << std::endl;
    particle.template setField<descriptors::FORCING,descriptors::FORCE>( force );
  }
}

/// Apply boundary force provided by force functor to the particle center as torque and force
/// - allows for additional specification of PCONDITION (as opposed to applySerializedParticleForce)
template<typename T, typename PARTICLETYPE, typename FORCEFUNCTOR, typename PCONDITION=conditions::valid_particles>
void applyLocalParticleForce( FORCEFUNCTOR& forceF, ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t iP0=0 )
{
  //Iterate over particles and apply functor's evaluate() directly
  for (std::size_t iP=iP0; iP!=particleSystem.size(); iP++) {
    auto particle = particleSystem.get(iP);
    //Execute F when particle meets condition
    doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
      [&](auto& particle){
      //Evaluate force functor
      T output[1];    //dummy output
      forceF.evaluate(output, particle, iP);
    });
  }
}

/// Initialize all fields in particle (necessary for clang)
template<typename T, typename PARTICLETYPE>
void initializeParticle( DynamicFieldGroupsD<T, typename PARTICLETYPE::fields_t>& dynamicFieldGroups, std::size_t iP )
{
  //Define init lambda expression
  typedef std::function<bool(const std::type_info&,int,std::string)> FunctionType;
  FunctionType initFunction = [](const std::type_info& typeInfo, int fieldSize, std::string fieldContentStr) {
    return true; //resetField=true
  };
  //Call recursive field traversal function with lambda expression
  descriptors::access_field_content<FunctionType,T,PARTICLETYPE,typename PARTICLETYPE::fields_t>::fieldsL2(
    initFunction, dynamicFieldGroups, iP );
}



////////////// Particle Lattice coupling Functions ////////////
//containing those, not included in particle tasks

/// Couple particle to lattice and detect contacts of resolved particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE,
  typename F=decltype(defaults::periodicity<DESCRIPTOR::d>)>
void coupleResolvedParticlesToLattice(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& contactContainer,
  const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
  SuperLattice<T,DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR> const& converter,
  std::vector<SolidBoundary<T,DESCRIPTOR::d>>& solidBoundaries,
  F getSetupPeriodicity = defaults::periodicity<DESCRIPTOR::d>)
{
  static_assert(DESCRIPTOR::template provides<descriptors::CONTACT_DETECTION>(),
                "The field CONTACT_DETECTION must be provided.");
  constexpr unsigned D = DESCRIPTOR::d;
  using namespace descriptors;

  contactContainer.cleanContacts();

  const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
  const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);

  //Loop over particles
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    //Write particle field
    setSuperParticleField( sGeometry, min, max, sLattice, converter,
                           particleSystem, contactContainer, iP, particle,
                           solidBoundaries, getSetupPeriodicity );

  }

  contact::communicateContacts<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>(contactContainer);
}

/// Couple particle to lattice and detect contacts of resolved particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE,
  typename F=decltype(defaults::periodicity<DESCRIPTOR::d>)>
void coupleResolvedParticlesToLattice(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& contactContainer,
  const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
  SuperLattice<T,DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR> const& converter,
  F getSetupPeriodicity = defaults::periodicity<DESCRIPTOR::d>)
{
  std::vector<SolidBoundary<T,DESCRIPTOR::d>> solidBoundaries = std::vector<SolidBoundary<T,DESCRIPTOR::d>>();
  coupleResolvedParticlesToLattice(particleSystem, contactContainer, sGeometry, sLattice, converter, solidBoundaries, getSetupPeriodicity);
}


/// Couple particle to lattice and detect contacts of resolved particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE,
  typename F=decltype(defaults::periodicity<DESCRIPTOR::d>)>
void coupleResolvedParticlesToLattice(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& contactContainer,
  const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
  SuperLattice<T,DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR> const& converter,
  std::vector<SolidBoundary<T,DESCRIPTOR::d>>& solidBoundaries,
  F getSetupPeriodicity = defaults::periodicity<DESCRIPTOR::d>)
{
  static_assert(DESCRIPTOR::template provides<descriptors::CONTACT_DETECTION>(),
                "The field CONTACT_DETECTION must be provided.");
  constexpr unsigned D = DESCRIPTOR::d;
  using namespace descriptors;

  const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
  const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);

  //Loop over particles
  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE, conditions::valid_particles>(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();
        //Write particle field
        setSuperParticleField(sGeometry, min, max, sLattice, converter,
                              particleSystem, contactContainer,
                              globalParticleID, particle, solidBoundaries,
                              getSetupPeriodicity, globiC );
      });
}

/// Couple particle to lattice and detect contacts of resolved particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE,
  typename F=decltype(defaults::periodicity<DESCRIPTOR::d>)>
void coupleResolvedParticlesToLattice(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& contactContainer,
  const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
  SuperLattice<T,DESCRIPTOR>& sLattice,
  UnitConverter<T,DESCRIPTOR> const& converter,
  F getSetupPeriodicity = defaults::periodicity<DESCRIPTOR::d>)
{
  std::vector<SolidBoundary<T,DESCRIPTOR::d>> solidBoundaries = std::vector<SolidBoundary<T,DESCRIPTOR::d>>();
  coupleResolvedParticlesToLattice(sParticleSystem, contactContainer, sGeometry, sLattice, converter, solidBoundaries, getSetupPeriodicity);
}


template<typename T, typename PARTICLETYPE>
T calcKineticEnergy( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  static_assert(D==3, "ERROR: 2D version of calcKineticEnergy not implemented yet!");

  T pMass = particle.template getField<PHYSPROPERTIES, MASS>();
  Vector<T,D> vel(particle.template getField<MOBILITY, VELOCITY>());
  Vector<T,Drot> pMofi(particle.template getField<PHYSPROPERTIES, MOFI>());
  Vector<T,Drot> angVel(particle.template getField<MOBILITY, ANG_VELOCITY>());

  T eTrans = T{0.5} * pMass * util::normSqr<T,D>(vel);
  T eRot = .5*(pMofi[0]*angVel[0]*angVel[0]
              +pMofi[1]*angVel[1]*angVel[1]
              +pMofi[2]*angVel[2]*angVel[2]);
  T eKin = eTrans+eRot;

  return eKin;
}


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
