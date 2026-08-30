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


/* This file containes particle tasks that can be passed to the particle manager.
 * Those need to provide an execute method, and two booleans specifying the coupling
 * and the necessity for beeing looped over all particles.
*/

#ifndef PARTICLE_TASKS_H
#define PARTICLE_TASKS_H


#include <cassert>

namespace olb {

namespace particles {

/// Couple lattice to particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
         typename FORCEFUNCTOR=SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>>
struct couple_lattice_to_particles_single_cuboid{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<bool,DESCRIPTOR::d> periodicity =
      Vector<bool,DESCRIPTOR::d>(false),
    std::size_t iP0=0 )
  {
    //Momentum Exchange/Loss
    FORCEFUNCTOR forceFunctor( sLattice, sGeometry,
                               particleSystem, converter, periodicity, iP0 );
    //Store boundary force in particles
    dynamics::applySerializableParticleForce( forceFunctor, particleSystem, iP0 );
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = false;
};

/// Apply external acceleration (e.g. for apply gravity)
template<typename T, typename PARTICLETYPE>
struct apply_external_acceleration_single_cuboid{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& ParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d> externalAcceleration,
    T timeStepSize, int globiC=0)
  {
    using namespace descriptors;
    //Apply acceleration
    Vector<T,PARTICLETYPE::d> force = access::getForce(particle);
    const T mass = access::getMass(particle);
    force += externalAcceleration * mass;
    particle.template setField<FORCING,FORCE>( force );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};

/// Process particle dynamics
template<typename T, typename PARTICLETYPE>
struct process_dynamics_single_cuboid{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& externalAcceleration,
    T timeStepSize, int globiC=0 )
  {
    //Call process on particle dynamics
    if constexpr (!access::providesDynamicsID<PARTICLETYPE>()){
      particle.process(timeStepSize);
    } else {
      unsigned short dynamicsID = access::getDynamicsID( particle );
      particle.process(timeStepSize, dynamicsID);
    }
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};

/// Couple particles to lattice
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
struct couple_particles_to_lattice_single_cuboid{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    Particle<T,PARTICLETYPE>& particle,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    int globiC = 0,
    Vector<bool,DESCRIPTOR::d> periodicity =
      Vector<bool,DESCRIPTOR::d>(false) )
  {
    using namespace descriptors;
    //Write particle field
    setSuperParticleField( sGeometry, sLattice, converter, particle, periodicity );
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = true;
};







#ifdef PARALLEL_MODE_MPI



/// Process particle dynamics
template<typename T, typename PARTICLETYPE>
struct process_dynamics_parallel{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& externalAcceleration,
    T timeStepSize, int globiC )
  {
    using PCONDITION = std::conditional_t<
      access::providesSurface<PARTICLETYPE>(),
      conditions::valid_particle_centres,       //only consider centre for resolved
      conditions::valid_particles>;             //only consider valid
    //Call process on particle dynamics when meeting PCONDITION
    doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
      [&](Particle<T,PARTICLETYPE> particle){
      particle.process(timeStepSize);
    }, globiC );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};


/// Couple lattice to parallel particles
/// - Resolved default: MomentumExchange, serialized
/// - Subgrid default:  StokesDrag, non-serialized
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
         typename FORCEFUNCTOR=BlockLatticeMomentumExchangeForce<T,DESCRIPTOR,PARTICLETYPE>>
struct couple_lattice_to_parallel_particles{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<bool,DESCRIPTOR::d> periodicity =
      Vector<bool,DESCRIPTOR::d>(false),
    std::size_t iP0=0 )
  {
    constexpr unsigned D = DESCRIPTOR::d;
    const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
    const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);

    //Test output
    communication::forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

      auto& blockLattice = sLattice.getBlock(iC);
      auto& blockGeometry = sGeometry.getBlockGeometry(iC);

      //Momentum Exchange/Loss
      FORCEFUNCTOR forceFunctor( blockLattice, blockGeometry,
                                 particleSystem, converter, min, max, periodicity, iP0 );

      //Store boundary force in particles
      // - if not serialized, allows for additional filtering (e.g. active_particles)
      if constexpr (FORCEFUNCTOR::serializeForce){
        dynamics::applySerializableParticleForce( forceFunctor, particleSystem, iP0 );
      } else {
        using PCONDITION = conditions::active_particles; //TODO: remove hardcoded active_particle limitation
        dynamics::applyLocalParticleForce<T,PARTICLETYPE,FORCEFUNCTOR,PCONDITION>(
          forceFunctor, particleSystem, iP0 );
      }
    });
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = false;
};



/// Couple lattice to parallel particles
/// - Resolved default: MomentumExchange, serialized
/// - Subgrid default:  StokesDrag, non-serialized
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
         typename FORCEFUNCTOR=BlockLatticeMomentumExchangeForce<T,DESCRIPTOR,PARTICLETYPE>,
         typename FORCEFUNCTOR2=BlockLatticeMomentumExchangeForce<T,DESCRIPTOR,PARTICLETYPE>>
struct couple_lattice_to_parallel_particles_two_forces{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<bool,DESCRIPTOR::d> periodicity =
      Vector<bool,DESCRIPTOR::d>(false),
    std::size_t iP0=0 )
  {
    constexpr unsigned D = DESCRIPTOR::d;
    const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
    const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);

    //Test output
    communication::forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

      auto& blockLattice = sLattice.getBlock(iC);
      auto& blockGeometry = sGeometry.getBlockGeometry(iC);

      //Momentum Exchange/Loss
      FORCEFUNCTOR forceFunctor( blockLattice, blockGeometry,
                                 particleSystem, converter, min, max, periodicity, iP0 );
      FORCEFUNCTOR2 forceFunctor2( blockLattice, blockGeometry,
                                 particleSystem, converter, min, max, periodicity, iP0 );


      //Store boundary force in particles
      // - if not serialized, allows for additional filtering (e.g. active_particles)
      if constexpr (FORCEFUNCTOR::serializeForce){
        dynamics::applySerializableParticleForce( forceFunctor, particleSystem, iP0 );
        dynamics::applySerializableParticleForce( forceFunctor2, particleSystem, iP0 );
      } else {
        using PCONDITION = conditions::active_particles; //TODO: remove hardcoded active_particle limitation
        dynamics::applyLocalParticleForce<T,PARTICLETYPE,FORCEFUNCTOR,PCONDITION>(
          forceFunctor, particleSystem, iP0 );
           dynamics::applyLocalParticleForce<T,PARTICLETYPE,FORCEFUNCTOR2,PCONDITION>(
          forceFunctor2, particleSystem, iP0 );
      }
    });
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = false;
};




/// Update particle core distribution of parallel particles
template<typename T, typename PARTICLETYPE>
struct update_particle_core_distribution{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    const T physDeltaX,
    const communication::ParticleCommunicator& communicator,
    const Vector<bool, PARTICLETYPE::d>& periodicity
    )
  {
    //Update particle distribution
    communication::updateParticleCuboidDistribution( sParticleSystem, physDeltaX,
//#ifdef PARALLEL_MODE_MPI
        communicator.particleDistribution,
//#endif
        periodicity
      );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = false;
};

/// Communicate surface force of parallel particles
template<typename T, typename PARTICLETYPE>
struct communicate_surface_force{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    const communication::ParticleCommunicator& communicator
    )
  {
    using namespace descriptors;

    //Retrieve dimension of force data
    constexpr unsigned D = PARTICLETYPE::d;
    constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
    constexpr unsigned forceDataSize = D+Drot;

    //Create globalIdDataMap
    std::map<std::size_t, Vector<T,forceDataSize>> globalIdDataMap;

    //Communicate surface force and add to globalIdDatamap
    communication::communicateSurfaceForce( sParticleSystem, globalIdDataMap
//#ifdef PARALLEL_MODE_MPI
        , communicator.surfaceForceComm
//#endif
      );

    //Assign data in globalIdData Map to correct particle
    communication::assignSurfaceForce( sParticleSystem, globalIdDataMap );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = false;
};


/// Couple particles to lattice
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
struct couple_parallel_particles_to_lattice{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    int globiC,
    Vector<bool,DESCRIPTOR::d> periodicity =
      Vector<bool,DESCRIPTOR::d>(false) )
  {
    constexpr unsigned D = DESCRIPTOR::d;
    using namespace descriptors;

    //Retrieve load balancer and local iC
    auto& superStructure = sParticleSystem.getSuperStructure();
    auto& loadBalancer = superStructure.getLoadBalancer();
    int iC = loadBalancer.loc(globiC);

    //Retrieve blockLattice and blockGeometry
    auto& blockLattice = sLattice.getBlock(iC);
    auto& blockGeometry = sGeometry.getBlockGeometry(iC);

    //Write particle to field
    if(isPeriodic(periodicity)) {
      const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
      const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);

      setBlockParticleField( blockGeometry, blockLattice, converter, min, max,
                           particle, periodicity);
    } else {
      setBlockParticleField( blockGeometry, blockLattice, converter, particle );
    }
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = true;
};



/// Apply external acceleration (e.g. for apply gravity)
template<typename T, typename PARTICLETYPE>
struct apply_external_acceleration_parallel{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d> externalAcceleration,
    T timeStepSize, int globiC=0)
  {
    using namespace descriptors;
    //Apply acceleration
    Vector<T,PARTICLETYPE::d> force = access::getForce(particle);
    const T mass = access::getMass(particle);
    force += externalAcceleration * mass;
    particle.template setField<FORCING,FORCE>( force );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};


/// Clear force and torque values, apply when using more forces
template<typename T, typename PARTICLETYPE>
struct clear_force_and_torque{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d> externalAcceleration,
    T timeStepSize, int globiC=0)
  {
    using namespace olb::descriptors;
    Vector<T,3> zero (0.,0.,0.);
    particle.template setField<FORCING,FORCE>(zero);
    static_assert(PARTICLETYPE::template providesNested<FORCING,TORQUE>(), "Field FORCING,TORQUE has to be provided");
    particle.template setField<FORCING,TORQUE>(zero);
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};

///Euler Rotation -computes the orientation of the main axis of the spheroid to the streamline
template<typename T, typename PARTICLETYPE>
struct compute_orientation_index{
  auto static execute(
    SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
    Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d> externalAcceleration,
    T timeStepSize, int globiC=0)
  {

    using namespace olb::descriptors;
    using namespace olb::util;
    static_assert(PARTICLETYPE::template providesNested<NUMERICPROPERTIES,ORIENTATION>(), "Field NUMERICPROPERTIES,ORIENTATION has to be provided");
    static_assert(PARTICLETYPE::template providesNested<EULER_ROTATION,ORIENTATION_ANGLE>(), "Field EULER_ROTATION,ORIENTATION_INDEX has to be provided");
    Vector<T,3> zero (0.,0.,0.);
    Vector<T,3> fvel (particle.template getField<MOBILITY,FLUIDVEL>());
    Vector<T,3> porient (particle.template getField<NUMERICPROPERTIES,ORIENTATION>());
    particle.template setField<EULER_ROTATION,ORIENTATION_ANGLE>(acos(fabs(dotProduct3D(fvel, porient)/norm(fvel)/norm(porient)))/3.141592653*180.);
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};






/// Aliases
//- for use in examples
//- differentiating between resolved and subgrid particles

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
using couple_lattice_to_particles = std::conditional_t<
  access::providesSurface<PARTICLETYPE>(),
  std::conditional_t<
    access::providesParallelization<PARTICLETYPE>(),
    couple_lattice_to_parallel_particles<T,DESCRIPTOR,PARTICLETYPE,
     BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>>,
    couple_lattice_to_particles_single_cuboid<T,DESCRIPTOR,PARTICLETYPE,
     SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>>>,
  couple_lattice_to_parallel_particles<T,DESCRIPTOR,PARTICLETYPE,
    BlockLatticeStokesDragForce<T,DESCRIPTOR,PARTICLETYPE,false>>
>;

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
using couple_particles_to_lattice = std::conditional_t<
  access::providesParallelization<PARTICLETYPE>(),
  couple_parallel_particles_to_lattice<T,DESCRIPTOR,PARTICLETYPE>,
  couple_particles_to_lattice_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>
>;

template<typename T, typename PARTICLETYPE>
using process_dynamics = std::conditional_t<
  access::providesParallelization<PARTICLETYPE>(),
  process_dynamics_parallel<T,PARTICLETYPE>,
  process_dynamics_single_cuboid<T,PARTICLETYPE>
>;

template<typename T, typename PARTICLETYPE>
using apply_gravity = std::conditional_t<
  access::providesParallelization<PARTICLETYPE>(),
  apply_external_acceleration_parallel<T,PARTICLETYPE>,
  apply_external_acceleration_single_cuboid<T,PARTICLETYPE>
>;

template<typename T, typename PARTICLETYPE>
using communicate_parallel_surface_force = communicate_surface_force<T,PARTICLETYPE>;

#else
//NON-MPI support for resolved particles

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
using couple_lattice_to_particles =
  couple_lattice_to_particles_single_cuboid<T,DESCRIPTOR,PARTICLETYPE,
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>>;

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
using couple_particles_to_lattice =
  couple_particles_to_lattice_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>;

template<typename T, typename PARTICLETYPE>
using process_dynamics =
  process_dynamics_single_cuboid<T,PARTICLETYPE>;

template<typename T, typename PARTICLETYPE>
using apply_gravity =
  apply_external_acceleration_single_cuboid<T,PARTICLETYPE>;


#endif




} //namespace particles

} //namespace olb


#endif
