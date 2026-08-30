/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_MANAGER_HH
#define PARTICLE_MANAGER_HH


#include "resolved/blockLatticeInteraction.h"

//#define VERBOSE_PARTICLEMANAGER

namespace olb {

namespace particles {

namespace dynamics {

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::ParticleManager(
   XParticleSystem<T,PARTICLETYPE>& xParticleSystem,
   SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
   SuperLattice<T,DESCRIPTOR>& sLattice,
   UnitConverter<T,DESCRIPTOR> const& converter,
   Vector<T,PARTICLETYPE::d> externalAcceleration,
   Vector<bool,PARTICLETYPE::d> periodic )
   : _xParticleSystem(xParticleSystem), _sGeometry(sGeometry),
     _sLattice(sLattice), _converter(converter), _externalAcceleration(externalAcceleration),
     _periodic(periodic)
{}

//Unpack tasks requiring loop
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename taskList, typename ISEQ>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::unpackTasksLooped(
  Particle<T,PARTICLETYPE>& particle, T timeStepSize, ISEQ indexSequence, int globiC)
{
  //Define function for task execution
  auto executeTask = [&](auto task){
    //Derive type of task
    using TASK = typename decltype(task)::type;
#ifdef VERBOSE_PARTICLEMANAGER
      OstreamManager clout( std::cout,"ParticleManager" );
      if constexpr (access::providesParallelization<PARTICLETYPE>()){
        clout.setMultiOutput(true);
      }
      std::string taskStr = typeid(TASK).name();
      taskStr.resize(50);
      clout << "-TaskLooped(globiC=" << globiC
            <<")-Particle " << particle.getId() << ": "
                << taskStr << std::endl;
      clout.setMultiOutput(false);
#endif
    //Evaluate coupling method
    if constexpr(TASK::latticeCoupling){
      TASK::execute( _xParticleSystem, particle, _sGeometry, _sLattice, _converter, globiC, _periodic );
    } else {
      TASK::execute( _xParticleSystem, particle, _externalAcceleration, timeStepSize, globiC );
    }
  };
  //Execute single tasks
  meta::list_for_each_index<taskList>(executeTask,indexSequence);
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename ...TASKLIST>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::execute(T timeStepSize){
#ifdef VERBOSE_PARTICLEMANAGER
      OstreamManager clout( std::cout,"ParticleManager" );
      clout << "TASKLIST:" << std::endl;
      singleton::mpi().synchronizeIO();
#endif
  //Create meta::list of TASKLIST parameter pack
  using taskList = typename meta::list<TASKLIST...>;
  //Retrieve index sequence for tasks to be executed without loop
  auto indexSequenceNoLoop =
    meta::filter_index_sequence<taskList,requires_no_loop>();
  //Define function for tasks not requiring loop
  auto executeNoLoopTask = [&](auto task){
    using TASK = typename decltype(task)::type;
#ifdef VERBOSE_PARTICLEMANAGER
      singleton::mpi().synchronizeIO();
      OstreamManager clout( std::cout,"ParticleManager" );
      if constexpr (access::providesParallelization<PARTICLETYPE>()){
        clout.setMultiOutput(true);
      }
      std::string taskStr = typeid(TASK).name();
      taskStr.resize(50);
      clout << "-Task: " << taskStr << std::endl;
      clout.setMultiOutput(false);
      singleton::mpi().synchronizeIO();
#endif
    if constexpr(TASK::latticeCoupling){
      TASK::execute( _xParticleSystem, _sGeometry, _sLattice, _converter, _periodic );
    } else {
      if constexpr (std::is_invocable_v<decltype(TASK::execute),
        XParticleSystem<T,PARTICLETYPE>&, const T,
        const communication::ParticleCommunicator&,
        const Vector<bool,PARTICLETYPE::d>&>)
      {
        TASK::execute( _xParticleSystem, _converter.getPhysDeltaX(),
            _communicator, _periodic);
      }
      else {
        TASK::execute( _xParticleSystem, _communicator);
      }
    }
#ifdef VERBOSE_PARTICLEMANAGER
      singleton::mpi().synchronizeIO(1000);
#endif
  };
  //Define function for sequence of tasks requiring loop
  auto executeLoopedTasks = [&](auto indexSequence){
    if constexpr( access::providesParallelization<PARTICLETYPE>() ){
      //Loop over valid particles
      communication::forParticlesInSuperParticleSystem<T,PARTICLETYPE,
        conditions::valid_particles>( _xParticleSystem,
        [&](Particle<T,PARTICLETYPE>& particle,
        ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC){
        //Unpack tasks to be looped
        unpackTasksLooped<taskList>( particle, timeStepSize, indexSequence, globiC );
      });
    } else {
      //Loop over all particles
      forParticlesInParticleSystem<T,PARTICLETYPE,conditions::all_particles>( _xParticleSystem,
        [&](Particle<T,PARTICLETYPE>& particle){
        //Unpack tasks to be looped
        unpackTasksLooped<taskList>( particle, timeStepSize, indexSequence, 0 );
      });
    }
   };
   //Evaluate index sequence and evaluate tasks or tasks sequences
   //- Example execution:
   //   Direct: T1
   //   Looped: [T2,T3,T4] (iterate taks for same particle)
   //   Direct: T5
   //   Direct: T6
   //   Looped: [T7,T8]    (iterate taks for same particle)
   meta::index_sequence_for_subsequence<taskList>(
    executeNoLoopTask, executeLoopedTasks, indexSequenceNoLoop);
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename ...TASKLIST>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::execute(){ execute<TASKLIST...>(_converter.getPhysDeltaT()); }

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
const communication::ParticleCommunicator& ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::getParticleCommunicator(){
  return _communicator;
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::setExternalAcceleration(const Vector<T,PARTICLETYPE::d>& externalAcceleration)
{
  _externalAcceleration = externalAcceleration;
}


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
