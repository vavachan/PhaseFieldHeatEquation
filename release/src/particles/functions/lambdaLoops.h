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


#ifndef LAMBDA_LOOPS_H
#define LAMBDA_LOOPS_H


namespace olb {

namespace particles {



//Do for particle
template<typename T, typename PARTICLETYPE, typename F>
void doForParticle(
  Particle<T,PARTICLETYPE>& particle, F f)
{
  if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&>){
    f( particle );
  } else {
    f();
  }
}


//Do for particle satisfying the provided particle condition
//- Provides dynamic and static support for
//  - PCONDITION::value(particle)
//  - f(particle)
//  - f()
//- Note: "conditions::all_particles" is faster, when applicable
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void doWhenMeetingCondition(
  Particle<T,PARTICLETYPE>& particle,
  F f )
{
  //Check whether dynamic (constexpression value not possible)
  if constexpr (PCONDITION::dynamic){
    if ( PCONDITION::template value<T,PARTICLETYPE>(particle) ){
      doForParticle( particle, f);
    }
  }
  //Use static constexpression value
  //- for now, implies no template use
  else {
    if constexpr (PCONDITION::value){
      doForParticle( particle, f);
    }
  }
}

//Do for particle satisfying the provided particle condition
//- Provides dynamic and static support for
//  - PCONDITION::value(particle)
//  - PCONDITION::value(particle,globiC)
//  - f(particle)
//  - f()
//- Note: "conditions::all_particles" is faster, when applicable
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void doWhenMeetingCondition(
  Particle<T,PARTICLETYPE>& particle,
  F f, int globiC )
{
  //Check whether dynamic (constexpression value not possible)
  if constexpr (PCONDITION::dynamic){
    //Check whether globiC has to be passed
    if constexpr (std::is_invocable_v<decltype(PCONDITION::template value<T,PARTICLETYPE>),
                  Particle<T,PARTICLETYPE>&,int>){
      if ( PCONDITION::template value<T,PARTICLETYPE>(particle,globiC) ){
        if constexpr (std::is_invocable_v<F,Particle<T,PARTICLETYPE>&>){
          f( particle );
        } else {
          f();
        }
      }
    } else {
      //Call without globiC
      if ( PCONDITION::template value<T,PARTICLETYPE>(particle) ){
        doForParticle( particle, f);
      }
    }
  }
  //Use static constexpression value
  else {
    //Call without globiC
    if constexpr (PCONDITION::value){
      doForParticle( particle, f);
    }
  }
}



//Iterate over particles in particle system satisfying the provided particle condition
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void forParticlesInParticleSystem(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  F f, int globiC )
{
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    //Execute F when particle meets condition
    doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
      [&](Particle<T,PARTICLETYPE> particle){
      f( particle );
    }, globiC );
  }
}

//Iterate over particles in particle system satisfying the provided particle condition
//- no globiC required
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void forParticlesInParticleSystem(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  F f )
{
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    //Execute F when particle meets condition
    doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
      [&](Particle<T,PARTICLETYPE> particle){
      f( particle );
    });
  }
}


namespace communication {

//Iterator over super particle system
template<typename T, typename PARTICLETYPE, typename F>
void forSystemsInSuperParticleSystem(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  F f )
{
  //Retrieve load balancer
  auto& loadBalancer = sParticleSystem.getSuperStructure().getLoadBalancer();
  //Loop over iCs
  for (int iC=0; iC<loadBalancer.size(); ++iC){
    int globiC = loadBalancer.glob(iC);
    //Retrieve container
    auto bParticleSystems = sParticleSystem.getBlockParticleSystems();
    auto& particleSystem = *bParticleSystems[iC];
    f( particleSystem, iC, globiC );
  }
}


//Iterate over particles in super particle system
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void forParticlesInSuperParticleSystem(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  F f )
{
  //Iterate over particle systems
  forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Iterate over particles
    forParticlesInParticleSystem<T,PARTICLETYPE,PCONDITION>( particleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){
      if constexpr (std::is_invocable_v<F,
        Particle<T,PARTICLETYPE>&,ParticleSystem<T,PARTICLETYPE>&,int>)
      {
        f( particle, particleSystem, globiC );
      } else {
        f( particle );
      }
    },globiC);
  });
}

//Do for particle in superParticleSystem
//- enables retrieval via particle locator, hence circumvents
// having to search for a speficic particle again
template<typename T, typename PARTICLETYPE, typename F>
void doForParticle(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  ParallelParticleLocator& locator, F f)
{
  //Retrieve rank and only proceed for correct one
  int rank = singleton::mpi().getRank();
  auto& loadBalancer = sParticleSystem.getSuperStructure().getLoadBalancer();
  if (rank==loadBalancer.rank(locator.globiC)){
    //Retrieve ParticleSystem
    auto bParticleSystems = sParticleSystem.getBlockParticleSystems();
    auto& particleSystem = *bParticleSystems[loadBalancer.loc(locator.globiC)];
    //Retrieve particle
    auto particle = particleSystem.get(locator.localID);
    //Call function f for particle
    doForParticle( particle, f );
  }
}



} //namespace communication


/// Iterate over particles in x particle system
/// - works both for ParticleSystem and SuperParticleSystem
/// - can both be run for
///   - [](Particle<T,PARTICLETYPE>,ParticleSystem<T,PARTICLETYPE>,int){}
///   - [](Particle<T,PARTICLETYPE>){}
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles, typename F>
void forParticlesInXParticleSystem(
  XParticleSystem<T,PARTICLETYPE>& xParticleSystem,
  F f )
{
  if constexpr( access::providesParallelization<PARTICLETYPE>() ){
    //Iterate over particle systems
    communication::forSystemsInSuperParticleSystem( xParticleSystem,
      [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
      //Iterate over particles
      forParticlesInParticleSystem<T,PARTICLETYPE,PCONDITION>( particleSystem,
        [&](Particle<T,PARTICLETYPE>& particle){
        if constexpr (std::is_invocable_v<F,
          Particle<T,PARTICLETYPE>&,ParticleSystem<T,PARTICLETYPE>&,int>)
        {
          f( particle, particleSystem, globiC );
        } else {
          f( particle );
        }
      },globiC);
    });
  } else {
    //Loop over all particles
    forParticlesInParticleSystem<T,PARTICLETYPE,PCONDITION>( xParticleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){
        if constexpr (std::is_invocable_v<F,
          Particle<T,PARTICLETYPE>&,ParticleSystem<T,PARTICLETYPE>&,int>)
        {
          f( particle, xParticleSystem, 0 );
        } else {
          f( particle );
        }
    });
  }
}





} //namespace particles

} //namespace olb


#endif
