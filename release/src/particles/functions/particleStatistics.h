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



#ifndef PARTICLE_STATISTICS_H
#define PARTICLE_STATISTICS_H


namespace olb {

namespace particles {

namespace statistics {

/// Evaluate particle statistics for each cuboid (each ParticleSystem)
template<typename T, typename S, typename PARTICLETYPE, typename F>
void evaluateParticleSystemStatistics(
  SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, std::vector<S>& data, int noSampledQuantities, F f, int root=0 )
{
  //Set number of quantities (including globiC at first position)
  int noQ = noSampledQuantities + 1;
  //Retrieve load balancer and set up samples vector
  auto& loadBalancer = sParticleSystem.getSuperStructure().getLoadBalancer();
  std::vector<S> samples(loadBalancer.size()*noQ);
  //Iterate over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Set index for samples array
    int idx = iC*noQ;
    samples[idx] = globiC;
    //Perform requested sampling (assuming globiC at first pos)
    f ( particleSystem, globiC, idx+1, samples );
  });

  //Execute gathering
  const int numTasks = singleton::mpi().getSize();
  const int recvBufSize = samples.size()*numTasks;
  data.resize(recvBufSize);
  int sendCount = samples.size();
  int recvCounts[numTasks];
  int displs[numTasks];
  std::fill(recvCounts,&recvCounts[numTasks],sendCount);
  std::generate(displs,&displs[numTasks],[&sendCount,n = 0] () mutable { return sendCount*n++; });
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().gatherv<S>(samples.data(), sendCount, data.data(),
                recvCounts, displs, root);
#endif
}

/// Gather number of active particles and total number for each cuboid
template<typename T, typename PARTICLETYPE, typename PCONDITION = conditions::valid_particles>
std::vector<std::size_t> gatherActivity(
  SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem )
{
  //Set sampling properties
  unsigned noSampledQuantities = 2;
  std::vector<std::size_t> dataV;
  //Evaluate particle system statistics
  evaluateParticleSystemStatistics( sParticleSystem, dataV, noSampledQuantities,
    [](ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC, int idx, auto& samples){
    //Introduce couters for sampled quantities
    std::size_t noP = 0;
    std::size_t noActive = 0;
    //Iterate over particles
    forParticlesInParticleSystem<T,PARTICLETYPE,PCONDITION>( particleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){
      //Evaluate counter increment
      if ( access::isActive(particle) ){ ++noActive; }
      ++noP;
    },globiC); //forParticlesInParticleSystem<T,PARTICLETYPE,PCONDITION>
   //Update samples vector with counters
   samples[idx+0] = noActive;
   samples[idx+1] = noP;
  });
  //Return gathered dataV
  return dataV;
}


/// Print activity gathered a main processor
template<typename T, typename PARTICLETYPE>
void printActivityGathered( SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem ){
  OstreamManager clout ( std::cout, "Activity" );
  auto activityStatistics = gatherActivity( sParticleSystem );
  unsigned noQ = 3;
  for (int idx=0; idx<activityStatistics.size(); idx+=noQ ){
    clout << "ParticleSystem " << activityStatistics[idx] << ":"
          << " active=" << activityStatistics[idx+1]
          << "/" << activityStatistics[idx+2]
          << std::endl;
  }
}





} //namespace statistics

} //namespace particles

} //namespace olb


#endif
