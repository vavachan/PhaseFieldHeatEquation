/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Nicolas Hafen, Mathias J. Krause
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



#ifndef SUBGRID_UTILITIES_H
#define SUBGRID_UTILITIES_H

#include "particles/functions/particleIoFunctions.h"

namespace olb {

namespace particles {

namespace subgrid {


// Initialize particle velocities
template<typename T, typename PARTICLETYPE, typename DESCRIPTOR>
void initializeParticleVelocity(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  UnitConverter<T,DESCRIPTOR> const& converter,
  SuperParticleSystem<T,PARTICLETYPE>& supParticleSystem )
{
  using namespace descriptors;
  constexpr unsigned D = DESCRIPTOR::d;
  // Iterate over individual particle systems
  communication::forSystemsInSuperParticleSystem( supParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    // Retrieve block quantities
    auto& blockLattice = sLattice.getBlock(iC);
    auto& blockGeometry = superGeometry.getBlockGeometry(iC);
    const auto& cuboid = blockGeometry.getCuboid();
    // Create block velocity interpolation functor
    BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> blockInterpPhysVelF(
      blockLattice, converter, cuboid);
    // Loop over all particles in particle system
    forParticlesInParticleSystem<T,PARTICLETYPE>( particleSystem,
      [&](Particle<T,PARTICLETYPE>& particle){
      // Retrieve particle position
      Vector<T,D> position = particle.template getField<GENERAL,POSITION>();
      // Calculate interpolated velocity
      Vector<T,D> fluidVel(0.);
      blockInterpPhysVelF(fluidVel.data(), position.data());
      //Apply velocity to particle
      particle.template setField<MOBILITY,VELOCITY>( fluidVel );
    });
  });
}


// Add particles in indicator domain
template<typename T, typename PARTICLETYPE>
void addParticles( SuperParticleSystem<T,PARTICLETYPE>& supParticleSystem,
  IndicatorF3D<T>& ind, T partRho, T radius,
  const std::size_t noOfParticles, util::Randomizer<T>& randomizer )
{
  Vector<T,3> pos(0.);
  std::size_t noP = 0;
  while (noP<noOfParticles){
    //Generate randomized position
       pos[0] = ind.getMin()[0]
         + randomizer.generate() * (ind.getMax()[0] - ind.getMin()[0]);
       pos[1] = ind.getMin()[1]
         + randomizer.generate() * (ind.getMax()[1] - ind.getMin()[1]);
       pos[2] = ind.getMin()[2]
         + randomizer.generate() * (ind.getMax()[2] - ind.getMin()[2]);
    //Call sub-grid creator if in indicator domain
    bool inside[1] = { false };
    ind(inside, pos.data());
    if (inside[0]){
      creators::addSubgridSphere3D( supParticleSystem, pos, radius, partRho );
      noP++;
    }
  }
}


#ifdef PARALLEL_MODE_MPI

// Calculate capture and escape rate
template<typename T, typename PARTICLETYPE, bool verbose=true>
void captureStatistics( SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicatorOutput,
  std::size_t& noActive, std::size_t reprParticleID=0 )
{
  noActive = 0;
  std::size_t noFound = 0;
  std::size_t noExcaped = 0;
  communication::forParticlesInSuperParticleSystem(sParticleSystem,
    [&](Particle<T,PARTICLETYPE>& particle,
      ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC){
    //Check activity
    bool active = access::isActive(particle);
    if (active){
      noActive++;
      //Check whether representative particle and potentially print
      std::size_t globalId = access::getGlobalID(particle);
      if constexpr (verbose){
        if (globalId==reprParticleID){
          particles::io::printSubgridParticleInfo( particle );
        }
      }
    }
    //Check individual material vicinities
    bool outputVicinity = boundaries::checkMaterialVicinity(
      materialIndicatorOutput, particle );
    if (outputVicinity){ noExcaped++; }
    //Increase particle counter
    noFound++;
  });
  singleton::mpi().reduceAndBcast(noActive, MPI_SUM);
  singleton::mpi().reduceAndBcast(noFound, MPI_SUM);
  singleton::mpi().reduceAndBcast(noExcaped, MPI_SUM);

  //Calculate escape and capture rate
  T escapeRate = T(noExcaped)/T(noFound);
  T captureRate = 1.-escapeRate;

  //Output
  if constexpr (verbose){
    OstreamManager clout( std::cout, "captureStatistics" );
    clout << "globalNumOfParticles=" << noFound;
    clout << "; activeParticles=" << noActive;
    clout << "; escapedParticles=" << noExcaped << std::endl;
    clout << "captureRate=" << captureRate << "; escapeRate=" << escapeRate;
    clout << std::endl;
  }
}

#endif


} //namespace subgrid

} //namespace particles

} //namespace olb


#endif
