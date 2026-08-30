/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-2023 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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



#ifndef PARTICLE_UTILITIES_H
#define PARTICLE_UTILITIES_H

#include <cassert>
#include <unordered_set>

#include "particles/communication/utilities.h"

namespace olb {

namespace particles {




//Get points on spherical hull
//TODO: has to be reworked thoroughly
// - separate position
// - pre-check necessary check by bounding box first
// - ordered by nested looping: x,y,z
template<typename T>
std::array<Vector<T,3>,26> discretePointsOnSphericalHull(
  Vector<T,3> position, T radius )
{
  //TODO: check correctness
  T distA = radius;
  T distB = (1/std::sqrt(2))*distA;
  T distC = (1/std::sqrt(2))*distB;

  //Point on spherical hull
  std::array<Vector<T,3>,26> pointsOnSphericalHull = {
    Vector<T,3>(position[0]-distC, position[1]-distC, position[2]-distC),
    Vector<T,3>(position[0],       position[1]-distB, position[2]-distB),
    Vector<T,3>(position[0]+distC, position[1]-distC, position[2]-distC),

    Vector<T,3>(position[0]-distB, position[1],       position[2]-distB),
    Vector<T,3>(position[0],       position[1],       position[2]-distA), //bottom
    Vector<T,3>(position[0]+distB, position[1],       position[2]-distB),

    Vector<T,3>(position[0]-distC, position[1]+distC, position[2]-distC),
    Vector<T,3>(position[0],       position[1]+distB, position[2]-distB),
    Vector<T,3>(position[0]+distC, position[1]+distC, position[2]-distC),


    Vector<T,3>(position[0]-distB, position[1]-distB, position[2]      ),
    Vector<T,3>(position[0],       position[1]-distA, position[2]      ), //right
    Vector<T,3>(position[0]+distB, position[1]-distB, position[2]      ),

    Vector<T,3>(position[0]-distA, position[1],       position[2]      ), //back
    // Centre excluded here!
    Vector<T,3>(position[0]+distA, position[1],       position[2]      ), //front

    Vector<T,3>(position[0]-distB, position[1]+distB, position[2]      ),
    Vector<T,3>(position[0],       position[1]+distA, position[2]      ), //left
    Vector<T,3>(position[0]+distB, position[1]+distB, position[2]      ),


    Vector<T,3>(position[0]-distC, position[1]-distC, position[2]+distC),
    Vector<T,3>(position[0],       position[1]-distB, position[2]+distB),
    Vector<T,3>(position[0]+distC, position[1]-distC, position[2]+distC),

    Vector<T,3>(position[0]-distB, position[1],       position[2]+distB),
    Vector<T,3>(position[0],       position[1],       position[2]+distA), //top
    Vector<T,3>(position[0]+distB, position[1],       position[2]+distB),

    Vector<T,3>(position[0]-distC, position[1]+distC, position[2]+distC),
    Vector<T,3>(position[0],       position[1]+distB, position[2]+distB),
    Vector<T,3>(position[0]+distC, position[1]+distC, position[2]+distC)
  };

  return pointsOnSphericalHull;
}

template<typename T>
std::array<Vector<T,2>,8> discretePointsOnSphericalHull(
  Vector<T,2> position, T radius )
{
  //TODO: check correctness
  T distA = radius;
  T distB = (1/std::sqrt(2))*distA;

  //Point on spherical hull
  std::array<Vector<T,2>,8> pointsOnSphericalHull = {
    Vector<T,2>(position[0]-distB, position[1]-distB),
    Vector<T,2>(position[0],       position[1]-distA),  //right
    Vector<T,2>(position[0]+distB, position[1]-distB),

    Vector<T,2>(position[0]-distA, position[1]),        //back
    // Centre excluded here!
    Vector<T,2>(position[0]+distA, position[1]),        //front

    Vector<T,2>(position[0]-distB, position[1]+distB),
    Vector<T,2>(position[0],       position[1]+distA),  //left
    Vector<T,2>(position[0]+distB, position[1]+distB)
  };

  return pointsOnSphericalHull;
}

struct discrete_points_on_hull{
  template<typename T, unsigned D>
  static constexpr auto calculate( Vector<T,D> position, T radius )
  {
    if constexpr(D==3 || D==2){
      return discretePointsOnSphericalHull( position, radius );
    } else {
      std::cerr << "ERROR: Only 2D and 3D supported!" << std::endl;
      return false;
    }
  }
  template<typename T, typename PARTICLETYPE>
  static constexpr auto calculate( Particle<T,PARTICLETYPE>& particle )
  {
    auto position = access::getPosition( particle );
    T radius = access::getRadius( particle );
    return calculate( position, radius );
  }
};

template<typename T, typename PARTICLETYPE>
void purgeInvalidParticles(
  XParticleSystem<T,PARTICLETYPE>& xParticleSystem )
{
  using PCONDITION = conditions::invalid_particles;
  if constexpr (access::providesParallelization<PARTICLETYPE>()){
    auto& sParticleSystem = xParticleSystem;
    //Iterate over particle systems
    communication::forSystemsInSuperParticleSystem( sParticleSystem,
      [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
      //Iterate over particles
      std::size_t iP=0;
      while(iP < particleSystem.size()){
        auto particle = particleSystem.get(iP);
        ++iP;
        //Execute F when particle meets condition
        doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
          [&](Particle<T,PARTICLETYPE> particle){
          //Erase particle
          auto& container = particleSystem.get();
          container.erase( particle.getId() );
          --iP; //reset counter
        },globiC);
      }
    });
  } else {
    auto& particleSystem = xParticleSystem;
    std::size_t iP=0;
    while(iP < particleSystem.size()){
      auto particle = particleSystem.get(iP);
      ++iP;
      //Execute F when particle meets condition
      doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
        [&](Particle<T,PARTICLETYPE> particle){
        //Erase particle
        auto& container = particleSystem.get();
        container.erase( particle.getId() );
        --iP; //reset counter
      });
    }
  }
}

template<typename T, typename PARTICLETYPE, std::size_t selectedID, typename F>
void doForParticleMatchingID( XParticleSystem<T,PARTICLETYPE>& xParticleSystem, F f ){
  const bool hasFieldValid = access::providesValid<PARTICLETYPE>();
  using PCOND = std::conditional_t<
    hasFieldValid,
    conditions::template valid_particle_matching_ID<selectedID>,
    conditions::template particle_matching_ID<selectedID>
  >;
  //Loop over particle and find the one matching iD
  forParticlesInXParticleSystem<T,PARTICLETYPE,PCOND>( xParticleSystem, f );
}

//Search particle in ParticleSystem by globalParticleID and return localParticleID if found
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particles>
bool searchParticleLocally( ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t globalIDrequested, std::size_t& localParticleID )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<PARALLELIZATION,ID>(), "Field PARALLELIZATION:ID has to be provided");
  //Initialize return quantities
  localParticleID = 0;
  bool found = false;
  //Loop over particles
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    if (!found){
      auto particle = particleSystem.get(iP);
      //Execute F when particle meets condition
      doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
        [&](Particle<T,PARTICLETYPE> particle){
        auto globalID = particle.template getField<PARALLELIZATION,ID>();
        if (globalID==globalIDrequested){
          localParticleID = particle.getId();
          found=true;
        }
      });
    } else {
      break;
    }
  }
  return found;
}



namespace communication {



/// Search particle in SuperParticleSystem by globalParticleID and return globiC if found
/// - WARNING: Sync is enabled by default, to avoid errors. As it might usually not be necessary though,
///   communication can be disables to increase performance.
template<typename T, typename PARTICLETYPE, typename PCONDITION=conditions::valid_particle_centres, bool sync=true>
bool searchParticleGlobally( SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::size_t globalIDrequested, std::size_t& localParticleID, int& globiC )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<PARALLELIZATION,ID>(), "Field PARALLELIZATION:ID has to be provided");
  //Retrieve load balancer
  auto& loadBalancer = sParticleSystem.getSuperStructure().getLoadBalancer();
  //Initialize return quantities
  localParticleID = 0;
  globiC = 0;
  bool found = false;
  //Loop over particle systems
  for (int iC=0; iC<loadBalancer.size(); ++iC){
    if (!found){
      int globiC = loadBalancer.glob(iC);
      //Retrieve container
      auto bParticleSystems = sParticleSystem.getBlockParticleSystems();
      auto& particleSystem = *bParticleSystems[iC];
      //Search particle on particleSystem
      //NOTE: Can not be substituted with searchParticle( ParticleSystem ..) due to globiC
      for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
        if (!found){
          auto particle = particleSystem.get(iP);
          //Execute F when particle meets condition
          doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
            [&](Particle<T,PARTICLETYPE> particle){
            auto globalID = particle.template getField<PARALLELIZATION,ID>();
            if (globalID==globalIDrequested){
              localParticleID = particle.getId();
              found=true;
            }
          }, globiC);
        } else {
          break;
        }
      } //for (std::size_t iP=0; iP<particleSystem.size(); ++iP)
    } else {
      break;
    }
  } //for (int iC=0; iC<loadBalancer.size(); ++iC)
#ifdef PARALLEL_MODE_MPI
  if constexpr(sync){
    singleton::mpi().reduceAndBcast(found, MPI_LOR);
    if (found){
      singleton::mpi().reduceAndBcast(globiC, MPI_MAX);
      singleton::mpi().reduceAndBcast(localParticleID, MPI_MAX);
    }
  }
#endif
  return found;
}


/**
 * Get a set of surface touching iCs (that are not globiC)
 * Allows to run an optional function per unique globiC
 */
template<typename T, unsigned D, typename F=std::function<void(int)>,
  bool domainWarning=false, bool checkDiscretePoints=false>
std::set<int> getSurfaceTouchingICs(
  CuboidDecomposition<T,D>& cuboidDecomposition,
  Vector<T,D> position, T circumRadius,
  const Vector<bool, D>& periodicity,
  int globiC, F f=[](int){} )
{
  // This option is prone to errors
  if constexpr(checkDiscretePoints) {
    //Set of destination iCs preventing sending to same iC multiple times
    std::set<int> iCs;
    //Retrieve points on hull
    auto pointsOnHull = discrete_points_on_hull::calculate( position, circumRadius );
    //Iterate over points on hull
    for (const PhysR<T,D>& point : pointsOnHull){
      //Retrieve residence iC
      int globiConHull=0;
      const bool cuboidFound = getCuboid(cuboidDecomposition,
          periodicity, point, globiConHull);
      if constexpr(domainWarning){
        if (!cuboidFound){
          std::cerr << "Warning [during getSurfaceTouchingICs]: Cuboid not found for point "
                    << point << std::endl;
        }
      }
      //Check whether point has left responsible ic
      if (cuboidFound && globiConHull!=globiC){
        //Add destination ic to list (to avoid resending)
        if(iCs.insert(globiConHull).second) {
          f(globiConHull);
        }
      } //if (globiConHull!=globiC)
    } //for (auto point : pointsOnHull)
    return iCs;
  }
  else {
    const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR();
    std::vector<int> tmp = communication::getNeighborCuboids<T,D,false>(
          cuboidDecomposition, util::ceil(circumRadius/physDeltaX), globiC);
    std::set<int> iCs;
    for(const int entry : tmp) {
      if( globiC != entry ) {
        iCs.insert(entry);
        f(entry);
      }
    }
    return iCs;
  }
  __builtin_unreachable();
}

template<typename T, typename PARTICLETYPE, typename F=std::function<void(int)>>
std::set<int> getSurfaceTouchingICs(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  Vector<T,PARTICLETYPE::d> position, T circumRadius,
  const Vector<bool, PARTICLETYPE::d>& periodicity,
  int globiC, F f=[](int){} )
{
  const auto& neighbors = sParticleSystem.getCuboidNeighborhood()[globiC];
  for(int neighbor: neighbors) {
    f(neighbor);
  }
  return neighbors;
}


//Get list of surface touching iCs (that are not globiC)
template<typename T, unsigned D, typename F=std::function<void(int)>,
  bool domainWarning=false>
std::vector<int> getVectorOfSurfaceTouchingICs(
  CuboidDecomposition<T,D>& cuboidDecomposition,
  Vector<T,D> position, T circumRadius,
  const Vector<bool, D>& periodicity,
  int globiC, F f=[](int){} )
{
  const std::set<int> touchedICs{getSurfaceTouchingICs<T,D,F,domainWarning>(
      cuboidDecomposition, position, circumRadius, periodicity, globiC, f)};
  std::vector<int> listOfICs;
  listOfICs.insert(listOfICs.end(), touchedICs.begin(), touchedICs.end());
  return listOfICs;
}


} //namespace communication


namespace sorting {


//Sort particleSystem by provided nested fields
//TODO: include derived field function
//TODO: for now only boolean, numeric comparison should be provided as well
//TODO: possibly generalize to work with all containers or objects, such as blockLattice
//       - this should be possible at least, when particleSystem and Lattice have the same
//         structure/base type
template<typename T, typename PARTICLETYPE, typename ...NESTED_FIELDS>
std::size_t partitionParticleSystem( ParticleSystem<T,PARTICLETYPE>& particleSystem ){
  using namespace descriptors;
  //Find first occurence of element belonging to partition B
  std::size_t iPfirstB=particleSystem.size(); //Fallback (simple to check)
  bool valB = true;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    auto particle = particleSystem.get(iP);
    bool active = particle.template getField<NESTED_FIELDS...>();
    if (active==valB){
      iPfirstB=iP;
      break;
    }
  }
  //Find succeeding elements belonging to partition A and move them to othe beginning
  for (std::size_t iP=iPfirstB; iP<particleSystem.size(); ++iP){
    auto particle = particleSystem.get(iP);
    bool active = particle.template getField<NESTED_FIELDS...>();
    if (active!=valB){
      particleSystem.swapParticles(iP,iPfirstB);
      ++iPfirstB;
    }
  }
  //Return splitpoint between partitions
  return iPfirstB;
}

} //namespace sorting

} //namespace particles

} //namespace olb


#endif
