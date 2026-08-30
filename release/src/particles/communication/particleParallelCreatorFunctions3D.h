/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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


#ifndef PARTICLE_PARALLEL_CREATOR_FUNCTIONS_3D_H
#define PARTICLE_PARALLEL_CREATOR_FUNCTIONS_3D_H

#include <algorithm>

#include "particles/functions/lambdaLoops.h"
#include "particles/functions/particleCreatorHelperFunctions.h"

namespace olb {

namespace particles {

namespace creators {





//// RESOLVED

template<typename T, typename PARTICLETYPE>
std::set<int> prepareParallelResolvedParticle(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  Vector<T,PARTICLETYPE::d> position, T circumRadius,
  const T epsilon, int& globiCcentre, const T physDeltaX,
  const Vector<bool, PARTICLETYPE::d>& periodicity )
{
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Find globiC of particle centre
  const bool inDomain = particles::communication::getCuboid(
      cuboidDecomposition, periodicity, position, globiCcentre);
  //TODO: add error handling, when !inDomain
  if (!inDomain){
    std::cerr << "ERROR: Particle added outside domain!" << std::endl;
  }
  if constexpr (access::providesContactMaterial<PARTICLETYPE>()) {
    // without particle enlargement for contact, since it is 0 at the time of creation
    const T detectionDistance = T{0.5} * util::sqrt(PARTICLETYPE::d) * physDeltaX;
    const T halfEpsilon = T{0.5} * epsilon;
    circumRadius = circumRadius - halfEpsilon + util::max(halfEpsilon, detectionDistance);
  }
  sParticleSystem.updateOffsetFromCircumRadius(circumRadius);
  //Find globiCs touching the surface hull
  std::set<int> setOfICs = communication::getSurfaceTouchingICs( sParticleSystem,
    position, circumRadius, periodicity, globiCcentre );
  //Add globiCcentre to list
  setOfICs.insert(globiCcentre);
  //Return list of iCs
  return setOfICs;
}

//Assign particle to calling cuboid and return particle object
template<typename T, typename PARTICLETYPE>
Particle<T,PARTICLETYPE> assignParallelResolvedParticleToiC(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxSurface, std::size_t globID, int globiCcentre,
  const Vector<T,PARTICLETYPE::d>& position, T density=0.,
  const Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& angleInDegree
    = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> (0.),
  const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.))
{
  using namespace descriptors;
  //Add resolved object
  auto particle = creators::addResolvedObject( particleSystem, idxSurface,
    position, density, angleInDegree, velocity );
  //Add additional data due to parallization
  particle.template setField<PARALLELIZATION,ID>( globID );
  particle.template setField<PARALLELIZATION,IC>( globiCcentre );
  particle.template setField<SURFACE,SURFACE_ID>( idxSurface );
  //Return particle object
  return particle;
}

//Assigne particle to cuboid, if globiC listed in setOfICs and return, if assigned
//- effectively resulting in assigning to all cuboids touching
template<typename T, typename PARTICLETYPE>
bool assignParallelResolvedParticle(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::set<int>& setOfICs, int globiC,
  std::size_t idxSurface, std::size_t globID, std::size_t& localID, int globiCcentre,
  const Vector<T,PARTICLETYPE::d>& position, T density=0.,
  const Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& angleInDegree
    = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> (0.),
  const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.))
{
  //Check whether asignment needed for globiC
  bool assigned=false;
  if( setOfICs.find(globiC) != setOfICs.end() ) {
    //Attach particle to cuboid with globiC (and retrieve new particle object)
    auto particle = assignParallelResolvedParticleToiC( particleSystem, idxSurface,
      globID, globiCcentre, position, density, angleInDegree, velocity );
    //Declare assignement found
    assigned=true;
    //Retrieve localID
    localID = particle.getId();
  } //if( std::find() )
  //Return whether assigned
  return assigned;
}


//Add particle on residing iC with surface on all iCs (and return ParallelParticleLocator containing globiC, globID, localID)
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedObject(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::size_t idxSurface,
  const Vector<T,PARTICLETYPE::d>& position, T density=0.,
  const Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& angleInDegree
    = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> (0.),
  const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.),
  const Vector<bool, PARTICLETYPE::d>& periodicity = Vector<T,PARTICLETYPE::d>(false))
{
  //Retrieve circumRadius from smoothIndicator
  constexpr unsigned D = PARTICLETYPE::d;
  typedef SmoothIndicatorF<T,T,D,true> SIndicatorBaseType;
  auto bParticleSystems = sParticleSystem.getBlockParticleSystems();
  auto& particleSystem = *bParticleSystems[0]; //Use first randomly
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  T circumRadius = vectorOfIndicators[idxSurface]->getCircumRadius();
  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, vectorOfIndicators[idxSurface].get()->getEpsilon(), globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

/// Adds surface (SmoothIndicator) to all block particle systems
template<typename T, typename PARTICLETYPE, typename F>
std::size_t addSurface( SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
                        F createUniqueIndicatorPtr )
{
  typedef SmoothIndicatorF<T, T, PARTICLETYPE::d, true> SIndicatorBaseType;
  std::size_t idxSurface = 0;

  std::vector<olb::particles::ParticleSystem<T, PARTICLETYPE>*>&
      blockParticleSystems = sParticleSystem.getBlockParticleSystems();
  // Iteration over each element independent from cuboids to ensure that surfaces are added everywhere
  // as some blocks may not have a corresponding cuboid
  for (unsigned i = 0; i < blockParticleSystems.size(); ++i) {
    auto& particleSystem     = *blockParticleSystems[i];
    auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                               std::vector<std::unique_ptr<SIndicatorBaseType>>>();
    if(i > 0 && idxSurface != vectorOfIndicators.size()) {
      std::cerr << "ERROR: Number of particle surfaces don't match." << std::endl;
      assert(false);
    }
    idxSurface   = vectorOfIndicators.size();
    vectorOfIndicators.push_back( createUniqueIndicatorPtr() );
  }

  return idxSurface;
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedSphere3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T radius, T epsilon, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const Vector<bool, 3>& periodicity = Vector<bool,3>(false))
{
  //Get circumRadius and set angleInDegree
  T circumRadius = SmoothIndicatorSphere3D<T,T,true>(
    Vector<T,3> (0.), radius, epsilon).getCircumRadius();
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>
    angleInDegree(0.);

  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();

  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, epsilon, globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  std::size_t idxSurface = addSurface(sParticleSystem,
      [&](){
        return std::make_unique<SmoothIndicatorSphere3D<T, T, true>>(Vector<T,3> (0.), radius, epsilon );
      });

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Do for all iCs

    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedCuboid3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, const Vector<T,3>& extend,
  T epsilon, T density=0.,
  const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const Vector<bool, 3>& periodicity = Vector<bool,3>(false))
{
  //Get circumRadius
  T circumRadius = SmoothIndicatorCuboid3D<T,T,true>(
    extend[0], extend[1], extend[2], Vector<T,3>(0.), epsilon ).getCircumRadius();

  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();

  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, epsilon, globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  std::size_t idxSurface = addSurface(sParticleSystem,
      [&](){
        return std::make_unique<SmoothIndicatorCuboid3D<T, T, true>>(
            extend[0], extend[1], extend[2], Vector<T,3>(0.), epsilon );
      });

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Do for all iCs

    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedCylinder3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, const Vector<T,3>& normal,
  T height, T radius, T epsilon, T density=0.,
  const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const Vector<bool, 3>& periodicity = Vector<bool,3>(false))
{
  //Get circumRadius and set angleInDegree
  T circumRadius = SmoothIndicatorCylinder3D<T,T,true>(Vector<T,3>(0.), normal, radius, height, epsilon ).getCircumRadius();

  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();

  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, epsilon, globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  std::size_t idxSurface = addSurface(sParticleSystem,
      [&](){
        return std::make_unique<SmoothIndicatorCylinder3D<T, T, true>>(
            Vector<T,3>(0.), normal, radius, height, epsilon );
      });

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Do for all iCs

    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedArbitraryShape3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T latticeSpacing,
  std::shared_ptr<IndicatorF3D<T>> indPtr,
  T epsilon, T density=0.,
  const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const Vector<bool, 3>& periodicity = Vector<bool,3>(false))
{
  //Get circumRadius and set angleInDegree
  T circumRadius = SmoothIndicatorCustom3D<T,T,true>(
      latticeSpacing, indPtr, PhysR<T, 3>(T {}), epsilon, angleInDegree).getCircumRadius();

  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();

  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, epsilon, globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  std::size_t idxSurface = addSurface(sParticleSystem,
      [&](){
        return std::make_unique<SmoothIndicatorCustom3D<T, T, true>>(
          latticeSpacing, indPtr, PhysR<T, 3>(T {}), epsilon, angleInDegree);
      });

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Do for all iCs

    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addResolvedCircle2D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,2>& position, T radius, T epsilon, T density=0.,
  const Vector<T,2>& velocity = Vector<T,2> (0.),
  const Vector<bool, 2>& periodicity = Vector<bool,2>(false))
{
  //Get circumRadius and set angleInDegree
  T circumRadius = SmoothIndicatorCircle2D<T,T,true>(
    Vector<T,2> (0.), radius, epsilon).getCircumRadius();
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>
    angleInDegree(0.);

  //Get new global particle id
  std::size_t globID = sParticleSystem.getGlobID();

  //Prepare particle and get list of touching iCs
  int globiCcentre;
  std::set<int> setOfICs = prepareParallelResolvedParticle(sParticleSystem,
    position, circumRadius, epsilon, globiCcentre,
    sParticleSystem.getSuperStructure().getCuboidDecomposition().getMotherCuboid().getDeltaR(),
    periodicity);

  std::size_t idxSurface = addSurface(sParticleSystem,
      [&](){
        return std::make_unique<SmoothIndicatorCircle2D<T, T, true>>(Vector<T,2> (0.), radius, epsilon );
      });

  //Iterate over particle systems
  std::size_t localIDcentre=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){
    //Do for all iCs

    //Assigne resolved particle to iCs specified in setOfICs
    std::size_t localID=0;
    assignParallelResolvedParticle( particleSystem, setOfICs, globiC,
      idxSurface, globID, localID, globiCcentre,
      position, density, angleInDegree, velocity );
    //Set localID if centre
    if (globiC==globiCcentre){ localIDcentre=localID; }
  });

  return ParallelParticleLocator(globiCcentre, globID, localIDcentre);
}

//// SUBGRID

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addSubgridSphere3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T radius, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.) )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  if (auto iC = cuboidDecomposition.getC(position)) {
    globiCcentre = *iC;
  } else {
    std::cerr << "ERROR: Particle added outside domain!" << std::endl;
  }

  //Iterate over particle systems
  std::size_t localID=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      addSubgridObject( particleSystem,
        position, radius, density, velocity );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Set localID
      localID=particle.getId();
    }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localID);
}

//Add particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
ParallelParticleLocator addSubgridSphereWithSpecies3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T radius, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.),
  const int species=0 )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  if (auto iC = cuboidDecomposition.getC(position)) {
    globiCcentre = *iC;
  } else {
    std::cerr << "ERROR: Particle added outside domain!" << std::endl;
  }

  //Iterate over particle systems
  std::size_t localID=0; //Return localID of particle inside particleSystem holding centre
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      addSubgridObject( particleSystem,
        position, radius, density, velocity );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Add species type to particle
      if constexpr( access::providesSpecies<PARTICLETYPE>() ){
        particle.template setField<PHYSPROPERTIES,SPECIES>( species );
      }
      //Set localID
      localID=particle.getId();
    }
  });

  //Return globiC, globID, localiD
  return ParallelParticleLocator(globiCcentre, globID, localID);
}






} //namespace creators

} //namespace particles

} //namespace olb

#endif
