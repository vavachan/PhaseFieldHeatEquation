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


#ifndef SURFACE_FORCE_COMMUNICATION_H
#define SURFACE_FORCE_COMMUNICATION_H

//Increase verbosity for further development
//#define VERBOSE_COMMUNICATION

namespace olb {

namespace particles {

namespace communication {



#ifdef PARALLEL_MODE_MPI

//Communicate particle surface force to cores containing their centres
template<typename T, typename PARTICLETYPE, unsigned forceDataSize>
void communicateSurfaceForce(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::map<std::size_t, Vector<T,forceDataSize>>& globalIdDataMap
#ifdef PARALLEL_MODE_MPI
  , MPI_Comm surfaceForceComm
#endif
  )
{
  //TODO: add intra and inter core differentiation
  using namespace descriptors;
  //Retrieve dimensions
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Retrieve loadBalancer
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer = superStructure.getLoadBalancer();

  //Create FieldArrayD for globalID, force and torque
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,descriptors::FORCE> fieldForce(1);
  using FORCING_EVAL = typename PARTICLETYPE
                     ::template derivedField<descriptors::FORCING>;
  using TORQUE_EVAL = typename FORCING_EVAL
                     ::template derivedField<descriptors::TORQUE>;
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,TORQUE_EVAL> fieldTorque(1);
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,descriptors::ID> fieldID(1);

  //Create communicatables
  auto communicatableForce = ConcreteCommunicatable(fieldForce);
  auto communicatableTorque = ConcreteCommunicatable(fieldTorque);
  auto communicatableID = ConcreteCommunicatable(fieldID);
  //Retrieve serial size
  const std::vector<unsigned int> indices{0};
  auto serialSize = communicatableForce.size(indices)
                  + communicatableTorque.size(indices)
                  + communicatableID.size(indices);
  //Create rankDataMap
  std::multimap<int,std::unique_ptr<std::uint8_t[]>> rankDataMap;

  //Loop over valid particles in particle system
  communication::forParticlesInSuperParticleSystem<T,PARTICLETYPE,
    conditions::valid_particle_surfaces>( sParticleSystem,
    [&](Particle<T,PARTICLETYPE>& particle,
    ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC){
    //Retrieve globiC of particle center
    int globiCcentre = particle.template getField<PARALLELIZATION,IC>();
    int rankDest = loadBalancer.rank(globiCcentre);
    //Retrieve force and torque and set FieldArray
    auto force = particle.template getField<FORCING,FORCE>();
    auto torque = particle.template getField<FORCING,TORQUE>();
    auto globalID = particle.template getField<PARALLELIZATION,ID>();
#ifdef VERBOSE_COMMUNICATION
    std::cout << "EVALUATING: force=" << force << ", torque=" << torque << ", globalID=" << globalID << std::endl;
#endif
    fieldForce.setField(0, force);
    fieldTorque.setField(0, torque);
    fieldID.setField(0, globalID);
    //Create buffer
    std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSize]{ });
    std::uint8_t* bufferRaw = buffer.get();
    //Serialize force
    std::size_t serialIdx = communicatableForce.serialize(indices, bufferRaw);
    //Serialice torque with offset serialIdx
    serialIdx += communicatableTorque.serialize(indices, &bufferRaw[serialIdx]);
    //Serialice globalID with offset serialIdx
    communicatableID.serialize(indices, &bufferRaw[serialIdx]);
    //Add serialized data to rankDataMap
    rankDataMap.insert(std::make_pair(rankDest, std::move(buffer)));
  }); //forParticlesInSuperParticleSystem()

  //Retrieve rank of direct neighbours
  auto& listNeighbourRanks = sParticleSystem.getNeighbourRanks();

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  //Create ranDataMapSorted (WARNING: has to be existent until all data is received! c.f. #290)
  std::map<int, std::vector<std::uint8_t> > rankDataMapSorted;

  //Fill send buffer
  fillSendBuffer( rankDataMap, rankDataMapSorted, serialSize );

  //Send mapped data
  communication::sendMappedData( rankDataMapSorted,
      listNeighbourRanks, serialSize, surfaceForceComm, mpiNbHelper );

  //Receive data and iterate buffer
  communication::receiveAndExecuteForData( listNeighbourRanks, serialSize,
      surfaceForceComm, mpiNbHelper,
    [&](int rankOrig, std::uint8_t* bufferRaw){
    //Deserialize force
    std::size_t serialIdx = communicatableForce.deserialize(indices, bufferRaw);
    //Deserialize torque with offset serialIdx
    serialIdx += communicatableTorque.deserialize(indices, &bufferRaw[serialIdx]);
    //Deserialize globalID with offset serialIdx
    communicatableID.deserialize(indices, &bufferRaw[serialIdx]);
    //Retrieve vectors
    auto force = fieldForce.getField(0);
    auto torque = fieldTorque.getField(0);
    auto globalID = fieldID.getField(0);

#ifdef VERBOSE_COMMUNICATION
    int rank = singleton::mpi().getRank();
    std::cout << "Received (on rank=" << rank
              << " from rank=" << rankOrig
              << "): force=" << force
              << std::endl
              << "                  "
              << ", torque=" << torque
              << ", globalID=" << globalID
              << std::endl;
#endif

    //Create forceData Vector
    Vector<T,forceDataSize> forceData;
    for (unsigned iDim=0; iDim<D; ++iDim) {
      forceData[iDim] = force[iDim];
    }

    if constexpr (utilities::dimensions::convert<PARTICLETYPE::d>::rotation == 1 ) {
      forceData[D] = torque;
    }
    else {
      for (unsigned iDim=0; iDim<Drot; ++iDim) {
        forceData[iDim+D] = torque[iDim];
      }
    }

    //Add id-data-pair to map OR only add data
    auto mapPair = globalIdDataMap.insert( std::pair<std::size_t,
      Vector<T,forceDataSize>>(globalID, forceData) );
    if (mapPair.second==false){ //Checks previous existance
      globalIdDataMap[globalID] += forceData;
    }
  }); //receiveAndExecuteForData()
}

#endif


//Assign (communicated) surface force to particle
//- by matchin its globalID
template<typename T, typename PARTICLETYPE, unsigned forceDataSize>
void assignSurfaceForce(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::map<std::size_t, Vector<T,forceDataSize>>& globalIdDataMap )
{
  using namespace descriptors;
  //Retrieve dimension
  constexpr unsigned D = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Loop over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC,
      int globiC){
    //Loop over globalIdDataPairs
    for ( auto globalIdDataPair : globalIdDataMap ){
      //Loop over particles
      // -! Lambda-expression unavailable due to break command
      for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
        auto particle = particleSystem.get(iP);
        //Only consider valid particle centres
        if (conditions::valid_particle_centres::value(particle, globiC)){
          //Retrieve globalIDs
          std::size_t globalID = globalIdDataPair.first;
          std::size_t globalIDiP = particle.template getField<PARALLELIZATION,ID>();
          //Check whether globalID matches
          if (globalIDiP==globalID){
            //Retrieve force data
            auto forceData = globalIdDataPair.second;
            auto force = particle.template getField<FORCING,FORCE>();
            auto torque = particle.template getField<FORCING,TORQUE>();
            //Add force to particle
            for (unsigned iDim=0; iDim<D; ++iDim) {
              force[iDim] += forceData[iDim];
            }
            //Add torque to particle
            if constexpr (utilities::dimensions::convert<PARTICLETYPE::d>::rotation == 1 ) {
              torque += forceData[D];
            }
            else {
              for (unsigned iDim=0; iDim<Drot; ++iDim) {
                torque[iDim] += forceData[iDim+D];
              }
            }

#ifdef VERBOSE_COMMUNICATION
    int rank = singleton::mpi().getRank();
    std::cout << "Assigned (on rank=" << rank
              << "): force=" << force
              << std::endl
              << "                  "
              << ", torque=" << torque
              << ", globalID=" << globalID
              << std::endl;
#endif
            //Update particle
            particle.template setField<FORCING,FORCE>( force );
            particle.template setField<FORCING,TORQUE>( torque );
            //Break, since particle already found
            break;
          }
        }
      }
    }
  });
}


} //namespace communication

} //namespace particles

} //namespace olb


#endif
