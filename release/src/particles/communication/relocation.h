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


#ifndef PARTICLE_RELOCATION_H
#define PARTICLE_RELOCATION_H

#include "particles/contact/contactFunctions.h"
#include "particles/functions/particleDynamicsFunctions.h"

//Increase verbosity for further development
//#define VERBOSE_COMMUNICATION

namespace olb {

namespace particles {

namespace communication {

//Update surface pointer by considering SURFACE_ID
template<typename T, typename PARTICLETYPE>
void updateSurfacePtr(
  ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t iP )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<SURFACE,SURFACE_ID>(), "Field SURFACE_ID has to be provided");
  //Retrieve particle
  auto particle = particleSystem.get(iP);
  //Retrieve surface id
  std::size_t idxSurface = particle.template getField<SURFACE,SURFACE_ID>();
  //Retrieve surface ptr
  using SIndicatorBaseType = SmoothIndicatorF<T,T,PARTICLETYPE::d,true>;
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
    std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  auto sIndicatorPtr = vectorOfIndicators.at( idxSurface ).get();
  //Update particle
  particle.template setField<SURFACE,SINDICATOR>( sIndicatorPtr );
}

//Attach serialized particle
//TODO: Change dynamics treatement, as those cannot be serialized
template<typename T, typename PARTICLETYPE>
std::size_t appendSerializedParticle(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::uint8_t rawBuffer[])
{
  //Extend particle container
  auto& particleContainer = particleSystem.get();
  particleContainer.push_back();
  //Retrieve new particle
  std::size_t idxParticleNew = particleSystem.size()-1;
  auto particleNew = particleSystem.get(idxParticleNew);
  //Deserialise particle
  int globiCdest;
  deserializeIcDestAndParticle(globiCdest,particleNew,rawBuffer);
  //Update surface pointer
  if constexpr( access::providesSurface<PARTICLETYPE>() ){
    updateSurfacePtr( particleSystem, idxParticleNew );
  }
  return particleNew.getId();
}

//Insert serialized particle at specific idxParticle
template<typename T, typename PARTICLETYPE>
std::size_t insertSerializedParticle(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  std::uint8_t rawBuffer[])
{
  using namespace descriptors;
  //Retrieve particle
  auto particle = particleSystem.get(idxParticle);

#ifdef VERBOSE_COMMUNICATION
  auto globiCcentreLast = particle.template getField<PARALLELIZATION,IC>();
#endif

  //Deserialise particle
  int globiCdest;
  deserializeIcDestAndParticle(globiCdest,particle,rawBuffer);

#ifdef VERBOSE_COMMUNICATION
  auto globiCcentre = particle.template getField<PARALLELIZATION,IC>();

  if(globiCcentre!=globiCcentreLast){
    int rank = singleton::mpi().getRank();
    std::cout << std::endl << "[Particle] "
              << "CENTRE SWITCH("
              << globiCcentreLast << "->" << globiCcentre << ")"
              << " at receiving at iC=" << "X"
              << "(rank=" << rank << ")" << std::endl;
  }
#endif


  //Update surface pointer
  if constexpr( access::providesSurface<PARTICLETYPE>() ){
    updateSurfacePtr( particleSystem, idxParticle );
  }
  return particle.getId();
}




//Attach serializes particle to new ParticleSystem
template<typename T, typename PARTICLETYPE>
std::size_t attachSerializedParticle(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::uint8_t* rawBuffer )
{
  using namespace descriptors;

  std::size_t particleID;
  //Check whether resolved (or subgrid)
  //- subgrid has to be appended always, as valid particle data will not
  //  be present on new core due to previous direct invalidation
  //- TODO: Check consistency, as for now, invalidated slots are considered
  //  for resolved, but not for sub grid particles
  if constexpr( access::providesSurface<PARTICLETYPE>() ){
    //Attach serialized particle (either at end or at existing index)
    if (particleSystem.size()==0){
      particleID = appendSerializedParticle( particleSystem, rawBuffer );
    } else {
      //Create temporary particle system to retrieve global id
      ParticleSystem<T,PARTICLETYPE> particleSystemTmp(1);
      auto particleTmp = particleSystemTmp.get(0);
      int globiCdummy;
      deserializeIcDestAndParticle(globiCdummy,particleTmp,rawBuffer);
      std::size_t globalID = particleTmp.template getField<PARALLELIZATION,ID>();
      //Check whether particle already existing
      // -! Lambda-expression unavailable due to break command
      bool found=false;
      //Loop over particles
      //- including invalid slots, as they are revalidated
      //  anyways when inserting particle data)
      for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
        auto particle = particleSystem.get(iP);
        auto globalIDiP = particle.template getField<PARALLELIZATION,ID>();
        if (globalIDiP==globalID){
          particleID = insertSerializedParticle( particleSystem, iP, rawBuffer );
          found=true;
          break;
        }
      }
      if (!found){
        particleID = appendSerializedParticle( particleSystem, rawBuffer );
      }
    }
  } else { //if constexpr ( access::providesSurface<PARTICLETYPE>() )
    //Attach serialized particle
    particleID = appendSerializedParticle( particleSystem, rawBuffer );
  }
  return particleID;
}


//Prepare inter relocation or intra relocation
//-inter: separate iCs on different cores
//-intra: separate iCs on same core
template<typename T, typename PARTICLETYPE>
void prepareRelocation(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  Particle<T,PARTICLETYPE>& particle,
  int globiC, int globiCdest,
  std::multimap<int,std::unique_ptr<std::uint8_t[]>>& rankDataMapRelocationInter,
  std::vector<std::unique_ptr<std::uint8_t[]>>& dataListRelocationIntra,
  std::size_t serialSize )
{
  //Retrieve loadBalancer
  auto& loadBalancer = sParticleSystem.getSuperStructure().getLoadBalancer();
  //Create buffer pointer
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSize]{ });

  //Serialize particle content into buffer
  serializeIcDestAndParticle( globiCdest, particle, buffer.get() );

  //Check, whether iCs are on same core
  int rank = singleton::mpi().getRank();
  int rankDest = loadBalancer.rank(globiCdest);
  if (rank != rankDest){
    //Add rank-buffer pair to rankDataMapRelocationInter
    rankDataMapRelocationInter.insert(std::make_pair(rankDest, std::move(buffer)));
  } else {
    //If on same core location can be done without mpi
    dataListRelocationIntra.push_back(std::move(buffer));
  } //else (rank != rankDest)

#ifdef VERBOSE_COMMUNICATION
  auto globalID = particle.template getField<
    descriptors::PARALLELIZATION,descriptors::ID>();
  std::cout << std::endl << "[Particle] SENDING";
  if (rank != rankDest){ std::cout << "(inter)"; } else { std::cout << "(intra)"; }
  std::cout << " Particle"
            << "(globalID=" << globalID << ") "
            << "from iC" << globiC << "(rank=" << rank << ") to iC"
            << globiCdest << "(rank=" << rankDest << ")" << std::endl;
  particle.template print<true>();
#endif
}


//Check relocation of discrete surface elements (resolved only)
template<typename T, typename PARTICLETYPE>
void checkSurfaceRelocation(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const T physDeltaX,
  Particle<T,PARTICLETYPE>& particle, int globiC,
  std::multimap<int,std::unique_ptr<std::uint8_t[]>>& rankDataMapRelocationInter,
  std::vector<std::unique_ptr<std::uint8_t[]>>& dataListRelocationIntra,
  std::size_t serialSize, const Vector<bool, PARTICLETYPE::d>& periodicity )
{
  //Retrieve cuboid geometry
  /* auto& superStructure = sParticleSystem.getSuperStructure(); */
  /* auto& cuboidDecomposition = superStructure.getCuboidDecomposition(); */
  //Retrieve position and circumRadius
  auto position = access::getPosition( particle );
  const T circumRadius = contact::evalCircumRadius(particle, physDeltaX);
  const int globiCcentre = particle.template getField<descriptors::PARALLELIZATION,descriptors::IC>();
  //Get unique iCs on surface and prepare relocation
  getSurfaceTouchingICs(sParticleSystem, position, circumRadius, periodicity, globiCcentre,
      [&](const int globiCdest){
        prepareRelocation( sParticleSystem, particle, globiC, globiCdest,
          rankDataMapRelocationInter, dataListRelocationIntra, serialSize);
      });
}

//Check whether relocation necessary
//Also differentiate between inter- and intra-core relocation
//Also check domain exit
template<typename T, typename PARTICLETYPE>
void checkRelocation(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const T physDeltaX,
  std::multimap<int,std::unique_ptr<std::uint8_t[]>>& rankDataMapRelocationInter,
  std::vector<std::unique_ptr<std::uint8_t[]>>& dataListRelocationIntra,
  std::size_t serialSize, const Vector<bool, PARTICLETYPE::d>& periodicity )
{
  using namespace descriptors;
  using namespace particles::access;
  constexpr unsigned D = PARTICLETYPE::d;

  //Retrieve load balancer and cuboid geometry
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& cuboidDecomposition = superStructure.getCuboidDecomposition();

  //Loop over valid particles
  forParticlesInSuperParticleSystem( sParticleSystem,
    [&](Particle<T,PARTICLETYPE>& particle,
    ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC){

    //Retrieve centre iC in previous timestep.
    //- used to detect iC change
    //- unly updated for resolved particles
    int globiCcentreLast= -1;
    if constexpr( providesSurface<PARTICLETYPE>() ){
      globiCcentreLast = particle.template getField<PARALLELIZATION,IC>();
    }
    // Retrieve and update particle position
    // if a periodic setup is used and the particle is out of bounds
    const PhysR<T,D> position = applyPeriodicityToPosition(cuboidDecomposition,
        periodicity, particle.template getField<GENERAL,POSITION>());
    //Find particle system by position
    int globiCcentre = -1;
    // getC doesn't treat periodicity! We can only use it because we updated it before
    auto iC = cuboidDecomposition.getC(position);
    const bool particleCentreInDomain = iC.operator bool();
    if (particleCentreInDomain) {
      globiCcentre = *iC;
    }
    /*
    // Retrieve particle position
    PhysR<T,D> position = particle.template getField<GENERAL,POSITION>();
    //Find particle system by position
    int globiCcentre = -1;
    const bool particleCentreInDomain = getCuboid(
        cuboidDecomposition, periodicity, position, globiCcentre);
    */
    // Replace the position with the updated position
    // if no periodicity is used or the particle isn't out of bounds nothing happens
    // (must be here, because it can be that the particle remains on the same block after crossing the boundary)
    // TODO: Add constexpr check if the setup is periodic or not before continuously overwriting the position here
    particle.template setField<GENERAL, POSITION>(position);
    /*if(isPeriodic(periodicity)) {
      const PhysR<T,D> min = getCuboidMin<T,D>(cuboidDecomposition);
      const PhysR<T,D> max = getCuboidMax<T,D>(cuboidDecomposition, min);

      for(unsigned iD=0; iD<D; ++iD) {
        position[iD] = applyPeriodicityToPosition(
            periodicity[iD], position[iD], max[iD], min[iD]);
      }
      particle.template setField<GENERAL, POSITION>(position);
    }*/

    //Update iC, if resolved
    if constexpr( providesSurface<PARTICLETYPE>() ){
      particle.template setField<PARALLELIZATION,IC>(globiCcentre);
    }

    //Check whether particle has lost the global domain
    if (particleCentreInDomain){
      //Check whether particle centre resides in current iC
      //If yes:
      //- subgrid: no data to be sent
      //- resolved: check surface
      if (globiCcentre == globiC){
        //Check, whether surface parts have to be sent (resolved only)
        if constexpr( providesSurface<PARTICLETYPE>() ){
          checkSurfaceRelocation( sParticleSystem, physDeltaX, particle, globiC,
            rankDataMapRelocationInter, dataListRelocationIntra, serialSize, periodicity );
        }
      //If particle centre does not reside in iC
      } else { // if (globiCcentre == globiC)
        //Differentiate between resolved and subgrid treatment
        if constexpr( providesSurface<PARTICLETYPE>() ){

          //Detect center switch
          //- to send data one last time
          //- necessary, since new responsible iC can only send at next timestep
          if(globiCcentre!=globiCcentreLast){
#ifdef VERBOSE_COMMUNICATION
            int rank = singleton::mpi().getRank();
            std::cout << std::endl << "[Particle] "
                      << "CENTRE SWITCH("
                      << globiCcentreLast << "->" << globiCcentre << ")"
                      << " at sending from iC=" << globiC
                      << "(rank=" << rank << ")" << std::endl;
#endif
            //Prepare inter relocation or intra relocation
            //-here: to new responsible ic
            prepareRelocation( sParticleSystem,
              particle, globiC, globiCcentre,
              rankDataMapRelocationInter, dataListRelocationIntra, serialSize);

            if constexpr( providesSurface<PARTICLETYPE>() ){
              //Sending to self, to avoid invalidation on receiving
              //- necessary, as responsibility changed
              prepareRelocation( sParticleSystem,
                particle, globiCcentre, globiC,
                rankDataMapRelocationInter, dataListRelocationIntra, serialSize);

              //Check, whether surface parts have to be sent (resolved only)
              checkSurfaceRelocation( sParticleSystem, physDeltaX, particle, globiC,
                rankDataMapRelocationInter, dataListRelocationIntra, serialSize, periodicity );
            }

          } //if(globiCcentre!=globiCcentreLast)
        } else { //if constexpr( providesSurface<PARTICLETYPE>() )
          //Prepare inter relocation or intra relocation
          prepareRelocation( sParticleSystem,
            particle, globiC, globiCcentre,
            rankDataMapRelocationInter, dataListRelocationIntra, serialSize);
          //Immediate invalidation (subgrid only)
          particle.template setField<GENERAL,INVALID>(true);
        }
      } // else (globiCcentre == globiC)
    } else { //particleCentreInDomain
      //Invalidate particle
      particle.template setField<GENERAL,INVALID>(true);
      //Output warning
      int rank = singleton::mpi().getRank();
      auto globalID = particle.template getField<PARALLELIZATION,ID>();
      std::cout << ">>>> WARNING: Particle " << globalID
           << " lost domain at iC="<< globiC
           << " (rank=" << rank << ") <<<<" << std::endl;
    }
  }); // forParticlesInSuperParticleSystem
}

//Fill send buffer with data in rankDataMap (particle agnostic)
void fillSendBuffer(
  std::multimap<int,std::unique_ptr<std::uint8_t[]>>& rankDataMap,
  std::map<int, std::vector<std::uint8_t>>& rankDataMapSorted,
  std::size_t serialSize )
{
  //Iterate over rank data pairs
  for (auto& rankDataPair : rankDataMap) {
    //Decompose pair
    int rankDest = rankDataPair.first;
    std::uint8_t* buffer = rankDataPair.second.get();
    //Write data to sorted map
    // - creates vector via 3-arg initializer list
    rankDataMapSorted[rankDest].insert(
      rankDataMapSorted[rankDest].end(),
      buffer,
      buffer + serialSize );
  }
}

#ifdef PARALLEL_MODE_MPI
//Send data for inter-core relocation (particle agnostic)
void sendMappedData(
  std::map<int, std::vector<std::uint8_t> >& rankDataMapSorted,
  const std::set<int>& availableRanks,
  const std::size_t serialSize,
  MPI_Comm Comm,
  singleton::MpiNonBlockingHelper& mpiNbHelper )
{
  //Create map with vector holding individual data (of type std::uint8_t)
  //Iterate over available ranks and retrieve total number of ranks taking part in communication
  int noRelevantRanks = availableRanks.size();
  //Check whether any rank is taking part in communication
  //Allocate mpi requests for relevant ranks
  mpiNbHelper.allocate(noRelevantRanks);
  //Iterate over ranks with index iRank
  int iRank = 0;
  for (int rankEval : availableRanks) {
    //Send data via mpi
    singleton::mpi().iSend<std::uint8_t>(
      rankDataMapSorted[rankEval].data(),          //buffer: pointer to rank section in rankDataMapSorted
      rankDataMapSorted[rankEval].size(),          //count: number of bytes (per rank)
      rankEval,                                    //dest: destination rank
      &mpiNbHelper.get_mpiRequest()[iRank],      //request: request for iRank (iterator for relevant ranks)
      1,                                           //tag
      Comm);                                       //MPI_Comm scope
    iRank += 1;
  }
}
#endif

//Check wheter invalidation needed on receival
//General Idea:
// - Invalidate surface parts, when no information was received from
//   responsible iC
template<typename T, typename PARTICLETYPE>
void checkInvalidationOnReceival(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::vector<std::vector<std::size_t>>& receivedLocalIDs,
  const Vector<bool,PARTICLETYPE::d>& periodicity )
{
  using namespace descriptors;

  //Retrieve load balancer
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer = superStructure.getLoadBalancer();

  //Loop over valid surfaces and check, whether any data was received
  forParticlesInSuperParticleSystem<T,PARTICLETYPE,conditions::valid_particle_surfaces>(
    sParticleSystem,
    [&](Particle<T,PARTICLETYPE>& particle,
      ParticleSystem<T,PARTICLETYPE>& particleSystem, int globiC){

    //Retrieve local iC
    int lociC = loadBalancer.loc(globiC);

    //Retrieve local particle iD
    std::size_t locId = particle.getId();

#ifdef VERBOSE_COMMUNICATION
    std::cout << "ReceivedLocalIDs:";
    for (int i=0; i<receivedLocalIDs[lociC].size(); ++i){
      auto idPrint = receivedLocalIDs[lociC][i];
      std::cout << idPrint;
      if (i<receivedLocalIDs.size()){ std::cout << ","; }
    }
    std::cout << std::endl;
#endif

    //Check wether particle belongs to those having received data (here, if not found)
    if( std::find(receivedLocalIDs[lociC].begin(),receivedLocalIDs[lociC].end(),
      locId) == receivedLocalIDs[lociC].end()){
      //Invalidate particle
      particle.template setField<GENERAL,INVALID>(true);

#ifdef VERBOSE_COMMUNICATION
      int rank = singleton::mpi().getRank();
      std::size_t globalID = particle.template getField<PARALLELIZATION,ID>();
      std::cout << std::endl << "[Particle] "
                << "INVALIDATING Particle"
                << "(globalID=" << globalID << ") "
                << "on iC" << globiC << "(rank=" << rank << ")" << std::endl;
    particle.template print<true>();
#endif

    } // if( std::find() )
    /* else { */
    /*   // Check if the current block still touches the surface, if not then invalidate the particle */
    /*   auto& cuboidDecomposition = superStructure.getCuboidDecomposition(); */
    /*   const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR(); */
    /*   const T circumRadius = contact::evalCircumRadius(particle, physDeltaX); */
    /*   const int globiCcenter = particle.template getField<PARALLELIZATION, IC>(); */
    /*   bool touchesSurface = false; */
    /*   communication::getSurfaceTouchingICs( */
    /*       sParticleSystem, access::getPosition(particle), circumRadius, periodicity, globiCcenter, */
    /*       [&](int surfaceTouchingGlobiC){ */
    /*         if(surfaceTouchingGlobiC == globiC) { */
    /*           touchesSurface = true; */
    /*         } */
    /*       }); */
    /*   if(!touchesSurface) { */
    /*     particle.template setField<GENERAL,INVALID>(true); */
    /*   } */
    /* } */
  }); //forParticlesInSuperParticleSystem<T,PARTICLETYPE,valid_particle_surfaces>
}


//Assign particle to destination iC
//- Both inter and intra relocation
template<typename T, typename PARTICLETYPE>
void assignParticleToIC(
  SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
  std::vector<std::vector<std::size_t>>& receivedLocalIDs,
  int rankOrig,
  std::uint8_t* rawBuffer )
{
  using namespace descriptors;

  //Deserialize globiCdest
  int globiCdest;
  deserializeIcDest<T,PARTICLETYPE>( globiCdest,rawBuffer );

  //Loop over ParticleSystems in SuperParticleSystem
  forSystemsInSuperParticleSystem( sParticleSystem,
    [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC,
      int globiC){

    if (globiC == globiCdest){

      //Attach serialized particle to new particleSystem
      [[maybe_unused]] std::size_t particleID =
        attachSerializedParticle(
        particleSystem, rawBuffer);

      //Add local particle iD to list (resolved only)
      if constexpr( access::providesSurface<PARTICLETYPE>() ){
        receivedLocalIDs[iC].push_back(particleID);
      }

#ifdef VERBOSE_COMMUNICATION
      int rank = singleton::mpi().getRank();
      auto particle = particleSystem.get(particleID);
      std::size_t globalID = particle.template getField<PARALLELIZATION,ID>();
      std::cout << std::endl << "[Particle] RECEIVING";
      if (rankOrig != rank){ std::cout << "(inter)"; } else { std::cout << "(intra)"; }
      std::cout << " Particle"
                << "(globalID=" << globalID << ") "
                << "from iC" << "X" << "(rank=" << rankOrig << ") to iC"
                << globiCdest << "(rank=" << rank << ")" << std::endl;
      particle.template print<true>();
#endif

    } //if (globiC == globiCdest)
  }); //forSystemsInSuperParticleSystem
}




#ifdef PARALLEL_MODE_MPI


//For received data (particle agnostic)
//WARNING: Expects a message (potentially empty) from all availableRanks
template<typename F>
void receiveAndExecuteForData(
  const std::set<int>& availableRanks,
  std::size_t serialSize,
  MPI_Comm Comm,
  singleton::MpiNonBlockingHelper& mpiNbHelper,
  F f )
{
  //Iterate available ranks
  for (int rankOrig : availableRanks) {
    //Probe whether incomming data on rank
    std::size_t noBytesReceived = singleton::mpi().probeReceiveSize(rankOrig, MPI_BYTE, 1, Comm);
    //Create receive buffer as vector
    std::vector<std::uint8_t> recv_buffer(noBytesReceived);
#ifdef VERBOSE_COMMUNICATION
      int rank = singleton::mpi().getRank();
      std::cout << "Serial(rank=" << rank << "): "
                << "noBytesToBeReceived=" << noBytesReceived << "; " << std::endl;
#endif
    //Receive data via mpi
    singleton::mpi().receive<std::uint8_t>(
      recv_buffer.data(),
      noBytesReceived,
      rankOrig,
      1,
      Comm);
    //Loop over bytes in chunks of serialSize
    for (unsigned iByte=0; iByte < noBytesReceived; iByte += serialSize) {
      //Define raw buffer
      std::uint8_t* rawBuffer = &recv_buffer[iByte];
      //Execute passed function
      f( rankOrig, rawBuffer );
#ifdef VERBOSE_COMMUNICATION
      int rank = singleton::mpi().getRank();
      std::cout << "Serial(rank=" << rank << "): "
                << "noBytesReceived=" << noBytesReceived << "; "
                << "iByte=" << iByte << std::endl;
#endif
      //assignParticleToIC( sParticleSystem, receivedLocalIDs, rankOrig, rawBuffer );
    } // for (int iByte=0; iByte<noBytesReceived; iByte+=serialSize)
  } // for (auto rankOrig : availableRanks)
  singleton::mpi().waitAll(mpiNbHelper);
  //Free mpi helper
  mpiNbHelper.free();
}



//Recieve particles
template<typename T, typename PARTICLETYPE>
void receiveParticles(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  std::vector<std::unique_ptr<std::uint8_t[]>>& dataListRelocationIntra,
  std::size_t serialSize,
  MPI_Comm particleDistributionComm,
  singleton::MpiNonBlockingHelper& mpiNbHelper,
  const Vector<bool,PARTICLETYPE::d>& periodicity )
{
  using namespace particles::access;

  //Retrieve load balancer and cuboid geometry
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer = superStructure.getLoadBalancer();

  //Create list of local particle iDs to check afterwards,
  //  whether data was received (resolved only)
  std::vector<std::vector<std::size_t>> receivedLocalIDs(
    loadBalancer.size());

  //Intra core relocations (without MPI)
  for ( auto& buffer : dataListRelocationIntra ){
    //Define raw buffer
    std::uint8_t* rawBuffer = buffer.get();
    int rankOrig = singleton::mpi().getRank();
    assignParticleToIC( sParticleSystem, receivedLocalIDs, rankOrig, rawBuffer );
  }

  //Retrieve neighbour ranks
  auto& listNeighbourRanks = sParticleSystem.getExtendedNeighbourRanks();
  //Iterate over inter core received particle data
  receiveAndExecuteForData( listNeighbourRanks, serialSize,
      particleDistributionComm, mpiNbHelper,
    [&](int rankOrig, std::uint8_t* rawBuffer){
    //Assign particle to destination iC
    assignParticleToIC( sParticleSystem, receivedLocalIDs, rankOrig, rawBuffer );
  });

  //Check, whether invalidation nccessary
  if constexpr( providesSurface<PARTICLETYPE>() ){
    checkInvalidationOnReceival( sParticleSystem, receivedLocalIDs, periodicity );
  }
}
#endif

//Update particle cuboid distribution
template<typename T, typename PARTICLETYPE>
void updateParticleCuboidDistribution(
  SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
  const T physDeltaX,
#ifdef PARALLEL_MODE_MPI
  MPI_Comm particleDistributionComm,
#endif
  const Vector<bool,PARTICLETYPE::d>& periodicity
  )
{
  //Retrieve serial size
  std::size_t serialSize = sParticleSystem.getSerialSize();

  //Increase serial size by Field holding the globiC of the destination
  // - for that, using FieldArrayD with single entry
  //TODO: Optimization possible here:
  // - Necessary for resolved particles
  // - Not necessary for subgrid: This could be exchanged by on-rank retrieval
  //   of the ic depending on particle position. This would lead to a reduction
  //   in buffer size, which will most likely speed up communication significantly.
  //   A separate treatment will however be very error prone.
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,descriptors::IC> fieldiC(1);
  const ConstSpan<unsigned> index(0, 1);
  // Workaround for Clang argument deduction bug
  auto communicatable = ConcreteCommunicatable(fieldiC);
  std::size_t serialSizeField = communicatable.size(index);
  serialSize += serialSizeField;

  //Check rankDataMapRelocationInter necessity and create relocation map TODO: adapt description
  std::multimap<int,std::unique_ptr<std::uint8_t[]>> rankDataMapRelocationInter;
  std::vector<std::unique_ptr<std::uint8_t[]>> dataListRelocationIntra;
  checkRelocation( sParticleSystem, physDeltaX, rankDataMapRelocationInter,
      dataListRelocationIntra, serialSize, periodicity );

#ifdef PARALLEL_MODE_MPI

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  //Retrieve rank of direct neighbours
  auto& listNeighbourRanks = sParticleSystem.getExtendedNeighbourRanks();

  //Create ranDataMapSorted (WARNING: has to be existent until all data is received! c.f. #290)
  std::map<int, std::vector<std::uint8_t> > rankDataMapSorted;

  //Fill send buffer
  fillSendBuffer( rankDataMapRelocationInter, rankDataMapSorted, serialSize );

  //Send particles in rankDataMapRelocationInter map via mpi().iSend and return global number of bytes
  sendMappedData( rankDataMapSorted, listNeighbourRanks,
    serialSize, particleDistributionComm, mpiNbHelper );

  //Receive particles via mpi().receive
  receiveParticles( sParticleSystem, dataListRelocationIntra,
    serialSize, particleDistributionComm, mpiNbHelper, periodicity );

#else
  std::cerr
    << "ERROR: Invalidation treatement for sequential runs not implemented yet!"
    << std::endl;
  throw;
#endif

}

#ifdef PARALLEL_MODE_MPI

template<typename DATA>
void collectDataAndAppendVector( std::vector<DATA>& dataVector,
  std::set<int>& availableRanks,
  singleton::MpiNonBlockingHelper& mpiNbHelper,
  MPI_Comm commGroup = MPI_COMM_WORLD)
{
  DATA data;
  auto communicatable = ConcreteCommunicatable(data);
  const int serialSize = communicatable.size();
  //Receive data
  receiveAndExecuteForData( availableRanks,
    serialSize, commGroup, mpiNbHelper,
    [&](int rankOrig, std::uint8_t* rawBuffer){
      //Deserialize buffer
      communicatable.deserialize( rawBuffer );
      //Append data to std::vector
      dataVector.push_back(data);
  });
}

#endif


} //namespace communication

} //namespace particles

} //namespace olb


#endif
