/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Jan E. Marquardt, Mathias J. Krause
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

/*
 * This file contains functions used for the particle-wall and particle-particle communication.
*/

#ifndef COMMUNICATION_EQUATIONS_OF_MOTION_RESULTS_H
#define COMMUNICATION_EQUATIONS_OF_MOTION_RESULTS_H

#include <unordered_set>

#include "particles/particles.h"

namespace olb {

namespace particles {

namespace communication {

/// evaluates ranks that need the new particle data
template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
std::set<int>
evalDestRanks(contact::ContactContainer<T, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>& contactContainer,
              const std::size_t                           globalParticleID)
{
  std::set<int> destRanks;

  for (const PARTICLECONTACTTYPE& contact : contactContainer.particleContacts) {
    if (
        //!contact.isEmpty() &&
        (globalParticleID == contact.getIDs()[0] ||
         globalParticleID == contact.getIDs()[1])) {
      destRanks.insert(contact.getResponsibleRank());
    }
  }
  for (const WALLCONTACTTYPE& contact : contactContainer.wallContacts) {
    if (
        //!contact.isEmpty() &&
        globalParticleID == contact.getParticleID()) {
      destRanks.insert(contact.getResponsibleRank());
    }
  }

  return destRanks;
}

// Communicate results of equations of motion to blocks on same rank
template <typename T, typename PARTICLETYPE>
void communicateEquationsOfMotionResultsIntra(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    const std::size_t& pID, const Vector<T, PARTICLETYPE::d>& position,
    const Vector<T, PARTICLETYPE::d>& velocity,
    const Vector<T, utilities::dimensions::convert<PARTICLETYPE::d>::rotation>&
        angle,
    const Vector<T, utilities::dimensions::convert<PARTICLETYPE::d>::rotation>&
        angVelocity)
{
  using namespace descriptors;

  particles::communication::forParticlesInSuperParticleSystem(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t currGlobalParticleID =
            particle.template getField<PARALLELIZATION, ID>();

        if (currGlobalParticleID == pID) {
          particle.template setField<GENERAL, POSITION>(position);
          particle.template setField<MOBILITY, VELOCITY>(velocity);
          particle.template setField<SURFACE, ANGLE>(
              utilities::dimensions::convert<
                  PARTICLETYPE::d>::serialize_rotation(angle));
          particle.template setField<MOBILITY, ANG_VELOCITY>(
              utilities::dimensions::convert<
                  PARTICLETYPE::d>::serialize_rotation(angVelocity));
        }
      });
}

// Communication of results of equation of motion
#ifdef PARALLEL_MODE_MPI
template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE>
std::size_t evalEquationsOfMotionResults(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    contact::ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
                                                         contactContainer,
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap)
{
  using namespace descriptors;

  //Create FieldArrayD for globalID, force and torque
  using GENERAL_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::GENERAL>;
  using MOBILITY_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::MOBILITY>;
  using SURFACE_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::SURFACE>;
  using POSITION_EVAL =
      typename GENERAL_EVAL ::template derivedField<descriptors::POSITION>;
  using VELOCITY_EVAL =
      typename MOBILITY_EVAL ::template derivedField<descriptors::VELOCITY>;
  using ANGLE_EVAL =
      typename SURFACE_EVAL ::template derivedField<descriptors::ANGLE>;
  using ANGVELOCITY_EVAL =
      typename MOBILITY_EVAL ::template derivedField<descriptors::ANG_VELOCITY>;
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, descriptors::ID> fieldID(1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, POSITION_EVAL> fieldPosition(
      1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, VELOCITY_EVAL> fieldVelocity(
      1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, ANGLE_EVAL> fieldAngle(1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, ANGVELOCITY_EVAL>
      fieldAngVelocity(1);

  //Create communicatables
  auto communicatableID          = ConcreteCommunicatable(fieldID);
  auto communicatablePosition    = ConcreteCommunicatable(fieldPosition);
  auto communicatableVelocity    = ConcreteCommunicatable(fieldVelocity);
  auto communicatableAngle       = ConcreteCommunicatable(fieldAngle);
  auto communicatableAngVelocity = ConcreteCommunicatable(fieldAngVelocity);

  //Retrieve serial size
  const std::vector<unsigned int> indices {0};
  const std::size_t               serialSize =
      communicatableID.size(indices) + communicatablePosition.size(indices) +
      communicatableVelocity.size(indices) + communicatableAngle.size(indices) +
      communicatableAngVelocity.size(indices);

  // Fill map with data
  communication::forParticlesInSuperParticleSystem(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        using PCONDITION = std::conditional_t<
            access::providesSurface<PARTICLETYPE>(),
            conditions::
                valid_particle_centres,   // only consider center for resolved
            conditions::valid_particles>; // only consider valid

        if (PCONDITION::value(particle, globiC)) {
          // Get quantities to package
          fieldID.setField(0,
                           particle.template getField<PARALLELIZATION, ID>());
          fieldPosition.setField(0, access::getPosition(particle));
          fieldVelocity.setField(0, access::getVelocity(particle));
          fieldAngle.setField(0, access::getAngle(particle));
          fieldAngVelocity.setField(0, access::getAngularVelocity(particle));

          communicateEquationsOfMotionResultsIntra(
              sParticleSystem,
              particle.template getField<PARALLELIZATION, ID>(),
              access::getPosition(particle), access::getVelocity(particle),
              access::getAngle(particle), access::getAngularVelocity(particle));

          // Evaluate which ranks need to know the data
          std::set<int> destRanks {
              evalDestRanks(contactContainer, fieldID.getField(0))};

          for (const int destRank : destRanks) {
            std::unique_ptr<std::uint8_t[]> buffer(
                new std::uint8_t[serialSize] {});
            std::uint8_t* bufferRaw = buffer.get();
            std::size_t   serialIdx =
                communicatableID.serialize(indices, bufferRaw);
            serialIdx += communicatablePosition.serialize(
                indices, &bufferRaw[serialIdx]);
            serialIdx += communicatableVelocity.serialize(
                indices, &bufferRaw[serialIdx]);
            serialIdx +=
                communicatableAngle.serialize(indices, &bufferRaw[serialIdx]);
            serialIdx += communicatableAngVelocity.serialize(
                indices, &bufferRaw[serialIdx]);
            dataMap.insert(std::make_pair(destRank, std::move(buffer)));
          }

#ifdef VERBOSE_CONTACT_COMMUNICATION
          std::cout << "Rank " << singleton::mpi().getRank()
                    << " sends particle data of ID " << fieldID.getField(0)
                    << " to ";
          for (const int destRank : destRanks) {
            std::cout << destRank << " ";
          }
          std::cout << std::endl;
#endif
        }
      });

  return serialSize;
}
#endif

#ifdef PARALLEL_MODE_MPI
template <typename T, typename PARTICLETYPE>
void receiveEquationsOfMotionResults(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    MPI_Comm                              equationsOfMotionComm,
    singleton::MpiNonBlockingHelper&      mpiNbHelper)
{
  using namespace descriptors;

  // Retrieve rank of direct neighbours
  auto& listNeighbourRanks = sParticleSystem.getNeighbourRanks();

  //Create FieldArrayD for globalID, force and torque
  using GENERAL_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::GENERAL>;
  using MOBILITY_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::MOBILITY>;
  using SURFACE_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::SURFACE>;
  using POSITION_EVAL =
      typename GENERAL_EVAL ::template derivedField<descriptors::POSITION>;
  using VELOCITY_EVAL =
      typename MOBILITY_EVAL ::template derivedField<descriptors::VELOCITY>;
  using ANGLE_EVAL =
      typename SURFACE_EVAL ::template derivedField<descriptors::ANGLE>;
  using ANGVELOCITY_EVAL =
      typename MOBILITY_EVAL ::template derivedField<descriptors::ANG_VELOCITY>;
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, descriptors::ID> fieldID(1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, POSITION_EVAL> fieldPosition(
      1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, VELOCITY_EVAL> fieldVelocity(
      1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, ANGLE_EVAL> fieldAngle(1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, ANGVELOCITY_EVAL>
      fieldAngVelocity(1);

  //Create communicatables
  auto communicatableID          = ConcreteCommunicatable(fieldID);
  auto communicatablePosition    = ConcreteCommunicatable(fieldPosition);
  auto communicatableVelocity    = ConcreteCommunicatable(fieldVelocity);
  auto communicatableAngle       = ConcreteCommunicatable(fieldAngle);
  auto communicatableAngVelocity = ConcreteCommunicatable(fieldAngVelocity);

  // Retrieve serial size
  const std::vector<unsigned int> indices {0};
  const std::size_t               serialSize =
      communicatableID.size(indices) + communicatablePosition.size(indices) +
      communicatableVelocity.size(indices) + communicatableAngle.size(indices) +
      communicatableAngVelocity.size(indices);

  // Receive data and iterate buffer
  // TODO: Improve processing of data, do not iterate through the all particles every time (maybe also add main iC to the communicated data to reduce processed particles)
  communication::receiveAndExecuteForData(
      listNeighbourRanks, serialSize, equationsOfMotionComm, mpiNbHelper,
      [&](int rankOrig, std::uint8_t* buffer) {
        std::size_t serialIdx = communicatableID.deserialize(indices, buffer);
        serialIdx +=
            communicatablePosition.deserialize(indices, &buffer[serialIdx]);
        serialIdx +=
            communicatableVelocity.deserialize(indices, &buffer[serialIdx]);
        serialIdx +=
            communicatableAngle.deserialize(indices, &buffer[serialIdx]);
        serialIdx +=
            communicatableAngVelocity.deserialize(indices, &buffer[serialIdx]);

        communication::forParticlesInSuperParticleSystem(
            sParticleSystem,
            [&](Particle<T, PARTICLETYPE>&       particle,
                ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
              const std::size_t currentGlobalParticleID =
                  particle.template getField<PARALLELIZATION, ID>();
              if (fieldID.getField(0) == currentGlobalParticleID) {
                // Override particle data
                particle.template setField<GENERAL, POSITION>(
                    fieldPosition.getField(0));
                particle.template setField<MOBILITY, VELOCITY>(
                    fieldVelocity.getField(0));
                particle.template setField<SURFACE, ANGLE>(
                    fieldAngle.getField(0));
                particle.template setField<MOBILITY, ANG_VELOCITY>(
                    fieldAngVelocity.getField(0));
              }
            });
      });
}
#endif

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE>
void communicateEquationsOfMotionResults(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    contact::ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
        contactContainer
#ifdef PARALLEL_MODE_MPI
    ,
    MPI_Comm equationsOfMotionComm
#endif
)
{
  // TODO: Test if it is faster to only send & receive between processors that must exchange data (that way we don't need to synchronize all processors)

#ifdef PARALLEL_MODE_MPI
  // Create map for data
  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  std::size_t serialSize =
      evalEquationsOfMotionResults(sParticleSystem, contactContainer, dataMap);

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  std::map<int, std::vector<std::uint8_t>> rankDataMapSorted;
  communication::fillSendBuffer(dataMap, rankDataMapSorted, serialSize);

  // Send mapped data
  communication::sendMappedData(rankDataMapSorted,
                                sParticleSystem.getNeighbourRanks(), serialSize,
                                equationsOfMotionComm, mpiNbHelper);

  receiveEquationsOfMotionResults(sParticleSystem, equationsOfMotionComm,
                                  mpiNbHelper);
#endif
}

} // namespace communication

} // namespace particles

} // namespace olb

#endif
