/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt, Mathias J. Krause
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

#ifndef PARTICLE_CONTACT_COMMUNICATION_FUNCTIONS_H
#define PARTICLE_CONTACT_COMMUNICATION_FUNCTIONS_H

namespace olb {

namespace particles {

namespace contact {

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void communicateContacts(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "contactCommunication");

  if (singleton::mpi().getSize() > 1) {
    int sumParticleContact = 0;
    int sumWallContact     = 0;
    int countParticleContactsPerRank[singleton::mpi().getSize()];
    int countWallContactsPerRank[singleton::mpi().getSize()];
    int localCountParticleContacts = 0;
    int localCountWallContacts     = 0;

    //clout << "communicate number of contacts ..." << std::endl;
    if (singleton::mpi().isMainProcessor()) {
      countParticleContactsPerRank[singleton::mpi().bossId()] =
          contactContainer.particleContacts.size();
      sumParticleContact +=
          countParticleContactsPerRank[singleton::mpi().bossId()];
      countWallContactsPerRank[singleton::mpi().bossId()] =
          contactContainer.wallContacts.size();
      sumWallContact += countWallContactsPerRank[singleton::mpi().bossId()];
      for (int rank = 0; rank < singleton::mpi().getSize(); ++rank) {
        if (rank != singleton::mpi().bossId()) {
          singleton::mpi().receive(&countParticleContactsPerRank[rank], 1,
                                   rank);
          sumParticleContact += countParticleContactsPerRank[rank];
          singleton::mpi().receive(&countWallContactsPerRank[rank], 1, rank);
          sumWallContact += countWallContactsPerRank[rank];
        }
      }
    }
    else {
      localCountParticleContacts =
          contactContainer.particleContacts
              .size(); // size of std::vector (called contacts in container)
      singleton::mpi().send(&localCountParticleContacts, 1,
                            singleton::mpi().bossId()); // 0 is rank of master
      localCountWallContacts = contactContainer.wallContacts.size();
      singleton::mpi().send(&localCountWallContacts, 1,
                            singleton::mpi().bossId());
    }
    //clout << "communicate number of contacts ... OK" << std::endl;

    //clout << "broadcasting number of contacts ..." << std::endl;
    singleton::mpi().bCast(&sumParticleContact, 1, singleton::mpi().bossId());
    singleton::mpi().bCast(&sumWallContact, 1, singleton::mpi().bossId());
    //std::cout << sumParticleContact << std::endl;
    //singleton::mpi().barrier();
    //clout << "broadcasting number of contacts ... OK" << std::endl;

    if (sumParticleContact > 0 || sumWallContact > 0) {
      ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> container(
          sumParticleContact, sumWallContact);
      int particleContactObjSize = sizeof(PARTICLECONTACTTYPE);
      int wallContactObjSize     = sizeof(WALLCONTACTTYPE);

      //clout << "Bundling contact objects ..." << std::endl;
      if (singleton::mpi().isMainProcessor()) {
        PARTICLECONTACTTYPE* currParticle = container.particleContacts.data();
        WALLCONTACTTYPE*     currWall     = container.wallContacts.data();
        for (int rank = 0; rank < singleton::mpi().getSize(); ++rank) {
          std::uint8_t* bytesParticle = (std::uint8_t*)currParticle;
          std::uint8_t* bytesWall     = (std::uint8_t*)currWall;
          if (rank != singleton::mpi().bossId()) {
            singleton::mpi().receive(bytesParticle,
                                     countParticleContactsPerRank[rank] *
                                         particleContactObjSize,
                                     rank);
            singleton::mpi().receive(
                bytesWall, countWallContactsPerRank[rank] * wallContactObjSize,
                rank);
          }
          else {
            std::memcpy(
                bytesParticle,
                (std::uint8_t*)contactContainer.particleContacts.data(),
                countParticleContactsPerRank[singleton::mpi().bossId()] *
                    particleContactObjSize);
            std::memcpy(bytesWall,
                        (std::uint8_t*)contactContainer.wallContacts.data(),
                        countWallContactsPerRank[singleton::mpi().bossId()] *
                            wallContactObjSize);
          }
          currParticle += countParticleContactsPerRank[rank];
          currWall += countWallContactsPerRank[rank];
        }
      }
      else {
        PARTICLECONTACTTYPE* currParticle =
            contactContainer.particleContacts.data();
        WALLCONTACTTYPE* currWall      = contactContainer.wallContacts.data();
        std::uint8_t*    bytesParticle = (std::uint8_t*)currParticle;
        std::uint8_t*    bytesWall     = (std::uint8_t*)currWall;
        singleton::mpi().send(
            bytesParticle, localCountParticleContacts * particleContactObjSize,
            singleton::mpi().bossId());
        singleton::mpi().send(bytesWall,
                              localCountWallContacts * wallContactObjSize,
                              singleton::mpi().bossId());
      }
      //clout << "Bundling contact objects ... OK" << std::endl;

      if (singleton::mpi().isMainProcessor()) {
        container.combineContacts();
      }

      //clout << "Broadcasting contact objects ..." << std::endl;
      std::uint8_t* containerParticlesBytes =
          (std::uint8_t*)container.particleContacts.data();
      singleton::mpi().bCast(containerParticlesBytes,
                             sumParticleContact * particleContactObjSize,
                             singleton::mpi().bossId());
      std::uint8_t* containerWallBytes =
          (std::uint8_t*)container.wallContacts.data();
      singleton::mpi().bCast(containerWallBytes,
                             sumWallContact * wallContactObjSize,
                             singleton::mpi().bossId());
      //singleton::mpi().barrier();
      //clout << "Broadcasting contact objects ... OK" << std::endl;
      /*
            clout << "#######" << std::endl;
            clout << "id: " << container.wallContacts[3].id << std::endl;
            clout << "#######" << std::endl;
      */
      //std::cout << "size before: " << container.particleContacts.size() << std::endl;
      container.cleanContacts();
      //std::cout << "size after: " << container.particleContacts.size() << std::endl;
      contactContainer = container;
    }
  }
#endif
}

/// update positions in contact to account for periodic boundaries
template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, typename F>
void accountForPeriodicParticleBoundary(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    contact::ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
                                             contactContainer,
    const SuperGeometry<T, PARTICLETYPE::d>& sGeometry, F getSetupPeriodicity)
{
  using namespace descriptors;

  // if a periodic setup is used, we have to update min and max in rare occasions
  // i.e. when the first/only particle "jumps" from the end to the beginning
  // because then the min and max from the step before isn't valid anymore
  const T physDeltaX =
      sGeometry.getCuboidDecomposition().getMotherCuboid().getDeltaR();
  const PhysR<T, PARTICLETYPE::d> cellMin =
      particles::communication::getCuboidMin<T, PARTICLETYPE::d>(
          sGeometry.getCuboidDecomposition());
  const PhysR<T, PARTICLETYPE::d> cellMax =
      particles::communication::getCuboidMax<T, PARTICLETYPE::d>(
          sGeometry.getCuboidDecomposition(), cellMin);

  communication::forParticlesInSuperParticleSystem(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t currentGlobalParticleID =
            particle.template getField<PARALLELIZATION, ID>();

        for (WALLCONTACTTYPE& wallContact : contactContainer.wallContacts) {
          if (!wallContact.isEmpty() && !wallContact.isNew() &&
              wallContact.getParticleID() == currentGlobalParticleID &&
              wallContact.getResponsibleRank() == singleton::mpi().getRank()) {
            const PhysR<T, PARTICLETYPE::d> min = wallContact.getMin();
            const PhysR<T, PARTICLETYPE::d> max = wallContact.getMax();
            wallContact.resetMinMax();
            const PhysR<T, PARTICLETYPE::d> unifiedPos = contact::unifyPosition(
                particle, cellMin, cellMax, getSetupPeriodicity, physDeltaX);
            wallContact.setParticlePosition(unifiedPos);
            const PhysR<T, PARTICLETYPE::d> newMin =
                contact::evalContactPosition(
                    particle, wallContact.getParticlePosition(), min, cellMin,
                    cellMax, getSetupPeriodicity, physDeltaX);
            const PhysR<T, PARTICLETYPE::d> newMax =
                contact::evalContactPosition(
                    particle, wallContact.getParticlePosition(), max, cellMin,
                    cellMax, getSetupPeriodicity, physDeltaX);
            wallContact.updateMinMax(newMin);
            wallContact.updateMinMax(newMax);
          }
        }
        for (PARTICLECONTACTTYPE& particleContact :
             contactContainer.particleContacts) {
          if (!particleContact.isEmpty() && !particleContact.isNew() &&
              particleContact.getIDs()[0] == currentGlobalParticleID &&
              particleContact.getResponsibleRank() ==
                  singleton::mpi().getRank()) {
            const PhysR<T, PARTICLETYPE::d> min = particleContact.getMin();
            const PhysR<T, PARTICLETYPE::d> max = particleContact.getMax();
            particleContact.resetMinMax();
            std::array<PhysR<T, PARTICLETYPE::d>, 2> unifiedPos;
            auto particleB = sParticleSystem.get(particleContact.getIDs()[1]);
            contact::unifyPositions(particle, particleB, cellMin, cellMax,
                                    getSetupPeriodicity, unifiedPos,
                                    physDeltaX);
            particleContact.setParticlePositions(unifiedPos);
            const PhysR<T, PARTICLETYPE::d> newMin =
                contact::evalContactPosition(
                    particle, particleContact.getParticlePositions()[0], min,
                    cellMin, cellMax, getSetupPeriodicity, physDeltaX);
            const PhysR<T, PARTICLETYPE::d> newMax =
                contact::evalContactPosition(
                    particle, particleContact.getParticlePositions()[0], max,
                    cellMin, cellMax, getSetupPeriodicity, physDeltaX);
            particleContact.updateMinMax(newMin);
            particleContact.updateMinMax(newMax);
          }
        }
      });
}

/// evaluate responsible ranks for particle-particle contact
template <typename T, typename PARTICLETYPE, bool CONVEX>
std::set<int> evalDestRanksForDetectionCommunication(
    ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                                          contact,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
  using namespace descriptors;
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer   = superStructure.getLoadBalancer();
  /* auto& cuboidDecomposition = superStructure.getCuboidDecomposition(); */

  // Set instead of unordered_set so that it's sorted
  // Therefore, it is easier to find the the rank with smallest ID
  std::set<int>                intersectingRanks {};
  std::array<int, 2>           responsibleRank = {-1, -1};
  std::array<std::set<int>, 2> touchedRanks {};
  std::set<int>      destRanks {};

  // Function that is called on failure
  const auto errorTreatment = [&contact, &destRanks]() {
    contact.resetMinMax();
    resetResponsibleRank(contact);
    destRanks.clear();
  };

  // Iterate over all particles to populate responsible and touched ranks of each particle in the contact
  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE,
                                                   conditions::valid_particles>(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();

        for (unsigned short index = 0; index < 2; ++index) {
          if (globalParticleID == contact.getIDs()[index]) {
            const int globiCcenter =
                particle.template getField<PARALLELIZATION, IC>();
            responsibleRank[index] = loadBalancer.rank(globiCcenter);
            touchedRanks[index].insert(responsibleRank[index]);

            const PhysR<T, PARTICLETYPE::d> position =
                access::getPosition(particle);
            const T circumRadius = evalCircumRadius(particle, deltaX);
            communication::getSurfaceTouchingICs(
                sParticleSystem, position, circumRadius, periodicity,
                globiCcenter, [&](const int iC) {
                  touchedRanks[index].insert(loadBalancer.rank(iC));
                });
          }
        }
      });

  // If at least one particle has no touched cuboids then cancel function call
  // (this should be caused by an invalid particle that somehow still has an existing contact)
  if (touchedRanks[0].empty() || touchedRanks[1].empty() ||
      responsibleRank[0] < 0 || responsibleRank[1] < 0) {
    errorTreatment();
  }
  // If both particles have the same responsible rank then let that rank take care of the contact
  else if (responsibleRank[0] == responsibleRank[1]) {
    destRanks.insert(responsibleRank[0]);
    contact.setResponsibleRank(responsibleRank[0]);
  }
  // If not, then determine treatment rank from intersection of ranks that touch the particles
  else {
    std::set_intersection(
        touchedRanks[0].begin(), touchedRanks[0].end(), touchedRanks[1].begin(),
        touchedRanks[1].end(),
        std::inserter(intersectingRanks, intersectingRanks.begin()));

    /*
     * It's possible that no intersection is found, if a rank loses responsibility
     * for the contact, but still holds the old information. At the same time,
     * the contact of both particles breaks, and they are also on different ranks,
     * without any processor that knows both. This is very unlikely, but it can
     * happen.
     * However, the following simple if-condition shouldn't be expensive as it is
     * almost always true.
     */
    if (!intersectingRanks.empty()) {
      destRanks.insert(responsibleRank[0]);
      destRanks.insert(responsibleRank[1]);

      // Sort responsible ranks so that the following is consistent
      if (responsibleRank[0] > responsibleRank[1]) {
        std::swap(responsibleRank[0], responsibleRank[1]);
      }

      // First use one of the responsible ranks if they know both particles
      if (intersectingRanks.find(responsibleRank[0]) !=
          intersectingRanks.end()) {
        contact.setResponsibleRank(responsibleRank[0]);
      }
      else if (intersectingRanks.find(responsibleRank[1]) !=
               intersectingRanks.end()) {
        contact.setResponsibleRank(responsibleRank[1]);
      }
      // Otherwise use smallest rank of intersection
      else {
        const int responsibleRankContact = *intersectingRanks.begin();
        contact.setResponsibleRank(responsibleRankContact);
        destRanks.insert(responsibleRankContact);
      }
    }
    else {
      errorTreatment();
    }
  }

#ifdef VERBOSE_CONTACT_COMMUNICATION
  std::cout << "Contact of ids " << contact.getIDs()[0] << " & "
            << contact.getIDs()[1] << " goes to " << destRanks << " and "
            << contact.getResponsibleRank()
            << " is responsible for the calculation (Rank "
            << singleton::mpi().getRank() << ")"
            << " the intersection was " << intersectingRanks
            << " from touched ranks " << touchedRanks[0] << " and "
            << touchedRanks[1] << std::endl;
#endif

  return destRanks;
}

/// evaluate responsible rank for solid boundary contact
template <typename T, typename PARTICLETYPE, bool CONVEX>
std::set<int> evalDestRanksForDetectionCommunication(
    WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>& contact,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
    [[maybe_unused]] const Vector<bool, PARTICLETYPE::d>& periodicity)
{
  using namespace descriptors;
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer   = superStructure.getLoadBalancer();

  bool                    isSuccessful = false;
  std::set<int> destRanks {};
  int                     destRank;
  // Find responsible rank
  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE,
                                                   conditions::valid_particles>(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();
        if (globalParticleID == contact.getParticleID()) {
          std::size_t globiCcenter =
              particle.template getField<PARALLELIZATION, IC>();
          destRank     = loadBalancer.rank(globiCcenter);
          isSuccessful = true;
        }
      });

  if (isSuccessful) {
    contact.setResponsibleRank(destRank);
    if (destRank != singleton::mpi().getRank()) {
      destRanks.insert(destRank);
    }

#ifdef VERBOSE_CONTACT_COMMUNICATION
    std::cout << "Contact of id " << contact.getParticleID()
              << " & solid boundary " << contact.getWallID() << " goes to "
              << destRanks << std::endl;
#endif
  }

  return destRanks;
}

/// evaluate ranks that touch both particles of a particle-particle contact
template <typename T, typename PARTICLETYPE, bool CONVEX>
std::set<int> evalDestRanksForPostContactTreatmentCommunication(
    ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                                          contact,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
  using namespace descriptors;
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer   = superStructure.getLoadBalancer();
  auto& cuboidDecomposition = superStructure.getCuboidDecomposition();

  std::set<int>      destRanks;
  std::array<std::set<int>, 2> touchedRanks;

  // Iterate over all particles to populate touched ranks of each particle in the contact
  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE,
                                                   conditions::valid_particles>(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();

        for (unsigned short index = 0; index < 2; ++index) {
          if (globalParticleID == contact.getIDs()[index]) {
            const PhysR<T, PARTICLETYPE::d> position =
                access::getPosition(particle);
            // the particle moved, therefore, the iC may have changed
            int        globiCcenter;
            const bool cuboidFound = communication::getCuboid(
                cuboidDecomposition, periodicity, position, globiCcenter);
            if (cuboidFound) {
              touchedRanks[index].insert(loadBalancer.rank(globiCcenter));

              const T circumRadius = evalCircumRadius(particle, deltaX);
              communication::getSurfaceTouchingICs(
                  sParticleSystem, position, circumRadius, periodicity,
                  globiCcenter, [&](const int iC) {
                    touchedRanks[index].insert(loadBalancer.rank(iC));
                  });
            }
          }
        }
      });

  std::set_intersection(touchedRanks[0].begin(), touchedRanks[0].end(),
                        touchedRanks[1].begin(), touchedRanks[1].end(),
                        std::inserter(destRanks, destRanks.begin()));

#ifdef VERBOSE_CONTACT_COMMUNICATION
  std::cout << "Contact of ids " << contact.getIDs()[0] << " & "
            << contact.getIDs()[1] << " goes to " << destRanks << " and "
            << contact.getResponsibleRank()
            << " is responsible for the calculation (Rank "
            << singleton::mpi().getRank() << ")" << '\n'
            << "Touched ranks of 1. particle: " << touchedRanks[0]
            << " and touched ranks of 2. particle: " << touchedRanks[1]
            << std::endl;
#endif

  return destRanks;
}

/// evaluate ranks that touch the particle of a particle-wall contact
template <typename T, typename PARTICLETYPE, bool CONVEX>
std::set<int> evalDestRanksForPostContactTreatmentCommunication(
    WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>& contact,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
  using namespace descriptors;
  auto& superStructure = sParticleSystem.getSuperStructure();
  auto& loadBalancer   = superStructure.getLoadBalancer();
  auto& cuboidDecomposition = superStructure.getCuboidDecomposition();

  std::set<int> destRanks;

  // Find responsible rank
  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE,
                                                   conditions::valid_particles>(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();
        if (globalParticleID == contact.getParticleID()) {
          const PhysR<T, PARTICLETYPE::d> position =
              access::getPosition(particle);
          // the particle moved, therefore, the iC may have changed
          int        globiCcenter;
          const bool cuboidFound = communication::getCuboid(
              cuboidDecomposition, periodicity, position, globiCcenter);
          if (cuboidFound) {
            destRanks.insert(loadBalancer.rank(globiCcenter));

            const T circumRadius {evalCircumRadius(particle, deltaX)};
            communication::getSurfaceTouchingICs(
                sParticleSystem, position, circumRadius, periodicity,
                globiCcenter, [&](const int iC) {
                  destRanks.insert(loadBalancer.rank(iC));
                });
          }
        }
      });

#ifdef VERBOSE_CONTACT_COMMUNICATION
  std::cout << "Contact of id " << contact.getParticleID()
            << " & solid boundary " << contact.getWallID() << " goes to "
            << destRanks << " and rank " << contact.getResponsibleRank()
            << " is responsible for the calculation" << std::endl;
#endif

  return destRanks;
}

template <typename T, typename PARTICLETYPE, typename CONTACTTYPE, bool CONVEX,
          bool SANITYCHECK = false>
void receiveContact(
    ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                                          newcontact,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    std::vector<CONTACTTYPE>& contacts, [[maybe_unused]] int rankOrig)
{
  using namespace descriptors;

  bool addNewContact = true;

  // Optional sanity check
  if constexpr (SANITYCHECK) {
    if (newcontact.getResponsibleRank() == singleton::mpi().getRank()) {
      bool isParticleDataAvailable[2] = {false, false};
      // Iterate over all particles to check for errors
      communication::forParticlesInSuperParticleSystem(
          sParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            for (unsigned short i = 0; i < 2; ++i) {
              if (particle.template getField<PARALLELIZATION, ID>() ==
                  newcontact.getIDs()[i]) {
                isParticleDataAvailable[i] = true;
              }
            }
          });
      if (!isParticleDataAvailable[0] || !isParticleDataAvailable[1]) {
        addNewContact = false;
        std::cout << "WARNING: rank " << singleton::mpi().getRank()
                  << " received a contact of particles with the global ids "
                  << newcontact.getIDs()[0] << " and " << newcontact.getIDs()[1]
                  << " from rank " << rankOrig
                  << " for which no particle data are available." << std::endl;
      }
    }
  }

  if (addNewContact) {
    const std::array<std::size_t, 2> particleIDs = newcontact.getIDs();

    const std::function<bool(const CONTACTTYPE&)> condition =
        [&particleIDs](const CONTACTTYPE& contact) {
          return particleContactConsistsOfIDs(contact, particleIDs);
        };

    auto contactIt = getContactIterator(contacts, condition);

    if (contactIt != std::end(contacts)) {
      contactIt->combineWith(newcontact);
    }
    else {
      contacts.push_back(newcontact);
    }

#ifdef VERBOSE_CONTACT_COMMUNICATION
    std::cout << "Rank " << singleton::mpi().getRank()
              << " received contact of ids " << newcontact.getIDs()[0] << " & "
              << newcontact.getIDs()[1] << ". Rank "
              << newcontact.getResponsibleRank()
              << " is responsible for the calculation." << std::endl;
#endif
  }
}

template <typename T, typename PARTICLETYPE, typename CONTACTTYPE, bool CONVEX,
          bool SANITYCHECK = false>
void receiveContact(WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d,
                                                          CONVEX>& newcontact,
                    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
                    std::vector<CONTACTTYPE>&             contacts,
                    [[maybe_unused]] int                  rankOrig)
{
  using namespace descriptors;

  bool addNewContact = true;

  if constexpr (SANITYCHECK) {
    if (newcontact.getResponsibleRank() == singleton::mpi().getRank()) {
      bool isParticleDataAvailable = false;
      // Iterate over all particles to check for errors
      communication::forParticlesInSuperParticleSystem(
          sParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            if (particle.template getField<PARALLELIZATION, ID>() ==
                newcontact.getParticleID()) {
              isParticleDataAvailable = true;
            }
          });
      if (!isParticleDataAvailable) {
        addNewContact = false;
        std::cout << "WARNING: rank " << singleton::mpi().getRank()
                  << " received a contact of particle with the global id "
                  << newcontact.getParticleID() << " from rank " << rankOrig
                  << " for which no particle data is available." << std::endl;
      }
    }
  }

  if (addNewContact) {
    const std::size_t particleID            = newcontact.getParticleID();
    const std::size_t solidBoundaryMaterial = newcontact.getWallID();

    const std::function<bool(const CONTACTTYPE&)> condition =
        [&particleID, &solidBoundaryMaterial](const CONTACTTYPE& contact) {
          return particleID == contact.getParticleID() &&
                 solidBoundaryMaterial == contact.getWallID();
        };

    auto contactIt = getContactIterator(contacts, condition);

    if (contactIt != std::end(contacts)) {
      contactIt->combineWith(newcontact);
    }
    else {
      contacts.push_back(newcontact);
    }
  }
}

template <typename T, typename PARTICLETYPE, typename CONTACTTYPE>
void communicateParallelContacts(
    std::vector<CONTACTTYPE>&             contacts,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
#ifdef PARALLEL_MODE_MPI
    MPI_Comm contactDetectionComm,
#endif
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "contactCommunication");

  // Retrieve rank of direct neighbours
  auto& listNeighbourRanks         = sParticleSystem.getNeighbourRanks();
  constexpr std::size_t serialSize = CONTACTTYPE::getSerialSize();

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  for (CONTACTTYPE& contact : contacts) {
    // Reset responsible rank because we evaluate the new one in the next step
    resetResponsibleRank(contact);

    if (!contact.isEmpty()) {
      std::set<int> destRanks {evalDestRanksForDetectionCommunication(
          contact, sParticleSystem, deltaX, periodicity)};

      // Only new contacts need communication because they originate from the on-lattice
      // contact detection. Known contacts are not updated during that step.
      if (particles::contact::isResponsibleRankSet(contact)
          // not new contacts need to be communicated as the responsible processors for each particle must know the responsible rank for contact treatment, which is evaluated above
          && (contact.isNew() ||
              contact.getResponsibleRank() == singleton::mpi().getRank())) {
        for (const int destRank : destRanks) {
          if (destRank != singleton::mpi().getRank()) {
            std::unique_ptr<std::uint8_t[]> buffer(
                new std::uint8_t[serialSize] {});
            contact.serialize(buffer.get());
            dataMap.insert(std::make_pair(destRank, std::move(buffer)));
          }
        }
      }

      // Reset contact on other ranks to avoid unnecessary communication
      //if (destRanks.find(singleton::mpi().getRank()) == destRanks.end()) {
      //  contact.resetMinMax();
      //}
    }
  }

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  std::map<int, std::vector<std::uint8_t>> rankDataMapSorted;
  communication::fillSendBuffer(dataMap, rankDataMapSorted, serialSize);

  //Send mapped data
  communication::sendMappedData(rankDataMapSorted, listNeighbourRanks,
                                serialSize, contactDetectionComm, mpiNbHelper);

  //Receive data and iterate buffer
  communication::receiveAndExecuteForData(
      listNeighbourRanks, serialSize, contactDetectionComm, mpiNbHelper,
      [&](int rankOrig, std::uint8_t* buffer) {
        CONTACTTYPE newcontact;
        newcontact.deserialize(buffer);

#ifdef VERBOSE_CONTACT_COMMUNICATION
        std::cout << "From Rank " << rankOrig << ": ";
#endif

        receiveContact(newcontact, sParticleSystem, contacts, rankOrig);
      });
#endif
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE>
void communicateContactsParallel(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
#ifdef PARALLEL_MODE_MPI
    MPI_Comm particleContactDetectionComm, MPI_Comm wallContactDetectionComm,
#endif
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "contactCommunication");

  communicateParallelContacts(contactContainer.wallContacts, sParticleSystem,
                              deltaX, wallContactDetectionComm, periodicity);

  communicateParallelContacts(contactContainer.particleContacts,
                              sParticleSystem, deltaX,
                              particleContactDetectionComm, periodicity);
#endif
}

template <typename T, typename PARTICLETYPE, typename CONTACTTYPE>
void communicatePostContactTreatmentContacts(
    std::vector<CONTACTTYPE>&             contacts,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
#ifdef PARALLEL_MODE_MPI
    MPI_Comm postContactTreatmentComm,
#endif
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "postContactTreatmentCommunication");

  // Retrieve rank of direct neighbours
  auto& listNeighbourRanks         = sParticleSystem.getNeighbourRanks();
  constexpr std::size_t serialSize = CONTACTTYPE::getSerialSize();

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  for (CONTACTTYPE& contact : contacts) {
    if (contact.getResponsibleRank() == singleton::mpi().getRank()) {
      if (!contact.isEmpty()) {
        std::set<int> destRanks {
            evalDestRanksForPostContactTreatmentCommunication(
                contact, sParticleSystem, deltaX, periodicity)};
        resetResponsibleRank(contact);
        for (const int destRank : destRanks) {
          if (destRank != singleton::mpi().getRank()) {
            std::unique_ptr<std::uint8_t[]> buffer(
                new std::uint8_t[serialSize] {});
            contact.serialize(buffer.get());
            dataMap.insert(std::make_pair(destRank, std::move(buffer)));
          }
        }
      }
    }
    // Otherwise reset min and max as they are not accurate anymore
    else {
      contact.resetMinMax();
    }
  }

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  std::map<int, std::vector<std::uint8_t>> rankDataMapSorted;
  communication::fillSendBuffer(dataMap, rankDataMapSorted, serialSize);

  //Send mapped data
  communication::sendMappedData(rankDataMapSorted, listNeighbourRanks,
                                serialSize, postContactTreatmentComm,
                                mpiNbHelper);

  // Remove all empty contacts (= all contacts that were empty and all contacts where the process wasn't responsible) to ensure that 100% of the data send by the responsible rank is used
  cleanContacts(contacts);

  //Receive data and iterate buffer
  communication::receiveAndExecuteForData(
      listNeighbourRanks, serialSize, postContactTreatmentComm, mpiNbHelper,
      [&](int rankOrig, std::uint8_t* buffer) {
        CONTACTTYPE newcontact;
        newcontact.deserialize(buffer);
        receiveContact(newcontact, sParticleSystem, contacts, rankOrig);
      });
#endif
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE>
/// It is necessary to communicate contacts so that no information is lost
void communicatePostContactTreatmentContacts(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem, const T deltaX,
#ifdef PARALLEL_MODE_MPI
    MPI_Comm particlePostContactTreatmentComm,
    MPI_Comm wallPostContactTreatmentComm,
#endif
    const Vector<bool, PARTICLETYPE::d>& periodicity)
{
#ifdef PARALLEL_MODE_MPI
  OstreamManager clout(std::cout, "postContactTreatmentCommunication");

  communicatePostContactTreatmentContacts(
      contactContainer.wallContacts, sParticleSystem, deltaX,
      wallPostContactTreatmentComm, periodicity);

  communicatePostContactTreatmentContacts(
      contactContainer.particleContacts, sParticleSystem, deltaX,
      particlePostContactTreatmentComm, periodicity);
#endif
}

// Contact force communication
#ifdef PARALLEL_MODE_MPI
template <typename T, unsigned D>
int getContactTreatmentResultsSerialSize()
{
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;

  // Create communicatables
  std::array<std::size_t, 1> commGlobalParticleID {1};
  ConcreteCommunicatable     communicatableID =
      ConcreteCommunicatable(commGlobalParticleID);
  Vector<T, D>           forceDummy(1);
  ConcreteCommunicatable communicatableForce =
      ConcreteCommunicatable(forceDummy);
  Vector<T, Drot>        torqueDummy(1);
  ConcreteCommunicatable communicatableTorque =
      ConcreteCommunicatable(torqueDummy);

  // Setting dimension dependent indices
  const std::vector<unsigned> indicesID {0};
  std::vector<unsigned>       indicesForce;
  if constexpr (D == 3) {
    indicesForce = std::vector<unsigned> {0, 1, 2};
  }
  else {
    indicesForce = std::vector<unsigned> {0, 1};
  }
  std::vector<unsigned> indicesTorque;
  if constexpr (D == 3) {
    indicesTorque = indicesForce;
  }
  else {
    indicesTorque = std::vector<unsigned> {0};
  }

  // Retrieve serial size
  return communicatableID.size(indicesID) +
         communicatableForce.size(indicesForce) +
         communicatableTorque.size(indicesTorque);
}
#endif

template <typename T, unsigned D>
void extendContactTreatmentResultsDataMap(
    const std::size_t& globalParticleID, Vector<T, D>& force,
    Vector<T, utilities::dimensions::convert<D>::rotation>& torque,
    const int                                               destRank,
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>&    dataMap)
{
#ifdef PARALLEL_MODE_MPI
  // Create communicatables
  std::array<std::size_t, 1> commGlobalParticleID {globalParticleID};
  ConcreteCommunicatable     communicatableID =
      ConcreteCommunicatable(commGlobalParticleID);
  ConcreteCommunicatable communicatableForce  = ConcreteCommunicatable(force);
  ConcreteCommunicatable communicatableTorque = ConcreteCommunicatable(torque);

  // Setting dimension dependent indices
  const std::vector<unsigned> indicesID {0};
  std::vector<unsigned>       indicesForce;
  if constexpr (D == 3) {
    indicesForce = std::vector<unsigned> {0, 1, 2};
  }
  else {
    indicesForce = std::vector<unsigned> {0, 1};
  }
  std::vector<unsigned> indicesTorque;
  if constexpr (D == 3) {
    indicesTorque = indicesForce;
  }
  else {
    indicesTorque = std::vector<unsigned> {0};
  }

  // Retrieve serial size
  const std::size_t serialSize = communicatableID.size(indicesID) +
                                 communicatableForce.size(indicesForce) +
                                 communicatableTorque.size(indicesTorque);

  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSize] {});
  std::uint8_t*                   bufferRaw = buffer.get();
  std::size_t serialIdx = communicatableID.serialize(indicesID, bufferRaw);
  serialIdx +=
      communicatableForce.serialize(indicesForce, &bufferRaw[serialIdx]);
  serialIdx +=
      communicatableTorque.serialize(indicesTorque, &bufferRaw[serialIdx]);

  dataMap.insert(std::make_pair(destRank, std::move(buffer)));
#ifdef VERBOSE_CONTACT_COMMUNICATION
  std::cout << "Rank " << singleton::mpi().getRank()
            << " sends contact treatment results of particle "
            << globalParticleID << " to " << destRank << std::endl;
#endif
#endif
}

#ifdef PARALLEL_MODE_MPI
template <typename T, typename PARTICLETYPE>
void receiveContactTreatmentResults(
    SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem,
    MPI_Comm contactTreatmentComm, singleton::MpiNonBlockingHelper& mpiNbHelper)
{
  using namespace descriptors;

  // Retrieve dimensions
  constexpr unsigned D    = PARTICLETYPE::d;
  constexpr unsigned Drot = utilities::dimensions::convert<D>::rotation;

  // Retrieve rank of direct neighbours
  auto& listNeighbourRanks = sParticleSystem.getNeighbourRanks();

#ifdef VERBOSE_CONTACT_COMMUNICATION
  std::cout << "Rank " << singleton::mpi().getRank() << " has the neighbours: ";
  for (auto& neighbourRank : listNeighbourRanks) {
    std::cout << neighbourRank << " ";
  }
  std::cout << std::endl;
#endif

  // Setting dimension dependent indices
  const std::vector<unsigned> indicesID {0};
  std::vector<unsigned>       indicesForce;
  if constexpr (D == 3) {
    indicesForce = std::vector<unsigned> {0, 1, 2};
  }
  else {
    indicesForce = std::vector<unsigned> {0, 1};
  }
  std::vector<unsigned> indicesTorque;
  if constexpr (D == 3) {
    indicesTorque = indicesForce;
  }
  else {
    indicesTorque = std::vector<unsigned> {0};
  }

  std::array<std::size_t, 1> commGlobalParticleID;
  Vector<T, D>               contactForce;
  Vector<T, Drot>            contactTorque;

  // Create communicatables
  ConcreteCommunicatable communicatableID =
      ConcreteCommunicatable(commGlobalParticleID);
  ConcreteCommunicatable communicatableForce =
      ConcreteCommunicatable(contactForce);
  ConcreteCommunicatable communicatableTorque =
      ConcreteCommunicatable(contactTorque);

  // Retrieve serial size
  const std::size_t serialSize = communicatableID.size(indicesID) +
                                 communicatableForce.size(indicesForce) +
                                 communicatableTorque.size(indicesTorque);

  // Receive data and iterate buffer
  // TODO: Improve processing of data, do not iterate through the all particles every time (maybe also add main iC to the communicated data to reduce processed particles)
  communication::receiveAndExecuteForData(
      listNeighbourRanks, serialSize, contactTreatmentComm, mpiNbHelper,
      [&](int rankOrig, std::uint8_t* buffer) {
        std::size_t serialIdx = communicatableID.deserialize(indicesID, buffer);
        serialIdx +=
            communicatableForce.deserialize(indicesForce, &buffer[serialIdx]);
        serialIdx +=
            communicatableTorque.deserialize(indicesTorque, &buffer[serialIdx]);

#ifdef VERBOSE_CONTACT_COMMUNICATION
        std::cout << "Rank " << singleton::mpi().getRank()
                  << " received contact treatment results for particle "
                  << commGlobalParticleID[0] << " (force = " << contactForce
                  << "; torque = " << contactTorque << ") from rank "
                  << rankOrig << std::endl;
#endif

        communication::forParticlesInSuperParticleSystem<
            T, PARTICLETYPE, conditions::valid_particles>(
            sParticleSystem,
            [&](Particle<T, PARTICLETYPE>&       particle,
                ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
              const std::size_t currentGlobalParticleID =
                  particle.template getField<PARALLELIZATION, ID>();
              if (commGlobalParticleID[0] == currentGlobalParticleID) {
                // Process contact force
                const Vector<T, D> force =
                    particle.template getField<FORCING, FORCE>();
                particle.template setField<FORCING, FORCE>(force +
                                                           contactForce);
                // Process torque from contact
                const Vector<T, Drot> torque(
                    particle.template getField<FORCING, TORQUE>());
                particle.template setField<FORCING, TORQUE>(
                    utilities::dimensions::convert<D>::serialize_rotation(
                        torque + contactTorque));
              }
            });
      });
}
#endif

template <typename T, typename PARTICLETYPE>
void communicateContactTreatmentResults(
    SuperParticleSystem<T, PARTICLETYPE>&                particleSystem,
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap
#ifdef PARALLEL_MODE_MPI
    ,
    MPI_Comm contactTreatmentComm
#endif
)
{
#ifdef PARALLEL_MODE_MPI
  if constexpr (access::providesParallelization<PARTICLETYPE>()) {

    std::size_t serialSize =
        getContactTreatmentResultsSerialSize<T, PARTICLETYPE::d>();

    //Create non blocking mpi helper
    singleton::MpiNonBlockingHelper mpiNbHelper;

    std::map<int, std::vector<std::uint8_t>> rankDataMapSorted;
    communication::fillSendBuffer(dataMap, rankDataMapSorted, serialSize);

    // Send mapped data
    communication::sendMappedData(
        rankDataMapSorted, particleSystem.getNeighbourRanks(), serialSize,
        contactTreatmentComm, mpiNbHelper);

    receiveContactTreatmentResults(particleSystem, contactTreatmentComm,
                                   mpiNbHelper);
  }
#endif
}

} //namespace contact

} //namespace particles

} //namespace olb

#endif
