/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2020 Mathias Krause
 *                2020      Adrian Kummerlaender
 *                2021      Nicolas Hafen
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

#ifndef BLOCK_LATTICE_INTERACTION_HH
#define BLOCK_LATTICE_INTERACTION_HH

#include "blockLatticeInteraction.h"

namespace olb {

namespace particles {

template <typename T, unsigned D>
bool getBlockParticleIntersection(const BlockGeometry<T, D>& blockGeometry,
                                  T invDeltaX, LatticeR<D>& start,
                                  LatticeR<D>& end, PhysR<T, D> position,
                                  T circumRadius)
{
  /// Setting block bounds excluding(!) padding
  LatticeR<D> blockMin(0);
  LatticeR<D> blockMax(blockGeometry.getExtent() - 1);

  /// Calculate relative start/end in block domain
  PhysR<T, D> particleMin(position - circumRadius);
  PhysR<T, D> relStart(particleMin - blockGeometry.getOrigin());
  PhysR<T, D> particleMax(position + circumRadius);
  PhysR<T, D> relEnd(particleMax - blockGeometry.getOrigin());

  /// Set latticeR start/end and check validity
  bool validRange = true;
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    start[iDim] =
        util::max(util::floor(invDeltaX * relStart[iDim]), blockMin[iDim]);
    end[iDim] = util::min(util::ceil(invDeltaX * relEnd[iDim]), blockMax[iDim]);
    int intersectionRange = end[iDim] - start[iDim];
    validRange            = validRange && (intersectionRange >= 0.);
  }

  return validRange;
}

template <typename T, unsigned D>
void checkSmoothIndicatorOutOfGeometry(
    bool& surfaceOutOfGeometry, PhysR<T, D>& ghostPos,
    const PhysR<T, D>& cellMin, const PhysR<T, D>& cellMax,
    const PhysR<T, D>& position, T circumRadius,
    const Vector<bool, D>& periodic)
{
  ghostPos             = position;
  surfaceOutOfGeometry = false;
  for (unsigned i = 0; i < D; ++i) {
    T const particleMin = position[i] - circumRadius;
    T const particleMax = position[i] + circumRadius;
    if (particleMin < cellMin[i] && periodic[i]) {
      surfaceOutOfGeometry = true;
      ghostPos[i] = communication::movePositionToEnd(
          position[i], cellMax[i], cellMin[i]);
    }
    else if (particleMax > cellMax[i] && periodic[i]) {
      surfaceOutOfGeometry = true;
      ghostPos[i] = communication::movePositionToStart(
          position[i], cellMax[i], cellMin[i]);
    }
  }
}

////TODO: periodic treatment not tested with new logic yet!
////TODO: rework logic, as it can most likely be simplified due to the new pos,rot,circumRadius passing
template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setBlockParticleField(const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
                           BlockLattice<T, DESCRIPTOR>&           blockLattice,
                           UnitConverter<T, DESCRIPTOR> const&    converter,
                           const PhysR<T, DESCRIPTOR::d>&         cellMin,
                           const PhysR<T, DESCRIPTOR::d>&         cellMax,
                           Particle<T, PARTICLETYPE>&             particle,
                           const Vector<bool, DESCRIPTOR::d>&     periodic)
{
  constexpr unsigned D = DESCRIPTOR::d;
  using namespace descriptors;
  const T           circumRadius = particles::access::getRadius(particle);
  const PhysR<T, D> originalPosition(particles::access::getPosition(particle));
  // Checking if the particle's surface leaves the domain
  bool        surfaceOutOfGeometry;
  PhysR<T, D> ghostPos;
  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, originalPosition, circumRadius,
                                    periodic);

  //Do the normal routine if the particle is in the geometry
  if (!surfaceOutOfGeometry) {
    setBlockParticleField(blockGeometry, blockLattice, converter, particle);
  }
  else {
    //sets the Particle to ghost position on the other side of the domain and sets the field
    particle.template setField<GENERAL, POSITION>(ghostPos);
    setBlockParticleField(blockGeometry, blockLattice, converter, particle);
    //Reverting Particle to its Previous position and setting the field
    particle.template setField<GENERAL, POSITION>(originalPosition);
    setBlockParticleField(blockGeometry, blockLattice, converter, particle);

    // For parallized particles, the update of the position is processed during relocation
    if constexpr(!particles::access::providesParallelization<PARTICLETYPE>()) {
      // Check if center of mass left the simulation domain
      bool centerOfMassOutOfDomain = false;
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        centerOfMassOutOfDomain = centerOfMassOutOfDomain ||
                                  originalPosition[iDim] < cellMin[iDim] ||
                                  originalPosition[iDim] > cellMax[iDim];
      }

      // Replace the position of the center of mass with the ghost position
      // (If the center of mass is outside the domain)
      if (centerOfMassOutOfDomain) {
        particle.template setField<GENERAL, POSITION>(ghostPos);
      }
    }
  }
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE,
          typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE, typename F>
void setBlockParticleField(
    const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T, DESCRIPTOR>&           blockLattice,
    UnitConverter<T, DESCRIPTOR> const&    converter,
    const PhysR<T, DESCRIPTOR::d>&         cellMin,
    const PhysR<T, DESCRIPTOR::d>&         cellMax,
    ParticleSystem<T, PARTICLETYPE>&       particleSystem,
    contact::ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
           contactContainer,
    size_t iP, Particle<T, PARTICLETYPE>& particle,
    std::vector<SolidBoundary<T, DESCRIPTOR::d>>& solidBoundaries,
    F                                                   getSetupPeriodicity)
{
  constexpr unsigned D = DESCRIPTOR::d;
  using namespace descriptors;
  const T circumRadius = particles::access::getRadius(particle);
  // Checking if the Particle leaves the Domain
  bool        surfaceOutOfGeometry;
  PhysR<T, D> originalPosition(particles::access::getPosition(particle));
  PhysR<T, D> ghostPos;
  checkSmoothIndicatorOutOfGeometry(
      surfaceOutOfGeometry, ghostPos, cellMin, cellMax, originalPosition,
      circumRadius, getSetupPeriodicity());

  //Do the normal routine if the particle is in the geometry
  if (!surfaceOutOfGeometry) {
    setBlockParticleField(blockGeometry, blockLattice, converter,
                          particleSystem, contactContainer, iP, particle,
                          solidBoundaries);
  }
  else {
    //sets the Particle to ghost position on the other side of the domain and sets the field
    particle.template setField<GENERAL, POSITION>(ghostPos);
    setBlockParticleField(blockGeometry, blockLattice, converter,
                          particleSystem, contactContainer, iP, particle,
                          solidBoundaries, cellMin, cellMax,
                          getSetupPeriodicity);
    //Reverting Particle to its Previous position and setting the field
    particle.template setField<GENERAL, POSITION>(originalPosition);
    setBlockParticleField(blockGeometry, blockLattice, converter,
                          particleSystem, contactContainer, iP, particle,
                          solidBoundaries, cellMin, cellMax,
                          getSetupPeriodicity);
  }
}

//Interation for block particle intersection for provided lambda function F
template <typename T, typename DESCRIPTOR, typename F>
void forSpatialLocationsInBlockParticleIntersection(
    const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T, DESCRIPTOR>& blockLattice, int padding,
    Vector<T, DESCRIPTOR::d> position, T circumRadius, F f)
{
  constexpr unsigned D = DESCRIPTOR::d;
  LatticeR<D>        start;
  LatticeR<D>        end;
  if (getBlockParticleIntersection(blockGeometry,
                                   T {1} / blockGeometry.getDeltaR(), start,
                                   end, position, circumRadius)) {
    //Extend range by padding if desired
    start -= LatticeR<D>(padding);
    end += LatticeR<D>(padding);
    //For all cells in subdomain defined by start/end
    blockLattice.forSpatialLocations(start, end,
                                     [&](const LatticeR<D>& latticeR) {
                                       f(latticeR);
                                     });
  }
}

//TODO: remove material 1 from hardcoded version
//- Can be done with unordered_set, analogous to other version below
template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setBlockParticleField(const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
                           BlockLattice<T, DESCRIPTOR>&           blockLattice,
                           UnitConverter<T, DESCRIPTOR> const&    converter,
                           Particle<T, PARTICLETYPE>&             particle)
{
  constexpr unsigned D = DESCRIPTOR::d;

  using namespace descriptors;
  auto circumRadius =
      particle.template getField<SURFACE, SINDICATOR>()->getCircumRadius();
  auto position = particle.template getField<GENERAL, POSITION>();

  //For all cells in block particle intersection
  int padding = blockGeometry.getPadding();
  forSpatialLocationsInBlockParticleIntersection(
      blockGeometry, blockLattice, padding, position, circumRadius,
      [&](const LatticeR<D>& latticeR) {
        //Get phys position
        T physR[D] = {};
        Vector<T, D> physRV;
        blockGeometry.getPhysR(physRV, latticeR);
        for (unsigned iD = 0; iD != D; ++iD) {
          physR[iD] = physRV[iD];
        }
        //Retrieve porosity
        FieldD<T, DESCRIPTOR, descriptors::POROSITY> porosity {};
        particles::resolved::evalSolidVolumeFraction(porosity.data(), physR, particle);

        //For cells containing particle bits
        if (!util::nearZero(porosity[0])) {
          auto material = blockGeometry.get(latticeR);
          //TODO: remove material 1 from hardcoded version
          if (material == 1) {
            //Retrieve local velocity
            FieldD<T, DESCRIPTOR, descriptors::VELOCITY> velocity {};
            velocity = particles::dynamics::calculateLocalVelocity(
                particle, PhysR<T, D>(physR));
            for (unsigned iDim = 0; iDim != D; ++iDim) {
              velocity[iDim] = converter.getLatticeVelocity(velocity[iDim]);
            }
            //Calculate solid volume fraction
            T solidVolumeFraction = 1. - porosity[0];
            //Write to fields
            {
              auto cell = blockLattice.get(latticeR);
              if constexpr (!DESCRIPTOR::template provides<descriptors::VELOCITY_SOLID>()) {
                cell.template getFieldPointer<
                    descriptors::VELOCITY_NUMERATOR>() += porosity[0] * velocity;
                cell.template getFieldPointer<
                    descriptors::VELOCITY_DENOMINATOR>() += porosity;
              }
              else {
                cell.template getFieldPointer<
                    descriptors::VELOCITY_SOLID>() += velocity;
              }
              cell.template getFieldPointer<descriptors::POROSITY>() *=
                  solidVolumeFraction;
            }
          } //if (material==1)
        }   //if (!util::nearZero(porosity[0]))
      });
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE,
          typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE, typename F>
void setBlockParticleField(
    const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T, DESCRIPTOR>&           blockLattice,
    UnitConverter<T, DESCRIPTOR> const&    converter,
    ParticleSystem<T, PARTICLETYPE>&       particleSystem,
    contact::ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
           contactContainer,
    size_t iP, Particle<T, PARTICLETYPE>& particle,
    std::vector<SolidBoundary<T, DESCRIPTOR::d>>& solidBoundaries,
    const PhysR<T, DESCRIPTOR::d>&                      cellMin,
    const PhysR<T, DESCRIPTOR::d>& cellMax, F getSetupPeriodicity)
{
  using namespace descriptors;
  constexpr unsigned D      = DESCRIPTOR::d;
  const T deltaX = blockGeometry.getDeltaR();
  auto    sIndicator   = particle.template getField<SURFACE, SINDICATOR>();
  const T detectionDistance = particles::contact::evalContactDetectionDistance(particle, deltaX);
  const T circumRadius = particles::contact::evalCircumRadius(
      detectionDistance, sIndicator->getCircumRadius(), sIndicator->getEpsilon());
  auto position = particle.template getField<GENERAL, POSITION>();
  //Retrieve wallMaterials
  std::set<int> wallMaterials = getLatticeMaterials( solidBoundaries );

  //For all cells in block particle intersection
  int padding = blockGeometry.getPadding();
  forSpatialLocationsInBlockParticleIntersection(
      blockGeometry, blockLattice, padding, position, circumRadius,
      [&](const LatticeR<D>& latticeR) {
        //Get phys position
        T physR[D] = {};
        Vector<T, D> physRV;
        blockGeometry.getPhysR(physRV, latticeR);
        for (unsigned iD = 0; iD != D; ++iD) {
          physR[iD] = physRV[iD];
        }
        T const signedDist = particles::resolved::signedDistanceToParticle(
            particle, PhysR<T, D>(physR));
        //Retrieve porosity
        FieldD<T, DESCRIPTOR, descriptors::POROSITY> porosity {};
        sdf::evalSolidVolumeFraction(porosity.data(), signedDist,
                          sIndicator->getEpsilon());

        //For cells containing particle bits
        if (!util::nearZero(porosity[0])) {
          //Retrieve material at cell and ensure location outside of boundary
          auto material = blockGeometry.getMaterial(latticeR);
          if (material != 0 && wallMaterials.find(material)==wallMaterials.end()) {
            //Retrieve local velocity
            FieldD<T, DESCRIPTOR, descriptors::VELOCITY> velocity {};
            velocity = particles::dynamics::calculateLocalVelocity(
                particle, PhysR<T, D>(physR));
            for (unsigned iDim = 0; iDim != D; ++iDim) {
                velocity[iDim] = converter.getLatticeVelocity(velocity[iDim]);
            }
            //Calculate solid volume fraction
            T solidVolumeFraction = 1. - porosity[0];
            //Write to fields
            {
              auto cell = blockLattice.get(latticeR);
              if constexpr (!DESCRIPTOR::template provides<descriptors::VELOCITY_SOLID>()) {
                cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>() +=
                  porosity[0] * velocity;
                cell.template getFieldPointer<
                  descriptors::VELOCITY_DENOMINATOR>() += porosity;
              }
              else {
                cell.template getFieldPointer<descriptors::VELOCITY_SOLID>() += velocity;
              }
              cell.template getFieldPointer<descriptors::POROSITY>() *=
                  solidVolumeFraction;
            }
          } //if (wallMaterials.find(material)==wallMaterials.end())
        } //if (!util::nearZero(porosity[0]))

        // For contacts we have another distance
        // Particle-particle detection (0.5 * deltaX layer for each equals deltaX which is then the same as for the wall detection)
        if (signedDist < detectionDistance) {
          // We must ignore cells on padding
          bool ignoreCell = false;
          for (unsigned i = 0; i < D; ++i) {
            if (latticeR[i] < 0 ||
                latticeR[i] > blockGeometry.getExtent()[i] - 1) {
              ignoreCell = true;
            }
          }

          // processing contact detection
          if (!ignoreCell) {
            auto contactHelper = blockLattice.get(latticeR)
                                       .template getFieldPointer<
                                           descriptors::CONTACT_DETECTION>();
            // this checks if another particle is already on this cell, if yes, the contacts are updated. If no, the ID (shifted by 1) is stored.
            std::size_t iP2 = contactHelper[0] - 1;
            if (contactHelper[0] > 0 && iP2 != iP) {
              std::size_t localiP2 = iP2;
              if constexpr (particles::access::providesParallelization<
                                PARTICLETYPE>()) {
                for (std::size_t localiP = 0; localiP < particleSystem.size();
                     ++localiP) {
                  const std::size_t globalParticleID =
                      particleSystem.get(localiP)
                          .template getField<PARALLELIZATION, ID>();
                  if (globalParticleID == iP2) {
                    localiP2 = localiP;
                    break;
                  }
                }
              }
              auto particle2 = particleSystem.get(localiP2);

              // Check that at least one particle should not ignore contacts
              bool detectParticleContacts = true;
              if constexpr (access::providesComputeContact<PARTICLETYPE>()) {
                detectParticleContacts =
                    access::isContactComputationEnabled(particle, particle2);
              }

              if (detectParticleContacts) {
                particles::contact::updateContacts(
                    contactContainer, std::array<std::size_t, 2>({iP2, iP}),
                    PhysR<T, D>(physR), particle2, particle, cellMin, cellMax,
                    getSetupPeriodicity, deltaX);
              }
            }
            else {
              contactHelper[0] = iP + 1;
            }

            // Check if contacts should be ignored
            bool detectWallContacts = true;
            if constexpr (access::providesComputeContact<PARTICLETYPE>()) {
              detectWallContacts =
                  access::isContactComputationEnabled(particle);
            }

            if (detectWallContacts) {
              for (unsigned short solidBoundaryID = 0; SolidBoundary<T, DESCRIPTOR::d>& solidBoundary : solidBoundaries){
                const T solidBoundaryDetectionDistance =
                    0.5 * util::sqrt(D) * deltaX +
                    solidBoundary.getEnlargementForContact();
                if (solidBoundary.getIndicator()->signedDistance(
                        PhysR<T, D>(physR)) < solidBoundaryDetectionDistance) {
                  particles::contact::updateContacts(
                      contactContainer, iP, solidBoundaryID, PhysR<T, D>(physR),
                      particle, cellMin, cellMax, getSetupPeriodicity, deltaX);
                }
                ++solidBoundaryID;
              }
            }
          }
        }
      });
}

template <typename T, typename DESCRIPTOR>
void resetBlockParticleField(
    const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T, DESCRIPTOR>&           blockLattice)
{
  constexpr unsigned D       = DESCRIPTOR::d;
  int                padding = blockGeometry.getPadding();
  LatticeR<D>        blockMin(0. - padding);
  LatticeR<D>        blockMax(blockGeometry.getExtent() - 1 + padding);

  using namespace descriptors;
  blockLattice.forSpatialLocations(
      blockMin, blockMax, [&](const LatticeR<D>& latticeR) {
        auto cell = blockLattice.get(latticeR);
        particles::resetAllParticleRelatedFields<DESCRIPTOR, decltype(cell), typename decltype(cell)::value_t>(cell);
      });
}

//TODO: HOTFIX ONLY
template <typename T, typename DESCRIPTOR>
void resetBlockContactField(
    const BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T, DESCRIPTOR>&           blockLattice)
{
  constexpr unsigned D       = DESCRIPTOR::d;
  int                padding = blockGeometry.getPadding();
  LatticeR<D>        blockMin(0. - padding);
  LatticeR<D>        blockMax(blockGeometry.getExtent() - 1 + padding);

  using namespace descriptors;
  blockLattice.forSpatialLocations(
      blockMin, blockMax, [&](LatticeR<D> latticeR) {
        auto cell = blockLattice.get(latticeR);
        particles::resetParticleContactRelatedFields<DESCRIPTOR, decltype(cell), typename decltype(cell)::value_t>(cell);
      });
}

} //namespace particles

} // namespace olb

#endif
