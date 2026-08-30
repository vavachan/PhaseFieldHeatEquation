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

#ifndef WALL_CONTACT_HH
#define WALL_CONTACT_HH

#include "wallContact.h"

namespace olb {
namespace particles {
namespace contact {

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::WallContactArbitraryFromOverlapVolume()
    : WallContactArbitraryFromOverlapVolume<T, D, CONVEX>(std::size_t(),
                                                          unsigned())
{}

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::WallContactArbitraryFromOverlapVolume(std::size_t particleID,
                                                         unsigned    wallID)
    : particleID(std::array<std::size_t, 1>({particleID}))
    , wallID(std::array<unsigned, 1>({wallID}))

{
  resetMinMax();
}

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    WallContactArbitraryFromOverlapVolume(
        const WallContactArbitraryFromOverlapVolume<T, D, CONVEX>& contact)
{
  min                        = contact.min;
  max                        = contact.max;
  particleID[0]              = contact.particleID[0];
  wallID[0]                  = contact.wallID[0];
  particlePosition           = contact.particlePosition;
  particlePositionUpdated[0] = contact.isParticlePositionUpdated();
  dampingFactor[0]           = contact.getDampingFactor();
  newContact[0]              = contact.isNew();
  responsibleRank[0]         = contact.getResponsibleRank();
}

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    WallContactArbitraryFromOverlapVolume(
        WallContactArbitraryFromOverlapVolume<T, D, CONVEX>&& contact)
{
  min                     = std::move(contact.min);
  max                     = std::move(contact.max);
  particleID              = std::move(contact.particleID);
  wallID                  = std::move(contact.wallID);
  particlePosition        = std::move(contact.particlePosition);
  particlePositionUpdated = std::move(contact.particlePositionUpdated);
  dampingFactor           = std::move(contact.dampingFactor);
  newContact              = std::move(contact.newContact);
  responsibleRank         = std::move(contact.responsibleRank);
}

template <typename T, unsigned D, bool CONVEX>
constexpr const std::size_t&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getParticleID() const
{
  return this->particleID[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr unsigned
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getWallID() const
{
  return this->wallID[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getParticlePosition() const
{
  return particlePosition;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::setParticlePosition(
    const PhysR<T, D>& position)
{
  particlePosition           = position;
  particlePositionUpdated[0] = true;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::setResponsibleRank(
    const int& rank)
{
  responsibleRank[0] = rank;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const int&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getResponsibleRank() const
{
  return responsibleRank[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getMin() const
{
  return this->min;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getMax() const
{
  return this->max;
}

template <typename T, unsigned D, bool CONVEX>
constexpr T
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::getDampingFactor() const
{
  return dampingFactor[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::setDampingFactor(
    const T newDampingFactor)
{
  dampingFactor[0] = newDampingFactor;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    setDampingFactorFromInitialVelocity(
        const T coefficientOfRestitution,
        const T initialRelativeVelocityMagnitude)
{
  setDampingFactor(particles::contact::evalDampingFactor(
      coefficientOfRestitution, initialRelativeVelocityMagnitude));
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::resetMinMax()
{
  for (unsigned iD = 0; iD < D; ++iD) {
    min[iD] = std::numeric_limits<olb::BaseType<T>>::max();
    max[iD] = -std::numeric_limits<olb::BaseType<T>>::max();
  }
  particlePositionUpdated[0] = false;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::updateMinMax(
    const PhysR<T, D>& positionInsideTheContact)
{
  particles::contact::updateMinMax(this->min, this->max,
                                   positionInsideTheContact);
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::increaseMinMax(
    const Vector<T, D>& increaseBy)
{
  this->max += increaseBy;
  this->min -= increaseBy;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::combineWith(
    WallContactArbitraryFromOverlapVolume<T, D, CONVEX>& contact)
{
  if (particleID[0] == contact.particleID[0] &&
      wallID[0] == contact.wallID[0]) {
    if(contact.isNew() == newContact[0]) {
      // Determine bounding box of combined overlap area
      for (unsigned iD = 0; iD < D; ++iD) {
        min[iD] = util::min(min[iD], contact.min[iD]);
        max[iD] = util::max(max[iD], contact.max[iD]);
      }
      // The damping factors should be either the same or one is -1 and the other has the correct value which is > 0
      dampingFactor[0] = util::max(dampingFactor[0], contact.getDampingFactor());
      responsibleRank[0] =
          util::min(responsibleRank[0], contact.getResponsibleRank());
      particlePositionUpdated[0] =
          particlePositionUpdated[0] || contact.isParticlePositionUpdated();
      if (!particlePositionUpdated[0] && contact.isParticlePositionUpdated()) {
        particlePosition = contact.getParticlePosition();
      }
    }
    else {
      if(!contact.isNew()) {
        min = contact.getMin();
        max = contact.getMax();
        dampingFactor[0] = contact.getDampingFactor();
        particlePositionUpdated[0] = contact.isParticlePositionUpdated();
        particlePosition = contact.getParticlePosition();
      }
      responsibleRank[0] =
          util::min(responsibleRank[0], contact.getResponsibleRank());
      newContact[0] = newContact[0] && contact.isNew();
    }

    // Ignore second contact by resetting it
    contact.resetMinMax();
  }
}

template <typename T, unsigned D, bool CONVEX>
constexpr bool
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::isEmpty() const
{
  for (unsigned iD = 0; iD < D; ++iD) {
    if (min[iD] > max[iD]) {
      return true;
    }
  }
  return false;
}

template <typename T, unsigned D, bool CONVEX>
constexpr bool
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::isNew() const
{
  return newContact[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::isNew(
    const bool newContact)
{
  this->newContact[0] = newContact;
}

template <typename T, unsigned D, bool CONVEX>
constexpr bool
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::isParticlePositionUpdated()
    const
{
  return particlePositionUpdated[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::setParticlePositionUpdated(
    bool updated)
{
  particlePositionUpdated[0] = updated;
}

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::operator=(
    const WallContactArbitraryFromOverlapVolume<T, D, CONVEX>& contact)
{
  min                        = contact.min;
  max                        = contact.max;
  particleID[0]              = contact.particleID[0];
  wallID[0]                  = contact.wallID[0];
  particlePosition           = contact.particlePosition;
  particlePositionUpdated[0] = contact.isParticlePositionUpdated();
  dampingFactor[0]           = contact.getDampingFactor();
  newContact[0]              = contact.isNew();
  responsibleRank[0]         = contact.getResponsibleRank();

  return *this;
}

template <typename T, unsigned D, bool CONVEX>
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>&
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::operator=(
    WallContactArbitraryFromOverlapVolume<T, D, CONVEX>&& contact)
{
  min                     = std::move(contact.min);
  max                     = std::move(contact.max);
  particleID              = std::move(contact.particleID);
  wallID                  = std::move(contact.wallID);
  particlePosition        = std::move(contact.particlePosition);
  particlePositionUpdated = std::move(contact.particlePositionUpdated);
  dampingFactor           = std::move(contact.dampingFactor);
  newContact              = std::move(contact.newContact);
  responsibleRank         = std::move(contact.responsibleRank);

  return *this;
}

template <typename T, unsigned D, bool CONVEX>
template <typename F>
std::size_t
WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::processWithCommunicatables(
    F f)
{
  //Create communicatables
  auto communicatablePosition = ConcreteCommunicatable(particlePosition);
  auto communicatableMin      = ConcreteCommunicatable(min);
  auto communicatableMax      = ConcreteCommunicatable(max);
  auto communicatableID       = ConcreteCommunicatable(particleID);
  auto communicatableMaterial = ConcreteCommunicatable(wallID);
  auto communicatableDamping  = ConcreteCommunicatable(dampingFactor);
  auto communicatableParticlePositionUpdated =
      ConcreteCommunicatable(particlePositionUpdated);
  auto communicatableIsNew = ConcreteCommunicatable(newContact);
  auto communicatableRank  = ConcreteCommunicatable(responsibleRank);

  return f(communicatablePosition, communicatableMin, communicatableMax,
           communicatableID, communicatableMaterial, communicatableDamping,
           communicatableParticlePositionUpdated, communicatableIsNew,
           communicatableRank);
}

template <typename T, unsigned D, bool CONVEX>
std::size_t WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::serialize(
    std::uint8_t* buffer)
{
  return processWithCommunicatables(
      [&](auto& communicatablePosition, auto& communicatableMin,
          auto& communicatableMax, auto& communicatableID,
          auto& communicatableMaterial, auto& communicatableDamping,
          auto& communicatableParticlePositionUpdated,
          auto& communicatableIsNew, auto& communicatableRank) {
        std::size_t serialIdx =
            communicatablePosition.serialize(this->indicesDim, buffer);
        serialIdx +=
            communicatableMin.serialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableMax.serialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableID.serialize(this->indicesSingle, &buffer[serialIdx]);
        serialIdx += communicatableMaterial.serialize(this->indicesSingle,
                                                      &buffer[serialIdx]);
        serialIdx += communicatableDamping.serialize(this->indicesSingle,
                                                     &buffer[serialIdx]);
        serialIdx += communicatableParticlePositionUpdated.serialize(
            this->indicesSingle, &buffer[serialIdx]);
        serialIdx += communicatableIsNew.serialize(this->indicesSingle,
                                                   &buffer[serialIdx]);
        serialIdx += communicatableRank.serialize(this->indicesSingle,
                                                  &buffer[serialIdx]);

        return serialIdx;
      });
}

template <typename T, unsigned D, bool CONVEX>
std::size_t WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::deserialize(
    std::uint8_t* buffer)
{
  return processWithCommunicatables(
      [&](auto& communicatablePosition, auto& communicatableMin,
          auto& communicatableMax, auto& communicatableID,
          auto& communicatableMaterial, auto& communicatableDamping,
          auto& communicatableParticlePositionUpdated,
          auto& communicatableIsNew, auto& communicatableRank) {
        std::size_t serialIdx =
            communicatablePosition.deserialize(this->indicesDim, buffer);
        serialIdx +=
            communicatableMin.deserialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableMax.deserialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx += communicatableID.deserialize(this->indicesSingle,
                                                  &buffer[serialIdx]);
        serialIdx += communicatableMaterial.deserialize(this->indicesSingle,
                                                        &buffer[serialIdx]);
        serialIdx += communicatableDamping.deserialize(this->indicesSingle,
                                                       &buffer[serialIdx]);
        serialIdx += communicatableParticlePositionUpdated.deserialize(
            this->indicesSingle, &buffer[serialIdx]);
        serialIdx += communicatableIsNew.deserialize(this->indicesSingle,
                                                     &buffer[serialIdx]);
        serialIdx += communicatableRank.deserialize(this->indicesSingle,
                                                    &buffer[serialIdx]);

        return serialIdx;
      });
}

template <typename T, unsigned D, bool CONVEX>
void WallContactArbitraryFromOverlapVolume<T, D, CONVEX>::print()
{
  OstreamManager clout(std::cout, "WallContact");
  clout.setMultiOutput(true);
#ifdef PARALLEL_MODE_MPI
  int rank = singleton::mpi().getRank();
  if (rank == responsibleRank[0]) {
#endif
    clout << "Min=" << this->min << ", Max=" << this->max << std::endl;
    clout << "particle ID=" << this->particleID[0]
          << ", wall ID=" << this->wallID[0]
          << ", DampingFactor=" << dampingFactor[0] << std::endl;
    clout << "Position=" << particlePosition << std::endl;
#ifdef PARALLEL_MODE_MPI
  }
#endif
  clout.setMultiOutput(false);
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
