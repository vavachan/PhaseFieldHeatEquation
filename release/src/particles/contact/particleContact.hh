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

#ifndef PARTICLE_CONTACT_HH
#define PARTICLE_CONTACT_HH

#include "particleContact.h"
#include "particles/functions/particleContactDetectionFunctions.h"

namespace olb {
namespace particles {
namespace contact {

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::ParticleContactArbitraryFromOverlapVolume()
    : ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>(
          std::array<std::size_t, 2>())
{}

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    ParticleContactArbitraryFromOverlapVolume(
        const std::array<std::size_t, 2>& particleIDs)
    : ids(sortParticleIDs(particleIDs))
{
  resetMinMax();
}

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    ParticleContactArbitraryFromOverlapVolume(
        const ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc)
{
  min                        = pc.getMin();
  max                        = pc.getMax();
  ids                        = pc.getIDs();
  particlePositions          = pc.getParticlePositions();
  particlePositionUpdated[0] = pc.isParticlePositionUpdated();
  dampingFactor[0]           = pc.getDampingFactor();
  newContact[0]              = pc.isNew();
  responsibleRank[0]         = pc.getResponsibleRank();
}

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    ParticleContactArbitraryFromOverlapVolume(
        ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&& pc)
{
  min                     = std::move(pc.min);
  max                     = std::move(pc.max);
  ids                     = std::move(pc.ids);
  particlePositions       = std::move(pc.particlePositions);
  particlePositionUpdated = std::move(pc.particlePositionUpdated);
  dampingFactor           = std::move(pc.dampingFactor);
  newContact              = std::move(pc.newContact);
  responsibleRank         = std::move(pc.responsibleRank);
}

template <typename T, unsigned D, bool CONVEX>
constexpr const std::array<std::size_t, 2>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getIDs() const
{
  return this->ids;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const std::array<PhysR<T, D>, 2>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getParticlePositions()
    const
{
  return particlePositions;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::setParticlePositions(
    const std::array<PhysR<T, D>, 2>& positions)
{
  particlePositions          = positions;
  particlePositionUpdated[0] = true;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getParticlePosition(
    const std::size_t& id) const
{
  if (id == ids[0]) {
    return particlePositions[0];
  }
  return particlePositions[1];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::setParticlePosition(
    const PhysR<T, D>& position, const std::size_t& id)
{
  particlePositions[id]      = position;
  particlePositionUpdated[0] = true;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::setResponsibleRank(
    const int& rank)
{
  responsibleRank[0] = rank;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const int&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getResponsibleRank()
    const
{
  return responsibleRank[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getMin() const
{
  return this->min;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const PhysR<T, D>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getMax() const
{
  return this->max;
}

template <typename T, unsigned D, bool CONVEX>
constexpr const T
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::getDampingFactor()
    const
{
  return dampingFactor[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::setDampingFactor(
    const T newDampingFactor)
{
  dampingFactor[0] = newDampingFactor;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::
    setDampingFactorFromInitialVelocity(
        const T coefficientOfRestitution,
        const T initialRelativeVelocityMagnitude)
{
  setDampingFactor(particles::contact::evalDampingFactor(
      coefficientOfRestitution, initialRelativeVelocityMagnitude));
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::resetMinMax()
{
  for (unsigned iD = 0; iD < D; ++iD) {
    min[iD] = std::numeric_limits<olb::BaseType<T>>::max();
    max[iD] = -std::numeric_limits<olb::BaseType<T>>::max();
  }
  particlePositionUpdated[0] = false;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::updateMinMax(
    const PhysR<T, D>& positionInsideTheContact)
{
  particles::contact::updateMinMax(this->min, this->max,
                                   positionInsideTheContact);
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::increaseMinMax(
    const Vector<T, D>& increaseBy)
{
  this->max += increaseBy;
  this->min -= increaseBy;
}

template <typename T, unsigned D, bool CONVEX>
constexpr void
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::combineWith(
    ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc)
{
  if (particleContactConsistsOfIDs<
          ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>, true>(pc,
                                                                         ids)) {
    if(newContact[0] == pc.isNew()) {
      newContact[0] = newContact[0] && pc.isNew();

      // Determine bounding box of combined overlap area
      for (unsigned iD = 0; iD < D; ++iD) {
        min[iD] = util::min(min[iD], pc.getMin()[iD]);
        max[iD] = util::max(max[iD], pc.getMax()[iD]);
      }
      // The damping factors should be either the same or one is -1 and the other has the correct value which is > 0
      dampingFactor[0]   = util::max(dampingFactor[0], pc.getDampingFactor());
      responsibleRank[0] = util::min(responsibleRank[0], pc.getResponsibleRank());
      particlePositionUpdated[0] =
          particlePositionUpdated[0] || pc.isParticlePositionUpdated();
      if (!particlePositionUpdated[0] && pc.isParticlePositionUpdated()) {
        particlePositions = pc.getParticlePositions();
      }
    }
    else {
      if(!pc.isNew()) {
        min = pc.getMin();
        max = pc.getMax();
        dampingFactor[0] = pc.getDampingFactor();
        particlePositionUpdated[0] = pc.isParticlePositionUpdated();
        particlePositions = pc.getParticlePositions();
      }
      responsibleRank[0] = util::min(responsibleRank[0], pc.getResponsibleRank());
      newContact[0] = newContact[0] && pc.isNew();
    }

    // Ignore second contact by resetting it
    pc.resetMinMax();
  }
}

template <typename T, unsigned D, bool CONVEX>
constexpr bool
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::isEmpty() const
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
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::isNew() const
{
  return newContact[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::isNew(
    const bool newContact)
{
  this->newContact[0] = newContact;
}

template <typename T, unsigned D, bool CONVEX>
constexpr bool ParticleContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::isParticlePositionUpdated() const
{
  return particlePositionUpdated[0];
}

template <typename T, unsigned D, bool CONVEX>
constexpr void ParticleContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::setParticlePositionUpdated(bool updated)
{
  particlePositionUpdated[0] = updated;
}

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::operator=(
    const ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc)
{
  min                        = pc.getMin();
  max                        = pc.getMax();
  ids                        = pc.getIDs();
  particlePositions          = pc.getParticlePositions();
  particlePositionUpdated[0] = pc.isParticlePositionUpdated();
  dampingFactor[0]           = pc.getDampingFactor();
  newContact[0]              = pc.isNew();
  responsibleRank[0]         = pc.getResponsibleRank();

  return *this;
}

template <typename T, unsigned D, bool CONVEX>
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::operator=(
    ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&& pc)
{
  min                     = std::move(pc.min);
  max                     = std::move(pc.max);
  ids                     = std::move(pc.ids);
  particlePositions       = std::move(pc.particlePositions);
  particlePositionUpdated = std::move(pc.particlePositionUpdated);
  dampingFactor           = std::move(pc.dampingFactor);
  newContact              = std::move(pc.newContact);
  responsibleRank         = std::move(pc.responsibleRank);

  return *this;
}

template <typename T, unsigned D, bool CONVEX>
template <typename F>
std::size_t ParticleContactArbitraryFromOverlapVolume<
    T, D, CONVEX>::processWithCommunicatables(F f)
{
  //Create communicatables
  auto communicatablePositions = ConcreteCommunicatable(particlePositions);
  auto communicatableMin       = ConcreteCommunicatable(min);
  auto communicatableMax       = ConcreteCommunicatable(max);
  auto communicatableIDs       = ConcreteCommunicatable(ids);
  auto communicatableDamping   = ConcreteCommunicatable(dampingFactor);
  auto communicatableParticlePositionUpdated =
      ConcreteCommunicatable(particlePositionUpdated);
  auto communicatableIsNew = ConcreteCommunicatable(newContact);
  auto communicatableRank  = ConcreteCommunicatable(responsibleRank);

  return f(communicatablePositions, communicatableMin, communicatableMax,
           communicatableIDs, communicatableDamping,
           communicatableParticlePositionUpdated, communicatableIsNew,
           communicatableRank);
}

template <typename T, unsigned D, bool CONVEX>
std::size_t ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::serialize(
    std::uint8_t* buffer)
{
  return processWithCommunicatables(
      [&](auto& communicatablePositions, auto& communicatableMin,
          auto& communicatableMax, auto& communicatableIDs,
          auto& communicatableDamping,
          auto& communicatableParticlePositionUpdated,
          auto& communicatableIsNew, auto& communicatableRank) {
        std::size_t serialIdx =
            communicatablePositions.serialize(this->indicesPart, buffer);
        serialIdx +=
            communicatableMin.serialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableMax.serialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableIDs.serialize(this->indicesPart, &buffer[serialIdx]);
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
std::size_t
ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::deserialize(
    std::uint8_t* buffer)
{
  return processWithCommunicatables(
      [&](auto& communicatablePositions, auto& communicatableMin,
          auto& communicatableMax, auto& communicatableIDs,
          auto& communicatableDamping,
          auto& communicatableParticlePositionUpdated,
          auto& communicatableIsNew, auto& communicatableRank) {
        std::size_t serialIdx =
            communicatablePositions.deserialize(this->indicesPart, buffer);
        serialIdx +=
            communicatableMin.deserialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx +=
            communicatableMax.deserialize(this->indicesDim, &buffer[serialIdx]);
        serialIdx += communicatableIDs.deserialize(this->indicesPart,
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
void ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>::print()
{
  OstreamManager clout(std::cout, "ParticleContact");
  clout.setMultiOutput(true);
  int rank = singleton::mpi().getRank();
  if (rank == responsibleRank[0]) {
    clout << "Min=" << this->min << ", Max=" << this->max << std::endl;
    clout << "IDs=" << this->ids << ", DampingFactor=" << dampingFactor[0]
          << std::endl;
    clout << "Positions=" << particlePositions << std::endl;
  }
  clout.setMultiOutput(false);
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
