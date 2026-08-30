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

#ifndef PARTICLE_CONTACT_H
#define PARTICLE_CONTACT_H

#include "contactFunctions.h"
#include "contactHelpers.h"
#include "core/blockStructure.h"
#include "core/util.h"
#include <array>
#include <limits>
#include <utility>

namespace olb {
namespace particles {
namespace contact {

template <unsigned D>
struct ParticleContact {
  static_assert(D == 2 || D == 3, "Only D=2 and D=3 are supported");
  static const std::vector<unsigned int> indicesDim;
  static const std::vector<unsigned int> indicesPart;
  static const std::vector<unsigned int> indicesSingle;
};

template <>
const std::vector<unsigned int> ParticleContact<2>::indicesDim({0, 1});
template <>
const std::vector<unsigned int> ParticleContact<3>::indicesDim({0, 1, 2});
template <unsigned D>
const std::vector<unsigned int> ParticleContact<D>::indicesPart({0, 1});
template <unsigned D>
const std::vector<unsigned int> ParticleContact<D>::indicesSingle({0});

// Contact type to be used, if particle contact should be ignored (e.g., if only wall contact is of interest)
// TODO: struct IgnoreParticleContact : ParticleContact {};

/// An object holding data for a contact which is described analog to Nassauer and Kuna (2013)
template <typename T, unsigned D, bool CONVEX = true>
struct ParticleContactArbitraryFromOverlapVolume : ParticleContact<D> {
private:
  /// position of the particle (necessary for periodic setups)
  std::array<PhysR<T, D>, 2> particlePositions;
  /// minimal coordinates of the overlap area
  PhysR<T, D> min;
  /// maximal coordinates of the overlap area
  PhysR<T, D> max;
  /// IDs of particles in contact
  std::array<std::size_t, 2> ids;
  /// Damping factor calculated from the magnitude of the initial relative impact velocity in direction of the contact normal and the wanted coefficient of restitution
  std::array<T, 1> dampingFactor = {-1};
  /// Helper variable to identify if the particle positions are up to date (should be resetted in the contact force calculation)
  std::array<bool, 1> particlePositionUpdated = {false};
  /// indicates if the contact was just found
  std::array<bool, 1> newContact = {true};
  /// Rank that is responsible for contact treatment
  std::array<int, 1> responsibleRank = {std::numeric_limits<int>::max()};

  /// Size of contained data
  static const std::size_t serialSize =
      sizeof(particlePositionUpdated) + sizeof(newContact) +
      sizeof(particlePositions) + sizeof(min) + sizeof(max) + sizeof(ids) +
      sizeof(responsibleRank) + sizeof(dampingFactor);

  template <typename F>
  std::size_t processWithCommunicatables(F f);

public:
  static_assert(D == 2 || D == 3, "Only D=2 and D=3 are supported");
  /// Constructor
  ParticleContactArbitraryFromOverlapVolume();
  /// Constructor
  ParticleContactArbitraryFromOverlapVolume(
      const std::array<std::size_t, 2>& particleIDs);
  /// Copy Constructor
  ParticleContactArbitraryFromOverlapVolume(
      const ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc);
  /// Move Constructor
  ParticleContactArbitraryFromOverlapVolume(
      ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&& pc);

  /// Read access to particle positions
  constexpr const std::array<PhysR<T, D>, 2>& getParticlePositions() const;
  /// Return particle position
  constexpr const PhysR<T, D>& getParticlePosition(const std::size_t& id) const;
  /// Read access to min
  constexpr const PhysR<T, D>& getMin() const;
  /// Read access to max
  constexpr const PhysR<T, D>& getMax() const;
  /// Read access to particle IDs
  constexpr const std::array<std::size_t, 2>& getIDs() const;
  /// Read access to damping factor
  constexpr const T getDampingFactor() const;
  /// Read access to the responsible rank
  constexpr const int& getResponsibleRank() const;

  /// Set particle positions
  constexpr void
  setParticlePositions(const std::array<PhysR<T, D>, 2>& positions);
  /// Set position of specific particle
  constexpr void setParticlePosition(const PhysR<T, D>& position,
                                     const std::size_t& id);
  /// Set damping factor for contact
  constexpr void setDampingFactor(const T dampingFactor);
  /// Set damping factor from the magnitude of the initial relative impact velocity in direction of contact normal and the wanted coefficient of restitution
  constexpr void
  setDampingFactorFromInitialVelocity(const T coefficientOfRestitution,
                                      const T initialRelativeVelocityMagnitude);
  /// Set processor that is responsible for contact treatment
  constexpr void setResponsibleRank(const int& rank);

  /// Reset min and max to default values
  constexpr void resetMinMax();
  /// Update min and max with given position inside the contact
  constexpr void updateMinMax(const PhysR<T, D>& positionInsideTheContact);
  /// Increase bounding box size
  constexpr void increaseMinMax(const Vector<T, D>& increaseBy);

  /// Combining two contacts, if the particle ids are the same
  constexpr void
  combineWith(ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc);
  /// Returns if contact holds data
  constexpr bool isEmpty() const;
  /// Returns if the contact is a new contact
  constexpr bool isNew() const;
  /// Sets 'isNew'
  constexpr void isNew(const bool newContact);
  /// Returns if the particle position is up-to-date
  constexpr bool isParticlePositionUpdated() const;
  // Set if the particle position is up-to-date
  constexpr void setParticlePositionUpdated(bool updated);
  /// Copy assignment
  ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&
  operator=(const ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>& pc);
  /// Move assignment
  ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&
  operator=(ParticleContactArbitraryFromOverlapVolume<T, D, CONVEX>&& pc);

  /// Serialize contact data
  std::size_t serialize(std::uint8_t* buffer);
  /// Deserialize contact data and save in object
  std::size_t deserialize(std::uint8_t* buffer);
  /// Get serial size
  constexpr static std::size_t getSerialSize() { return serialSize; };

  /// Print relevant quantities
  void print();
};

} // namespace contact
} // namespace particles
} // namespace olb
#endif
