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
 * This file contains functions used for the particle-wall and particle-particle contact force calculation.
*/

#ifndef CONTACT_FUNCTIONS_H
#define CONTACT_FUNCTIONS_H

#include <numbers>

namespace olb {

namespace particles {

namespace contact {

template <typename T>
constexpr T evalEffectiveYoungModulus(T E1, T E2, T nu1, T nu2)
{
  const T denominator = (1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2;
  return T {1} / denominator;
}

/// Calculates the damping factor according to Carvalho & Martins (2019) (10.1016/j.mechmachtheory.2019.03.028)
template <typename T>
constexpr T evalDampingFactor(const T coefficientOfRestitution,
                              const T initialRelativeVelocityMagnitude)
{
  // This case should never happen, but could for example due to initial conditions
  if (initialRelativeVelocityMagnitude <= 0) {
    return 0;
  }
  return 1.5 * (1 - coefficientOfRestitution) *
         (11 - coefficientOfRestitution) / (1 + 9 * coefficientOfRestitution) /
         initialRelativeVelocityMagnitude;
}

template <typename T, typename PARTICLETYPE>
T evalContactDetectionDistance(Particle<T, PARTICLETYPE>& particle,
                               T const                    physDeltaX)
{
  constexpr unsigned D = PARTICLETYPE::d;

  if constexpr (particles::access::providesContactMaterial<PARTICLETYPE>()) {
    constexpr T factor = []() {
      static_assert(D == 2 || D == 3, "Only D=2 and D=3 are supported");
      return T{0.5} * (D == 3 ? std::numbers::sqrt3_v<T> : std::numbers::sqrt2_v<T>);
    }();
    return factor * physDeltaX +
           particles::access::getEnlargementForContact(particle);
  }
  else {
    return T {0};
  }

  __builtin_unreachable();
}

template <typename T>
T evalCircumRadius(T const contactDetectionDistance, T const circumRadius,
                   T const epsilon)
{
  return circumRadius - epsilon + util::max(epsilon, contactDetectionDistance);
}

template <typename T, typename PARTICLETYPE>
T evalCircumRadius(Particle<T, PARTICLETYPE>& particle, T const physDeltaX,
                   T const circumRadius, T const epsilon)
{
  const T contactDetectionDistance =
      evalContactDetectionDistance(particle, physDeltaX);
  return evalCircumRadius(contactDetectionDistance, circumRadius, epsilon);
}

template <typename T, typename PARTICLETYPE>
T evalCircumRadius(Particle<T, PARTICLETYPE>& particle, T const physDeltaX)
{
  using namespace descriptors;
  auto sIndicator = particle.template getField<SURFACE, SINDICATOR>();
  return evalCircumRadius(particle, physDeltaX, sIndicator->getCircumRadius(),
                          sIndicator->getEpsilon());
}

template <typename CONTACTTYPE>
void cleanContacts(std::vector<CONTACTTYPE>& contacts)
{
  contacts.erase(std::remove_if(contacts.begin(), contacts.end(),
                                [](const CONTACTTYPE& contact) {
                                  return contact.isEmpty();
                                }),
                 contacts.end());
}

template <typename CONTACTTYPE>
void resetResponsibleRank(CONTACTTYPE& contact)
{
  contact.setResponsibleRank(std::numeric_limits<int>::max());
}

template <typename CONTACTTYPE>
bool isResponsibleRankSet(CONTACTTYPE& contact)
{
  return contact.getResponsibleRank() < std::numeric_limits<int>::max();
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
