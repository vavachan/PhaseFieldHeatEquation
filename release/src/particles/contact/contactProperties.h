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

#ifndef CONTACT_PROPERTIES_H
#define CONTACT_PROPERTIES_H

namespace olb {
namespace particles {
namespace contact {

template <typename T>
struct ContactProperty {
  constexpr ContactProperty() = default;

  constexpr ContactProperty(const T _effectiveYoungsModulus,
                            const T _coefficientOfRestitution,
                            const T _coefficientKineticFriction,
                            const T _coefficientStaticFriction,
                            const T _staticKineticTransitionVelocity)
      : effectiveYoungsModulus(_effectiveYoungsModulus)
      , coefficientOfRestitution(_coefficientOfRestitution)
      , coefficientOfKineticFriction(_coefficientKineticFriction)
      , coefficientOfStaticFriction(_coefficientStaticFriction)
      , staticKineticTransitionVelocity(_staticKineticTransitionVelocity) {};

  T effectiveYoungsModulus {};
  T coefficientOfRestitution {};
  T coefficientOfKineticFriction {};
  T coefficientOfStaticFriction {};
  T staticKineticTransitionVelocity {};
};

/**
 * Object that stores properties which are necessary for the computation of contact forces
 * `N` = number of different materials
 * The `material` here is an identifier of a solid material with certain (mechanical) properties
 * This `material` is something completely different from the lattice's material number, which is used to assign boundary conditions and dynamics.
 */
template <typename T, unsigned N, bool ENABLE_RANGE_CHECK = false>
struct ContactProperties {
private:
  std::array<ContactProperty<T>, N*(N + 1) / 2> data;

  constexpr inline unsigned getIndex(unsigned materialA,
                                     unsigned materialB) const;

public:
  /// Constructor
  constexpr ContactProperties() = default;

  /// Set contact properties
  constexpr void set(const unsigned materialA, const unsigned materialB,
                     const T effectiveYoungsModulus,
                     const T coefficientOfRestitution,
                     const T coefficientKineticFriction,
                     const T coefficientStaticFriction,
                     const T staticKineticTransitionVelocity = T {0.01});

  /// Get effective modulus of elasticity
  constexpr T getEffectiveYoungsModulus(const unsigned materialA,
                                        const unsigned materialB) const;
  /// Get damping constant
  constexpr T getCoefficientOfRestitution(const unsigned materialA,
                                          const unsigned materialB) const;
  /// Get coefficient of kinetic friction
  constexpr T getKineticFrictionCoefficient(const unsigned materialA,
                                            const unsigned materialB) const;
  /// Get coefficient of static friction
  constexpr T getStaticFrictionCoefficient(const unsigned materialA,
                                           const unsigned materialB) const;
  /// Get transition velocity (static to kinetic)
  constexpr T
  getStaticKineticTransitionVelocity(const unsigned materialA,
                                     const unsigned materialB) const;
};

template <typename T, unsigned N, typename F>
constexpr ContactProperties<T, N> createContactProperties(F f)
{
  ContactProperties<T, N> contactProperties {};
  f(contactProperties);
  return contactProperties;
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
