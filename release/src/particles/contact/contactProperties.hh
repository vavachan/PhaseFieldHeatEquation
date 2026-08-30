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

#ifndef CONTACT_PROPERTIES_HH
#define CONTACT_PROPERTIES_HH

#include "contactProperties.h"
#include "core/olbDebug.h"
#include "utilities/oalgorithm.h"

namespace olb {
namespace particles {
namespace contact {

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr inline unsigned
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getIndex(unsigned materialA,
                                                      unsigned materialB) const
{
  const unsigned minMat {util::min(materialA, materialB)};
  const unsigned maxMat {util::max(materialA, materialB)};
  const unsigned index {(maxMat + 1) * maxMat / 2 + minMat};

  // Optional range check - for debugging
  if constexpr (ENABLE_RANGE_CHECK) {
    OLB_ASSERT(index >= 0 && index <= N,
               "ERROR: ContactProperties material is not defined");
  }

  return index;
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr void ContactProperties<T, N, ENABLE_RANGE_CHECK>::set(
    const unsigned materialA, const unsigned materialB,
    const T effectiveYoungsModulus, const T coefficientOfRestitution,
    const T coefficientKineticFriction, const T coefficientStaticFriction,
    const T staticKineticTransitionVelocity)
{
  if (materialA >= N || materialB >= N) {
    std::cerr << "WARNING: At least one of the provided materials exceeds the "
                 "total number of materials."
              << std::endl;
  }

  data[getIndex(materialA, materialB)] =
      ContactProperty(effectiveYoungsModulus, coefficientOfRestitution,
                      coefficientKineticFriction, coefficientStaticFriction,
                      staticKineticTransitionVelocity);
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr T
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getEffectiveYoungsModulus(
    const unsigned materialA, const unsigned materialB) const
{
  return data[getIndex(materialA, materialB)].effectiveYoungsModulus;
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr T
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getCoefficientOfRestitution(
    const unsigned materialA, const unsigned materialB) const
{
  return data[getIndex(materialA, materialB)].coefficientOfRestitution;
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr T
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getKineticFrictionCoefficient(
    const unsigned materialA, const unsigned materialB) const
{
  return data[getIndex(materialA, materialB)].coefficientOfKineticFriction;
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr T
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getStaticFrictionCoefficient(
    const unsigned materialA, const unsigned materialB) const
{
  return data[getIndex(materialA, materialB)].coefficientOfStaticFriction;
}

template <typename T, unsigned N, bool ENABLE_RANGE_CHECK>
constexpr T
ContactProperties<T, N, ENABLE_RANGE_CHECK>::getStaticKineticTransitionVelocity(
    const unsigned materialA, const unsigned materialB) const
{
  return data[getIndex(materialA, materialB)].staticKineticTransitionVelocity;
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
