/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

/* This file contains functions used for the calculation of particle related dynamics.
 *
*/

#ifndef SMOOTHINDICATOR_INTERACTION_H
#define SMOOTHINDICATOR_INTERACTION_H

#include "core/blockStructure.h"
#include "functors/analytical/indicator/indicatorBase.h"
#include "functors/primitive/sdf.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {

namespace particles {

namespace resolved {

template <typename T, typename S, typename PARTICLETYPE>
PhysR<S, PARTICLETYPE::d> transformInput(Particle<T, PARTICLETYPE>& particle,
                                         const PhysR<S, PARTICLETYPE::d>& input)
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;

  const PhysR<T, D> position = particle.template getField<GENERAL, POSITION>();
  if constexpr (PARTICLETYPE::template providesNested<SURFACE, ROT_MATRIX>()) {
    const Vector<T, utilities::dimensions::convert<D>::matrix> rotationMatrix =
        util::invertRotationMatrix<T, D>(
            particle.template getField<SURFACE, ROT_MATRIX>());
    return util::executeRotation<T, D, true>(input, rotationMatrix, position);
  }
  else {
    return input - position;
  }

  __builtin_unreachable();
}

template <typename T, typename S, typename PARTICLETYPE>
PhysR<S, PARTICLETYPE::d>
transformDirection(Particle<T, PARTICLETYPE>&       particle,
                   const PhysR<S, PARTICLETYPE::d>& direction)
{
  constexpr unsigned D = PARTICLETYPE::d;
  using namespace descriptors;

  if constexpr (PARTICLETYPE::template providesNested<SURFACE, ROT_MATRIX>()) {
    const Vector<T, utilities::dimensions::convert<D>::matrix> rotationMatrix =
        util::invertRotationMatrix<T, D>(
            particle.template getField<SURFACE, ROT_MATRIX>());
    return util::executeRotation<T, D>(direction, rotationMatrix);
  }
  else {
    return direction;
  }

  __builtin_unreachable();
}

template <typename T, typename S, typename PARTICLETYPE>
bool isInsideCircumRadius(Particle<T, PARTICLETYPE>&       particle,
                          const PhysR<S, PARTICLETYPE::d>& input)
{
  using namespace descriptors;
  auto sIndicator = access::getSmoothIndicatorPtr(particle);
  return norm(access::getPosition(particle) - input) <=
         sIndicator->getCircumRadius();
}

template <typename T, typename S, typename PARTICLETYPE>
bool isInsideParticle(Particle<T, PARTICLETYPE>&       particle,
                      const PhysR<S, PARTICLETYPE::d>& input)
{
  using namespace descriptors;
  auto sIndicator = access::getSmoothIndicatorPtr(particle);

  if (!isInsideCircumRadius(particle, input)) {
    return false;
  }

  PhysR<S, PARTICLETYPE::d> newInput = transformInput(particle, input);
  return sIndicator->signedDistance(newInput) <= 0;
}

template <typename T, typename S, typename PARTICLETYPE>
const S signedDistanceToParticle(Particle<T, PARTICLETYPE>&       particle,
                                 const PhysR<S, PARTICLETYPE::d>& input)
{
  using namespace descriptors;
  auto sIndicator = particle.template getField<SURFACE, SINDICATOR>();
  PhysR<S, PARTICLETYPE::d> newInput = transformInput(particle, input);

  //Check whether elongation provided
  if constexpr (access::providesElongation<PARTICLETYPE>()) {
    const unsigned D = PARTICLETYPE::d;
    //Retrieve original sdf
    std::function<S(const Vector<S, D>&)> sdf =
        [&sIndicator](const Vector<S, D>& input) {
          return sIndicator->signedDistance(input);
        };
    //Retrieve elongation
    auto elongation = access::getElongation(particle);
    //Return elongated whapper
    constexpr bool symmetryCheck = false; //User responsibility here.
    return sdf::elongation<T, symmetryCheck>(sdf, newInput, elongation);
  }
  else {
    return sIndicator->signedDistance(newInput);
  }
}

template <typename T, typename S, typename PARTICLETYPE>
const S distanceToParticle(Particle<T, PARTICLETYPE>&        particle,
                           const Vector<S, PARTICLETYPE::d>& origin,
                           const Vector<S, PARTICLETYPE::d>& direction,
                           S precision, S pitch)
{
  using namespace descriptors;
  auto sIndicator = particle.template getField<SURFACE, SINDICATOR>();
  PhysR<S, PARTICLETYPE::d> input = transformInput(particle, origin);
  PhysR<S, PARTICLETYPE::d> dir   = transformDirection(particle, direction);

  T distance(T(0.));
  if (util::distance(
          distance, input, dir, precision, pitch,
          [&](bool output[1], const T input[3]) {
            output[0] = sIndicator->signedDistance(
                            PhysR<T, PARTICLETYPE::d>(input)) <= 0;
            return output[0];
          },
          [&](const Vector<T, 3>& pos) {
            // bigger bounding box is usually necessary
            return norm(pos) <= 1.5 * sIndicator->getCircumRadius();
          })) {
    return distance;
  }
  else {
#ifdef OLB_DEBUG
    OstreamManager clout(std::cout, "distanceCalculation");
    clout << "WARNING: Distance to particle in direction "
          << normalize(direction) << " couldn't be calculated." << std::endl;
#endif
    return T(0.);
  }

  __builtin_unreachable();
}

template <typename T, typename S, typename PARTICLETYPE>
const S distanceToParticle(Particle<T, PARTICLETYPE>&        particle,
                           const Vector<S, PARTICLETYPE::d>& origin,
                           const Vector<S, PARTICLETYPE::d>& direction,
                           S                                 precision)
{
  using namespace descriptors;
  auto sIndicator = particle.template getField<SURFACE, SINDICATOR>();
  PhysR<S, PARTICLETYPE::d> input = transformInput(particle, origin);
  PhysR<S, PARTICLETYPE::d> dir   = transformDirection(particle, direction);

  T distance(T(0.));
  if (util::distance(
          distance, input, dir, precision,
          [&](const Vector<T, PARTICLETYPE::d>& input) {
            return sIndicator->signedDistance(input);
          },
          [&](const Vector<T, PARTICLETYPE::d>& pos) {
            // bigger bounding box is usually necessary
            return norm(pos) <= 1.5 * sIndicator->getCircumRadius();
          })) {
    return distance;
  }
  else {
#ifdef OLB_DEBUG
    OstreamManager clout(std::cout, "distanceCalculation");
    clout << "WARNING: Distance to particle in direction "
          << normalize(direction) << " couldn't be calculated." << std::endl;
#endif
    return T(0.);
  }

  __builtin_unreachable();
}

template <typename T, typename S, typename PARTICLETYPE>
Vector<S, PARTICLETYPE::d>
normalOnParticleSurface(Particle<T, PARTICLETYPE>&        particle,
                        const Vector<S, PARTICLETYPE::d>& pos, const T meshSize)
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  auto sIndicator      = particle.template getField<SURFACE, SINDICATOR>();

  // the following gives us the option to overload the calculation of the surface normal for different indicators
  return sIndicator->surfaceNormal(pos, meshSize, [&](const Vector<S, D>& pos) {
    return transformInput(particle, pos);
  });
}

template <typename T, typename S, typename PARTICLETYPE>
bool evalSolidVolumeFraction(T output[], const S input[],
                             Particle<T, PARTICLETYPE>& particle)
{
  using namespace descriptors;
  auto    sIndicator = particle.template getField<SURFACE, SINDICATOR>();
  T const signedDist =
      signedDistanceToParticle(particle, PhysR<S, PARTICLETYPE::d>(input));
  return sdf::evalSolidVolumeFraction(output, signedDist,
                                      sIndicator->getEpsilon());
}

} //namespace resolved

} //namespace particles

} //namespace olb

#endif
