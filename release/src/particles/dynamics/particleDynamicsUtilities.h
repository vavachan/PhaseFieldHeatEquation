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

#ifndef PARTICLE_DYNAMICS_UTILITIES_H
#define PARTICLE_DYNAMICS_UTILITIES_H

#include <cassert>

namespace olb {

namespace particles {

//Forward declaration of particle
template <typename T, typename DESCRIPTOR> class Particle;

namespace dynamics {

//TODO: should be moved to more accessible location
template<typename T, typename PARTICLETYPE>
void updateRotationMatrix( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  auto rotationMatrix = util::calculateRotationMatrix<T,PARTICLETYPE::d>(
                          particle.template getField<SURFACE,ANGLE>() );
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>() ) {
    particle.template setField<SURFACE,ROT_MATRIX>( rotationMatrix );
  } else {
    std::cerr << "ERROR: The field ROT_MATRIX must be provided." << std::endl;
  }
}


/// Helper functions

//TODO: Check whether any are actually used anymore after particle dynamics adaption

//TODO: Check, whether obsolet with new doAtParticleWallContact
template<typename T, typename PARTICLETYPE, typename F>
void doAtCubicBoundPenetration(
  Particle<T,PARTICLETYPE>& particle,
  Vector<T,PARTICLETYPE::d> domainMin,
  Vector<T,PARTICLETYPE::d> domainMax,
  F boundTreatment )
{
  using namespace particles::access;
  auto radius = getRadius<T,PARTICLETYPE>( particle );
  auto position = getPosition<T,PARTICLETYPE>( particle );
  for (int iDim=0; iDim<PARTICLETYPE::d; ++iDim) {
    T distToLowerBound = (position[iDim] - radius) - domainMin[iDim];
    T distToUpperBound = domainMax[iDim] - (position[iDim] + radius);
    Vector<T,PARTICLETYPE::d> normal(0.);
    if ( distToLowerBound <= 0. ) {
      normal[iDim] = 1.;
      boundTreatment(iDim, normal, distToLowerBound);
    }
    if ( distToUpperBound <= 0. ) {
      normal[iDim] = -1.;
      boundTreatment(iDim, normal, distToUpperBound);
    }
  }
}

template<typename T, typename PARTICLETYPE>
void resetDirection( Particle<T,PARTICLETYPE>& particle,
                     Vector<T,PARTICLETYPE::d> positionPre, int iDir )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  position[iDir] = positionPre[iDir];
  velocity[iDir] = 0.;
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
}

template<typename T, typename PARTICLETYPE>
void resetMovement( Particle<T,PARTICLETYPE>& particle,
                    Vector<T,PARTICLETYPE::d> positionPre,
                    Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> anglePre
                      = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.) )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  particle.template setField<GENERAL,POSITION>( positionPre );
  particle.template setField<MOBILITY,VELOCITY>( Vector<T,D>(0.) );
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ANGLE>() ) {
    particle.template setField<SURFACE,ANGLE>( anglePre );
    particle.template setField<MOBILITY,ANG_VELOCITY>( Vector<T,Drot>(0.) );
  }
}

//TODO: remove in favour of ParametersD
template<typename T, typename PARTICLETYPE>
struct ParticleDynamicStateNoAngle{
  Vector<T,PARTICLETYPE::d> position;
  static const bool hasAngle = false;

  ParticleDynamicStateNoAngle()
  {}

  ParticleDynamicStateNoAngle(
    Vector<T,PARTICLETYPE::d> position )
    : position(position)
  {}
};

template<typename T, typename PARTICLETYPE>
struct ParticleDynamicStateAngle{
  Vector<T,PARTICLETYPE::d> position;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle;
  static const bool hasAngle = true;

  ParticleDynamicStateAngle()
  {}

  ParticleDynamicStateAngle(
    Vector<T,PARTICLETYPE::d> position )
    : position(position)
  {}

  ParticleDynamicStateAngle(
    Vector<T,PARTICLETYPE::d> position,
    Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle )
    : position(position), angle(angle)
  {}
};

template <typename T, typename PARTICLETYPE>
using DynState = std::conditional_t<
  PARTICLETYPE::template providesNested<descriptors::SURFACE,descriptors::ANGLE>(),
  ParticleDynamicStateAngle<T,PARTICLETYPE>,
  ParticleDynamicStateNoAngle<T,PARTICLETYPE>
>;

} //namespace dynamics

} //namespace particles

} //namespace olb

#endif
