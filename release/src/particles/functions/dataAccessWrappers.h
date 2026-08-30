/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Jan E. Marquardt, Martin Sadric, Mathias J. Krause
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


/* Wrappers for simplified access of particle data.
 * - features automatic differentiation for different particle types (e.g. radius)
 * - includes asserts for respective field access
 *
*/

#ifndef PARTICLE_DATA_ACCESS_WRAPPERS_H
#define PARTICLE_DATA_ACCESS_WRAPPERS_H

#include <cassert>

namespace olb {

namespace particles {

namespace access {

//Calculate angular acceleration when considering 3 dimensions
//TODO: remove to separate location
template<typename T>
Vector<T,3> calcAngAcceleration3D( Vector<T,3>& torque, Vector<T,3>& momentOfInertia,
                                   const Matrix<T, 3, 3>& rotationMatrix )
{
  T data[3][3] = {
    {momentOfInertia[0],    T {0},    T {0}},
    {   T {0}, momentOfInertia[1],    T {0}},
    {   T {0},    T {0}, momentOfInertia[2]}
  };
  Matrix<T, 3, 3>       inertiaTensor(data);
  inertiaTensor = (rotationMatrix * inertiaTensor) * (rotationMatrix.transpose());

  /*
   * In the following, we use the solution of a linear equation system
   * torque = inertia tensor * angular acceleration + angular velocity x (inertia tensor * angular velocity)
   * Here, we want to know the components of the angular acceleration
   * In the following, left = torque - angular velocity x (inertia tensor * angular velocity)
   */

  // const Vector<T,3> angularVelocity = getAngularVelocity(particle);
  // Using Euler's equation for rigid body dynamics
  // Currently not, as it leads to inaccurate results of some tests
  // Add the cross product to use Euler's equation for rigid body dynamics
  const Vector<T,3> left = torque;//- crossProduct3D(angularVelocity, (inertiaTensor * angularVelocity));

  const T denominator = util::pow(inertiaTensor[0][2], 2) * inertiaTensor[1][1]
    + inertiaTensor[0][0] * util::pow(inertiaTensor[1][2], 2)
    - 2 * inertiaTensor[0][1] * inertiaTensor[0][2] * inertiaTensor[1][2]
    + util::pow(inertiaTensor[0][1], 2) * inertiaTensor[2][2]
    - inertiaTensor[0][0] * inertiaTensor[1][1] * inertiaTensor[2][2];
  const T factor = T{1} / denominator;

  Vector<T,3> angularAcceleration;
  angularAcceleration[0] = (util::pow(inertiaTensor[1][2], 2) - inertiaTensor[1][1] * inertiaTensor[2][2]) * left[0]
    + (- inertiaTensor[0][2] * inertiaTensor[1][2] + inertiaTensor[0][1] * inertiaTensor[2][2]) * left[1]
    + (inertiaTensor[0][2] * inertiaTensor[1][1] - inertiaTensor[0][1] * inertiaTensor[1][2]) * left[2];
  angularAcceleration[1] = (- inertiaTensor[0][2] * inertiaTensor[1][2] + inertiaTensor[0][1] * inertiaTensor[2][2]) * left[0]
    + (util::pow(inertiaTensor[0][2], 2) - inertiaTensor[0][0] * inertiaTensor[2][2]) * left[1]
    + (- inertiaTensor[0][1] * inertiaTensor[0][2] + inertiaTensor[0][0] * inertiaTensor[1][2]) * left[2];
  angularAcceleration[2] = (inertiaTensor[0][2] * inertiaTensor[1][1] - inertiaTensor[0][1] * inertiaTensor[1][2]) * left[0]
    + (- inertiaTensor[0][1] * inertiaTensor[0][2] + inertiaTensor[0][0] * inertiaTensor[1][2]) * left[1]
    + (util::pow(inertiaTensor[0][1], 2) - inertiaTensor[0][0] * inertiaTensor[1][1]) * left[2];

  for(unsigned iDim = 0; iDim < 3; ++iDim) {
    angularAcceleration[iDim] *= factor;
  }

  return angularAcceleration;
}

//// Checks for specific particle fields and groups

/// Provides field ID
template<typename PARTICLETYPE>
constexpr bool providesID()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PARALLELIZATION,ID>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesID(Particle<T,PARTICLETYPE>& particle){
  return providesID<PARTICLETYPE>();
}

/// Provides field ID
template<typename PARTICLETYPE>
constexpr bool providesInvalid()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<GENERAL,INVALID>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesInvalid(Particle<T,PARTICLETYPE>& particle){
  return providesInvalid<PARTICLETYPE>();
}

/// Provides field POSITION
template<typename PARTICLETYPE>
constexpr bool providesPosition()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<GENERAL,POSITION>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesPosition(Particle<T,PARTICLETYPE>& particle){
  return providesPosition<PARTICLETYPE>();
}

/// Provides field RADIUS
template<typename PARTICLETYPE>
constexpr bool providesRadius()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PHYSPROPERTIES,RADIUS>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesRadius(Particle<T,PARTICLETYPE>& particle){
  return providesRadius<PARTICLETYPE>();
}

/// Provides field DENSITY
template<typename PARTICLETYPE>
constexpr bool providesDensity()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PHYSPROPERTIES,DENSITY>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesDensity(Particle<T,PARTICLETYPE>& particle){
  return providesDensity<PARTICLETYPE>();
}

/// Provides field MASS
template<typename PARTICLETYPE>
constexpr bool providesMass()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PHYSPROPERTIES,MASS>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesMass(Particle<T,PARTICLETYPE>& particle){
  return providesMass<PARTICLETYPE>();
}

/// Provides field MASS or DENSITY
template<typename PARTICLETYPE>
constexpr bool providesMassOrDensity()
{
  using namespace descriptors;
  return providesMass<PARTICLETYPE>() || providesDensity<PARTICLETYPE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesMassOrDensity(Particle<T,PARTICLETYPE>& particle){
  return providesMassOrDensity<PARTICLETYPE>();
}

/// Provides field ANGLE
template<typename PARTICLETYPE>
constexpr bool providesAngle()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,ANGLE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesAngle(Particle<T,PARTICLETYPE>& particle){
  return providesAngle<PARTICLETYPE>();
}

/// Provides field VELOCITY
template<typename PARTICLETYPE>
constexpr bool providesVelocity()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesVelocity(Particle<T,PARTICLETYPE>& particle){
  return providesVelocity<PARTICLETYPE>();
}

/// Provides field VELOCITY
template<typename PARTICLETYPE>
constexpr bool providesAngVelocity()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<MOBILITY,ANG_VELOCITY>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesAngVelocity(Particle<T,PARTICLETYPE>& particle){
  return providesAngVelocity<PARTICLETYPE>();
}

/// Provides field FORCE
template<typename PARTICLETYPE>
constexpr bool providesForce()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<FORCING,FORCE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesForce(Particle<T,PARTICLETYPE>& particle){
  return providesForce<PARTICLETYPE>();
}

/// Provides field TORQUE
template<typename PARTICLETYPE>
constexpr bool providesTorque()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<FORCING,TORQUE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesTorque(Particle<T,PARTICLETYPE>& particle){
  return providesTorque<PARTICLETYPE>();
}

/// Provides field ADHESION
template<typename PARTICLETYPE>
constexpr bool providesAdhesion()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<FORCING,ADHESION>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesAdhesion(Particle<T,PARTICLETYPE>& particle){
  return providesAdhesion<PARTICLETYPE>();
}

/// Provides field MOFI
template<typename PARTICLETYPE>
constexpr bool providesMomentOfInertia()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PHYSPROPERTIES,MOFI>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesMomentOfInertia(Particle<T,PARTICLETYPE>& particle){
  return providesMomentOfInertia<PARTICLETYPE>();
}

/// Provides field INVALID
template<typename PARTICLETYPE>
constexpr bool providesValid()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<GENERAL,INVALID>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesValid(Particle<T,PARTICLETYPE>& particle){
  return providesValid<PARTICLETYPE>();
}

/// Provides field DYNAMICS_ID
template<typename PARTICLETYPE>
constexpr bool providesDynamicsID()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DYNAMICS_ID>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesDynamicsID(Particle<T,PARTICLETYPE>& particle){
  return providesDynamicsID<PARTICLETYPE>();
}

/// Provides field ACTIVE
template<typename PARTICLETYPE>
constexpr bool providesActive()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<DYNBEHAVIOUR,ACTIVE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesActive(Particle<T,PARTICLETYPE>& particle){
  return providesActive<PARTICLETYPE>();
}

/// Provides field COMPUTE_MOTION
template<typename PARTICLETYPE>
constexpr bool providesComputeMotion()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<DYNBEHAVIOUR,COMPUTE_MOTION>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesComputeMotion(Particle<T,PARTICLETYPE>& particle){
  return providesComputeMotion<PARTICLETYPE>();
}

/// Provides field COMPUTE_CONTACT
template<typename PARTICLETYPE>
constexpr bool providesComputeContact()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<DYNBEHAVIOUR,COMPUTE_CONTACT>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesComputeContact(Particle<T,PARTICLETYPE>& particle){
  return providesComputeContact<PARTICLETYPE>();
}

/// Provides field ROT_MATRIX
template<typename PARTICLETYPE>
constexpr bool providesRotationMatrix()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesRotationMatrix(Particle<T,PARTICLETYPE>& particle){
  return providesRotationMatrix<PARTICLETYPE>();
}

/// Provides group SURFACE
template<typename PARTICLETYPE>
constexpr bool providesSurface()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesSurface(Particle<T,PARTICLETYPE>& particle){
  return providesSurface<PARTICLETYPE>();
}

// Provides field SINDICATOR
template<typename PARTICLETYPE>
constexpr bool providesSmoothIndicator()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,SINDICATOR>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesSmoothIndicator(Particle<T,PARTICLETYPE>& particle)
{
  return providesSmoothIndicator<PARTICLETYPE>();
}

/// Provides group PARALLELIZATION
template<typename PARTICLETYPE>
constexpr bool providesParallelization()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PARALLELIZATION>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesParallelization(Particle<T,PARTICLETYPE>& particle){
  return providesParallelization<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesSpecies()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<PHYSPROPERTIES,SPECIES>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesSpecies(Particle<T,PARTICLETYPE>& particle){
  return providesSpecies<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesDetaching()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DETACHING>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesDetaching(Particle<T,PARTICLETYPE>& particle){
  return providesDetaching<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesCORoffset()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesCORoffset(Particle<T,PARTICLETYPE>& particle){
  return providesCORoffset<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesElongation()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,ELONGATION>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesElongation(Particle<T,PARTICLETYPE>& particle){
  return providesElongation<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesAccelerationStrd()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<MOBILITY,ACCELERATION_STRD>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesAccelerationStrd(Particle<T,PARTICLETYPE>& particle){
  return providesAccelerationStrd<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesAngAccelerationStrd()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<MOBILITY,ANG_ACC_STRD>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesAngAccelerationStrd(Particle<T,PARTICLETYPE>& particle){
  return providesAngAccelerationStrd<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesEnlargementForContactTreatment()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<NUMERICPROPERTIES,ENLARGEMENT_FOR_CONTACT>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesEnlargementForContactTreatment(Particle<T,PARTICLETYPE>& particle){
  return providesEnlargementForContactTreatment<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesIsInContact()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<NUMERICPROPERTIES,IS_IN_CONTACT>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesIsInContact(Particle<T,PARTICLETYPE>& particle){
  return providesIsInContact<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool providesContactMaterial()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<MECHPROPERTIES,MATERIAL>();
}
template<typename T, typename PARTICLETYPE>
constexpr bool providesContactMaterial(Particle<T,PARTICLETYPE>& particle){
  return providesContactMaterial<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool is2D()
{
  return (PARTICLETYPE::d==2);
}
template<typename T, typename PARTICLETYPE>
constexpr bool is2D(Particle<T,PARTICLETYPE>& particle){
  return is2D<PARTICLETYPE>();
}

template<typename PARTICLETYPE>
constexpr bool is3D()
{
  return (PARTICLETYPE::d==3);
}
template<typename T, typename PARTICLETYPE>
constexpr bool is3D(Particle<T,PARTICLETYPE>& particle){
  return is3D<PARTICLETYPE>();
}

//// Getter for specific particle fields

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getPosition( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesPosition<PARTICLETYPE>(),
      "Field GENERAL:POSITION has to be provided");
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  return position;
}

template<bool ensureAngularBounds=false, typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngle( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  static_assert(providesAngle<PARTICLETYPE>(), "Field SURFACE:ANGLE has to be provided");
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() );
  if constexpr (ensureAngularBounds){
    for (unsigned iRot=0; iRot<Drot; ++iRot) {
      angle[iRot] = util::fmod( angle[iRot], 2.*M_PI );
    }
  }
  return angle;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::matrix> getRotationMatrix( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned DrotMat = utilities::dimensions::convert<PARTICLETYPE::d>::matrix;
  static_assert(providesRotationMatrix<PARTICLETYPE>(), "Field SURFACE:ROT_MATRIX has to be provided");
  Vector<T,DrotMat> rotationMatrix( particle.template getField<SURFACE,ROT_MATRIX>() );
  return rotationMatrix;
}

template<unsigned dir=2, typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getSurfaceNormal( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  auto rotationMatrix = getRotationMatrix( particle );
  //Note: Convention here, positively pointing towards dir (default z-dir)
  Vector<T,D> normal;
  for (unsigned iDim=0; iDim<D; ++iDim) {
    unsigned iMat = D*iDim+dir;
    normal[iDim] = rotationMatrix[iMat];
  }
  return normal;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getVelocity( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesVelocity<PARTICLETYPE>(), "Field MOBILITY:VELOCITY has to be provided");
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  return velocity;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngularVelocity( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  static_assert(providesAngVelocity<PARTICLETYPE>(), "Field MOBILITY:ANG_VELOCITY has to be provided");
  Vector<T,Drot> angVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() );
  return angVelocity;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getForce( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesForce<PARTICLETYPE>(), "Field FORCING:FORCE has to be provided");
  Vector<T,D> force( particle.template getField<FORCING,FORCE>() );
  return force;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getTorque( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  static_assert(providesTorque<PARTICLETYPE>(), "Field FORCING:TORQUE has to be provided");
  Vector<T,Drot> torque( particle.template getField<FORCING,TORQUE>() );
  return torque;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getMomentOfInertia( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  static_assert(providesMomentOfInertia<PARTICLETYPE>(), "Field PHYSPROPERTIES:MOFI has to be provided");
  Vector<T,Drot> mofi( particle.template getField<PHYSPROPERTIES,MOFI>() );
  return mofi;
}

// Consisting of normal and tangential adhesion component of particle.
template<typename T, typename PARTICLETYPE>
Vector<T,2> getAdhesion( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;

  static_assert(providesAdhesion<PARTICLETYPE>(), "Field FORCING:ADHESION has to be provided");
  Vector<T,2> adhesion( particle.template getField<FORCING,ADHESION>() );
  return adhesion;
}

template<typename T, typename PARTICLETYPE>
bool isValid( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  bool valid = true;
  if constexpr(providesValid<PARTICLETYPE>()) {
    valid = !particle.template getField<GENERAL,INVALID>();
  }
  return valid;
}

template<typename T, typename PARTICLETYPE>
bool isActive( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  bool active = true;
  if constexpr(providesActive<PARTICLETYPE>()) {
    active = particle.template getField<DYNBEHAVIOUR,ACTIVE>();
  }
  return active;
}

//Get smooth indicator pointer
template<typename T, typename PARTICLETYPE>
auto getSmoothIndicatorPtr( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  static_assert(providesSmoothIndicator<PARTICLETYPE>(), "Field FORCING:SINDICATOR has to be provided");
  auto sIndicatorPtr = particle.template getField<SURFACE,SINDICATOR>();
  return sIndicatorPtr;
}

//Get radius
template<typename T, typename PARTICLETYPE>
T getRadius( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  T radius;
  if constexpr (providesSmoothIndicator<PARTICLETYPE>() ) {
    radius = getSmoothIndicatorPtr(particle)->getCircumRadius();
  }
  else if constexpr ( providesRadius<PARTICLETYPE>() ) {
    radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
  }
  else {
    std::cerr << "ERROR: no Field found providing radius!" << std::endl;
  }
  return radius;
}

/// Returns the volume of a particle (for subgrid the volume of a sphere/circle is multiplied by the shapeFactor),
/// i.e., the shape factor is the volume ratio of the shape of interest to a sphere/circle
template<typename T, typename PARTICLETYPE>
T getVolume( Particle<T,PARTICLETYPE> particle, [[maybe_unused]] T shapeFactor = T{1} )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  T volume;
  if constexpr(providesSurface(particle)) {
    if constexpr (D == 3) {
      volume = getSmoothIndicatorPtr(particle)->getVolume();
    }
    else {
      volume = getSmoothIndicatorPtr(particle)->getArea();
    }
  }
  else {
    const T radius = getRadius(particle);
    if constexpr(D == 3) {
      volume = (4./3.) * M_PI * util::pow(radius, D);
    }
    else {
      volume = M_PI * util::pow(radius, D);
    }
    volume *= shapeFactor;
  }
  return volume;
}

template<typename T, typename PARTICLETYPE>
T getDensity( Particle<T,PARTICLETYPE> particle, [[maybe_unused]] T shapeFactor = T{1} )
{
  using namespace descriptors;
  static_assert(providesMassOrDensity(particle),
      "Field PHYSPROPERTIES:MASS or PHYSPROPERTIES:DENSITY has to be provided");
  T density;
  // Always rather use the density field than the mass field
  if constexpr(providesDensity(particle)) {
    density = particle.template getField<PHYSPROPERTIES, DENSITY>();
  }
  else {
    const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
    density = mass / getVolume(particle, shapeFactor);
  }
  return density;
}

template<typename T, typename PARTICLETYPE>
T getMass( Particle<T,PARTICLETYPE> particle, [[maybe_unused]] T shapeFactor = T{1} )
{
  using namespace descriptors;
  static_assert(providesMassOrDensity(particle),
      "Field PHYSPROPERTIES:MASS or PHYSPROPERTIES:DENSITY has to be provided");
  T mass;
  // Always rather use the mass field than the density field
  if constexpr(providesMass(particle)) {
    mass = particle.template getField<PHYSPROPERTIES, MASS>();
  }
  else {
    const T density = particle.template getField<PHYSPROPERTIES, DENSITY>();
    mass = density * getVolume(particle, shapeFactor);
  }
  return mass;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getAccelerationStrd( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;

  if(providesAccelerationStrd(particle)){
    return particle.template getField<MOBILITY,ACCELERATION_STRD>();
  }
}

template<typename T, typename PARTICLETYPE>
void setAccelerationStrd( Particle<T,PARTICLETYPE> particle, Vector<T,PARTICLETYPE::d> acceleration)
{
  using namespace descriptors;
  static_assert(providesAccelerationStrd<PARTICLETYPE>(), "Field MOBILITY:ACCELERATION has to be provided");
  if constexpr(providesAccelerationStrd(particle)) {
    particle.template setField<MOBILITY, ACCELERATION_STRD>(acceleration);
  }
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getAcceleration( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;

  static_assert(providesForce<PARTICLETYPE>(), "Field FORCING:FORCE has to be provided");
  static_assert(providesMass<PARTICLETYPE>(), "Field PHYSPROPERTIES:MASS has to be provided");

  Vector<T,D> acceleration;
  Vector<T,D> force = getForce(particle);
  T mass = getMass(particle);
  for (unsigned iDim=0; iDim<D; ++iDim) {
    acceleration[iDim] = force[iDim] / mass;
  }
  return acceleration;
}

// no need for setAcceleration (non-std) as it is calculated from force and mass

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngAcceleration(
  Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;

  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  static_assert(providesTorque<PARTICLETYPE>(), "Field FORCING:TORQUE has to be provided");
  static_assert(providesMomentOfInertia<PARTICLETYPE>(), "Field PHYSPROPERTIES:MOFI has to be provided");
  Vector<T,Drot> angularAcceleration;
  Vector<T,Drot> torque( getTorque(particle));
  Vector<T,Drot> momentOfInertia( getMomentOfInertia(particle) );
  if constexpr (PARTICLETYPE::d == 3 && providesRotationMatrix<PARTICLETYPE>()) {
    const Matrix<T,3,3> rotationMatrix(getRotationMatrix(particle));
    angularAcceleration = calcAngAcceleration3D( torque, momentOfInertia, rotationMatrix );
  } else {
    for (unsigned iRot=0; iRot<Drot; ++iRot) {
      angularAcceleration[iRot] = torque[iRot] / momentOfInertia[iRot];
    }
  }
  return angularAcceleration;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngAccelerationStrd( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;

  if(providesAngAccelerationStrd(particle)){
    return particle.template getField<MOBILITY,ANG_ACC_STRD>();
  }
}

//Get global id
template<typename T, typename PARTICLETYPE>
auto getGlobalID( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  static_assert(providesID<PARTICLETYPE>(), "Field PARALLELIZATION:ID has to be provided");
  auto globalID = particle.template getField<PARALLELIZATION,ID>();
  return globalID;
}

//Get globiC (e.g allowing for the determination, whether dealing with particle centre)
template<typename T, typename PARTICLETYPE>
auto getGlobalIC( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  static_assert(providesID<PARTICLETYPE>(), "Field PARALLELIZATION:ID has to be provided");
  auto globalIC = particle.template getField<PARALLELIZATION,IC>();
  return globalIC;
}

//Get dynamics id
template<typename T, typename PARTICLETYPE>
unsigned short getDynamicsID( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  static_assert(providesDynamicsID<PARTICLETYPE>(),
                "Field DYNBEHAVIOUR:DYNAMICS_ID has to be provided");
  unsigned short dynamicsID = particle.template getField<DYNBEHAVIOUR,DYNAMICS_ID>();
  return dynamicsID;
}

//Get detaching state
template<typename T, typename PARTICLETYPE>
bool isDetaching( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  static_assert(providesDetaching<PARTICLETYPE>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  bool detaching = particle.template getField<DYNBEHAVIOUR,DETACHING>();
    return detaching;
  }

//Get extent of cuboid surface
//- throws error during static_cast, when assumption of cuboid shape is false
template<typename T, typename PARTICLETYPE>
auto getCuboidSurfaceExtent( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesSmoothIndicator<PARTICLETYPE>(), "Field FORCING:SINDICATOR has to be provided");
  using SIndicatorType = std::conditional_t<
    D == 2,
    SmoothIndicatorCuboid2D<T,T,true>,
    SmoothIndicatorCuboid3D<T,T,true>
  >;
  //Retrieve surface indicator
  auto sIndicatorPtr = getSmoothIndicatorPtr(particle);
  auto sIndicatorCuboidPtr = static_cast<SIndicatorType*>(sIndicatorPtr);
  auto& indicatorCuboid = sIndicatorCuboidPtr->getIndicator();
  //Retrieve extent
  if constexpr(D==2){
    Vector<T,D> extent(
      indicatorCuboid.getxLength(),
      indicatorCuboid.getyLength()
    );
    return extent;
  } else {
    Vector<T,D> extent(
      indicatorCuboid.getxLength(),
      indicatorCuboid.getyLength(),
      indicatorCuboid.getzLength()
    );
    return extent;
  }
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getCORoffset( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesCORoffset<PARTICLETYPE>(), "Field SURFACE:COR_OFFSET has to be provided");
  Vector<T,D> offsetCOR( particle.template getField<SURFACE,COR_OFFSET>() );
  return offsetCOR;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getElongation( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(providesElongation<PARTICLETYPE>(), "Field SURFACE:ELONGATION has to be provided");
  Vector<T,D> elongation( particle.template getField<SURFACE,ELONGATION>() );
  return elongation;
}

// Get particle enlargement for contact treatment
template<typename T, typename PARTICLETYPE>
T getEnlargementForContact( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  if constexpr (access::providesEnlargementForContactTreatment<PARTICLETYPE>()) {
    return particle.template getField<NUMERICPROPERTIES,ENLARGEMENT_FOR_CONTACT>();
  }
  else {
    return T{0.};
  }

  __builtin_unreachable();
}

template<typename T, typename PARTICLETYPE>
unsigned getContactMaterial( Particle<T,PARTICLETYPE> particle)
{
  using namespace descriptors;
  static_assert(providesContactMaterial(particle),
                "Field MECHPROPERTIES::MATERIAL has to be provided");
  return particle.template getField<MECHPROPERTIES, MATERIAL>();
}

//// Setter for specific particle fields

template<typename T, typename PARTICLETYPE>
void setDensity( Particle<T,PARTICLETYPE> particle, T density, [[maybe_unused]] T shapeFactor = T{1} )
{
  using namespace descriptors;
  static_assert(providesMassOrDensity(particle),
      "Field PHYSPROPERTIES:MASS or PHYSPROPERTIES:DENSITY has to be provided");

  if constexpr(providesDensity(particle)) {
    particle.template setField<PHYSPROPERTIES, DENSITY>(density);
  }
  if constexpr(providesMass(particle)) {
    const T mass = density * getVolume(particle, shapeFactor);
    particle.template setField<PHYSPROPERTIES, MASS>(mass);
  }
}

template<typename T, typename PARTICLETYPE>
void setMass( Particle<T,PARTICLETYPE> particle, T mass, [[maybe_unused]] T shapeFactor = T{1} )
{
  using namespace descriptors;
  static_assert(providesMassOrDensity(particle),
      "Field PHYSPROPERTIES:MASS or PHYSPROPERTIES:DENSITY has to be provided");

  if constexpr(providesMass(particle)) {
    particle.template setField<PHYSPROPERTIES, MASS>(mass);
  }
  if constexpr(providesDensity(particle)){
    const T density = mass / getVolume(particle, shapeFactor);
    particle.template setField<PHYSPROPERTIES, DENSITY>(density);
  }
}

template<typename T, typename PARTICLETYPE>
void setPosition( Particle<T,PARTICLETYPE> particle, Vector<T,PARTICLETYPE::d> position)
{
  using namespace descriptors;
  static_assert(providesPosition(particle),
      "Field GENERAL:POSITION has to be provided");
  particle.template setField<GENERAL, POSITION>(position);
}

template<typename T, typename PARTICLETYPE>
void setContactMaterial( Particle<T,PARTICLETYPE> particle, unsigned material )
{
  using namespace descriptors;
  static_assert(providesContactMaterial(particle),
      "Field MECHPROPERTIES::MATERIAL has to be provided");
  particle.template setField<MECHPROPERTIES, MATERIAL>(material);
}

template<typename T, typename PARTICLETYPE>
void setAngle( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle )
{
  using namespace descriptors;
  static_assert(providesAngle<PARTICLETYPE>(), "Field SURFACE:ANGLE has to be provided");
  particle.template setField<SURFACE, ANGLE>(utilities::dimensions::convert<
                  PARTICLETYPE::d>::serialize_rotation(angle));
}

template<typename T, typename PARTICLETYPE>
void setRotationMatrix( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::matrix> rotMatrix )
{
  using namespace descriptors;
  static_assert(providesRotationMatrix<PARTICLETYPE>(), "Field SURFACE:ROT_MATRIX has to be provided");
  particle.template setField<SURFACE,ROT_MATRIX>(rotMatrix);
}

template<typename T, typename PARTICLETYPE>
void setVelocity( Particle<T,PARTICLETYPE> particle, Vector<T,PARTICLETYPE::d> velocity )
{
  using namespace descriptors;
  static_assert(providesVelocity<PARTICLETYPE>(), "Field MOBILITY:VELOCITY has to be provided");
  particle.template setField<MOBILITY,VELOCITY>(velocity);
}

template<typename T, typename PARTICLETYPE>
void setAngularVelocity( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angVelocity)
{
  using namespace descriptors;
  static_assert(providesAngVelocity<PARTICLETYPE>(), "Field MOBILITY:ANG_VELOCITY has to be provided");
  particle.template setField<MOBILITY,ANG_VELOCITY>(utilities::dimensions::convert<
                                             PARTICLETYPE::d>::serialize_rotation(angVelocity));
}

template<typename T, typename PARTICLETYPE>
void setAngAccelerationStrd( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angAcceleration)
{
  using namespace descriptors;
  static_assert(providesAngAccelerationStrd<PARTICLETYPE>(), "Field MOBILITY:ANG_ACC_STRD has to be provided");
  particle.template setField<MOBILITY,ANG_ACC_STRD>(utilities::dimensions::convert<
                                             PARTICLETYPE::d>::serialize_rotation(angAcceleration));
}

// no need for setAngAcceleration (non-strd) as it is calculated from torque and moment of inertia

template<typename T, typename PARTICLETYPE>
void setForce( Particle<T,PARTICLETYPE> particle, Vector<T,PARTICLETYPE::d> force )
{
  using namespace descriptors;
  static_assert(providesForce<PARTICLETYPE>(), "Field FORCING:FORCE has to be provided");
  particle.template setField<FORCING,FORCE>(force);
}

template<typename T, typename PARTICLETYPE>
void setTorque( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> torque )
{
  using namespace descriptors;
  static_assert(providesTorque<PARTICLETYPE>(), "Field FORCING:TORQUE has to be provided");
  particle.template setField<FORCING,TORQUE>(utilities::dimensions::convert<
                                             PARTICLETYPE::d>::serialize_rotation(torque));
}

template<typename T, typename PARTICLETYPE>
void setMomentOfInertia( Particle<T,PARTICLETYPE> particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> mofi)
{
  using namespace descriptors;
  static_assert(providesMomentOfInertia<PARTICLETYPE>(), "Field PHYSPROPERTIES:MOFI has to be provided");
  particle.template setField<PHYSPROPERTIES,MOFI>(utilities::dimensions::convert<
                                              PARTICLETYPE::d>::serialize_rotation(mofi));
}

// Consisting of normal and tangential adhesion component of particle.
template<typename T, typename PARTICLETYPE>
void setAdhesion( Particle<T,PARTICLETYPE>& particle, Vector<T,2> adhesion )
{
  using namespace descriptors;
  static_assert(providesAdhesion<PARTICLETYPE>(), "Field FORCING:ADHESION has to be provided");
  particle.template setField<FORCING,ADHESION>(adhesion);
}

template<typename T, typename PARTICLETYPE>
void setInvalid( Particle<T,PARTICLETYPE> particle, bool value=true )
{
  using namespace descriptors;
  static_assert(providesInvalid<PARTICLETYPE>(), "Field GENERAL:INVALID has to be provided");
  particle.template setField<GENERAL,INVALID>( value );
}

template<typename T, typename PARTICLETYPE>
void setValid( Particle<T,PARTICLETYPE> particle, bool value=true )
{
  setInvalid(particle, !value);
}

template<typename T, typename PARTICLETYPE>
void setActive( Particle<T,PARTICLETYPE> particle, bool value=true)
{
  using namespace descriptors;
  static_assert(providesActive<PARTICLETYPE>(), "Field DYNBEHAVIOUR:ACTIVE has to be provided");
  particle.template setField<DYNBEHAVIOUR,ACTIVE>( value );
}

template<typename T, typename PARTICLETYPE>
void setInactive( Particle<T,PARTICLETYPE> particle, bool value=true)
{
  setActive(particle, !value);
}

//Set smooth indicator pointer
template<typename T, typename PARTICLETYPE>
void setSmoothIndicatorPtr( Particle<T,PARTICLETYPE> particle,
    SmoothIndicatorF<T,T,PARTICLETYPE::d,true>* sindicator )
{
  using namespace descriptors;
  static_assert(providesSmoothIndicator<PARTICLETYPE>(), "Field FORCING:SINDICATOR has to be provided");
  particle.template setField<SURFACE,SINDICATOR>(sindicator);
}

//Set radius
template<typename T, typename PARTICLETYPE>
void setRadius( Particle<T,PARTICLETYPE>& particle, T radius )
{
  using namespace descriptors;
  static_assert(providesRadius<PARTICLETYPE>(), "Field PHYSPROPERTIES:RADIUS has to be provided");
  particle.template setField<PHYSPROPERTIES,RADIUS>(radius);
}

template<typename T, typename PARTICLETYPE>
void setGlobalID( Particle<T,PARTICLETYPE> particle, std::size_t id )
{
  using namespace descriptors;
  static_assert(providesID<PARTICLETYPE>(), "Field PARALLELIZATION:ID has to be provided");
  particle.template setField<PARALLELIZATION,ID>(id);
}

template<typename T, typename PARTICLETYPE>
void setGlobalIC( Particle<T,PARTICLETYPE> particle, int id )
{
  using namespace descriptors;
  static_assert(providesID<PARTICLETYPE>(), "Field PARALLELIZATION:ID has to be provided");
  particle.template setField<PARALLELIZATION,IC>(id);
}

template<typename T, typename PARTICLETYPE>
void setDynamicsID( Particle<T,PARTICLETYPE>& particle, unsigned short dynamicsID )
{
  using namespace descriptors;
  static_assert(providesDynamicsID<PARTICLETYPE>(),
                "Field DYNBEHAVIOUR:DYNAMICS_ID has to be provided");
  particle.template setField<DYNBEHAVIOUR,DYNAMICS_ID>(dynamicsID);
}

template<typename T, typename PARTICLETYPE>
void setDetaching( Particle<T,PARTICLETYPE>& particle, bool value)
{
  using namespace descriptors;
  static_assert(providesDetaching<PARTICLETYPE>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  particle.template setField<DYNBEHAVIOUR,DETACHING>(value);
}

template<typename T, typename PARTICLETYPE>
void setCORoffset( Particle<T,PARTICLETYPE> particle, Vector<T, PARTICLETYPE::d> offsetCOR )
{
  using namespace descriptors;
  static_assert(providesCORoffset<PARTICLETYPE>(), "Field SURFACE:COR_OFFSET has to be provided");
  particle.template setField<SURFACE,COR_OFFSET>(offsetCOR);
}

template<typename T, typename PARTICLETYPE>
void setElongation( Particle<T,PARTICLETYPE> particle, Vector<T,PARTICLETYPE::d> elongation )
{
  using namespace descriptors;
  static_assert(providesElongation<PARTICLETYPE>(), "Field SURFACE:ELONGATION has to be provided");
  particle.template setField<SURFACE,ELONGATION>(elongation);
}

template<typename T, typename PARTICLETYPE>
void setEnlargementForContact(Particle<T, PARTICLETYPE>& particle, T value)
{
  using namespace descriptors;
  static_assert(access::providesEnlargementForContactTreatment<PARTICLETYPE>(), "Field NUMERICPROPERTIES:ENLARGEMENT_FOR_CONTACT has to be provided");
  particle.template setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(value);
}

/// Check if motion is enabled
template<typename T, typename PARTICLETYPE>
bool isMotionComputationEnabled( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;

  // if the field is not provided, return true by default
  if constexpr(!providesComputeMotion<PARTICLETYPE>()) {
    return true;
  }
  else {
    return particle.template getField<DYNBEHAVIOUR,COMPUTE_MOTION>();
  }

  __builtin_unreachable();
}

/// Check if contact should be regarded (specification for a single particle)
template<typename T, typename PARTICLETYPE>
bool isContactComputationEnabled( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;

  // if the field is not provided, return true by default
  // (meaning that the contact forces are computed and applied to the force field)
  if constexpr(!providesComputeContact<PARTICLETYPE>()) {
    return true;
  }
  else {
    return particle.template getField<DYNBEHAVIOUR,COMPUTE_CONTACT>();
  }

  __builtin_unreachable();
}

/// Check if contact should be regarded (interaction of two known particles)
template<typename T, typename PARTICLETYPE>
bool isContactComputationEnabled(Particle<T,PARTICLETYPE>& particleA,
                      Particle<T,PARTICLETYPE>& particleB)
{
  return (isContactComputationEnabled(particleA) || isContactComputationEnabled(particleB));
}

template<typename T, typename PARTICLETYPE>
void enableMotionComputation(Particle<T, PARTICLETYPE>& particle, bool value = true)
{
  using namespace descriptors;
  static_assert(providesComputeMotion<PARTICLETYPE>(), "Field DYNBEHAVIOUR:COMPUTE_MOTION has to be provided");
  particle.template setField<DYNBEHAVIOUR, COMPUTE_MOTION>(value);
}

template<typename T, typename PARTICLETYPE>
void enableContactComputation(Particle<T, PARTICLETYPE>& particle, bool value = true)
{
  using namespace descriptors;
  static_assert(providesComputeContact<PARTICLETYPE>(), "Field DYNBEHAVIOUR:COMPUTE_CONTACT has to be provided");
  particle.template setField<DYNBEHAVIOUR, COMPUTE_CONTACT>(value);
}

template<typename T, typename PARTICLETYPE>
void disableMotionComputation(Particle<T, PARTICLETYPE>& particle, bool value = true)
{
  enableMotionComputation(particle, !value);
}

template<typename T, typename PARTICLETYPE>
void disableContactComputation(Particle<T, PARTICLETYPE>& particle, bool value = true)
{
  enableContactComputation(particle, !value);
}

template<typename T, typename PARTICLETYPE>
void setRestingParticle( Particle<T,PARTICLETYPE> particle)
{
  if constexpr(providesVelocity(particle)){
      setVelocity(particle, Vector<T,PARTICLETYPE::d>(0.));
  }

  if constexpr(providesAccelerationStrd(particle)){
      setAccelerationStrd(particle, Vector<T,PARTICLETYPE::d>(0.));
  }


  if constexpr(providesAngVelocity(particle)){
      setAngularVelocity(particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.));
  }

  if constexpr(providesAngAccelerationStrd(particle)){
      setAngAccelerationStrd(particle, Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.));
  }
}

} //namespace access

} //namespace particles

} //namespace olb




#endif
