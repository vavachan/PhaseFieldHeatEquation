/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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


/* Treatment of particle dymamics while detaching from surfaces
 * - mainly used for wall-flow applicatio:
 *     https://gitlab.com/openlb/olb/-/blob/private/nicolas/main4/apps/nicolas/A_WallFlow/wallFlow/README.md
 * - for now, mainly intended to be used in 3D
 *
*/

#ifndef PARTICLE_SURFACE_DETACHMENT_H
#define PARTICLE_SURFACE_DETACHMENT_H

//#define VERBOSE_ADHESION_FORCE

#include <cassert>

namespace olb {

namespace particles {

namespace interaction {


//Calculate 2D position on cylcic hull
template<typename T>
Vector<T,2> pos2DOnCyclicHull(Vector<T,2> position, T radius, T angle){
  Vector<T,2> relPos( std::cos(angle-M_PI)*radius,
                     -std::sin(angle-M_PI)*radius );
  return position + relPos;
}

//Eccentric position of position rotated around center of roateion (COR)
template<typename T>
Vector<T,3> eccentricPosition3D(
  Vector<T,3>& position, Vector<T,3> relPosCOR,
  T angle1DRad, unsigned axis, bool verbose=false)
{
  //Retrieve new absolute centre of rotation
  Vector<T,3> CORnew = position+relPosCOR;
  //Calculate radius of rotation
  T radiusRot = norm(relPosCOR);
  //Convert to 2D space
  Vector<T,2> CORnew2D;
  if (axis==0){
    CORnew2D[0] = CORnew[1];
    CORnew2D[1] = CORnew[2];
  } else if (axis==1){
    CORnew2D[0] = CORnew[0];
    CORnew2D[1] = CORnew[2];
  } else if (axis==2){
    CORnew2D[0] = CORnew[0];
    CORnew2D[1] = CORnew[1];
  } else {
    std::cerr << "ERROR: Unknown axis " << axis << "!" << std::endl;
  }
  //Retrieve position on cyclic hull
  Vector<T,2> posEccentric2D = pos2DOnCyclicHull(CORnew2D,radiusRot,angle1DRad);
  //Convert to 3D space
  Vector<T,3> posEccentric;
  if (axis==0){
    posEccentric[0] = position[0];
    posEccentric[1] = posEccentric2D[0];
    posEccentric[2] = posEccentric2D[1];
  } else if (axis==1){
    posEccentric[0] = posEccentric2D[0];
    posEccentric[1] = position[1];
    posEccentric[2] = posEccentric2D[1];
  } else if (axis==2){
    posEccentric[0] = posEccentric2D[0];
    posEccentric[1] = posEccentric2D[1];
    posEccentric[2] = position[2];
  } else {
    std::cerr << "ERROR: Unknown axis " << axis << "!" << std::endl;
  }
  //Output
  if (verbose){
    OstreamManager clout(std::cout, "EccentricRot");
    clout << "---Eccentric Position ---" << std::endl;
    clout << "Position= " << position << std::endl;
    clout << "RelPosCOR= " << relPosCOR << std::endl;
    clout << "CORnew= " << CORnew << std::endl;
    clout << "radiusRot= " << radiusRot << std::endl;
    clout << "angle1DRad (repr. rotation at initial position)= " << angle1DRad << std::endl;
    clout << "posEccentric2D= " << posEccentric2D << std::endl;
    clout << "-------------------------" << std::endl;
  }
  return posEccentric;
}

//Get constant cuboid diagonal angle
template<typename T>
T getCuboid3DDiagonalAngle1D( Vector<T,3>& extent ){
  return -std::atan2(-extent[2],extent[0]);
}




//Handle detachment
// - Corrects eccentric particle position, and updates
//   offset of centre of rotation
template<typename T, typename PARTICLETYPE>
void handleDetachment( SolidBoundary<T,PARTICLETYPE::d>& wall,
                       Vector<T,PARTICLETYPE::d>& mainFlowDirection,
                       Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  using namespace access;
  static_assert(PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>(),
                "Field SURFACE:COR_OFFSET has to be provided");
  //Retrieve surface normal of surface closest to particle
  auto surfaceNormal = boundaries::getNormalOnClosestSurface( wall, particle );
  //Retrieve rotation axis for detachmant
  auto rotationAxisNormal = crossProduct( surfaceNormal, mainFlowDirection );
  unsigned rotationAxis = util::maxElementAbsPos( rotationAxisNormal );
  //Calculate relative angle between particle and surface
  auto normalParticleSurface = getSurfaceNormal( particle );
  auto orientationRelative = util::angleBetweenVectors(
    normalParticleSurface, surfaceNormal );
  //Retrieve geometry specific angle offset
  auto extent = getCuboidSurfaceExtent( particle );
  //Retrieve current position and offset of center of rotation (COR)
  auto position = getPosition( particle );
  auto CORoffset = particle.template getField<SURFACE,COR_OFFSET>();
  //Retrieve eccentric position
  decltype(position) posEccentric;
  if constexpr (is3D<PARTICLETYPE>()){
    T constAngleOffset = getCuboid3DDiagonalAngle1D( extent );
    //Calculate effective angle around rotationAxis (considering const offset)
    T signedDirection = util::maxElementAbs(surfaceNormal);
    T angleEffective1D = orientationRelative[rotationAxis]
                       + constAngleOffset * signedDirection;
    //Calculate eccentric position of particle centre
    posEccentric =  eccentricPosition3D( position,
      CORoffset, angleEffective1D, rotationAxis );
  } else {
    std::cerr << "ERROR: 2D version not implemented yet!" << std::endl;
  }
  //Calculate difference between original and new position and update CORoffset
  auto posCubeDiffEcc = posEccentric-position;
  CORoffset -= posCubeDiffEcc;
  //Set new quantities to particle
  particle.template setField<GENERAL,POSITION>( posEccentric );
  particle.template setField<SURFACE,COR_OFFSET>( CORoffset );
}

template<typename T, unsigned D>
void getDetachmentAxes( Vector<T,D> mainFlowDirection, Vector<T,D> surfaceNormal,
  unsigned short& axisFlow, unsigned short& axisSurface, unsigned short& axisRot )
{
  axisFlow = util::maxElementAbsPos( mainFlowDirection );
  axisSurface = util::maxElementAbsPos( surfaceNormal );
  if (axisFlow==0){
    if      (axisSurface==1){ axisRot=2; }
    else if (axisSurface==2){ axisRot=1; }
    else { std::cerr << "Error: Axes not unique!" << std::endl; }
  } else if (axisFlow==1){
    if      (axisSurface==0){ axisRot=2; }
    else if (axisSurface==2){ axisRot=0; }
    else { std::cerr << "Error: Axes not unique!" << std::endl; }
  } else {
    if      (axisSurface==0){ axisRot=1; }
    else if (axisSurface==1){ axisRot=0; }
    else { std::cerr << "Error: Axes not unique!" << std::endl; }
  }
}



template<typename T, typename PARTICLETYPE>
void setCORcuboid3Dflush( SolidBoundary<T,PARTICLETYPE::d>& wall,
                            Vector<T,PARTICLETYPE::d>& mainFlowDirection,
                            Particle<T,PARTICLETYPE>& particle )
{
  using namespace particles::access;
  using namespace descriptors;
  //Retrieve cuboid extent
  auto extent = getCuboidSurfaceExtent( particle );
  //Retrieve surface normal of attached surface
  auto surfaceNormal = boundaries::getNormalOnClosestSurface( wall, particle );
  const unsigned short axisFlow = util::maxElementAbsPos( mainFlowDirection );
  const unsigned short axisSurface = util::maxElementAbsPos( surfaceNormal );
  const unsigned axisZunrotated = 2; //Assuming z-oriented flush placement
  //Calculate COR
  Vector<T,PARTICLETYPE::d> CORoffset(0.);
  CORoffset[axisFlow] =  .5*extent[axisFlow] * mainFlowDirection[axisFlow];
  CORoffset[axisSurface] = -.5*extent[axisZunrotated] * surfaceNormal[axisSurface];
  //Set COR offset
  particle.template setField<SURFACE,COR_OFFSET>(CORoffset);
}

template<typename T, typename PARTICLETYPE>
void initializeDetachment( SolidBoundary<T,PARTICLETYPE::d>& wall,
                           Particle<T,PARTICLETYPE>& particle,
                           Vector<T,PARTICLETYPE::d>& mainFlowDirection )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,ACTIVE>(),
                "Field DYNBEHAVIOUR:ACTIVE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DETACHING>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  //Set COR
  if constexpr (access::is3D<PARTICLETYPE>()){
    setCORcuboid3Dflush( wall, mainFlowDirection, particle );
  } else {
    std::cerr << "ERROR: 2D version not implemented yet!" << std::endl;
  }
  //Set detaching and active state
  particle.template setField<DYNBEHAVIOUR,DETACHING>(true);
  particle.template setField<DYNBEHAVIOUR,ACTIVE>( false );
  //Set contact state if provided
  if constexpr(access::providesComputeContact<PARTICLETYPE>()){
    particle.template setField<DYNBEHAVIOUR,COMPUTE_CONTACT>(false);
  }
}

//TODO: for now only for cuboid 3D
//INFO: Prototype version which leads to particle oscillation.
template<typename T, typename PARTICLETYPE>
void applyAdhesionForce( SolidBoundary<T,PARTICLETYPE::d>& wall,
                         Vector<T,PARTICLETYPE::d>& mainFlowDirection,
                         Particle<T,PARTICLETYPE>& particle )
{
  using namespace particles::access;
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DETACHING>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  static_assert(PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>(),
                "Field SURFACE:COR_OFFSET has to be provided");
  //Retrieve surface normal and distance of surface closest to particle
  auto surfaceNormal = boundaries::getNormalOnClosestSurface( wall, particle );
  //Retrieve cuboid surface extent
  auto extent = getCuboidSurfaceExtent( particle );
  const unsigned axisZunrotated = 2; //Assuming z-oriented flush placement
  T offsetCentreWall = .5*extent[axisZunrotated];
  //Retrieve current CORoffset for force lever
  auto CORoffset = getCORoffset( particle );
  //Evaluate axes
  const unsigned short axisFlow = util::maxElementAbsPos( mainFlowDirection );
  const unsigned short axisSurface = util::maxElementAbsPos( surfaceNormal );

  T surfaceDistance = -CORoffset[axisSurface]-offsetCentreWall;
  T lever = CORoffset[axisFlow];

  //Retrieve torque
  Vector<T,3> torque = getTorque( particle );
  Vector<T,3> force = getForce( particle );

  //TODO: adapt hardcoded axis
  unsigned short rotAxis = 1;

  //Retrieve normal and tangential adhesion component of particle.
  Vector<T,2> adhesion = getAdhesion( particle );
  T adhesionNormal = adhesion[0];
  T adhesionTangential = adhesion[1];

  //When distance greater than 500 Angstroem (for now random)
  //TODO: should be depending on Lennard-Jones Potential as well
  T fullForceDist = 500e-10;

  //Only consider force up fullForceDist
  if (surfaceDistance <= fullForceDist){

    //Ensure force to only act as adhesion and not as repulsion
    if (surfaceDistance > 0){

    //TODO: linear for now. Should be Lennard-Jones instead
    T distFactor = surfaceDistance/fullForceDist;

    //Calculate force and torque contributions due to adhesion
    T addhesionForceNormalContribution = -distFactor*adhesionNormal;
    T adhesionTorqueContribution = addhesionForceNormalContribution*lever;

    //Calculate new force and torque
    force[axisSurface] += addhesionForceNormalContribution * surfaceNormal[axisSurface];
    torque[rotAxis] += adhesionTorqueContribution;

    //Apply new force and torque
    particle.template setField<FORCING,FORCE>( force );
    particle.template setField<FORCING,TORQUE>( torque );

#ifdef VERBOSE_ADHESION_FORCE
      auto orient = getAngle(particle);
      std::cout << "Adhesion:" << std::endl
                << " - distance=" << surfaceDistance
                << " (factor=" << (distFactor) << ")"
                << ", angle=" << orient[1] << std::endl
                << " - forceNormalContrib=" << addhesionForceNormalContribution
                << ", torqueContrib=" << adhesionTorqueContribution << std::endl
                << " - resForce=" << force[axisSurface]
                << ", resTorque=" << torque[rotAxis] << std::endl
                << " - resForce3D=" << force
                << ", resTorque3D=" << torque
                << std::endl;
#endif

    }
  }
}

/// Calculation of rotation induced normal force
template<typename T, typename PARTICLETYPE>
T getRotationInducedNormalForce(Particle<T,PARTICLETYPE>& particle,
  Vector<T,PARTICLETYPE::d>& surfaceNormal, Vector<T,PARTICLETYPE::d>& mainFlowDirection)
{
  using namespace particles::access;
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>(),
                "Field SURFACE:COR_OFFSET has to be provided");

  //Calculate tilt axis
  Vector<T,3> tiltAxis = crossProduct( mainFlowDirection, surfaceNormal );
  unsigned axisRot = util::maxElementAbsPos(tiltAxis);
  unsigned axisFlow = util::maxElementAbsPos(mainFlowDirection);
  //Retrieve particle torque
  auto torque = getTorque(particle);
  //Calculate tilt torque
  T torqueTilt = torque[axisRot]*-tiltAxis[axisRot];
  //Retrieve lever by COR offset
  auto CORoffset = getCORoffset( particle );
  T lever = CORoffset[axisFlow];
  //Calculate rotation induced normal force (ensure no division by zero)
  T forceNormalRot = 0.;
  if (lever!=0){ forceNormalRot = torqueTilt / lever; }
  //Return rotation induced normal force
  return forceNormalRot;
}


/// Check adhesion and return true if still adhering
template<typename T, typename PARTICLETYPE>
bool checkAdhesion( SolidBoundary<T,PARTICLETYPE::d>& wall,
                    Vector<T,PARTICLETYPE::d>& mainFlowDirection,
                    Particle<T,PARTICLETYPE>& particle )
{
  using namespace particles::access;
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,ACTIVE>(),
                "Field DYNBEHAVIOUR:ACTIVE has to be provided");
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DETACHING>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  static_assert(PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>(),
                "Field SURFACE:COR_OFFSET has to be provided");

  //Retrieve surface normal of surface closest to particle
  auto surfaceNormal = boundaries::getNormalOnClosestSurface( wall, particle );
  //Retrieve current adhesion (normal component)
  Vector<T,2> forceAdhesion = getAdhesion( particle );
  T forceAdhesionN = forceAdhesion[0];
  //Calculate rotation induced normal force
  T forceNormalRot = getRotationInducedNormalForce( particle, surfaceNormal, mainFlowDirection );
  //Evaluate threshold
  if (forceNormalRot>forceAdhesionN){
    //Set active, when passing threshold once
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( true );
    return false;
  } else {
    return true;
  }
}



template<typename T, typename PARTICLETYPE>
void evaluateDetachmentState( SolidBoundary<T,PARTICLETYPE::d>& wall,
                              Particle<T,PARTICLETYPE>& particle,
                              T tiltThreshold = 0.3*M_PI )
{
  using namespace particles::access;
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,DETACHING>(),
                "Field DYNBEHAVIOUR:DETACHING has to be provided");
  static_assert(PARTICLETYPE::template providesNested<SURFACE,COR_OFFSET>(),
                "Field SURFACE:COR_OFFSET has to be provided");
  //Retrieve position and radius
  auto position = getPosition( particle );
  auto radius = getRadius( particle );
  //Retrieve relative orientation to surface
  auto normalParticleSurface = getSurfaceNormal( particle );
  auto orientationRelative = boundaries::getRelativeSurfaceOrientation(
    wall, position, normalParticleSurface, radius );
  //Evaluate tilt intensity
  T tiltIntensity = util::max_element(abs( orientationRelative ));
  if ( tiltIntensity > tiltThreshold){ //For now, random condition here
    //Set state (standard verlet integration with collision)
    particle.template setField<DYNBEHAVIOUR,DETACHING>( false );
    //TODO: Generalize call as detachment is not directly connected to contact
    if constexpr( providesComputeContact<PARTICLETYPE>() ){
      particle.template setField<DYNBEHAVIOUR,COMPUTE_CONTACT>( true );
    }
    //Reset centre of rotation
    particle.template setField<SURFACE,COR_OFFSET>( Vector<T,PARTICLETYPE::d>(0) );
  }
}


/// Check particle re-deposition
// - assuming that particle activity can be measured by its kinetic energy, as long
//   not undergoing the inflection point (eKin=0) during collision. In this case,
//   however, the particle experiences a significant force, which is therefore
//   evaluated as well. When both force and kinetic energy are negligible, particle
//   is deactivated after the predefined physical time.
template<typename T, typename PARTICLETYPE, typename DESCRIPTOR>
bool checkParticleReDeposition(
  Particle<T,PARTICLETYPE>& particle,
  UnitConverter<T,DESCRIPTOR> const& converter,
  T forceAbsoluteThreshold, T kinEnergyThreshold,
  T timeNoActiveThreshold, std::size_t iTinterval=1 )
{
  using namespace descriptors;
  //Retrieve force and kinetic energy
  auto force = access::getForce( particle );
  T forceAbs = norm(force);
  T kinEnergy = dynamics::calcKineticEnergy( particle );
  //Evaluate force threshold
  if (forceAbs<forceAbsoluteThreshold){
    //Evaluate kinetic energy threshold
    if (kinEnergy<kinEnergyThreshold){
      //Retrieve deactive counter
      std::size_t counterDeactive =
        particle.template getField<DYNBEHAVIOUR,COUNTER<ACTIVE>>();
      //Increment deactive counter and calculate no no active time
      counterDeactive+=iTinterval;
      T noActiveTime = converter.getPhysTime( counterDeactive );
      //If no active time not passed
      if (noActiveTime < timeNoActiveThreshold){
        //Write counter
        particle.template setField<DYNBEHAVIOUR,COUNTER<ACTIVE>>(counterDeactive);
      } else {
        //If time passed, deactivate particle
        particle.template setField<DYNBEHAVIOUR,ACTIVE>(false);
        if constexpr( access::providesComputeContact<PARTICLETYPE>() ){
          particle.template setField<DYNBEHAVIOUR,COMPUTE_CONTACT>(false);
        }
        return true;
      } //if (noActiveTime < noActiveTimeDeactivate)
    } else {
      //If kinetic energy threshold not passed, reset counter
      particle.template setField<DYNBEHAVIOUR,COUNTER<ACTIVE>>(0);
    } //if (kinEnergy<eKinThres)
  }
  return false;
}




} //namespace interaction

} //namespace particles

} //namespace olb


#endif
