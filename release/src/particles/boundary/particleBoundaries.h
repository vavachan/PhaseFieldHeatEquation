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



#ifndef PARTICLE_BOUNDARIES_H
#define PARTICLE_BOUNDARIES_H


namespace olb {

namespace particles {

namespace boundaries {


/// Velocity wall reflection
/// - applying resurfacing for current timeStep
/// - inverting normal velocity component at wall impact for consequtive timeStep
template<bool useCubicBounds=false, typename T, typename PARTICLETYPE>
void velocityWallReflection( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  T coefficientOfRestitution=1.0 )
{
  using namespace descriptors;
  //Execute wall treatment
  doAtParticleWallContact<useCubicBounds>( particle, solidBoundary,
    [&]( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal, T penetrationDepth ){
    //Retrieve velocity
    auto velocity = access::getVelocity( particle );
    auto position = access::getPosition( particle );
    //Modify velocity and set after-reflection position
    for (int iDim=0; iDim<PARTICLETYPE::d; ++iDim){
      //1. Resurface
      position[iDim] += normal[iDim] * 2.*penetrationDepth;
      //2. Reflection
      velocity[iDim] += std::abs(normal[iDim]) * -2.*velocity[iDim];
      velocity[iDim] *= coefficientOfRestitution;
    }
    //Apply velocity
    access::setVelocity( particle, velocity );
    access::setPosition( particle, position );
  });
}

/// Alias for cubic version of velocity wall reflection
/// - mainly used as legacy support
template<typename T, typename PARTICLETYPE>
void cuboidVelocityWallReflection( Particle<T,PARTICLETYPE>& particle,
  Vector<T,PARTICLETYPE::d> origin, Vector<T,PARTICLETYPE::d> end )
{
  //Create solid wall
  Vector<T,PARTICLETYPE::d> extent = end-origin;
  IndicatorCuboid3D<T> cuboid(extent, origin);
  std::unique_ptr<IndicatorF<T,PARTICLETYPE::d>> boundaryIndicator =
      std::make_unique<IndicInverse<T,PARTICLETYPE::d>>(cuboid);
  SolidBoundary<T,PARTICLETYPE::d> solidBoundary(std::move(boundaryIndicator), 2, 0);
  //Execute velocity wall reflection
  velocityWallReflection<true>( particle, solidBoundary );
}

/// Wall slip
/// - applying resurfacing for current timeStep
/// - setting normal velocity component to 0 at wall impact for consequtive timeStep
template<bool useCubicBounds=false, typename T, typename PARTICLETYPE>
void wallSlip( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
{
  using namespace descriptors;
  //Execute wall treatment
  doAtParticleWallContact<useCubicBounds>( particle, solidBoundary,
    [&]( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal, T penetrationDepth ){
    //Retrieve velocity
    auto velocity = access::getVelocity( particle );
    auto position = access::getPosition( particle );
    //Modify velocity and set after-reflection position
    for (int iDim=0; iDim<PARTICLETYPE::d; ++iDim){
      //1. Resurface
      position[iDim] += normal[iDim] * 2.*penetrationDepth;
      //2. Velocity swallow //TODO: include momentum transfer treatment
      velocity[iDim] = 0.;
    }
    //Apply velocity
    access::setVelocity( particle, velocity );
    access::setPosition( particle, position );
  });
}

/// Wall capture
/// - setting particle field ACTIVE=false at wall impact
/// - setting resting particle
template<bool useCubicBounds=false, typename T, typename PARTICLETYPE>
void wallCapture( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
{
  using namespace descriptors;
  //Execute wall treatment
    bool isActive = particle.template getField<DYNBEHAVIOUR,ACTIVE>();
  if (isActive){
  doAtParticleWallContact<useCubicBounds>( particle, solidBoundary,
    [&]( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal ){
    //Deactivate particle

    access::setInactive( particle );
    access::setRestingParticle( particle );
  });
  }

}

/// Wall capture based on material rather than SolidBoundary
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at material contact
/// - setting resting particle
template<typename T, typename PARTICLETYPE>
void materialCapture( Particle<T,PARTICLETYPE>& particle,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( materialIndicator, particle );
    if (vicinity){
      //Deactivate particle
      access::setInactive( particle );
      access::setRestingParticle( particle );
    }
  }
}

/// Wall capture with material awareness
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at wall impact
/// - setting resting particle
/// - represents combination of wallCapture and materialCapture
/// - intended to use aprior material check to speed up complex SolidBoudaries carying STLs
template<typename T, typename PARTICLETYPE>
void wallCaptureMaterialAware( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( materialIndicator, particle );
    if (vicinity){
      //Apply wall capture
      constexpr bool useCubicBounds=false;
      wallCapture<useCubicBounds>(particle, solidBoundary);
    }
  }
}


/// Escape
/// - setting particle field ACTIVE=false
template<bool useCubicBounds=false, typename T, typename PARTICLETYPE>
void escape ( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
{
  using namespace descriptors;
  //Execute wall treatment
  doAtParticleWallContact<useCubicBounds>( particle, solidBoundary,
    [&]( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal ){
    //Deactivate particle
    access::setInactive( particle );
  });
}

/// Escape boundary based on material rather than SolidBoundary
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at material contact
template<typename T, typename PARTICLETYPE>
void materialEscape( Particle<T,PARTICLETYPE>& particle,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( materialIndicator, particle );
    if (vicinity){
      //Deactivate particle
      access::setInactive( particle );
    }
  }
}

/// Escape boundary with material awareness
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at wall impact
/// - represents combination of wallCapture and materialCapture
/// - intended to use aprior material check to speed up complex SolidBoudaries carying STLs
template<typename T, typename PARTICLETYPE>
void escapeMaterialAware( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( materialIndicator, particle );
    if (vicinity){
      constexpr bool useCubicBounds=false;
      boundaries::escape<useCubicBounds>(particle, solidBoundary);
    }
  }
}


/// Escape and capture based on material rather than SolidBoundary
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at material contact
/// - setting resting particle (only capture)
template<typename T, typename PARTICLETYPE>
void materialCaptureAndEscape( Particle<T,PARTICLETYPE>& particle,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& captureMaterialIndicator,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& escapeMaterialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( captureMaterialIndicator, particle );
    if (vicinity){
      access::setInactive( particle );
      access::setRestingParticle( particle );
    }
    else {
      vicinity = checkMaterialVicinity( escapeMaterialIndicator, particle );
      if (vicinity){
        access::setInactive( particle );
      }
    }
  }
}

/// Escape boundary with material awareness
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at wall impact
/// - represents combination of wallCapture and materialCapture
/// - intended to use aprior material check to speed up complex SolidBoudaries carying STLs
template<typename T, typename PARTICLETYPE>
void wallCaptureAndEscapeMaterialAware( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& captureMaterialIndicator,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& escapeMaterialIndicator )
{
  using namespace descriptors;
  //Check whether still active
  bool isActive = access::isActive( particle );
  if (isActive){
    bool vicinity = checkMaterialVicinity( captureMaterialIndicator, particle );
    if (vicinity){
      constexpr bool useCubicBounds=false;
      boundaries::wallCapture<useCubicBounds>(particle, solidBoundary);
    }
    else {
      vicinity = checkMaterialVicinity( escapeMaterialIndicator, particle );
      if (vicinity){
        constexpr bool useCubicBounds=false;
        boundaries::escape<useCubicBounds>(particle, solidBoundary);
      }
    }
  }
}


 template<bool useCubicBounds=false, typename T, typename PARTICLETYPE>
void wallCaptureSpheroidSubgrid( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,ACTIVE>(), "Field DYNBEHAVIOUR:ACTIVE has to be provided");

  doAtParticleWallContactSpheroid<useCubicBounds>( particle, solidBoundary,
    [&]( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal ){
    //Deactivate particle
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( false );
 });
}


  /// Wall capture with material awareness
/// - implies necessity of SuperGeometry, so no DEM only possible
/// - setting particle field ACTIVE=false at wall impact
/// - represents combination of wallCapture and materialCapture
/// - intended to use aprior material check to speed up complex SolidBoudaries carying STLs
template<typename T, typename PARTICLETYPE>
void wallCaptureMaterialAwareSpheroidSubgrid( Particle<T,PARTICLETYPE>& particle,
  SolidBoundary<T,PARTICLETYPE::d>& solidBoundary,
  SuperIndicatorMaterial<T,PARTICLETYPE::d>& materialIndicator )
{
  using namespace descriptors;
  static_assert(PARTICLETYPE::template providesNested<DYNBEHAVIOUR,ACTIVE>(), "Field DYNBEHAVIOUR:ACTIVE has to be provided");
  //Check whether still active
  bool isActive = particle.template getField<DYNBEHAVIOUR,ACTIVE>();
  if (isActive){
    bool vicinity = checkMaterialVicinity( materialIndicator, particle );
    if (vicinity){
      //Apply wall capture
      //constexpr bool useCubicBounds=false;
  doAtParticleWallContactSpheroid(particle, solidBoundary,
  [&] ( Particle<T,PARTICLETYPE>& particle,
    Vector<T,PARTICLETYPE::d>& normal ) {
    //Deactivate particle
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( false );
 }
 );
    }
  }
}


} //namespace boundaries

} //namespace particles

} //namespace olb


#endif
