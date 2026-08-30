/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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

//TODO: WARNING: For now, there is a lot of code duplication, which has to be sorted out in the future.

/* This file contains particle creator functions.
 * Those are generally meant for users to be able to create particles
 * in a comfortable way. As this section might be subject to frequent changes,
 * those should, however, be avoided inside header files. Here, direct field access
 * as shown below should be preferred!
 *
*/


#ifndef PARTICLE_CREATOR_FUNCTIONS_3D_H
#define PARTICLE_CREATOR_FUNCTIONS_3D_H


#include "particles/particleSystem.h"
#include "particles/functions/particleCreatorHelperFunctions.h"
#include "particles/functions/eulerRotation.h"


namespace olb {

namespace particles {

namespace creators {

//// SPHERE 3D

/// Set resolved sphere for existing particle but new surface
template<typename T, typename PARTICLETYPE, bool IGNORE_WARNINGS=false>
void setResolvedSphere3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, T radius, T epsilon, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorSphere3D<T, T, true>(Vector<T,3> (0.), radius, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  // Safety mechanism for wrong particle type - It is better not to use rotation matrix!
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>()
                 && !IGNORE_WARNINGS ) {
    OstreamManager clout(std::cout, "creatorSphere3D");
    clout << "WARNING: A rotation matrix is provided but is not necessary for a sphere." << std::endl;
  }

  /// Set resolved sphere 3D with given suface
  setResolvedObject<T,PARTICLETYPE,true>(particleSystem, idxParticle, idxSurface,
                                         position, density, Vector<T,3>(0.), velocity);
}

/// Add resolved sphere as new particle with new surface
template<typename T, typename PARTICLETYPE, bool IGNORE_WARNINGS=false>
void addResolvedSphere3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          const Vector<T,3>& position, T radius, T epsilon, T density=0.,
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved sphere 3D at given index
  setResolvedSphere3D<T,PARTICLETYPE,IGNORE_WARNINGS>(
      particleSystem, idxParticle, position, radius, epsilon, density, velocity );
}


//// CYLINDER 3D

/// Set resolved cylinder for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCylinder3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, const Vector<T,3>& normal,
  T height, T radius, T epsilon, T density=0.,
  const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
  const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCylinder3D<T, T, true>(Vector<T,3>(0.), normal, radius, height, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  /// Set resolved cylinder 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angleInDegree, velocity);
}

/// Add resolved cylinder as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCylinder3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            const Vector<T,3>& position, const Vector<T,3>& normal,
                            T height, T radius, T epsilon, T density=0.,
                            const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                            const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cylinder 3D at given index
  setResolvedCylinder3D( particleSystem, idxParticle, position,
                         normal, height, radius, epsilon, density, angleInDegree, velocity );
}


//// CUBOID 3D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCuboid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          std::size_t idxParticle,
                          const Vector<T,3>& position, const Vector<T,3>& extend, T epsilon, T density=0.,
                          const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCuboid3D<T, T, true>(extend[0], extend[1], extend[2], Vector<T,3>(0.), epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  /// Set resolved cuboid 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angleInDegree, velocity);
}


/// Add resolved cuboid as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCuboid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                          const Vector<T,3>& position, const Vector<T,3>& extend, T epsilon, T density=0.,
                          const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                          const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cuboid 3D at given index
  setResolvedCuboid3D( particleSystem, idxParticle, position,
                       extend, epsilon, density, angleInDegree, velocity );
}


//// CONE 3D

/// Set resolved cone for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedCone3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                        std::size_t idxParticle,
                        const Vector<T,3>& position1, const Vector<T,3>& position2,
                        T radius1, T radius2, T epsilon, T density=0,
                        const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                        const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCone3D<T, T, true>(position1, position2, radius1, radius2, epsilon ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  const Vector<T,3> position = sIndicatorPtr->calcCenterOfMass();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  /// Set resolved cone 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angleInDegree, velocity);
}


/// Add resolved cone as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedCone3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                        const Vector<T,3>& position1, const Vector<T,3>& position2,
                        T radius1, T radius2, T epsilon, T density=0,
                        const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                        const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved cone 3D at given index
  setResolvedCone3D( particleSystem, idxParticle,
                     position1, position2, radius1, radius2, epsilon, density, angleInDegree, velocity );
}


//// Ellipsoid 3D

/// Set resolved ellipsoid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedEllipsoid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                             std::size_t idxParticle,
                             const Vector<T,3>& position,
                             const Vector<T,3>& radius,
                             T epsilon, T density=0.,
                             const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                             const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorEllipsoid3D<T, T, true>(Vector<T,3>(0.), radius, epsilon, Vector<T,3>(0.) ));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  /// Set resolved ellipsoid 3D with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angleInDegree, velocity);
}


/// Add resolved ellipsoid as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedEllipsoid3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                             const Vector<T,3>& position,
                             const Vector<T,3>& radius,
                             T epsilon, T density=0.,
                             const Vector<T,3>& angleInDegree = Vector<T,3> (0.),
                             const Vector<T,3>& velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved ellipsoid 3D at given index
  setResolvedEllipsoid3D( particleSystem, idxParticle, position,
                          radius, epsilon, density, angleInDegree, velocity );
}


//// ARBITRARY SHAPE 3D

/// Set resolved cuboid for existing particle but new surface
template<typename T, typename PARTICLETYPE>
void setResolvedArbitraryShape3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  std::size_t idxParticle,
                                  Vector<T,3> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF3D<T>> indPtr,
                                  T epsilon, T density=0.,
                                  Vector<T,3> angleInDegree = Vector<T,3> (0.),
                                  Vector<T,3> velocity = Vector<T,3> (0.))
{
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T, T, true> SIndicatorBaseType;

  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCustom3D<T, T, true>(latticeSpacing, indPtr, Vector<T,3> (0.), epsilon, Vector<T,3> (0.)));

  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  std::size_t idxSurface = vectorOfIndicators.size();
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  /// Set resolved arbitrary shape with given suface
  setResolvedObject(particleSystem, idxParticle, idxSurface,
                    position, density, angleInDegree, velocity);
}

/// Add resolved arbitrary shape as new particle with new surface
template<typename T, typename PARTICLETYPE>
void addResolvedArbitraryShape3D( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                  Vector<T,3> position,
                                  T latticeSpacing,
                                  std::shared_ptr<IndicatorF3D<T>> indPtr,
                                  T epsilon, T density=0.,
                                  Vector<T,3> angleInDegree = Vector<T,3> (0.),
                                  Vector<T,3> velocity = Vector<T,3> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved arbitrary shape 3D at given index
  setResolvedArbitraryShape3D( particleSystem, idxParticle, position,
      latticeSpacing, indPtr, epsilon, density, angleInDegree, velocity );
}


//// SUBGRID

//Add subgrid particle
template<typename T, typename PARTICLETYPE>
void addSubgridSphere3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const Vector<T,3>& position, T radius, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.) )
{
  addSubgridObject( particleSystem,
    position, radius, density, velocity );
}

//// SUBGRID HAIDER LEVENSPIEL

/// Set new particle
template<typename T, typename PARTICLETYPE>
void setSubgridHaiderLevenspielObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,PARTICLETYPE::d>& position,
  T radius, T sphericity, T mass, T density,
  const Vector<T,PARTICLETYPE::d>& velocity, T lengthP = 1., T diameterP = 1.)
{
  using namespace descriptors;
  using namespace access;

  //Initialize fields (Seems to be necessary for clang-1000.10.44.4 but not for gcc)
  particles::dynamics::initializeParticle<T,PARTICLETYPE>(
    particleSystem.get().data(), idxParticle);

  //Set values
  auto particle = particleSystem.get( idxParticle );
  using namespace Haider_Levenspiel;
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<PHYSPROPERTIES,MASS>( mass );
  particle.template setField<PHYSPROPERTIES,DENSITY>( density );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<PHYSPROPERTIES,RADIUS>( radius );
  particle.template setField<NUMERICPROPERTIES, A> (util::exp(2.3288-6.4581*sphericity + 2.4486*sphericity*sphericity));
  particle.template setField<NUMERICPROPERTIES, B> (0.0964 + 0.5565*sphericity);
  particle.template setField<NUMERICPROPERTIES, C> (util::exp(4.905-13.8944*sphericity + 18.4222*sphericity*sphericity-10.2599*sphericity*sphericity*sphericity));
  particle.template setField<NUMERICPROPERTIES, D> (util::exp(1.4681+12.2584*sphericity-20.7322*sphericity*sphericity+15.8855*sphericity*sphericity*sphericity));
  if constexpr ( providesActive(particle) ){
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( true );
  }
  T beta = lengthP/diameterP;
  T a = diameterP/2.;
  Vector<T,3> scaling(beta*a,a,a);
  particle.template setField<NUMERICPROPERTIES,SCALING>( scaling); //used for visualisation in Paraview

}

/// Add subgrid object as new particle
template<typename T, typename PARTICLETYPE>
void addSubgridHaiderLevenspielObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const Vector<T,PARTICLETYPE::d>& position, T radius, T sphericity, T mass=0., T density=1.,
  const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.), T lengthP=1., T diameterP=1.)
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set subgrid object 3D at given index
  setSubgridHaiderLevenspielObject( particleSystem, idxParticle, position, radius, sphericity, mass, density, velocity, lengthP, diameterP );
}


//Add subgrid particle
template<typename T, typename PARTICLETYPE>
void addSubgridHaiderLevenspiel3D(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const Vector<T,3>& position, T surface, T volume, T density=0.,
  const Vector<T,3>& velocity = Vector<T,3> (0.), T lengthP = 1., T diameterP = 1.)
{
  //volume equivalent radius used
  T radius=util::pow(3*volume/(4*3.141592653),0.33333333333333);
  T sphericity = (4*M_PI*std::pow(3*volume/(4*M_PI),0.66666666666))/surface;
  T mass = density*volume;
  addSubgridHaiderLevenspielObject( particleSystem,
    position, radius, sphericity, mass, velocity, lengthP, diameterP );
}


//Add Haider Levenspiel particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
void addSubgridHaiderLevenspiel3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T surface, T volume, T density=1.,
  const Vector<T,3>& velocity = Vector<T,3> (0.) , T lengthP = 1., T diameterP = 1.)
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = false;
  if (auto iC = cuboidDecomposition.getC(position)) {
    inDomain = true;
    globiCcentre = *iC;
  }
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }
  //Iterate over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

  //Do only on residence iC
  if (globiC == globiCcentre){
  //Run creator on particleSystem
  T radius=util::pow(3*volume/(4*3.141592653),0.33333333333333);//volume equivalent
  T sphericity = (4*M_PI*std::pow(3*volume/(4*M_PI),0.66666666666))/surface;
  T mass = density*volume;
  addSubgridHaiderLevenspielObject( particleSystem,
  position, radius, sphericity, mass, density, velocity, lengthP, diameterP );
  //Add additional data due to parallization
  auto particle = particleSystem.get(particleSystem.size()-1);
  particle.template setField<PARALLELIZATION,ID>( globID );
  //Reinforcing default values
  particle.template setField<GENERAL,INVALID>( false );
    }
  });
}


////SUBGRID TRAN-CONG

/// Set new particle Tran Cong
template<typename T, typename PARTICLETYPE>
void setSubgridTranCongObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,PARTICLETYPE::d>& position,
  T radius, T c, T A_p, T dA_dnn, T mass, T density,
  const Vector<T,PARTICLETYPE::d>& velocity)
{
  using namespace descriptors;
  using namespace access;

  //Initialize fields (Seems to be necessary for clang-1000.10.44.4 but not for gcc)
  particles::dynamics::initializeParticle<T,PARTICLETYPE>(
    particleSystem.get().data(), idxParticle);
  //Set values
  auto particle = particleSystem.get( idxParticle );
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<PHYSPROPERTIES,MASS>( mass );
  particle.template setField<PHYSPROPERTIES,DENSITY>( density );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<PHYSPROPERTIES,RADIUS>( radius );
 particle.template setField <NUMERICPROPERTIES,CIRCULARITY> (c);
 particle.template setField <NUMERICPROPERTIES,PROJECTED_SURFACE> (A_p);
 particle.template setField <NUMERICPROPERTIES,DA_DN> (dA_dnn);

  if constexpr ( providesActive(particle) ){
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( true );
  }
}


template<typename T, typename PARTICLETYPE>
void addSubgridTranCongObject(ParticleSystem<T,PARTICLETYPE>& particleSystem,
    const Vector<T,PARTICLETYPE::d>& position, T radius, T c, T A_p , T dA_dn , T mass =0., T density=1., const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.) )
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address

  particleSystem.extend();
  setSubgridTranCongObject( particleSystem, idxParticle,
    position, radius, c, A_p , dA_dn, mass, density, velocity );
}



//Add Tran Cong particle on residing iC with surface on all iCs
template<typename T, typename PARTICLETYPE>
void addSubgridTranCongFromCylinder3D(
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T length, T diameter, T _angle = 0., T density=1.,
  const Vector<T,3>& velocity = Vector<T,3> (0.) )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = false;
  if (auto iC = cuboidDecomposition.getC(position)) {
    inDomain = true;
    globiCcentre = *iC;
  }
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }
  //Iterate over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      T r = diameter/2.;
      T angle = _angle/360.*(2*M_PI);
      T volume = M_PI*r*r*length;
      T A_p = util::abs(M_PI*diameter*diameter/4.*olb::util::cos(angle)) + util::abs(diameter*length*olb::util::sin(angle));
      T P_p = 2*M_PI*r*olb::util::cos(angle) + (2*r+length)*2*olb::util::sin(angle);
      T d_n=util::pow(6*volume/M_PI, 0.33333333333); //volume equivalent diameter
      T d_A = util::sqrt(4*A_p/M_PI);//projected surface equivalent diameter7
      T c = (M_PI*d_A)/P_p;//circularity
      T dA_dn = d_A/d_n;
      T radius = r*olb::util::cos(angle) + length/2.*olb::util::sin(angle);//calculated for the deposition function
      //T radius = r;//calculated for the deposition function


  T mass = density*volume;
  addSubgridTranCongObject( particleSystem,
    position, radius, c, A_p,  dA_dn, mass, density, velocity );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      T beta = length/diameter;
      T a = r;
        Vector<T,3> scaling(beta*a,a,a);
      particle.template setField<NUMERICPROPERTIES,SCALING>( scaling); //used for visualisation in Paraview
      //Reinforcing default values
      particle.template setField<GENERAL,INVALID>( false );
    }
  });
}


///EULER ROTATION SPHEROIDS



/// Set new particle from quaternion
template<typename T, typename PARTICLETYPE>
void setSubgridEulerRotationSpheroid( //assume the rotational axis in z-axis
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
   Vector<T,4> quaternion = Vector<T,4> (0.,0.,0.,1.),
  Vector<T,3> velocity = Vector<T,3> (0.,0.,0.),
Vector<T,3> ang_velocity = Vector<T,3> (0.,0.,0.) /*arround principle axis*/, int label = 0  )
{
  using namespace descriptors;
  using namespace access;


  //Initialize fields (Seems to be necessary for clang-1000.10.44.4 but not for gcc)
  particles::dynamics::initializeParticle<T,PARTICLETYPE>(
    particleSystem.get().data(), idxParticle);
  //Set values
  T beta = length/diameter;
  T a = diameter/2.;
  T mass = 4./3.*M_PI*a*a*a*beta*density;//spheroid
  //T mass = 2*M_PI*a*a*a*beta*density;//cylinder
  Vector<T,3> mom_of_in (((1+beta*beta)*a*a)/5.*mass,((1+beta*beta)*a*a)/5.*mass, 2.*a*a/5.*mass);//spheroid
   //Vector<T,3> mom_of_in (((1./4.+1./6.*beta*beta)*a*a)*mass,((1./4.+1./6.*beta*beta)*a*a)*mass, a*a/2.*mass);//cylinder

  Vector<T,3> abg0 (beta*beta/(beta*beta-1)+beta/(2*util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))),
  beta*beta/(beta*beta-1)+beta/(2*util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))),
 (-2)/(beta*beta-1)-beta/(util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))));

  auto particle = particleSystem.get( idxParticle );
  particle.template setField<NUMERICPROPERTIES,ABG0> (abg0);
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<PHYSPROPERTIES,MASS>( mass );
  particle.template setField<NUMERICPROPERTIES,MOMENT_OF_INERTIA>( mom_of_in );
  particle.template setField<PHYSPROPERTIES,DENSITY>( density );
  eler::qnormalize(quaternion);
  particle.template setField<NUMERICPROPERTIES,ROT_MATRIX>(eler::computeTransformMatrixFromEulQuat<T>(quaternion));
  particle.template setField<NUMERICPROPERTIES,QUATERNION>( quaternion );
  particle.template setField<NUMERICPROPERTIES,PARAQUATERNION> (eler::getParaquaternion(particle.template getField<NUMERICPROPERTIES,QUATERNION>()));
  particle.template setField<NUMERICPROPERTIES,ANG_VELOCITY> (ang_velocity);
  particle.template setField<NUMERICPROPERTIES, ORIENTATION> (eler::rotate(eler::inverseTransformMatrix(particle.template getField<NUMERICPROPERTIES,ROT_MATRIX>()), Vector<T,3>(0.,0.,1.)));
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<PHYSPROPERTIES,RADIUS>( a );//equal to semi_minor_axis
  Vector<T,3> scaling(beta*a,a,a);
  particle.template setField<NUMERICPROPERTIES,SCALING>( scaling); //used for visualisation in Paraview
  Vector<T,9> dyadic(0.);
  dyadic[0] = (16.*(beta*beta-1))/((2*beta*beta-3)*olb::util::log(beta+olb::util::sqrt(beta*beta-1))/(olb::util::sqrt(beta*beta-1))+beta);
  dyadic[4] = dyadic[0];
  dyadic[8] = (8.*(beta*beta-1))/((2*beta*beta-1)*olb::util::log(beta+olb::util::sqrt(beta*beta-1))/(olb::util::sqrt(beta*beta-1))-beta);
  particle.template setField<NUMERICPROPERTIES,TRANSLATIONAL_DYADIC>( dyadic );
  particle.template setField<NUMERICPROPERTIES,BETA> (beta);
  particle.template setField<NUMERICPROPERTIES,LABEL> (label);

  if constexpr ( providesActive(particle) ){
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( true );
  }
}


template<typename T, typename PARTICLETYPE>
void setSubgridEulerRotationSpheroid( //assume the rotational axis in z-axis
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
   Vector<T,3> eul_ang = Vector<T,3> (0.,0.,0.),
  Vector<T,3> velocity = Vector<T,3> (0.,0.,0.),
Vector<T,3> ang_velocity = Vector<T,3> (0.,0.,0.) /*arround principle axis*/, int label = 0  )
{
  using namespace descriptors;
  using namespace access;

  //Initialize fields (Seems to be necessary for clang-1000.10.44.4 but not for gcc)
  particles::dynamics::initializeParticle<T,PARTICLETYPE>(
    particleSystem.get().data(), idxParticle);
  //Set values
  T beta = length/diameter;
  T a = diameter/2.;
  T mass = 4./3.*M_PI*a*a*a*beta*density;//
  //T mass = 2*M_PI*a*a*a*beta*density;//cylinder
  Vector<T,3> mom_of_in (((1+beta*beta)*a*a)/5.*mass,((1+beta*beta)*a*a)/5.*mass, 2.*a*a/5.*mass);
   //Vector<T,3> mom_of_in (((1./4.+1./6.*beta*beta)*a*a)*mass,((1./4.+1./6.*beta*beta)*a*a)*mass, a*a/2.*mass);//cylinder

  Vector<T,3> abg0 (beta*beta/(beta*beta-1)+beta/(2*util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))),
  beta*beta/(beta*beta-1)+beta/(2*util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))),
 (-2)/(beta*beta-1)-beta/(util::pow(beta*beta-1,1.5))*util::log((beta- util::sqrt(beta*beta-1))/(beta+util::sqrt(beta*beta-1))));

  auto particle = particleSystem.get( idxParticle );
    particle.template setField<NUMERICPROPERTIES,ABG0> (abg0);
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<PHYSPROPERTIES,MASS>( mass );
  particle.template setField<NUMERICPROPERTIES,MOMENT_OF_INERTIA>( mom_of_in );
  particle.template setField<PHYSPROPERTIES,DENSITY>( density );
  particle.template setField<NUMERICPROPERTIES,ROT_MATRIX>(eler::computeTransformMatrixFromEulAng<T>(eul_ang));
  particle.template setField<NUMERICPROPERTIES,QUATERNION>( eler::computeEulQuatFromEulAng<T>(eul_ang) );
  particle.template setField<NUMERICPROPERTIES,PARAQUATERNION> (eler::getParaquaternion(particle.template getField<NUMERICPROPERTIES,QUATERNION>()));
  particle.template setField<NUMERICPROPERTIES,ANG_VELOCITY> (ang_velocity);
  particle.template setField<NUMERICPROPERTIES, ORIENTATION> (eler::rotate(eler::inverseTransformMatrix(particle.template getField<NUMERICPROPERTIES,ROT_MATRIX>()), Vector<T,3>(0.,0.,1.)));
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<PHYSPROPERTIES,RADIUS>( a );//equal to semi_minor_axis
  Vector<T,3> scaling(beta*a,a,a);//semi axis saved, for whole dimensions multiply 2x
  particle.template setField<NUMERICPROPERTIES,SCALING>( scaling); //used for visualisation in Paraview
  Vector<T,9> dyadic(0.);
  dyadic[0] = (16.*(beta*beta-1))/((2*beta*beta-3)*olb::util::log(beta+olb::util::sqrt(beta*beta-1))/(olb::util::sqrt(beta*beta-1))+beta);
  dyadic[4] = dyadic[0];
  dyadic[8] = (8.*(beta*beta-1))/((2*beta*beta-1)*olb::util::log(beta+olb::util::sqrt(beta*beta-1))/(olb::util::sqrt(beta*beta-1))-beta);
  particle.template setField<NUMERICPROPERTIES,TRANSLATIONAL_DYADIC>( dyadic );
  particle.template setField<NUMERICPROPERTIES,BETA> (beta);
  particle.template setField<NUMERICPROPERTIES,LABEL> (label);

  if constexpr ( providesActive(particle) ){
    particle.template setField<DYNBEHAVIOUR,ACTIVE>( true );
  }
}


template<typename T, typename PARTICLETYPE>
void addSubgridEulerRotationSpheroid( //assume the rotational axis in z-axis
 ParticleSystem<T,PARTICLETYPE>& ParticleSystem,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
  Vector<T,4> quaternion = Vector<T,4> (0.,0.,0.,1.),
   Vector<T,3> velocity = Vector<T,3> (0.,0.,0.),
 Vector<T,3> ang_velocity = Vector<T,3> (0.,0.,0.) /*arround principle axis*/,
 int label = 0)
{
  //Retrieve new index
  std::size_t idxParticle = ParticleSystem.size();

  //Initialize particle address

  ParticleSystem.extend();
  setSubgridEulerRotationSpheroid( ParticleSystem, idxParticle,
    position, length, diameter, density , quaternion, velocity, ang_velocity, label );
}

template<typename T, typename PARTICLETYPE>
void addSubgridEulerRotationSpheroid( //assume the rotational axis in z-axis
 ParticleSystem<T,PARTICLETYPE>& ParticleSystem,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
  Vector<T,3> eul_ang = Vector<T,3> (0.,0.,0.),
   Vector<T,3> velocity = Vector<T,3> (0.,0.,0.),
 Vector<T,3> ang_velocity = Vector<T,3> (0.,0.,0.) /*arround principle axis*/,
 int label = 0)
{
  //Retrieve new index
  std::size_t idxParticle = ParticleSystem.size();

  //Initialize particle address

  ParticleSystem.extend();
  setSubgridEulerRotationSpheroid( ParticleSystem, idxParticle,
    position, length, diameter, density , eul_ang, velocity, ang_velocity, label );
}



template<typename T, typename PARTICLETYPE>
void addSubgridEulerRotationSpheroid3D( //assume the rotational axis in z-axis
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
  const Vector<T,3> velocity = Vector<T,3> (0.),
  const Vector<T,3> eul_ang = Vector<T,3> (0.,0.,0.),
const Vector<T,3> ang_velocity = Vector<T,3> (0.) /*arround principle axis*/, int label = 0 )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = false;
  if (auto iC = cuboidDecomposition.getC(position)) {
    inDomain = true;
    globiCcentre = *iC;
  }
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }
  //Iterate over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

  //Do only on residence iC
  if (globiC == globiCcentre){
    //Run creator on particleSystem

  addSubgridEulerRotationSpheroid( particleSystem,
    position, length, diameter, density, eul_ang, velocity, ang_velocity, label );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Reinforcing default values
      particle.template setField<GENERAL,INVALID>( false );
    }
  });
}


template<typename T, typename PARTICLETYPE>
void addSubgridEulerRotationSpheroid3D( //assume the rotational axis in z-axis
  SuperParticleSystem<T,PARTICLETYPE>& sParticleSystem,
  const Vector<T,3>& position, T length, T diameter, T density=1.,
  const Vector<T,3> velocity = Vector<T,3> (0.),
  const Vector<T,4> eul_quat = Vector<T,4> (0.,0.,0.,1.),
const Vector<T,3> ang_velocity = Vector<T,3> (0.) /*arround principle axis*/, int label = 0 )
{
  using namespace descriptors;
  //Retrieving cuboid geometry
  auto& cuboidDecomposition= sParticleSystem.getSuperStructure().getCuboidDecomposition();
  //Get global particle id
  std::size_t globID = sParticleSystem.getGlobID();
  //Find globiC the particle belongs to
  int globiCcentre;
  bool inDomain = false;
  if (auto iC = cuboidDecomposition.getC(position)) {
    inDomain = true;
    globiCcentre = *iC;
  }
  if (!inDomain){ std::cerr << "ERROR: Particle added outside domain!" << std::endl; }
  //Iterate over particle systems
  communication::forSystemsInSuperParticleSystem( sParticleSystem,
  [&](ParticleSystem<T,PARTICLETYPE>& particleSystem, int iC, int globiC){

    //Do only on residence iC
    if (globiC == globiCcentre){
      //Run creator on particleSystem
      T a = diameter/2.;//semi_minor_axis

  addSubgridEulerRotationSpheroid( particleSystem,
    position, length, diameter, density, eul_quat, velocity, ang_velocity, label );
      //Add additional data due to parallization
      auto particle = particleSystem.get(particleSystem.size()-1);
      particle.template setField<PARALLELIZATION,ID>( globID );
      //Reinforcing default values
      particle.template setField<GENERAL,INVALID>( false );
    }
  });
}












} //namespace creators

} //namespace particles

} //namespace olb

#endif
