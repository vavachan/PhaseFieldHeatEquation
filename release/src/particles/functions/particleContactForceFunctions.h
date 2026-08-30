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

#ifndef PARTICLE_CONTACT_FORCE_FUNCTIONS_H
#define PARTICLE_CONTACT_FORCE_FUNCTIONS_H

#include "core/vector.h"
#include "functors/analytical/indicator/indicatorBase.h"
#include "particles/particles.h"

namespace olb {

namespace particles {

namespace defaults {
template <typename T, unsigned D>
const auto processContactForce =
    [](const PhysR<T, D>& contactPosition, const Vector<T, D>& contactNormal,
       const Vector<T, D>& normalForce, const Vector<T, D>& tangentialForce,
       const std::array<std::size_t, 2>& ids, T overlapVolume, T indentation,
       bool isParticleWallContact) {};
}

namespace contact {

template <typename T, unsigned D, typename CONTACTTYPE>
constexpr T
evalCurrentDampingFactor(CONTACTTYPE& contact, const T coefficientOfRestitution,
                         const Vector<T, D>& initialRelativeNormalVelocity)
{
  if (contact.getDampingFactor() < 0) {
    contact.setDampingFactorFromInitialVelocity(
        coefficientOfRestitution, norm(initialRelativeNormalVelocity));
  }
  return contact.getDampingFactor();
}

/// Returns the relative velocity of two particles at a defined position
template <typename T, typename PARTICLETYPE>
constexpr Vector<T, PARTICLETYPE::d>
evalRelativeVelocity(Particle<T, PARTICLETYPE>&        particleA,
                     Particle<T, PARTICLETYPE>&        particleB,
                     const Vector<T, PARTICLETYPE::d>& position)
{
  return particles::dynamics::calculateLocalVelocity(particleA, position) -
         particles::dynamics::calculateLocalVelocity(particleB, position);
}

template <typename T, typename PARTICLETYPE>
constexpr Vector<T, PARTICLETYPE::d>
evalRelativeVelocity(Particle<T, PARTICLETYPE>&        particle,
                     const Vector<T, PARTICLETYPE::d>& position)
{
  return particles::dynamics::calculateLocalVelocity(particle, position);
}

template <typename T, unsigned D>
constexpr Vector<T, D>
evalRelativeNormalVelocity(const Vector<T, D>& contactNormal,
                           const Vector<T, D>& relVel)
{
  return Vector<T, D>((relVel * contactNormal) * contactNormal);
}

template <typename T, typename F>
constexpr void forEachPosition(Vector<T, 3> startPos, Vector<T, 3> endPos, F f)
{
  Vector<T, 3> pos;
  for (pos[0] = startPos[0]; pos[0] <= endPos[0]; ++pos[0]) {
    for (pos[1] = startPos[1]; pos[1] <= endPos[1]; ++pos[1]) {
      for (pos[2] = startPos[2]; pos[2] <= endPos[2]; ++pos[2]) {
        f(pos);
      }
    }
  }
}

template <typename T, typename F>
constexpr void forEachPosition(Vector<T, 2> startPos, Vector<T, 2> endPos, F f)
{
  Vector<T, 2> pos;
  for (pos[0] = startPos[0]; pos[0] <= endPos[0]; ++pos[0]) {
    for (pos[1] = startPos[1]; pos[1] <= endPos[1]; ++pos[1]) {
      f(pos);
    }
  }
}

template <typename T, typename F>
constexpr void forEachPositionWithBreak(Vector<T, 3> startPos,
                                        Vector<T, 3> endPos, F f)
{
  bool         breakLoops = false;
  Vector<T, 3> pos;
  for (pos[0] = startPos[0]; pos[0] <= endPos[0]; ++pos[0]) {
    for (pos[1] = startPos[1]; pos[1] <= endPos[1]; ++pos[1]) {
      for (pos[2] = startPos[2]; pos[2] <= endPos[2]; ++pos[2]) {
        f(pos, breakLoops);
        if (breakLoops) {
          return;
        }
      }
    }
  }
}

template <typename T, typename F>
constexpr void forEachPositionWithBreak(Vector<T, 2> startPos,
                                        Vector<T, 2> endPos, F f)
{
  bool         breakLoops = false;
  Vector<T, 2> pos;
  for (pos[0] = startPos[0]; pos[0] <= endPos[0]; ++pos[0]) {
    for (pos[1] = startPos[1]; pos[1] <= endPos[1]; ++pos[1]) {
      f(pos, breakLoops);
      if (breakLoops) {
        return;
      }
    }
  }
}
template <typename T, unsigned D>
constexpr Vector<T, D>
evalContactDeltaX(const Vector<T, D>& min, const Vector<T, D>& max,
                  const T        physDeltaX,
                  const unsigned contactBoxResolutionPerDirection)
{
  // Adaptiv resolution of found contact
  const Vector<T, D>  boxSize           = max - min;
  Vector<T, D>        contactPhysDeltaX = boxSize;
  std::array<bool, D> fallbackUsed {[]() {
    static_assert(D == 2 || D == 3, "Only D=2 and D=3 are supported");
    if constexpr (D == 3) {
      return std::array<bool, D> {false, false, false};
    }
    else {
      return std::array<bool, D> {false, false};
    }
  }()};

  for (unsigned iD = 0; iD < D; ++iD) {
    contactPhysDeltaX[iD] /= contactBoxResolutionPerDirection;
    if (util::nearZero(contactPhysDeltaX[iD])) {
      contactPhysDeltaX[iD] = physDeltaX / contactBoxResolutionPerDirection;
      fallbackUsed[iD]      = true;
    }
  }

  // If fallback was used, check if in another dimension the deltaX is known.
  // If yes, then use the smallest known deltaX. Otherwise, keep the fallback.
  // This is possible because any dimenison for which we cannot find a proper deltaX has probably even smaller lenghts than the others.
  for (unsigned iD = 0; iD < D; ++iD) {
    if (fallbackUsed[iD]) {
      for (unsigned iD2 = 0; iD2 < D; ++iD2) {
        if (iD != iD2) {
          contactPhysDeltaX[iD] =
              util::min(contactPhysDeltaX[iD], contactPhysDeltaX[iD2]);
        }
      }
    }
  }

  return contactPhysDeltaX;
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, bool CONVEX>
void prepareForceCalculation(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer)
{

  communicateContacts(contactContainer);

  if constexpr (!CONVEX) {
    OstreamManager clout(std::cout, "forceCalcPrep");
    // TODO: Split contacts if the points are not connected
    // How to remember properties then? Store old min and max (maybe with a little extra) and then check for to which contact a position belongs?
    clout << "WARNING: Concave particles aren't supported." << std::endl;
  }
}

/// Elastic and viscous force from overlap volume according to Nassauer & Kuna (2013)
template <typename T, unsigned D>
Vector<T, D> calcElasticAndViscousForce(T effectiveE, T k, T overlapVolume,
                                        T indentation, T dampingFactor,
                                        const Vector<T, D>& relVelNormal,
                                        const Vector<T, D>& contactNormal)
{
  // Sum is necessary to get the sign of the damping since velocity and normal have the exact same/opposite direction
  T sum = relVelNormal[0] * contactNormal[0];
  for (unsigned iD = 1; iD < D; ++iD) {
    sum += relVelNormal[iD] * contactNormal[iD];
  }

  const T forceMagnitude =
      effectiveE * k * util::sqrt(overlapVolume * indentation) *
      (1 + dampingFactor * util::sign(-sum) * norm(relVelNormal));
  return util::max(forceMagnitude, T {0}) * contactNormal;
}

/// Takes a normalized normal and the voxel sizes and approximates the corresponding surface
template <typename T, unsigned D>
T evalApproxSurface(const Vector<T, D>& normalizedNormal,
                    const PhysR<T, D>&  contactPhysDeltaX)
{
  T approxSurface;
  if constexpr (D == 3) {
    approxSurface =
        contactPhysDeltaX[1] * contactPhysDeltaX[2] *
        /* util::abs(normalizedNormal[0]); // * normalizedNormal[0]; */
        normalizedNormal[0] * normalizedNormal[0];
    approxSurface +=
        contactPhysDeltaX[0] * contactPhysDeltaX[2] *
        /* util::abs(normalizedNormal[1]); // * normalizedNormal[1]; */
        normalizedNormal[1] * normalizedNormal[1];
    approxSurface +=
        contactPhysDeltaX[0] * contactPhysDeltaX[1] *
        /* util::abs(normalizedNormal[2]); // * normalizedNormal[2]; */
        normalizedNormal[2] * normalizedNormal[2];
  }
  else {
    approxSurface =
        contactPhysDeltaX[0] *
        /* util::abs(normalizedNormal[1]); // * normalizedNormal[1]; */
        normalizedNormal[1] * normalizedNormal[1];
    approxSurface +=
        contactPhysDeltaX[1] *
        /* util::abs(normalizedNormal[0]); // * normalizedNormal[0]; */
        normalizedNormal[0] * normalizedNormal[0];
  }
  return approxSurface;
  /* return approxSurface / norm(contactPhysDeltaX); */
}

template <typename T, unsigned D>
Vector<T, D> getTangentialForceFromOverlapVolume(
    const Vector<T, D>& contactNormal, const T coefficientStaticFriction,
    const T coefficientKineticFriction, const T staticKineticTransitionVelocity,
    const T k, const Vector<T, D>& relVel, const Vector<T, D>& normalForce)
{
  // Only regard if static friction is > 0, because we obtain nan otherwise
  if (coefficientStaticFriction > T {0.}) {
    const Vector<T, D> relVelNormal =
        evalRelativeNormalVelocity(contactNormal, relVel);
    const Vector<T, D> relVelTangential(relVel - relVelNormal);
    T const            magnitudeTangentialVelocity = norm(relVelTangential);
    if (magnitudeTangentialVelocity != T {0}) {
      T const velocityQuotient =
          magnitudeTangentialVelocity / staticKineticTransitionVelocity;

      T magnitudeTangentialForce =
          2 * coefficientStaticFriction *
          (1 - T {0.09} * util::pow(coefficientKineticFriction /
                                        coefficientStaticFriction,
                                    4));
      magnitudeTangentialForce -= coefficientKineticFriction;
      magnitudeTangentialForce *= velocityQuotient * velocityQuotient /
                                  (util::pow(velocityQuotient, 4) + 1);
      magnitudeTangentialForce += coefficientKineticFriction;
      magnitudeTangentialForce -= coefficientKineticFriction /
                                  (velocityQuotient * velocityQuotient + 1);
      magnitudeTangentialForce *= norm(normalForce);

      return -magnitudeTangentialForce * relVelTangential /
             magnitudeTangentialVelocity;
    }
  }

  return Vector<T, D>(T {0.});
}

template <typename T, unsigned D>
Vector<T, D>
getNormalForceFromOverlapVolume(const Vector<T, D>& contactNormal,
                                const T indentation, const T overlapVolume,
                                const T effectiveE, const T dampingFactor,
                                const T k, const Vector<T, D>& relVel)
{
  const Vector<T, D> relVelNormal =
      evalRelativeNormalVelocity(contactNormal, relVel);
  return calcElasticAndViscousForce(effectiveE, k, overlapVolume, indentation,
                                    dampingFactor, relVelNormal, contactNormal);
}

template <typename T, typename PARTICLETYPE>
void applyForceOnParticle(
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
    Particle<T, PARTICLETYPE>&                           particle,
    XParticleSystem<T, PARTICLETYPE>&                    particleSystem,
    const Vector<T, PARTICLETYPE::d>&                    contactPoint,
    const Vector<T, PARTICLETYPE::d>&                    normalForce,
    const Vector<T, PARTICLETYPE::d>&                    tangentialForce)
{
  using namespace descriptors;
  constexpr unsigned D            = PARTICLETYPE::d;
  constexpr unsigned Drot         = utilities::dimensions::convert<D>::rotation;
  Vector<T, D>       contactForce = normalForce + tangentialForce;
  const PhysR<T, D>  position = particle.template getField<GENERAL, POSITION>();
  const PhysR<T, D>  contactPointDistance = contactPoint - position;
  Vector<T, Drot>    torqueFromContact(
      particles::dynamics::torque_from_force<D, T>::calculate(
          contactForce, contactPointDistance));

  if (access::isContactComputationEnabled(particle)) {
    bool applyNow = true;

    if constexpr (access::providesParallelization<PARTICLETYPE>()) {
      auto& loadBalancer = particleSystem.getSuperStructure().getLoadBalancer();
      const std::size_t globalParticleIC =
          particle.template getField<PARALLELIZATION, IC>();
      int responsibleRank = loadBalancer.rank(globalParticleIC);
      if (responsibleRank != singleton::mpi().getRank()) {
        applyNow = false;
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();
        extendContactTreatmentResultsDataMap(globalParticleID, contactForce,
                                             torqueFromContact, responsibleRank,
                                             dataMap);
      }
    }

    if (applyNow) {
      if constexpr (access::providesParallelization<PARTICLETYPE>()) {
        // Apply force only to the main particle
        const std::size_t globalParticleID =
            particle.template getField<PARALLELIZATION, ID>();
        particles::communication::forParticlesInSuperParticleSystem<
            T, PARTICLETYPE, conditions::valid_particles>(
            particleSystem,
            [&](Particle<T, PARTICLETYPE>&       particle,
                ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
              const int globalParticleIC =
                  particle.template getField<PARALLELIZATION, IC>();
              const std::size_t currGlobalParticleID =
                  particle.template getField<PARALLELIZATION, ID>();
              if (currGlobalParticleID == globalParticleID &&
                  globiC == globalParticleIC) {
                // Applying contact force on the particle
                const Vector<T, D> force =
                    particle.template getField<FORCING, FORCE>();
                particle.template setField<FORCING, FORCE>(force +
                                                           contactForce);

                // Torque calculation and application
                const Vector<T, Drot> torque(
                    particle.template getField<FORCING, TORQUE>());
                particle.template setField<FORCING, TORQUE>(
                    utilities::dimensions::convert<D>::serialize_rotation(
                        torque + torqueFromContact));
#ifdef VERBOSE_CONTACT_COMMUNICATION
                std::cout << "Rank " << singleton::mpi().getRank()
                          << " procsses data of particle " << globalParticleID
                          << " (contactForce = " << contactForce + force
                          << "; torque = " << torqueFromContact + torque << ")"
                          << std::endl;
#endif
              }
            });
      }
      else {
        // Applying contact force on the particle
        const Vector<T, D> force = particle.template getField<FORCING, FORCE>();
        particle.template setField<FORCING, FORCE>(force + contactForce);

        // Torque calculation and application
        const Vector<T, Drot> torque(
            particle.template getField<FORCING, TORQUE>());
        particle.template setField<FORCING, TORQUE>(
            utilities::dimensions::convert<D>::serialize_rotation(
                torque + torqueFromContact));
      }
    }
  }
}

template <
    typename T, typename PARTICLETYPE, typename CONTACTPROPERTIES,
    typename F = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
void applyForceFromOverlapVolume(
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
    Particle<T, PARTICLETYPE>& particleA, Particle<T, PARTICLETYPE>& particleB,
    XParticleSystem<T, PARTICLETYPE>& particleSystem,
    const PhysR<T, PARTICLETYPE::d>&  contactPoint,
    const Vector<T, PARTICLETYPE::d>& contactNormal, const T indentation,
    const T overlapVolume, const Vector<T, PARTICLETYPE::d>& relVel,
    const T dampingFactor, const CONTACTPROPERTIES& contactProperties, T k,
    const std::size_t& particleAID, const std::size_t& particleBID,
    F processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;

  const unsigned materialA =
      particleA.template getField<MECHPROPERTIES, MATERIAL>();
  const unsigned materialB =
      particleB.template getField<MECHPROPERTIES, MATERIAL>();

  // Calculation of contact force
  const Vector<T, D> normalForce = getNormalForceFromOverlapVolume(
      contactNormal, indentation, overlapVolume,
      contactProperties.getEffectiveYoungsModulus(materialA, materialB),
      dampingFactor, k, relVel);
  const Vector<T, D> tangentialForce = getTangentialForceFromOverlapVolume(
      contactNormal,
      contactProperties.getStaticFrictionCoefficient(materialA, materialB),
      contactProperties.getKineticFrictionCoefficient(materialA, materialB),
      contactProperties.getStaticFrictionCoefficient(materialA, materialB), k,
      relVel, normalForce);

  processContactForce(contactPoint, contactNormal, normalForce, tangentialForce,
                      std::array<std::size_t, 2>({particleAID, particleBID}),
                      overlapVolume, indentation, false);

  // Applying forces on the particle
  applyForceOnParticle(dataMap, particleA, particleSystem, contactPoint,
                       normalForce, tangentialForce);
  applyForceOnParticle(dataMap, particleB, particleSystem, contactPoint,
                       -1 * normalForce, -1 * tangentialForce);
}

template <
    typename T, typename PARTICLETYPE, typename CONTACTPROPERTIES,
    typename F = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
void applyForceFromOverlapVolume(
    std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
    Particle<T, PARTICLETYPE>&                           particle,
    XParticleSystem<T, PARTICLETYPE>&                    particleSystem,
    SolidBoundary<T, PARTICLETYPE::d>&                   solidBoundary,
    const PhysR<T, PARTICLETYPE::d>&                     contactPoint,
    const Vector<T, PARTICLETYPE::d>& contactNormal, const T indentation,
    const T overlapVolume, const Vector<T, PARTICLETYPE::d>& relVel,
    const T dampingFactor, const CONTACTPROPERTIES& contactProperties, T k,
    const std::size_t& particleID, const std::size_t& wallID,
    F processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;

  const unsigned particleMaterial =
      particle.template getField<MECHPROPERTIES, MATERIAL>();
  const unsigned wallMaterial = solidBoundary.getContactMaterial();

  // Calculation of contact force
  const Vector<T, D> normalForce = getNormalForceFromOverlapVolume(
      contactNormal, indentation, overlapVolume,
      contactProperties.getEffectiveYoungsModulus(particleMaterial,
                                                  wallMaterial),
      dampingFactor, k, relVel);
  const Vector<T, D> tangentialForce = getTangentialForceFromOverlapVolume(
      contactNormal,
      contactProperties.getStaticFrictionCoefficient(particleMaterial,
                                                     wallMaterial),
      contactProperties.getKineticFrictionCoefficient(particleMaterial,
                                                      wallMaterial),
      contactProperties.getStaticKineticTransitionVelocity(particleMaterial,
                                                           wallMaterial),
      k, relVel, normalForce);

  processContactForce(contactPoint, contactNormal, normalForce, tangentialForce,
                      std::array<std::size_t, 2>({particleID, wallID}),
                      overlapVolume, indentation, true);

  // Applying contact forces on the particle
  applyForceOnParticle(dataMap, particle, particleSystem, contactPoint,
                       normalForce, tangentialForce);
}

template <typename T, unsigned D>
void evalStartAndEndPos(Vector<long, D>& startPos, Vector<long, D>& endPos,
                        const PhysR<T, D>& max, const PhysR<T, D>& min,
                        const PhysR<T, D>& contactPhysDeltaX)
{
  for (unsigned iD = 0; iD < D; ++iD) {
    startPos[iD] = util::floor(min[iD] / contactPhysDeltaX[iD]);
    endPos[iD]   = util::ceil(max[iD] / contactPhysDeltaX[iD]);
  }
}

template <typename T, unsigned D>
T evalForceViaVolumeCoefficient(const Vector<T, D>& min,
                                const Vector<T, D>& max, T overlapVolume)
{
  const T a = util::log(2. / 3.) / util::log(2.);
  const T m = 2. / util::sqrt(M_PI) * util::pow(4. / M_PI, -a);

  Vector<T, D> length            = max - min;
  T            boundingBoxVolume = 1;
  for (unsigned iD = 0; iD < D; ++iD) {
    boundingBoxVolume *= length[iD];
  }

  const T k =
      m * util::pow(util::max(boundingBoxVolume / overlapVolume, T {1}), a);
  return k;
}

template <typename T>
constexpr T evalForceViaVolumeCoefficient(const T contactVolume,
                                          const T contactSurface,
                                          const T indentationDepth)
{
  constexpr T a = 0.6118229;
  constexpr T b = 1.900093;
  constexpr T c = 0.651237;

  const T x = indentationDepth * contactSurface / contactVolume;
  const T k = a + b * util::exp(-c * x);
  return k;
}

template <typename T, unsigned D, typename F1, typename F2, typename F3,
          typename F4, typename F5>
void processCell(const PhysR<T, D>& pos, const PhysR<T, D>& contactPhysDeltaX,
                 PhysR<T, D>& center, Vector<T, D>& contactNormal,
                 unsigned& volumeCellCount, T& surfaceCellCount, F1 isInside,
                 F2 isOnSurface, F3 calcNormalA, F4 calcNormalB,
                 F5 updateMinMax)
{
  // outside for normal calculation
  if (!isInside(pos)) {
    bool onSurfaceA = false, onSurfaceB = false;
    // calc normals
    Vector<T, D> normalA, normalB;
    // TODO: Improve mesh size
    T const meshSize = 1e-8 * util::average(contactPhysDeltaX);
    normalA          = calcNormalA(pos, meshSize);
    normalB          = calcNormalB(pos, meshSize);
    if (!util::nearZero(norm(normalA))) {
      normalA = normalize(normalA);
    }
    if (!util::nearZero(norm(normalB))) {
      normalB = normalize(normalB);
    }

    if (isOnSurface(pos, contactPhysDeltaX, onSurfaceA, onSurfaceB, normalA,
                    normalB)) {
      T approxSurface;
      if (!util::nearZero(norm(normalA)) && onSurfaceA) {
        approxSurface = evalApproxSurface(normalA, contactPhysDeltaX);
        surfaceCellCount += approxSurface;
        contactNormal += normalA * approxSurface;
      }
      if (!util::nearZero(norm(normalB)) && onSurfaceB) {
        approxSurface = evalApproxSurface(normalB, contactPhysDeltaX);
        surfaceCellCount += approxSurface;
        contactNormal -= normalB * approxSurface;
      }
    }
  }
  // inside volume calculation
  else {
    ++volumeCellCount;

    for (unsigned iD = 0; iD < D; ++iD) {
      center[iD] += pos[iD];
    }
    // Updating min and max for next sub-timestep
    updateMinMax(pos);
  }
}

template <typename T, unsigned D, typename F1, typename F2, typename F3,
          typename F4>
void processContactViaVolume(const Vector<T, D>& min, const Vector<T, D>& max,
                             T        physDeltaX,
                             unsigned contactBoxResolutionPerDirection,
                             F1 resetContactBox, F2 processCellWrapped,
                             F3 calculateIndentation, F4 applyForce)
{
  Vector<T, D> contactPhysDeltaX =
      evalContactDeltaX(min, max, physDeltaX, contactBoxResolutionPerDirection);

  // direction of force
  Vector<T, D> contactNormal(T(0.));
  // contact point
  Vector<T, D> center(T(0.));
  // number of cells in contact volume
  unsigned volumeCellCount = 0;
  // number of cells in contact surface
  T surfaceCellCount = 0.;
  // volume of cell describing the contact
  T infinitesimalVolume = contactPhysDeltaX[0] * contactPhysDeltaX[1];
  if constexpr (D == 3) {
    infinitesimalVolume *= contactPhysDeltaX[2];
  }

  Vector<long, D> startPos;
  Vector<long, D> endPos;
  evalStartAndEndPos(startPos, endPos, max, min, contactPhysDeltaX);
  // Increase the box 2*deltaX in each direction
  startPos -= Vector<long, D>(1);
  endPos += Vector<long, D>(1);

  resetContactBox();

  forEachPosition(startPos, endPos, [&](const Vector<long, D>& pos) {
    Vector<T, D> physPos;
    for (unsigned iD = 0; iD < D; ++iD) {
      physPos[iD] = pos[iD] * contactPhysDeltaX[iD];
    }
    processCellWrapped(contactNormal, center, volumeCellCount, surfaceCellCount,
                       physPos, contactPhysDeltaX);
  });

  if (!util::nearZero(norm(contactNormal))) {
    if (volumeCellCount > 0 && surfaceCellCount > 0) {
      // finishing center calculation
      center /= volumeCellCount;
      // finishing normal calculation
      contactNormal = normalize(contactNormal);
      // calculate precision depending on the size of the overlap
      T precision {0};
      for (unsigned iD = 0; iD < D; ++iD) {
        precision += util::pow(contactNormal[iD] * contactPhysDeltaX[iD], 2);
      }
      precision = T {1e-2} * util::sqrt(precision / D);
      const T indentation =
          calculateIndentation(center, contactNormal, precision);
      // currently, the normal points outside the particle (because we use SDF), but the force acts in the opposite way
      contactNormal *= -1;
      // Avoid indentation < 0
      if (indentation > 0) {
        applyForce(center, contactNormal, util::max(indentation, T {0}),
                   volumeCellCount * infinitesimalVolume);
      }
    }
  }
#ifdef OLB_DEBUG
  else {
    OstreamManager clout(std::cout, "contactProcessing");
    clout << "WARNING: The contact normal's norm is zero. Therefore, the "
             "detected contact is ignored."
          << std::endl;
  }
#endif
}

template <typename T, unsigned D, typename F1>
bool evalStartingPoint(const Vector<long, D>& startPos,
                       const Vector<long, D>& endPos,
                       const Vector<T, D>&    contactPhysDeltaX,
                       Vector<long, D>& currPos, F1 isInsideContact)
{
  bool pointFound = false;

  forEachPositionWithBreak(startPos, endPos,
                           [&](const Vector<long, D>& pos, bool& breakLoops) {
                             PhysR<T, D> physPos;
                             for (unsigned iD = 0; iD < D; ++iD) {
                               physPos[iD] = contactPhysDeltaX[iD] * pos[iD];
                             }
                             if (isInsideContact(physPos)) {
                               currPos    = pos;
                               pointFound = true;
                               breakLoops = true;
                             }
                           });

  return pointFound;
}

template <typename T, unsigned D, typename F1, typename F2>
void evalBoundingBoxRecursive(const Vector<long, D>&  startPos,
                              const Vector<long, D>&  endPos,
                              const Vector<T, D>&     contactPhysDeltaX,
                              const Vector<long, D>&  currPos,
                              const Vector<short, D>& dirIn, F1 isInsideContact,
                              F2 updateMinMax)
{
  PhysR<T, D> currPhysPos;
  for (unsigned iD = 0; iD < D; ++iD) {
    currPhysPos[iD] = contactPhysDeltaX[iD] * currPos[iD];
  }
  bool isInside   = isInsideContact(currPhysPos);
  bool isInsideBB = true;
  for (unsigned iD = 0; iD < D; ++iD) {
    if (currPos[iD] < startPos[iD] || currPos[iD] > endPos[iD]) {
      isInsideBB = false;
      break;
    }
  }

  // if point is inside, then we check the neighboring points too, otherwise they cannot be inside for convex shapes
  if (isInside || isInsideBB) {
    if (isInside) {
      updateMinMax(currPhysPos);
    }

    forEachPosition(Vector<unsigned short, D>(0),
                    Vector<unsigned short, D>(abs(dirIn)),
                    [&](const Vector<unsigned short, D>& dirOut) {
                      if (norm(dirOut) >= 1) {
                        Vector<short, D> signedDirOut;
                        for (unsigned iD = 0; iD < D; ++iD) {
                          signedDirOut[iD] = dirOut[iD] * dirIn[iD];
                        }
                        const Vector<long, D> nextPos = signedDirOut + currPos;
                        evalBoundingBoxRecursive(
                            startPos, endPos, contactPhysDeltaX, nextPos,
                            signedDirOut, isInsideContact, updateMinMax);
                      }
                    });
  }
}

template <typename T, unsigned D, typename F1, typename F2, typename F3>
void correctBoundingBox(const PhysR<T, D>& min, const PhysR<T, D>& max,
                        const T        physDeltaX,
                        const unsigned contactBoxResolutionPerDirection,
                        F1 isInsideContact, F2 resetContactBox, F3 updateMinMax)
{
  const Vector<T, D> contactPhysDeltaX =
      evalContactDeltaX(min, max, physDeltaX, contactBoxResolutionPerDirection);
  Vector<long, D> startPos;
  Vector<long, D> endPos;
  evalStartAndEndPos(startPos, endPos, max, min, contactPhysDeltaX);
  resetContactBox();

  Vector<long, D> origin;
  const bool      originFound = evalStartingPoint(
      startPos, endPos, contactPhysDeltaX, origin, isInsideContact);

  if (originFound) {
    PhysR<T, D> physOrigin;
    for (unsigned iD = 0; iD < D; ++iD) {
      physOrigin[iD] = contactPhysDeltaX[iD] * origin[iD];
    }
    updateMinMax(physOrigin);

    for (unsigned i = 0; i < utilities::dimensions::convert<D>::neighborsCount;
         ++i) {
      const Vector<short, D> direction =
          utilities::dimensions::convert<D>::neighborDirections[i];
      evalBoundingBoxRecursive(startPos, endPos, contactPhysDeltaX,
                               origin + direction, direction, isInsideContact,
                               updateMinMax);
    }
  }
}

template <typename T, unsigned D, typename F1, typename F2, typename F3>
void correctBoundingBoxNewContact(PhysR<T, D> min, PhysR<T, D> max,
                                  T        physDeltaX,
                                  unsigned contactBoxResolutionPerDirection,
                                  F1       signedDistanceToIntersection,
                                  F2 resetContactBox, F3 updateMinMax)
{
  Vector<T, D> contactPhysDeltaX =
      evalContactDeltaX(min, max, physDeltaX, contactBoxResolutionPerDirection);

  for (unsigned iD = 0; iD < D; ++iD) {
    if (util::nearZero(min[iD] - max[iD])) {
      min[iD] -= contactPhysDeltaX[iD];
      /* T {0.5} * contactBoxResolutionPerDirection * contactPhysDeltaX[iD]; */
      max[iD] += contactPhysDeltaX[iD];
      /* T {0.5} * contactBoxResolutionPerDirection * contactPhysDeltaX[iD]; */
    }
  }

  Vector<long, D> startPos;
  Vector<long, D> endPos;
  evalStartAndEndPos(startPos, endPos, max, min, contactPhysDeltaX);
  resetContactBox();

  forEachPosition(startPos, endPos, [&](const Vector<long, D>& pos) {
    bool isOnContactSurface = false;
    for (unsigned iD = 0; iD < D; ++iD) {
      if (pos[iD] == startPos[iD] || pos[iD] == endPos[iD]) {
        isOnContactSurface = true;
        break;
      }
    }

    if (isOnContactSurface) {
      // Determine physical position
      PhysR<T, D> physPos;
      for (unsigned iD = 0; iD < D; ++iD) {
        physPos[iD] = pos[iD] * contactPhysDeltaX[iD];
      }

      for (unsigned i = 0;
           i < utilities::dimensions::convert<D>::neighborsCount; ++i) {
        // do not regard direction if the neighbor is also on the surface
        bool regardDirection = true;
        //const Vector<long,D> neighborPos = pos + direction;
        const Vector<long, D> neighborPos =
            pos + Vector<long, D>(
                      utilities::dimensions::convert<D>::neighborDirections[i]);
        for (unsigned iD = 0; iD < D; ++iD) {
          if (neighborPos[iD] <= startPos[iD] ||
              neighborPos[iD] >= endPos[iD]) {
            regardDirection = false;
            break;
          }
        }

        if (regardDirection) {
          const PhysR<T, D> normalizedDirection =
              PhysR<T, D>(utilities::dimensions::convert<
                          D>::template normalizedNeighborDirections<T>[i]);
          T distance;

          // TODO: the sdf from intersections is not exact - replace with another option
          // Only change min & max if boundary was found
          /* if (util::distance( */
          /*         distance, physPos, direction, physDeltaX * T {1e-2}, */
          /*         [&](const Vector<T, D>& pos) { */
          /*           return signedDistanceToIntersection(pos); */
          /*         }, */
          /*         [&](const Vector<T, D>& pos) { */
          /*           return signedDistanceToIntersection(pos) < physDeltaX; */
          /*         })) { */
          if (util::distance<T, D, false>(
                  distance, physPos, normalizedDirection, physDeltaX * T {1e-2},
                  [&](const Vector<T, D>& pos) {
                    return signedDistanceToIntersection(pos);
                  },
                  physDeltaX)) {
            // Determine position of boundary in given direction
            const PhysR<T, D> posOnBoundary =
                physPos + distance * normalizedDirection;

            // Attempt to update min and max with found position on contact boundary
            updateMinMax(posOnBoundary);
          }
        }
      }
    }
  });
}

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, unsigned BBCORRECTIONMETHOD = 0,
          bool CONVEX = true, bool useSDF = true>
struct particle_particle;

template <typename T, typename PARTICLETYPE, typename WALLCONTACTTYPE,
          unsigned BBCORRECTIONMETHOD, bool CONVEX, bool useSDF>
struct particle_particle<
    T, PARTICLETYPE,
    ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>,
    WALLCONTACTTYPE, BBCORRECTIONMETHOD, CONVEX, useSDF> {
  static void resetContainer(
      ContactContainer<
          T,
          ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>,
          WALLCONTACTTYPE>& contactContainer)
  {
    for (auto& contact : contactContainer.particleContacts) {
      contact.resetMinMax();
    }
  }

  static void
  processCell(ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d,
                                                        CONVEX>& contact,
              Particle<T, PARTICLETYPE>&                         particleA,
              Particle<T, PARTICLETYPE>&                         particleB,
              Vector<T, PARTICLETYPE::d>&                        contactNormal,
              Vector<T, PARTICLETYPE::d>& center, unsigned& volumeCellCount,
              T& surfaceCellCount, const Vector<T, PARTICLETYPE::d>& pos,
              const Vector<T, PARTICLETYPE::d>& contactPhysDeltaX)
  {
    constexpr unsigned D = PARTICLETYPE::d;

    // Check if pos is inside contact
    const std::function<bool(const PhysR<T, D>&)> isInside =
        [&](const PhysR<T, D>& pos) {
          bool const isInsideParticleA =
              particles::resolved::signedDistanceToParticle(particleA, pos) <=
              particles::access::getEnlargementForContact(particleA);
          bool const isInsideParticleB =
              particles::resolved::signedDistanceToParticle(particleB, pos) <=
              particles::access::getEnlargementForContact(particleB);
          return isInsideParticleA && isInsideParticleB;
        };

    // Check if pos is on surface
    const std::function<bool(const PhysR<T, D>&, const PhysR<T, D>&, bool&,
                             bool&, const Vector<T, D>&, const Vector<T, D>&)>
        isOnSurface = [&](const PhysR<T, D>& pos,
                          const PhysR<T, D>& contactPhysDeltaX,
                          bool& onSurfaceA, bool& onSurfaceB,
                          const Vector<T, D>& normalA,
                          const Vector<T, D>& normalB) {
          Vector<T, D> neighbor;
          for (unsigned iD = 0; iD < D; ++iD) {
            neighbor[iD] = -normalA[iD] * contactPhysDeltaX[iD];
          }
          onSurfaceA = particles::resolved::signedDistanceToParticle(
                           particleA, pos + neighbor) <=
                       particles::access::getEnlargementForContact(particleA);
          for (unsigned iD = 0; iD < D; ++iD) {
            neighbor[iD] = -normalB[iD] * contactPhysDeltaX[iD];
          }
          onSurfaceB = particles::resolved::signedDistanceToParticle(
                           particleB, pos + neighbor) <=
                       particles::access::getEnlargementForContact(particleB);
          return onSurfaceA || onSurfaceB;
        };

    // Calculate normal to particle surface
    const std::function<Vector<T, D>(const PhysR<T, D>&, const T)> calcNormalA =
        [&](const PhysR<T, D>& pos, const T meshSize) {
          // also with an enlarged particle, the normal should still point the same way, since we apply an equal layer to the whole surface
          return particles::resolved::normalOnParticleSurface(particleA, pos,
                                                              meshSize);
        };
    // Calculate normal to indicator surface
    const std::function<Vector<T, D>(const PhysR<T, D>&, const T)> calcNormalB =
        [&](const PhysR<T, D>& pos, const T meshSize) {
          // also with an enlarged particle, the normal should still point the same way, since we apply an equal layer to the whole surface
          return particles::resolved::normalOnParticleSurface(particleB, pos,
                                                              meshSize);
        };
    // Wrapper for update min and max
    const std::function<void(const PhysR<T, D>&)> updateMinMax =
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        };

    particles::contact::processCell(pos, contactPhysDeltaX, center,
                                    contactNormal, volumeCellCount,
                                    surfaceCellCount, isInside, isOnSurface,
                                    calcNormalA, calcNormalB, updateMinMax);
  }

  template <
      typename CONTACTPROPERTIES,
      typename F = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
  static void apply(
      std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
      XParticleSystem<T, PARTICLETYPE>&                    particleSystem,
      ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                               contact,
      const CONTACTPROPERTIES& contactProperties, const T physDeltaX,
      const unsigned contactBoxResolutionPerDirection, const T k,
      F processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    // indentation depth
    T indentationDepth = T(0.);

    if (!contact.isEmpty()) {
      OstreamManager clout(std::cout, "forceApplication");
      // particles in contact
      auto particleA = particleSystem.get(contact.getIDs()[0]);
      auto particleB = particleSystem.get(contact.getIDs()[1]);

      // function to reset min & max of bounding box
      const std::function<void()> resetMinMax = [&contact]() {
        contact.resetMinMax();
      };
      // wrapper for cell processing to calculate contact point and normal
      const std::function<void(Vector<T, D>&, Vector<T, D>&, unsigned&, T&,
                               const Vector<T, D>&, const Vector<T, D>&)>
          processCellWrapped =
              [&](Vector<T, D>& contactNormal, Vector<T, D>& center,
                  unsigned& volumeCellCount, T& surfaceCellCount,
                  const Vector<T, D>& physPos,
                  const Vector<T, D>& contactPhysDeltaX) {
                processCell(contact, particleA, particleB, contactNormal,
                            center, volumeCellCount, surfaceCellCount, physPos,
                            contactPhysDeltaX);
              };
      // function to calculate indendation depth
      const std::function<T(const Vector<T, D>&, const Vector<T, D>&, const T)>
          calculateIndentation = [&](const Vector<T, D>& center,
                                     const Vector<T, D>& contactNormal,
                                     const T             distancePrecision) {
            const PhysR<T, D> originA =
                center -
                contactNormal *
                    particles::access::getEnlargementForContact(particleA);
            const PhysR<T, D> originB =
                center +
                contactNormal *
                    particles::access::getEnlargementForContact(particleB);
            T indentation =
                particles::resolved::distanceToParticle(
                    particleA, originA, contactNormal, distancePrecision) +
                particles::resolved::distanceToParticle(
                    particleB, originB, -1 * contactNormal, distancePrecision);
            return indentation;
          };
      // function to apply force on particles
      const std::function<void(const Vector<T, D>&, const Vector<T, D>&, T, T)>
          applyForce = [&](const Vector<T, D>& center,
                           const Vector<T, D>& contactNormal, T indentation,
                           T overlapVolume) {
            const unsigned materialA =
                particleA.template getField<MECHPROPERTIES, MATERIAL>();
            const unsigned materialB =
                particleB.template getField<MECHPROPERTIES, MATERIAL>();

            indentationDepth = indentation;
            const Vector<T, D> relVel =
                evalRelativeVelocity(particleA, particleB, center);
            const T dampingFactor = evalCurrentDampingFactor(
                contact,
                contactProperties.getCoefficientOfRestitution(materialA,
                                                              materialB),
                evalRelativeNormalVelocity(contactNormal, relVel));

            applyForceFromOverlapVolume(
                dataMap, particleA, particleB, particleSystem, center,
                contactNormal, indentation, overlapVolume, relVel,
                dampingFactor, contactProperties, k, contact.getIDs()[0],
                contact.getIDs()[1], processContactForce);
          };

      processContactViaVolume(contact.getMin(), contact.getMax(), physDeltaX,
                              contactBoxResolutionPerDirection, resetMinMax,
                              processCellWrapped, calculateIndentation,
                              applyForce);
    }
  }

  static void correctBoundingBoxNewContact(
      XParticleSystem<T, PARTICLETYPE>& particleSystem,
      ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
              contact,
      const T physDeltaX, const unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    auto particleA = particleSystem.get(contact.getIDs()[0]);
    auto particleB = particleSystem.get(contact.getIDs()[1]);

    particles::contact::correctBoundingBoxNewContact(
        contact.getMin(), contact.getMax(), physDeltaX,
        contactBoxResolutionPerDirection,
        [&](const Vector<T, D>& pos) {
          return sdf::intersection(
              sdf::rounding(
                  particles::resolved::signedDistanceToParticle(particleA, pos),
                  particles::access::getEnlargementForContact(particleA)),
              sdf::rounding(
                  particles::resolved::signedDistanceToParticle(particleB, pos),
                  particles::access::getEnlargementForContact(particleB)));
        },
        [&contact]() {
          contact.resetMinMax();
        },
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        });
  }

  static void correctBoundingBox(
      XParticleSystem<T, PARTICLETYPE>&               particleSystem,
      std::vector<SolidBoundary<T, PARTICLETYPE::d>>& solidBoundaries,
      WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
              contact,
      const T physDeltaX, const unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    Vector<T, D> contactPhysDeltaX =
        evalContactDeltaX(contact.getMin(), contact.getMax(), physDeltaX,
                          contactBoxResolutionPerDirection);
    contact.increaseMinMax(T {0.5} * util::sqrt(D) * contactPhysDeltaX);
    correctBoundingBoxNewContact(particleSystem, contact, physDeltaX,
                                 contactBoxResolutionPerDirection);
  }

  static void correctBoundingBoxExistingContact(
      XParticleSystem<T, PARTICLETYPE>& particleSystem,
      ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
          contact,
      T physDeltaX, unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    auto particleA = particleSystem.get(contact.getIDs()[0]);
    auto particleB = particleSystem.get(contact.getIDs()[1]);

    const std::function<bool(const PhysR<T, D>&)> isInsideContact =
        [&](const Vector<T, D>& pos) {
          return particles::resolved::isInsideParticle(particleA, pos) &&
                 particles::resolved::isInsideParticle(particleB, pos);
        };
    const std::function<void()> resetMinMax = [&contact]() {
      contact.resetMinMax();
    };
    const std::function<void(const PhysR<T, D>&)> updateMinMax =
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        };

    particles::contact::correctBoundingBox(
        contact.getMin(), contact.getMax(), physDeltaX,
        contactBoxResolutionPerDirection, isInsideContact, resetMinMax,
        updateMinMax);
  }

  template <
      typename CONTACTPROPERTIES,
      typename F1 = decltype(defaults::periodicity<PARTICLETYPE::d>),
      typename F2 = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
  static void apply(
      XParticleSystem<T, PARTICLETYPE>& particleSystem,
      ContactContainer<
          T,
          ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>,
          WALLCONTACTTYPE>&                                contactContainer,
      const CONTACTPROPERTIES&                             contactProperties,
      const SuperGeometry<T, PARTICLETYPE::d>&             sGeometry,
      std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
      const unsigned contactBoxResolutionPerDirection = 8,
      const T        k       = T {4. / (3 * util::sqrt(M_PI))},
      F1 getSetupPeriodicity = defaults::periodicity<PARTICLETYPE::d>,
      F2 processContactForce =
          defaults::processContactForce<T, PARTICLETYPE::d>)
  {
    using namespace descriptors;
    const T physDeltaX =
        sGeometry.getCuboidDecomposition().getMotherCuboid().getDeltaR();

    for (auto& particleContact : contactContainer.particleContacts) {
      if (!particleContact.isEmpty()) {
        bool isResponsible = true;

        if constexpr (access::providesParallelization<PARTICLETYPE>()) {
          isResponsible = singleton::mpi().getRank() ==
                          particleContact.getResponsibleRank();
        }

        if (isResponsible) {
          // Check if contact should be computed
          bool computeContact = true;
          if constexpr (access::providesComputeContact<PARTICLETYPE>()) {
            auto particleA = particleSystem.get(particleContact.getIDs()[0]);
            auto particleB = particleSystem.get(particleContact.getIDs()[1]);
            computeContact =
                access::isContactComputationEnabled(particleA, particleB);
          }

          if (computeContact) {
            std::array<PhysR<T, PARTICLETYPE::d>, 2> originalPos;
            if constexpr (isPeriodic(getSetupPeriodicity())) {
              for (unsigned i = 0; i < 2; ++i) {
                auto particle = particleSystem.get(particleContact.getIDs()[i]);
                originalPos[i] =
                    particle.template getField<GENERAL, POSITION>();
                particle.template setField<GENERAL, POSITION>(
                    particleContact.getParticlePositions()[i]);
              }
            }

            if (particleContact.isNew()) {
              correctBoundingBoxNewContact(particleSystem, particleContact,
                                           physDeltaX,
                                           contactBoxResolutionPerDirection);
              particleContact.isNew(particleContact.isEmpty());
            }
            else {
              if constexpr (BBCORRECTIONMETHOD == 1) {
                correctBoundingBoxExistingContact(
                    particleSystem, particleContact, physDeltaX,
                    contactBoxResolutionPerDirection);
              }
              else {
                correctBoundingBoxNewContact(particleSystem, particleContact,
                                             physDeltaX,
                                             contactBoxResolutionPerDirection);
              }
            }
            apply(dataMap, particleSystem, particleContact, contactProperties,
                  physDeltaX, contactBoxResolutionPerDirection, k,
                  processContactForce);

            if constexpr (isPeriodic(getSetupPeriodicity())) {
              for (unsigned i = 0; i < 2; ++i) {
                auto particle = particleSystem.get(particleContact.getIDs()[i]);
                particle.template setField<GENERAL, POSITION>(originalPos[i]);
              }
            }
            particleContact.setParticlePositionUpdated(false);
          }
          else {
            // if contact should be ignored, reset min and max so that it counts
            // as empty and can be deleted during the next cleanContact() call.
            particleContact.resetMinMax();
          }
        }
      }
    }
  }
};

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          typename WALLCONTACTTYPE, unsigned BBCORRECTIONMETHOD = 0,
          bool CONVEX = true, bool useSDF = true>
struct particle_wall;

template <typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
          unsigned BBCORRECTIONMETHOD, bool CONVEX, bool useSDF>
struct particle_wall<
    T, PARTICLETYPE, PARTICLECONTACTTYPE,
    WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>,
    BBCORRECTIONMETHOD, CONVEX, useSDF> {
  static void resetContainer(
      ContactContainer<
          T, PARTICLECONTACTTYPE,
          WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>>&
          contactContainer)
  {
    for (auto& contact : contactContainer.wallContacts) {
      contact.resetMinMax();
    }
  }

  static void
  processCell(WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                                                 contact,
              Particle<T, PARTICLETYPE>&         particle,
              SolidBoundary<T, PARTICLETYPE::d>& solidBoundary,
              Vector<T, PARTICLETYPE::d>&        contactNormal,
              Vector<T, PARTICLETYPE::d>& center, unsigned& volumeCellCount,
              T& surfaceCellCount, const Vector<T, PARTICLETYPE::d>& pos,
              const Vector<T, PARTICLETYPE::d>& contactPhysDeltaX)
  {
    /* OstreamManager     clout(std::cout, "cellProcessingForContactForce"); */
    constexpr unsigned D = PARTICLETYPE::d;

    // Check if pos is inside contact
    const std::function<bool(const PhysR<T, D>&)> isInside =
        [&](const PhysR<T, D>& pos) {
          bool const isInsideIndicator =
              solidBoundary.getIndicator()->signedDistanceExact(pos.data()) <=
              solidBoundary.getEnlargementForContact();
          bool const isInsideParticle =
              particles::resolved::signedDistanceToParticle(particle, pos) <=
              particles::access::getEnlargementForContact(particle);
          return isInsideParticle && isInsideIndicator;
        };

    // Check if pos is on surface
    const std::function<bool(const PhysR<T, D>&, const PhysR<T, D>&, bool&,
                             bool&, const Vector<T, D>&, const Vector<T, D>&)>
        isOnSurface = [&](const PhysR<T, D>& pos,
                          const PhysR<T, D>& contactPhysDeltaX,
                          bool& onSurfaceA, bool& onSurfaceB,
                          const Vector<T, D>& normalA,
                          const Vector<T, D>& normalB) {
          Vector<T, D> neighbor;
          for (unsigned iD = 0; iD < D; ++iD) {
            neighbor[iD] = -normalA[iD] * contactPhysDeltaX[iD];
          }
          onSurfaceA = particles::resolved::signedDistanceToParticle(
                           particle, pos + neighbor) <=
                       particles::access::getEnlargementForContact(particle);
          for (unsigned iD = 0; iD < D; ++iD) {
            neighbor[iD] = -normalB[iD] * contactPhysDeltaX[iD];
          }
          onSurfaceB = solidBoundary.getIndicator()->signedDistanceExact(
                           (pos + neighbor).data()) <=
                       solidBoundary.getEnlargementForContact();
          return onSurfaceA || onSurfaceB;
        };

    // Calculate normal to particle surface
    const std::function<Vector<T, D>(const PhysR<T, D>&, const T)> calcNormalA =
        [&](const PhysR<T, D>& pos, const T meshSize) {
          // also with an enlarged particle, the normal should still point the same way, since we apply an equal layer to the whole surface
          return particles::resolved::normalOnParticleSurface(particle, pos,
                                                              meshSize);
        };
    // Calculate normal to indicator surface
    const std::function<Vector<T, D>(const PhysR<T, D>&, const T)> calcNormalB =
        [&](const PhysR<T, D>& pos, const T meshSize) {
          return solidBoundary.getIndicator()->surfaceNormalExact(pos, meshSize);
        };
    // Wrapper for update of min and max
    const std::function<void(const PhysR<T, D>&)> updateMinMax =
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        };

    particles::contact::processCell(pos, contactPhysDeltaX, center,
                                    contactNormal, volumeCellCount,
                                    surfaceCellCount, isInside, isOnSurface,
                                    calcNormalA, calcNormalB, updateMinMax);
  }

  template <
      typename CONTACTPROPERTIES,
      typename F = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
  static void apply(
      std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
      XParticleSystem<T, PARTICLETYPE>&                    particleSystem,
      std::vector<SolidBoundary<T, PARTICLETYPE::d>>&      solidBoundaries,
      WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
                               contact,
      const CONTACTPROPERTIES& contactProperties, const T physDeltaX,
      const unsigned contactBoxResolutionPerDirection, const T k,
      F processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    // indentation depth
    T indentationDepth = T(0.);

    if (!contact.isEmpty()) {
      // particle and wall in contact
      auto particle = particleSystem.get(contact.getParticleID());
      SolidBoundary<T, D>& solidBoundary = solidBoundaries[contact.getWallID()];

      // function to reset min & max of bounding box
      const std::function<void()> resetMinMax = [&contact]() {
        contact.resetMinMax();
      };
      // wrapper for cell processing to calculate contact point and normal
      const std::function<void(Vector<T, D>&, Vector<T, D>&, unsigned&, T&,
                               const Vector<T, D>&, const Vector<T, D>&)>
          processCellWrapped =
              [&](Vector<T, D>& contactNormal, Vector<T, D>& center,
                  unsigned& volumeCellCount, T& surfaceCellCount,
                  const Vector<T, D>& physPos,
                  const Vector<T, D>& contactPhysDeltaX) {
                processCell(contact, particle, solidBoundary, contactNormal,
                            center, volumeCellCount, surfaceCellCount, physPos,
                            contactPhysDeltaX);
              };
      // function to calculate indendation depth
      const std::function<T(const Vector<T, D>&, const Vector<T, D>&, const T)>
          calculateIndentation = [&](const Vector<T, D>& center,
                                     const Vector<T, D>& contactNormal,
                                     const T             distancePrecision) {
            T indentation =
                particles::resolved::signedDistanceToParticle(particle,
                                                              center) +
                particles::access::getEnlargementForContact(particle);
            if (!util::nearZero(indentation)) {
              const PhysR<T, D> origin =
                  center -
                  contactNormal *
                      particles::access::getEnlargementForContact(particle);
              indentation = particles::resolved::distanceToParticle(
                  particle, origin, contactNormal, distancePrecision);
            }
            T distToWall;
            if (solidBoundary.getIndicator()->distance(distToWall, center,
                                                       distancePrecision,
                                                       -1 * contactNormal)) {
              indentation +=
                  distToWall + solidBoundary.getEnlargementForContact();
            }
#ifdef OLB_DEBUG
            else if (!util::nearZero(
                         solidBoundary.getIndicator()->signedDistanceExact(
                             center))) {
              OstreamManager clout(std::cout, "forceApplication");
              clout << "WARNING: No distance to wall determined." << std::endl;
            }
#endif
            return indentation;
          };
      // function to apply force on particles
      const std::function<void(const Vector<T, D>&, const Vector<T, D>&, T, T)>
          applyForce = [&](const Vector<T, D>& center,
                           const Vector<T, D>& contactNormal, T indentation,
                           T overlapVolume) {
            const unsigned particleMaterial =
                particle.template getField<MECHPROPERTIES, MATERIAL>();
            const unsigned wallMaterial = solidBoundary.getContactMaterial();

            indentationDepth          = indentation;
            const Vector<T, D> relVel = evalRelativeVelocity(particle, center);
            const T            dampingFactor = evalCurrentDampingFactor(
                contact,
                contactProperties.getCoefficientOfRestitution(particleMaterial,
                                                                         wallMaterial),
                evalRelativeNormalVelocity(contactNormal, relVel));

            applyForceFromOverlapVolume(
                dataMap, particle, particleSystem, solidBoundary, center,
                contactNormal, indentation, overlapVolume, relVel,
                dampingFactor, contactProperties, k, contact.getParticleID(),
                contact.getWallID(), processContactForce);
          };

      processContactViaVolume(contact.getMin(), contact.getMax(), physDeltaX,
                              contactBoxResolutionPerDirection, resetMinMax,
                              processCellWrapped, calculateIndentation,
                              applyForce);
    }
  }

  static void correctBoundingBoxNewContact(
      XParticleSystem<T, PARTICLETYPE>&               particleSystem,
      std::vector<SolidBoundary<T, PARTICLETYPE::d>>& solidBoundaries,
      WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
              contact,
      const T physDeltaX, const unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    Particle<T, PARTICLETYPE> particle =
        particleSystem.get(contact.getParticleID());
    SolidBoundary<T, D>& solidBoundary = solidBoundaries[contact.getWallID()];

    // TODO: Is there a better way to increase the particle size than sdf::rounding?
    particles::contact::correctBoundingBoxNewContact(
        contact.getMin(), contact.getMax(), physDeltaX,
        contactBoxResolutionPerDirection,
        [&](const PhysR<T, D>& pos) {
          return sdf::intersection(
              sdf::rounding(
                  particles::resolved::signedDistanceToParticle(particle, pos),
                  particles::access::getEnlargementForContact(particle)),
              sdf::rounding(solidBoundary.getIndicator()->signedDistanceExact(pos),
                            solidBoundary.getEnlargementForContact()));
        },
        [&contact]() {
          contact.resetMinMax();
        },
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        });
  }

  static void correctBoundingBox(
      XParticleSystem<T, PARTICLETYPE>&               particleSystem,
      std::vector<SolidBoundary<T, PARTICLETYPE::d>>& solidBoundaries,
      WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
              contact,
      const T physDeltaX, const unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    Vector<T, D> contactPhysDeltaX =
        evalContactDeltaX(contact.getMin(), contact.getMax(), physDeltaX,
                          contactBoxResolutionPerDirection);
    contact.increaseMinMax(T {0.5} * util::sqrt(D) * contactPhysDeltaX);
    correctBoundingBoxNewContact(particleSystem, solidBoundaries, contact,
                                 physDeltaX, contactBoxResolutionPerDirection);
  }

  static void correctBoundingBoxExistingContact(
      XParticleSystem<T, PARTICLETYPE>&               particleSystem,
      std::vector<SolidBoundary<T, PARTICLETYPE::d>>& solidBoundaries,
      WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>&
              contact,
      const T physDeltaX, const unsigned contactBoxResolutionPerDirection)
  {
    using namespace descriptors;
    constexpr unsigned D = PARTICLETYPE::d;

    auto                 particle = particleSystem.get(contact.getParticleID());
    SolidBoundary<T, D>& solidBoundary = solidBoundaries[contact.getWallID()];

    const std::function<bool(const PhysR<T, D>&)> isInsideContact =
        [&](const Vector<T, D>& pos) {
          bool                          isInsideWall;
          solidBoundary.getIndicator()->operator()(&isInsideWall, pos.data());
          return particles::resolved::isInsideParticle(particle, pos) &&
                 isInsideWall;
        };
    const std::function<void()> resetMinMax = [&contact]() {
      contact.resetMinMax();
    };
    const std::function<void(const PhysR<T, D>&)> updateMinMax =
        [&contact](const PhysR<T, D>& pos) {
          contact.updateMinMax(pos);
        };

    particles::contact::correctBoundingBox(
        contact.getMin(), contact.getMax(), physDeltaX,
        contactBoxResolutionPerDirection, isInsideContact, resetMinMax,
        updateMinMax);
  }

  template <
      typename CONTACTPROPERTIES,
      typename F1 = decltype(defaults::periodicity<PARTICLETYPE::d>),
      typename F2 = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
  static void
  apply(XParticleSystem<T, PARTICLETYPE>&               particleSystem,
        std::vector<SolidBoundary<T, PARTICLETYPE::d>>& solidBoundaries,
        ContactContainer<
            T, PARTICLECONTACTTYPE,
            WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, CONVEX>>&
                                                             contactContainer,
        const CONTACTPROPERTIES&                             contactProperties,
        const SuperGeometry<T, PARTICLETYPE::d>&             sGeometry,
        std::multimap<int, std::unique_ptr<std::uint8_t[]>>& dataMap,
        const unsigned contactBoxResolutionPerDirection = 8,
        const T        k       = T {4. / (3 * util::sqrt(M_PI))},
        F1 getSetupPeriodicity = defaults::periodicity<PARTICLETYPE::d>,
        F2 processContactForce =
            defaults::processContactForce<T, PARTICLETYPE::d>)
  {
    using namespace descriptors;
    const T physDeltaX =
        sGeometry.getCuboidDecomposition().getMotherCuboid().getDeltaR();

    for (auto& wallContact : contactContainer.wallContacts) {
      if (!wallContact.isEmpty()) {
        bool isResponsible = true;

        if constexpr (access::providesParallelization<PARTICLETYPE>()) {
          isResponsible =
              singleton::mpi().getRank() == wallContact.getResponsibleRank();
        }

        if (isResponsible) {
          // Check if contact should be computed
          bool computeContact = true;
          if constexpr (access::providesComputeContact<PARTICLETYPE>()) {
            auto particle  = particleSystem.get(wallContact.getParticleID());
            computeContact = access::isContactComputationEnabled(particle);
          }

          if (computeContact) {
            PhysR<T, PARTICLETYPE::d> originalPos;
            if constexpr (isPeriodic(getSetupPeriodicity())) {
              auto particle = particleSystem.get(wallContact.getParticleID());
              originalPos   = particle.template getField<GENERAL, POSITION>();
              particle.template setField<GENERAL, POSITION>(
                  wallContact.getParticlePosition());
            }

            if (wallContact.isNew()) {
              correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                           wallContact, physDeltaX,
                                           contactBoxResolutionPerDirection);
              wallContact.isNew(wallContact.isEmpty());
            }
            else {
              if constexpr (BBCORRECTIONMETHOD == 1) {
                correctBoundingBoxExistingContact(
                    particleSystem, solidBoundaries, wallContact, physDeltaX,
                    contactBoxResolutionPerDirection);
              }
              else {
                correctBoundingBox(particleSystem, solidBoundaries, wallContact,
                                   physDeltaX,
                                   contactBoxResolutionPerDirection);
              }
            }
            apply(dataMap, particleSystem, solidBoundaries, wallContact,
                  contactProperties, physDeltaX,
                  contactBoxResolutionPerDirection, k, processContactForce);

            if constexpr (isPeriodic(getSetupPeriodicity())) {
              particleSystem.get(wallContact.getParticleID())
                  .template setField<GENERAL, POSITION>(originalPos);
              // TODO: Set that unfied position is not set anymore
            }
            wallContact.setParticlePositionUpdated(false);
          }
          else {
            // if contact should be ignored, reset min and max so that it counts
            // as empty and can be deleted during the next cleanContact() call.
            wallContact.resetMinMax();
          }
        }
      }
    }
  }
};

/*
 * Combined contact processing
 * Adds an option to use different bounding box correction methods for existing contacts:
 * 0: Less fast, but more accurate option (this option is always used for new contacts)
 * 1: Faster, but inaccurate option (this option needs tests)
 */
template <
    typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
    typename WALLCONTACTTYPE, typename CONTACTPROPERTIES,
    unsigned BBCORRECTIONMETHOD = 0,
    typename F1 = decltype(defaults::periodicity<PARTICLETYPE::d>),
    typename F2 = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
void processContacts(
    SuperParticleSystem<T, PARTICLETYPE>&                      particleSystem,
    std::vector<SolidBoundary<T, PARTICLETYPE::d>>&            solidBoundaries,
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    const CONTACTPROPERTIES&                 contactProperties,
    const SuperGeometry<T, PARTICLETYPE::d>& sGeometry,
#ifdef PARALLEL_MODE_MPI
    MPI_Comm contactTreatmentComm,
#endif
    const unsigned contactBoxResolutionPerDirection = 8,
    const T        k = static_cast<T>(4. / (3 * util::sqrt(M_PI))),
    F1             getSetupPeriodicity = defaults::periodicity<PARTICLETYPE::d>,
    F2 processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)

{
  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE,
                    BBCORRECTIONMETHOD>::apply(particleSystem, contactContainer,
                                               contactProperties, sGeometry,
                                               dataMap,
                                               contactBoxResolutionPerDirection,
                                               k, getSetupPeriodicity,
                                               processContactForce);

  particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE,
                BBCORRECTIONMETHOD>::apply(particleSystem, solidBoundaries,
                                           contactContainer, contactProperties,
                                           sGeometry, dataMap,
                                           contactBoxResolutionPerDirection, k,
                                           getSetupPeriodicity,
                                           processContactForce);
  if constexpr (access::providesParallelization<PARTICLETYPE>()) {
    communicateContactTreatmentResults<T, PARTICLETYPE>(particleSystem, dataMap
#ifdef PARALLEL_MODE_MPI
                                                        ,
                                                        contactTreatmentComm
#endif
    );
  }
}

template <
    typename T, typename PARTICLETYPE, typename PARTICLECONTACTTYPE,
    typename WALLCONTACTTYPE, typename CONTACTPROPERTIES,
    unsigned BBCORRECTIONMETHOD = 0,
    typename F1 = decltype(defaults::periodicity<PARTICLETYPE::d>),
    typename F2 = decltype(defaults::processContactForce<T, PARTICLETYPE::d>)>
void processContacts(
    ParticleSystem<T, PARTICLETYPE>&                           particleSystem,
    std::vector<SolidBoundary<T, PARTICLETYPE::d>>&            solidBoundaries,
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& contactContainer,
    const CONTACTPROPERTIES&                 contactProperties,
    const SuperGeometry<T, PARTICLETYPE::d>& sGeometry,
    const unsigned contactBoxResolutionPerDirection = 8,
    const T        k = static_cast<T>(4. / (3 * util::sqrt(M_PI))),
    F1             getSetupPeriodicity = defaults::periodicity<PARTICLETYPE::d>,
    F2 processContactForce = defaults::processContactForce<T, PARTICLETYPE::d>)

{
  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE,
                    BBCORRECTIONMETHOD>::apply(particleSystem, contactContainer,
                                               contactProperties, sGeometry,
                                               dataMap,
                                               contactBoxResolutionPerDirection,
                                               k, getSetupPeriodicity,
                                               processContactForce);

  particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE,
                BBCORRECTIONMETHOD>::apply(particleSystem, solidBoundaries,
                                           contactContainer, contactProperties,
                                           sGeometry, dataMap,
                                           contactBoxResolutionPerDirection, k,
                                           getSetupPeriodicity,
                                           processContactForce);
}

} //namespace contact

} //namespace particles

} //namespace olb

#endif
