/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006 Jonas Latt
 *                2021 Julius Jessberger, Adrian Kummerländer
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

/** \file
 * Implementation of computation and definition of the moments density,
 * velocity, stress.
 */

#ifndef DYNAMICS_MOMENTA_ELEMENTS_H
#define DYNAMICS_MOMENTA_ELEMENTS_H

#include <type_traits>
#include <functional>

#include "descriptor/descriptor.h"
#include "dynamics/lbm.h"
#include "boundary/helper.h"

#include "core/util.h"
#include "core/cellD.h"


namespace olb {

template<typename T, typename DESCRIPTOR> class CellD;
template<typename DESCRIPTOR> struct lbHelpers;

namespace momenta {


// --------------------- EXTRACTED HELPER FUNCTIONS ---------------------------

template<int direction, int orientation, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
V velocityBMRho(CELL& cell, const U& u) any_platform
{
  constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
  constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

  V rhoOnWall{};
  for (auto e : onWallIndices) {
    rhoOnWall += cell[e];
  }

  V rhoNormal{};
  for (auto e : normalIndices) {
    rhoNormal += cell[e];
  }

  return (V{2}*rhoNormal+rhoOnWall+V{1}) / (V{1}+orientation*u[direction]);
}

template<int direction, int orientation, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
V velocityBMP(CELL& cell, const U& u) any_platform
{
  constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
  constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

  auto rho = cell.template getField<descriptors::RHO>();

  V pOnWall{};
  for (auto e : onWallIndices) {
    pOnWall += cell[e];
  }
  pOnWall -= cell[0];

  V pNormal{};
  for (auto e : normalIndices) {
    pNormal += cell[e];
  }

  return ((V{2}*pNormal+pOnWall-rho*orientation*u[direction])/descriptors::invCs2<V,DESCRIPTOR>()
        -0.5*rho*descriptors::t<V,DESCRIPTOR>(0)*u[direction]*u[direction]) / (V{1}-descriptors::t<V,DESCRIPTOR>(0));
}

template<int direction, int orientation, typename CELL, typename U, typename FLUX, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
V heatFluxBMRho(CELL& cell, const U& u, FLUX& flux) any_platform
{
  V rho = velocityBMRho<direction,orientation>(cell, u);
  rho -= orientation * flux / (V{1} + orientation * u[direction]);
  return rho;
}


// ---------------------- DENSITY STRUCTURES ----------------------------------

struct ZeroDensity {
  template <typename TYPE, typename CELL, typename RHO, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    rho = 0;
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "ZeroDensity";
  }
};

struct OneDensity {
  template <typename TYPE, typename CELL, typename RHO, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    rho = 1;
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "OneDensity";
  }
};

/// Standard computation for density in the bulk as zeroth moment of the population.
struct BulkDensity {
  template <typename TYPE, typename CELL, typename RHO, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    rho = lbm<DESCRIPTOR>::computeRho(cell);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "BulkDensity";
  }
};

/// Computation of hydrodynamic pressure as zeroth moment of the population.
template<typename MOMENTUM>
struct BulkPressure {
  template <typename TYPE, typename CELL, typename PRESSURE, typename DESCRIPTOR=typename CELL::descriptor_t, typename V=typename CELL::value_t>
  void compute(CELL& cell, PRESSURE& p) any_platform
  {
    V sum = 0;
    for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      sum += cell[iPop];
    }
    auto rho = cell.template getField<descriptors::RHO>();
    auto gradRho = cell.template getField<descriptors::NABLARHO>();
    Vector<V,DESCRIPTOR::d> u{};
    V uGradRho = 0;
    MOMENTUM().template computeU<TYPE>(cell, u);
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      uGradRho += u[iD] * gradRho[iD];
    }
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    const V s0 = -descriptors::t<V,DESCRIPTOR>(0) * descriptors::invCs2<V,DESCRIPTOR>()/V{2} * uSqr;
    p = 1/ descriptors::invCs2<V,DESCRIPTOR>() / (1-descriptors::t<V,DESCRIPTOR>(0)) * (sum + 0.5*uGradRho + rho*s0);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "BulkPressure";
  }
};

template<typename DENSITY>
struct SourcedDensity {
  template <typename TYPE, typename CELL, typename RHO>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    const auto source = cell.template getField<descriptors::SOURCE>();
    DENSITY().template compute<TYPE>(cell, rho);
    rho += 0.5 * source;
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform
  {
    const auto source = cell.template getField<descriptors::SOURCE>();
    rho -= 0.5 * source;
  }

  static std::string getName(){
    return "SourcedDensity<" + DENSITY().getName() + ">";
  }
};

/// The density is fixed and stored in the external field RHO.
struct FixedDensity {
  struct RHO : public descriptors::FIELD_BASE<1, 0, 0> { };

  template <typename TYPE, typename CELL, typename R>
  void compute(CELL& cell, R& rho) any_platform
  {
    rho = cell.template getField<RHO>();
  }

  template <typename TYPE, typename CELL, typename R>
  void define(CELL& cell, const R& rho) any_platform
  {
    cell.template setField<RHO>(rho);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t>
  void initialize(CELL& cell) any_platform
  {
    cell.template setField<RHO>(V{1});
  }

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "FixedDensity";
  }
};

/// The density is stored in descriptors::FORCE[0] (TODO: absurd, to be changed)
struct FreeEnergyInletOutletDensity {
  template <typename TYPE, typename CELL, typename RHO>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    rho = cell.template getField<descriptors::FORCE>()[0];
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform
  {
    cell.template getFieldPointer<descriptors::FORCE>()[0] = rho;
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t>
  void initialize(CELL& cell) any_platform
  {
    cell.template getFieldPointer<descriptors::FORCE>()[0] = V{1};
  }

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "FreeEnergyInletOutletDensity";
  }
};

/** For fixed heat flux, the density is computed from flux, velocity and
 * populations, similar to fixed velocity boundaries.
 */
template <int direction, int orientation>
struct HeatFluxBoundaryDensity {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    const auto uNS = cell.template getField<descriptors::VELOCITY>();

    V conduction[DESCRIPTOR::d];
    TYPE().computeU(cell, conduction);
    rho = heatFluxBMRho<direction,orientation>(cell, uNS.data(), conduction[direction]);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "HeatFluxBoundaryDensity";
  }
};

/// Density computation for fixed velocity boundary
template <int direction, int orientation>
struct VelocityBoundaryDensity {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    V u[DESCRIPTOR::d];
    TYPE().computeU(cell, u);
    rho = velocityBMRho<direction,orientation>(cell, u);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "VelocityBoundaryDensity<" + std::to_string(direction) + "," + std::to_string(orientation) + ">";
  }
};

/// Pressure computation for fixed velocity boundary
template <int direction, int orientation>
struct VelocityBoundaryPressure {
  template <typename TYPE, typename CELL, typename P, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, P& p) any_platform
  {
    V u[DESCRIPTOR::d];
    TYPE().computeU(cell, u);
    p = velocityBMP<direction,orientation>(cell, u);
  }

  template <typename TYPE, typename CELL, typename P>
  void define(CELL& cell, const P& p) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename P>
  void inverseShift(CELL& cell, P& p) any_platform {};

  static std::string getName(){
    return "VelocityBoundaryPressure<" + std::to_string(direction) + "," + std::to_string(orientation) + ">";
  }
};

template <int normalX, int normalY>
struct InnerCornerDensity2D {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    V u[DESCRIPTOR::d];
    TYPE().computeU(cell, u);
    const V rhoX = velocityBMRho<0,normalX>(cell, u);
    const V rhoY = velocityBMRho<1,normalY>(cell, u);
    rho = (rhoX + rhoY) / V{2};
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "InnerCornerDensity2D";
  }
};

template <int normalX, int normalY, int normalZ>
struct InnerCornerDensity3D {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    V u[DESCRIPTOR::d];
    TYPE().computeU(cell, u);
    const V rhoX = velocityBMRho<0,normalX>(cell, u);
    const V rhoY = velocityBMRho<1,normalY>(cell, u);
    const V rhoZ = velocityBMRho<2,normalZ>(cell, u);
    rho = (rhoX + rhoY + rhoZ) / V(3);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "InnerCornerDensity3D";
  }
};

template <int plane, int normal1, int normal2>
struct InnerEdgeDensity3D {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho) any_platform
  {
    V u[DESCRIPTOR::d];
    TYPE().computeU(cell, u);
    const V rho1 = velocityBMRho<(plane+1)%3, normal1>(cell, u);
    const V rho2 = velocityBMRho<(plane+2)%3, normal2>(cell, u);
    rho = (rho1 + rho2) / V(2);
  }

  template <typename TYPE, typename CELL, typename RHO>
  void define(CELL& cell, const RHO& rho) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "InnerEdgeDensity3D";
  }
};


// ---------------------- MOMENTUM STRUCTURES ---------------------------------

/** Standard computation for momentum in the bulk as first moment of the
 * population. Applies the methods from lbHelpers.
 */
struct BulkMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, u);
    const V rho = TYPE().computeRho(cell);
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "BulkMomentum";
  }
};

struct IncompressibleBulkMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, u);
    auto rho = cell.template getField<descriptors::RHO>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
      // TODO do this check relatively / in a more controlled manner
      if (util::fabs(u[iD]) < 1e-15) {
        u[iD] = 0;
      }
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "IncompressibleBulkMomentum";
  }
};

/// The velocity is fixed and stored in the external field U.
struct FixedVelocityMomentumGeneric {
  struct VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> { };

  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    const V rho = TYPE().computeRho(cell);
    auto uExt = cell.template getField<VELOCITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = uExt[iD] * rho;
    }
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    const auto uExt = cell.template getField<VELOCITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = uExt[iD];
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    Vector<V,DESCRIPTOR::d> u{};
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FixedVelocityMomentumGeneric";
  }
};

/** The velocity is stored in the external field U, except for the component
 * "direction", which is computed by means of the population and the pressure.
 */
template <int direction, int orientation>
struct FixedPressureMomentum {
  struct VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> { };

  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    computeU<TYPE>(cell, j);
    const V rho = TYPE().computeRho(cell);
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] *= rho;
    }
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    auto values = cell.template getField<VELOCITY>();
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = values[iD];
    }
    const V rho = TYPE().computeRho(cell);

    constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
    constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

    V rhoOnWall = V{};
    for (auto e : onWallIndices) {
      rhoOnWall += cell[e];
    }

    V rhoNormal = V{};
    for (auto e : normalIndices) {
      rhoNormal += cell[e];
    }

    u[direction] = orientation * ((V{2}*rhoNormal+rhoOnWall+V{1}) / rho-V{1});
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    V u[DESCRIPTOR::d] { V{} };
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FixedPressureMomentum";
  }
};

/** The velocity is stored in the external field U, except for the component
 * "direction", which is computed by means of the population and the pressure.
 * This is for incompressible multi-phase models.
 */
template <int direction, int orientation>
struct PressureBoundaryMomentum {
  struct VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> { };

  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    computeU<TYPE>(cell, j);
    auto rho = cell.template getField<descriptors::RHO>();
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] *= rho;
    }
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    auto values = cell.template getField<VELOCITY>();
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = values[iD];
    }
    auto rho = cell.template getField<descriptors::RHO>();
    auto gradRho = cell.template getField<descriptors::NABLARHO>();
    //V rho = 10.;
    //V gradRho[2] = {0.,0.};
    auto force = cell.template getField<descriptors::FORCE>();
    const V p = TYPE().computeRho(cell);

    constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
    constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

    V popOnWall = V{};
    for (auto e : onWallIndices) {
      popOnWall += cell[e];
    }
    popOnWall -= cell[0];

    V popNormal = V{};
    for (auto e : normalIndices) {
      popNormal += cell[e];
    }

    V w0 = descriptors::t<V,DESCRIPTOR>(0);
    V cs2 = 1./descriptors::invCs2<V,DESCRIPTOR>();
    V pops = V{2}*popNormal+popOnWall;
    V term = cs2*(V{2}*rho-gradRho[direction])/(V{2}*rho*w0);
    V otherU = 0;
    for (int iDim = 1; iDim < DESCRIPTOR::d; ++iDim) {
      u[(direction+iDim)%DESCRIPTOR::d] = 0.5/rho*force[(direction+iDim)%DESCRIPTOR::d];
      otherU += cs2*0.5*u[(direction+iDim)%DESCRIPTOR::d]*gradRho[(direction+iDim)%DESCRIPTOR::d]
                -0.5*w0*rho*u[(direction+iDim)%DESCRIPTOR::d]*u[(direction+iDim)%DESCRIPTOR::d];
    }
    u[direction] = orientation*(-term + util::sqrt(term*term + V{2}/rho/w0*(cs2*(pops+0.5*force[direction])-(V{1}-w0)*p+otherU)));
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    V u[DESCRIPTOR::d] { V{} };
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "PressureBoundaryMomentum";
  }
};

/// The velocity is stored in the external field descriptors::VELOCITY.
struct FixedVelocityMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    const V rho = TYPE().computeRho(cell);
    auto uExt = cell.template getField<descriptors::VELOCITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = uExt[iD] * rho;
    }
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    auto uExt = cell.template getField<descriptors::VELOCITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = uExt[iD];
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template setField<descriptors::VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    const V u[DESCRIPTOR::d] = { V{} };
    cell.template setField<descriptors::VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FixedVelocityMomentum";
  }
};

/** The conduction is computed from density and population.
 * Be careful that computeU(...) computes the conduction whereas compute(...)
 * computes the transport (= conduction + convective transport)
 */
template <int direction, int orientation>
struct FixedTemperatureMomentum {
  // compute (heat) transport
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    const V temp = TYPE().computeRho(cell);
    const auto uNS = cell.template getFieldPointer<descriptors::VELOCITY>();
    computeU<TYPE>(cell, j);
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] += temp * uNS[iD];
    }
  }

  // compute (heat) conduction
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
    constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

    const V temp = TYPE().computeRho(cell);
    const auto uNS = cell.template getField<descriptors::VELOCITY>();

    V uOnWall[DESCRIPTOR::d] = { V{} };
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
        uOnWall[iD] += (cell[onWallIndices[fIndex]]
                          - equilibrium<DESCRIPTOR>::firstOrder(onWallIndices[fIndex],temp,uNS.data()))
                          * descriptors::c<DESCRIPTOR>(onWallIndices[fIndex],iD);
      }
    }

    V uNormal[DESCRIPTOR::d] = { V{} };
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
        uNormal[iD] += (cell[normalIndices[fIndex]]
                          - equilibrium<DESCRIPTOR>::firstOrder(normalIndices[fIndex],temp,uNS.data()))
                          * descriptors::c<DESCRIPTOR>(normalIndices[fIndex],iD);
      }
    }

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = uOnWall[iD] + V(2) * uNormal[iD];
    }
  }

  // define the conduction
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {};

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FixedTemperatureMomentum";
  }
};


/** The first moment (the heat conduction) is fixed.
 * Implementation is identical to FixedVelocityMomentumGeneric, but the
 * compute(...) method is different: is computes the (heat) transport, similar
 * to FixedTemperatureMomentum.
 */
struct FixedVelocityMomentumAD {
  struct VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> { };

  // compute (heat) transport
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    const V temp = TYPE().computeRho(cell);
    const auto uNS = cell.template getFieldPointer<descriptors::VELOCITY>();
    computeU<TYPE>(cell, j);
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] += temp * uNS[iD];
    }
  }

  // compute (heat) conduction
  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    const auto uExt = cell.template getFieldPointer<VELOCITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = uExt[iD];
    }
  }

  // define the conduction
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    V u[DESCRIPTOR::d] = { V{} };
    cell.template setField<VELOCITY>(u);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FixedVelocityMomentumAD";
  }
};


template<int direction, int orientation>
struct FreeEnergyInletOutletMomentum {
  template <typename TYPE, typename CELL, typename J>
  void compute(CELL& cell, J& j) any_platform
  {
    FixedPressureMomentum<direction,orientation>().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = V(0);
    }
    u[direction] = cell.template getFieldPointer<descriptors::FORCE>()[1];
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    cell.template getFieldPointer<descriptors::FORCE>()[1] = u[direction];
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {
    FixedPressureMomentum<direction,orientation>().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform{};

  static std::string getName() {
    return "FreeEnergyInletOutletMomentum";
  }
};


struct FreeEnergyMomentum {
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    auto force = cell.template getField<descriptors::FORCE>();
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = force[iVel];
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {}

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {}

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "FreeEnergyMomentum";
  }
};

struct ForcedPSMMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    const V rho = TYPE().computeRho(cell);
    computeU<TYPE>(cell, j);
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      j[iVel] *= rho;
    }
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    V rho{};
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);
    V epsilon = 1. - cell.template getField<descriptors::POROSITY>();
    V omega = cell.template getField<descriptors::OMEGA>();
    V paramA = 1 / omega - 0.5;
    V paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    V paramC = (1. - paramB);
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * (u[iVel] + cell.template getFieldComponent<descriptors::FORCE>(iVel) * 0.5)
              + paramB * cell.template getFieldComponent<descriptors::VELOCITY_SOLID>(iVel);
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {}

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {}

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "ForcedPSMMomentum";
  }
};

struct GuoZhaoMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    lbm<DESCRIPTOR>::computeU(cell, u);

    const V epsilon = cell.template getField<descriptors::EPSILON>();
    const V nu      = cell.template getField<descriptors::NU>();
    const V k       = cell.template getField<descriptors::K>();

    auto bodyF = cell.template getFieldPointer<descriptors::BODY_FORCE>();

    for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      u[iDim] += 0.5*epsilon*bodyF[iDim];
    }

    const V uMag = util::sqrt( util::normSqr<V,DESCRIPTOR::d>(u) );
    const V Fe = 0.;//1.75/util::sqrt(150.*util::pow(epsilon,3));

    const V c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
    const V c_1 = 0.5*epsilon*Fe/util::sqrt(k);

    for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      u[iDim] /= (c_0 + util::sqrt(c_0*c_0 + c_1*uMag));
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform { }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform { }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "GuoZhaoMomentum";
  }
};


/// For offLattice boundary conditions
struct OffBoundaryMomentum {
  struct DISTANCES : public descriptors::FIELD_BASE<0, 0, 1> { };
  struct VELOCITY : public descriptors::FIELD_BASE<0, 0, 3> { };
  struct VELOCITY_COEFFICIENTS : public descriptors::FIELD_BASE<0, 0, 1> { };

  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] = V{};
    }
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    const auto distances = cell.template getFieldPointer<DISTANCES>();
    const auto velocities = cell.template getFieldPointer<VELOCITY>();

    for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
      u[iD] = V{};
    }
    unsigned counter = 0;
    for (int iPop = 0; iPop < DESCRIPTOR::q; iPop++) {
      if ( !util::nearZero(distances[iPop]+1) ) {
        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          u[iD] += velocities[3*iPop + iD];
        }
        counter++;
      }
    }
    if (counter!=0) {
      for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
        u[iD] /= counter;
      }
    }
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    const auto distances = cell.template getFieldPointer<DISTANCES>();
    auto velocities = cell.template getFieldPointer<VELOCITY>();
    auto velocityCoefficient = cell.template getFieldPointer<VELOCITY_COEFFICIENTS>();

    for (int iPop = 0; iPop < DESCRIPTOR::q; iPop++) {
      if ( !util::nearZero(distances[iPop]+1) ) {
        velocityCoefficient[iPop] = 0;
        // scalar product of c(iPop) and u
        for (int sum = 0; sum < DESCRIPTOR::d; sum++) { // +/- problem because of first stream than postprocess
          velocityCoefficient[iPop] -= descriptors::c<DESCRIPTOR>(iPop,sum)*u[sum];
        }
        // compute summand for boundary condition
        velocityCoefficient[iPop] *= 2*descriptors::invCs2<V,DESCRIPTOR>() * descriptors::t<V,DESCRIPTOR>(iPop);

        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          velocities[3 * iPop + iD] = u[iD];
        }
      }
    }
  }

  template <typename TYPE, typename CELL, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform
  {
    auto distances = cell.template getFieldPointer<DISTANCES>();
    for (int iPop = 0; iPop < DESCRIPTOR::q; iPop++) {
      distances[iPop] = -1;
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "OffBoundaryMomentum";
  }
};


/// Momentum computation for P1 dynamics.
struct P1Momentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    std::array<V,DESCRIPTOR::q> cellShifted;
    for (int iPop = 0; iPop <DESCRIPTOR::q; ++iPop) {
      cellShifted[iPop] = cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop);
    }
    std::array<V,DESCRIPTOR::d> moment1;
    moment1.fill( V(0) );
    // sum_j v_j f_j
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
        moment1[iDim] += descriptors::c<DESCRIPTOR>(iPop,iDim)*cellShifted[iPop];
      }
    }
    for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
      j[iDim] = moment1[iDim];
    }
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = V{};
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {}

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {}

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "P1Momentum";
  }
};


/// Momentum computation for Poisson dynamics.
struct PoissonMomentum {
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = V{};
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {}

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {}

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "PoissonMomentum";
  }
};


struct PorousGuoMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    lbm<DESCRIPTOR>::computeJ(cell, j);
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    lbm<DESCRIPTOR>::computeU(cell, u);

    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      u[iDim] *= porosity;
    }
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform { }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform { }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "PorousGuoMomentum";
  }
};


template<typename MOMENTUM>
struct PorousParticleMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    V rho{};
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);
    V u_tmp[3] = {0., 0., 0.};
    if (cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>()[0] > std::numeric_limits<V>::epsilon()) {
      for (int i=0; i<DESCRIPTOR::d; i++)  {
        u_tmp[i] = (V(1) - cell.template getField<descriptors::POROSITY>())
                  * ( cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[i]
                      / cell.template getField<descriptors::VELOCITY_DENOMINATOR>()
                      - u[i] );
        u[i] += V(0.5) * rho * u_tmp[i];
      }
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "PorousParticleMomentum<" + MOMENTUM().getName() + ">";
  }
};

/// Momentum of a moving porous medium
template <typename MOMENTUM>
struct MovingPorousMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    const V rho = TYPE().computeRho(cell);
    MOMENTUM().template computeU<TYPE>(cell, u);
    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] += V{0.5} * rho
             * (V{1} - porosity)
             * (cell.template getFieldComponent<descriptors::VELOCITY>(iD) - u[iD]);
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName(){
    return "MovingPorousMomentum<" + MOMENTUM().getName() + ">";
  }
};

template<typename STRESS>
struct MovingPorousStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    V uNew[DESCRIPTOR::d] { };
    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      uNew[iD] = u[iD] + V{0.5} * rho
                       * (V{1} - porosity)
                       * (cell.template getFieldComponent<descriptors::VELOCITY>(iD) - u[iD]);
    }

    STRESS().template compute<TYPE>(cell, rho, uNew, pi);
    V forceTensor[util::TensorVal<DESCRIPTOR>::n];
    // Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        forceTensor[iPi] = V{0.5} * rho * (uNew[iBeta] + uNew[iAlpha]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] += forceTensor[iPi];
    }
  }

  static std::string getName() {
    return "MovingPorousStress<" + STRESS().getName() + ">";
  }
};

template <typename MOMENTUM>
struct MovingPorousMomentumCombination {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    MOMENTUM().template computeU<TYPE>(cell, u);
    const V wmporosity = cell.template getField<descriptors::WMPOROSITY>();
    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = wmporosity*(porosity*u[iD] + (V{1} - porosity) * cell.template getFieldComponent<descriptors::VELOCITY>(iD)) + (V{1} - wmporosity)*cell.template getFieldComponent<descriptors::WMVELOCITY>(iD);
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName(){
    return "MovingPorousMomentumCombination<" + MOMENTUM().getName() + ">";
  }
};

template <typename MOMENTUM>
struct MovingPorousMomentumCombinationNoWM {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    MOMENTUM().template computeU<TYPE>(cell, u);
    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = porosity*u[iD] + (V{1} - porosity) * cell.template getFieldComponent<descriptors::VELOCITY>(iD);
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName(){
    return "MovingPorousMomentumCombinationNoWM<" + MOMENTUM().getName() + ">";
  }
};

/// Momentum is zero at solid material.
struct ZeroMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      j[iD] = V{};
    }
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      u[iD] = V{};
    }
  }

  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform {}

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform {}

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "ZeroMomentum";
  }
};


template<typename MOMENTUM>
struct ForcedMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    MOMENTUM().template computeU<TYPE>(cell, u);
    const auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * V(0.5);
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void inverseShift(CELL& cell, U& u) any_platform {
    const auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] -= force[iVel] * V(0.5);
    }
  }

  static std::string getName() {
    return "ForcedMomentum<" + MOMENTUM().getName() + ">";
  }
};

template<typename MOMENTUM>
struct GuoZhaoForcedMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    MOMENTUM().template computeU<TYPE>(cell, u);
    const V epsilon = cell.template getField<descriptors::EPSILON>();
    const V nu      = cell.template getField<descriptors::NU>();
    const V k       = cell.template getField<descriptors::K>();
    const auto bodyF = cell.template getFieldPointer<descriptors::BODY_FORCE>();

    const V uMag = util::sqrt( util::normSqr<V,DESCRIPTOR::d>(u) );
    const V Fe = 0.;//1.75/util::sqrt(150.*util::pow(epsilon,3));
    const V c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
    const V c_1 = 0.5*epsilon*Fe/util::sqrt(k);

    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += bodyF[iVel] * V(0.5) * epsilon;
      u[iVel] /= (c_0 + util::sqrt(c_0*c_0 + c_1*uMag));
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void inverseShift(CELL& cell, U& u) any_platform {
    const V epsilon = cell.template getField<descriptors::EPSILON>();
    const V nu      = cell.template getField<descriptors::NU>();
    const V k       = cell.template getField<descriptors::K>();
    const auto bodyF = cell.template getFieldPointer<descriptors::BODY_FORCE>();

    const V uMag = util::sqrt( util::normSqr<V,DESCRIPTOR::d>(u) );
    const V Fe = 0.;//1.75/util::sqrt(150.*util::pow(epsilon,3));
    const V c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
    const V c_1 = 0.5*epsilon*Fe/util::sqrt(k);

    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] *= (c_0 + util::sqrt(c_0*c_0 + c_1*uMag));
      u[iVel] -= bodyF[iVel] * V(0.5) * epsilon;
    }
  }

  static std::string getName() {
    return "GuoZhaoForcedMomentum<" + MOMENTUM().getName() + ">";
  }
};

template<typename MOMENTUM>
struct PorousMomentum {
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    MOMENTUM().template compute<TYPE>(cell, j);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    MOMENTUM().template computeU<TYPE>(cell, u);
    const V porosity = cell.template getField<descriptors::POROSITY>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] *= porosity;
    }
  }

  template <typename TYPE, typename CELL, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const U& u) any_platform
  {
    MOMENTUM().template define<TYPE>(cell, u);
  }

  template <typename TYPE, typename CELL>
  void initialize(CELL& cell) any_platform
  {
    MOMENTUM().template initialize<TYPE>(cell);
  }

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {}

  static std::string getName() {
    return "PorousMomentum<" + MOMENTUM().getName() + ">";
  }
};


// ---------------------- STRESS STRUCTURES -----------------------------------

/** Standard stress computation as second moment of the population.
 * Utilizes the implementation in lbHelpers.
 */
struct BulkStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    lbm<DESCRIPTOR>::computeStress(cell, rho, u, pi);
  }

  static std::string getName() {
    return "BulkStress";
  }
};

/** Standard stress computation as second moment of the population.
 * Utilizes the implementation in lbHelpers.
 */
struct IncompressibleBulkStress {
  template <typename TYPE, typename CELL, typename P, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const P& p, const U& u, PI& pi) any_platform
  {
    int iPi = 0;
    const auto rho = cell.template getField<descriptors::RHO>();
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        pi[iPi] = V();
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          pi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)*
                     descriptors::c<DESCRIPTOR>(iPop,iBeta) * cell[iPop];
        }
        // stripe off equilibrium contribution
        pi[iPi] -= rho*u[iAlpha]*u[iBeta];
        if (iAlpha==iBeta) {
          pi[iPi] -= p;
        }
        ++iPi;
      }
    }
  }

  static std::string getName() {
    return "IncompressibleBulkStress";
  }
};

/// Computation of the stress tensor for regularized boundary nodes
template <int direction, int orientation>
struct RegularizedBoundaryStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho, const U& u, PI& pi) any_platform
  {
    BoundaryHelpers<V,DESCRIPTOR,direction,orientation>::computeStress(
      cell, rho, u, pi);
  }

  static std::string getName() {
    return "RegularizedBoundaryStress";
  }
};

/// Computation of the stress tensor in an inner corner (2D case)
template <int normalX, int normalY>
struct InnerCornerStress2D {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, RHO& rho, const U& u, PI& pi) any_platform
  {
    PopulationCellD<V,DESCRIPTOR> newCell(
      cell.template getField<descriptors::POPULATION>());
    int v[DESCRIPTOR::d] = { -normalX, -normalY };
    int unknownF  = util::findVelocity<DESCRIPTOR >(v);

    if (unknownF != DESCRIPTOR::q) {
      int oppositeF = descriptors::opposite<DESCRIPTOR>(unknownF);

      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      newCell[unknownF] = newCell[oppositeF]
                          - equilibrium<DESCRIPTOR>::secondOrder(oppositeF, rho, u, uSqr)
                          + equilibrium<DESCRIPTOR>::secondOrder(unknownF, rho, u, uSqr);
    }

    lbm<DESCRIPTOR>::computeStress(newCell, rho, u, pi);
  }

  static std::string getName() {
    return "InnerCornerStress2D";
  }
};

/// Computation of the stress tensor in an inner corner (3D case)
template <int normalX, int normalY, int normalZ>
struct InnerCornerStress3D {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    PopulationCellD<V,DESCRIPTOR> newCell = cell.template getField<descriptors::POPULATION>();
    int v[DESCRIPTOR::d] = { -normalX, -normalY, -normalZ };
    int unknownF  = util::findVelocity<DESCRIPTOR >(v);

    if (unknownF != DESCRIPTOR::q) {
      int oppositeF = descriptors::opposite<DESCRIPTOR>(unknownF);

      V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

      newCell[unknownF] = newCell[oppositeF]
                          - equilibrium<DESCRIPTOR>::secondOrder(oppositeF, rho, u, uSqr)
                          + equilibrium<DESCRIPTOR>::secondOrder(unknownF, rho, u, uSqr);
    }

    lbm<DESCRIPTOR>::computeStress(newCell, rho, u, pi);
  }

  static std::string getName() {
    return "InnerCornerStress3D";
  }
};

/// Computation of the stress tensor in an inner edge
template <int plane, int normal1, int normal2>
struct InnerEdgeStress3D {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    PopulationCellD<V,DESCRIPTOR> newCell = cell.template getField<descriptors::POPULATION>();
    for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      if ( (descriptors::c<DESCRIPTOR>(iPop,(plane+1)%3) == -normal1) &&
          (descriptors::c<DESCRIPTOR>(iPop,(plane+2)%3) == -normal2) ) {
        int opp = descriptors::opposite<DESCRIPTOR>(iPop);
        newCell[iPop] = newCell[opp]
                        - equilibrium<DESCRIPTOR>::secondOrder(opp, rho, u, uSqr)
                        + equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      }
    }
    lbm<DESCRIPTOR>::computeStress(newCell, rho, u, pi);
  }

  static std::string getName() {
    return "InnerEdgeStress3D";
  }
};

/// Access to the stress computation is forbidden and raises an error.
struct NoStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    // TODO: Re-enable and fix downstream issue
    //throw std::bad_function_call();
  }

  static std::string getName() {
    return "NoStress";
  }
};

/// The stress is always zero.
struct ZeroStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] = V{};
    }
  }

  static std::string getName() {
    return "ZeroStress";
  }
};

/// The stress is fixed and stored in the external field PI.
struct FixedStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    auto piExt = cell.template getField<descriptors::TENSOR>();
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] = piExt[iPi];
    }
  }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    cell.template setField<descriptors::TENSOR>(pi);
  }

  static std::string getName() {
    return "FixedStress";
  }
};

template<typename STRESS>
struct ForcedStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    V uNew[DESCRIPTOR::d] { };
    const auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      uNew[iD] = u[iD] - V{0.5} * force[iD];
    }
    STRESS().template compute<TYPE>(cell, rho, uNew, pi);
    V forceTensor[util::TensorVal<DESCRIPTOR>::n];
    // Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        forceTensor[iPi] = V{0.5} * rho * (force[iAlpha]*uNew[iBeta] + uNew[iAlpha]*force[iBeta]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] += forceTensor[iPi];
    }
  }

  static std::string getName() {
    return "ForcedStress<" + STRESS().getName() + ">";
  }
};

template<typename STRESS>
struct ForcedIncompressibleStress {
  template <typename TYPE, typename CELL, typename P, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const P& p, const U& u, PI& pi) any_platform
  {
    V uNew[DESCRIPTOR::d] { };
    const auto force = cell.template getFieldPointer<descriptors::FORCE>();
    const auto partialRho = cell.template getFieldPointer<descriptors::NABLARHO>();
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      uNew[iD] = u[iD] - V{0.5} * force[iD];
    }
    STRESS().template compute<TYPE>(cell, p, uNew, pi);
    V forceTensor[util::TensorVal<DESCRIPTOR>::n];
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        forceTensor[iPi] = V{0.5} /descriptors::invCs2<V,DESCRIPTOR>() * (partialRho[iAlpha]*uNew[iBeta] + uNew[iAlpha]*partialRho[iBeta]);
        if (iAlpha == iBeta) forceTensor[iPi] += V{0.5} /descriptors::invCs2<V,DESCRIPTOR>() * partialRho[iAlpha]*uNew[iAlpha];
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] += forceTensor[iPi];
    }
  }

  static std::string getName() {
    return "ForcedIncompressibleStress<" + STRESS().getName() + ">";
  }
};

template<typename STRESS>
struct GuoZhaoForcedStress {
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    V uNew[DESCRIPTOR::d] { };
    const V epsilon = cell.template getField<descriptors::EPSILON>();
    const V nu      = cell.template getField<descriptors::NU>();
    const V k       = cell.template getField<descriptors::K>();
    const auto bodyF = cell.template getFieldPointer<descriptors::BODY_FORCE>();

    const V uMag = util::sqrt( util::normSqr<V,DESCRIPTOR::d>(u) );
    const V Fe = 0.;//1.75/util::sqrt(150.*util::pow(epsilon,3));
    const V c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
    const V c_1 = 0.5*epsilon*Fe/util::sqrt(k);

    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      uNew[iD] = u[iD] * (c_0 + util::sqrt(c_0*c_0 + c_1*uMag)) - V{0.5} * bodyF[iD] * epsilon;
    }
    STRESS().template compute<TYPE>(cell, rho, uNew, pi);
    V forceTensor[util::TensorVal<DESCRIPTOR>::n];
    // Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        forceTensor[iPi] = V{0.5} * rho * (bodyF[iAlpha]*uNew[iBeta] + uNew[iAlpha]*bodyF[iBeta]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
      pi[iPi] += forceTensor[iPi];
    }
  }

  static std::string getName() {
    return "GuoZhaoForcedStress<" + STRESS().getName() + ">";
  }
};

}  // namespace momenta

}  // namespace olb

#endif
