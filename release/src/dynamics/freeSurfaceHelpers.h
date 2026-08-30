/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
 *                2024-2025 Danial Khazaeipoul
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

#pragma once

#include "descriptor/fields.h"

#include "mrt.h"
#include "collisionLES.h"
#include "collisionMRT.h"

namespace olb {

namespace FreeSurface {

struct Stage0 {};
struct Stage1 {};
struct Stage2 {};
struct Stage3 {};
struct Stage4 {};
struct Stage5 {};

template<typename T>
constexpr T zeroThreshold() {
  if constexpr (std::is_same<T, float>::value) {
    return T{1e-6};
  } else {
    return T{1e-14};
  }
}

template <typename T>
inline bool isNearZero(T a, T threshold) {
  return util::abs(a) < threshold;
}

enum class SolverType {
  rowPivotingLU,
  completePivotingLU,
  columnPivotingQR
};
constexpr SolverType solver = SolverType::completePivotingLU;

enum class Type : std::uint8_t {
  Gas = 0,
  Interface = 1,
  Fluid = 2,
  Solid = 4
};

enum class Flags : std::uint8_t {
  None = 0,
  ToGas = 1,
  ToFluid = 2,
  NewInterface = 4
};

struct NeighbourInfo {
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;
  size_t interface_neighbours = 0;
};

// Global data structures storing: cell type, cell flags, mass, epsilon, velocity, and mass exchange
struct CELL_TYPE          : public descriptors::TYPED_FIELD_BASE<Type, 1> {};
struct CELL_FLAGS         : public descriptors::TYPED_FIELD_BASE<Flags, 1> {};
struct MASS               : public descriptors::FIELD_BASE<1> {};
struct EPSILON            : public descriptors::FIELD_BASE<1> {};
struct PREVIOUS_VELOCITY  : public descriptors::FIELD_BASE<0, 1, 0> {};
struct TEMP_MASS_EXCHANGE : public descriptors::FIELD_BASE<0, 0, 1> {};

// Parameters used in post processors
struct DROP_ISOLATED_CELLS       : public descriptors::FIELD_BASE<1> {};
struct TRANSITION                : public descriptors::FIELD_BASE<1> {};
struct LONELY_THRESHOLD          : public descriptors::FIELD_BASE<1> {};
struct HAS_SURFACE_TENSION       : public descriptors::FIELD_BASE<1> {};
struct SURFACE_TENSION_PARAMETER : public descriptors::FIELD_BASE<1> {};

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T, DESCRIPTOR>& sLattice);

template <typename CELL>
static bool isCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL>
static void setCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL>
static bool hasNeighbour(CELL& cell, const FreeSurface::Type& type) any_platform;

template <typename CELL>
static bool hasCellFlags(CELL& cell, const FreeSurface::Flags& type) any_platform;

template <typename CELL>
static void setCellFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

template <typename CELL>
static bool hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

template <typename CELL>
static NeighbourInfo getNeighbourInfo(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static bool isHealthyInterface(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static void computeMassExcessWeights
(
  CELL& cell, const Vector<V, CELL::descriptor_t::d>& normal,
  const bool& enableAllInterfaces, Vector<V, CELL::descriptor_t::q>& weights
) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V getClampedEpsilon(const CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V getSmoothEpsilon(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V getClampedSmoothEpsilon(const CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static Vector<V, CELL::descriptor_t::d> computeInterfaceNormal(CELL& cell) any_platform;

template<typename T, int N>
static void iterativeRefinement
(
  const std::array<T, N * N>& A, const std::array<T, N * N>& M,
  std::array<T, N>& x, const std::array<T, N>& b, const int Nsol,
  int maxIter = int(1)
) any_platform;

template<typename T, int N>
static void rowPivotingLUSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) any_platform;

template<typename T, int N>
static void completePivotingLUSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) any_platform;

template<typename T, int N>
static void columnPivotingQRSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) any_platform;

template<typename T>
static T plicCubeReduced(const T& volume, const Vector<T, 3>& n) any_platform;

template<typename T>
static T plicCube(const T& volume, const Vector<T, 3>& n) any_platform;

template<typename T>
static T plicCubeInverse(const T& d0, const Vector<T, 3>& n) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V computeCurvaturePLIC2D(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V computeCurvaturePLIC3D(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V computeCurvatureFDM3D(CELL& cell) any_platform;

template <typename CELL, typename V = typename CELL::value_t>
static V computeCurvature(CELL& cell) any_platform;

} // namespace FreeSurface

namespace forcing {

// Adapted from Guo et al.'s forcing scheme, with support for free surface dynamics.
template <template <typename> typename Forced = momenta::Forced>
struct FreeSurfaceGuo {
  static std::string getName() {
    return "FreeSurfaceGuoForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    constexpr static bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V epsilon = FreeSurface::getClampedEpsilon(cell);
      const auto force = epsilon * cell.template getField<descriptors::FORCE>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

// Adapted from Ladd and Verberg forcing scheme, with support for free surface dynamics.
struct FreeSurfaceLaddVerberg {
  static std::string getName() {
    return "FreeSurfaceLaddVerbergForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      const auto statistic = CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      // While this duplication can be resolved using CSE it should be extracted into a helper
      V rt[DESCRIPTOR::q] { }; // relaxation times vector.
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rt[iPop] = descriptors::s<V,DESCRIPTOR>(iPop);
      }
      for (int iPop=0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
        rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
      }
      V invM_S[DESCRIPTOR::q][DESCRIPTOR::q]; // relaxation times matrix
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
          invM_S[iPop][jPop] = V{};
          for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
            if (kPop == jPop) {
              invM_S[iPop][jPop] += descriptors::invM<V,DESCRIPTOR>(iPop,kPop) * rt[kPop];
            }
          }
        }
      }
      const V epsilon = FreeSurface::getClampedEpsilon(cell);
      const auto force = epsilon * cell.template getField<descriptors::FORCE>();
      mrt<DESCRIPTOR>::addExternalForce(cell, rho, u, invM_S, force);
      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

} // namespace forcing

// BGK collision step with Guo forcing scheme, with support for free surface dynamics.
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using FreeSurfaceForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::FreeSurfaceGuo<momenta::Forced>
>;

// Smagorinsky BGK collision step with Guo forcing scheme, with support for free surface dynamics.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyFreeSurfaceForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>,
  forcing::FreeSurfaceGuo<momenta::ForcedWithStress>
>;

// Smagorinsky MRT collision step with Ladd-Verberg forcing scheme, with support for free surface dynamics.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyFreeSurfaceForcedMRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::MRT>,
  forcing::FreeSurfaceLaddVerberg
>;

inline any_platform FreeSurface::Flags operator&(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) & static_cast<std::uint8_t>(rhs));
}

inline any_platform FreeSurface::Flags operator|(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}

} // namespace olb

#include "freeSurfaceHelpers.hh"
