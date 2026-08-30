/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcek, Davide Dapelo
 *                2022 Nando Suntoyo, Adrian Kummerlaender
 *                2024 Shota Ito
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
 * BGK Dynamics with adjusted omega -- header file.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 *
 * Carrau-Yasuda model similar to:  "https://doi.org/10.1016/j.camwa.2009.02.021"
 * Casson model similar to:         "https://doi.org/10.1063/1.2772250"
 *
 */

#ifndef NON_NEWTONIAN_BGK_DYNAMICS_H
#define NON_NEWTONIAN_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

// *INDENT-OFF*

namespace visco {

// Model parameters for the Carrau-Yasuda model
struct N : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct MU_INF : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct MU_ZERO : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct LAMBDA : public descriptors::FIELD_BASE<1> { };
struct A : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

// Model parameters for the Casson model
struct K_ZERO : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};
struct K_ONE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
 };


/// Compute and update cell-wise OMEGA using viscosity model MODEL
template <typename COLLISION, typename MODEL>
struct OmegaFromCell {
  using parameters = meta::merge<typename COLLISION::parameters, typename MODEL::parameters>;

  static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
    "Assumption violated: COLLISION needs to contain OMEGA");

  static constexpr bool is_vectorizable = false;

  static std::string getName()
  {
    return "visco::OmegaFromCell<" + COLLISION::getName() + ", " + MODEL::getName()  + ">";
  }

  template <typename CELL, typename PARAMETERS, typename OMEGA, typename RHO, typename PI, typename V=typename CELL::value_t>
  static V computeOmega(CELL& cell, PARAMETERS& parameters, OMEGA& omega0, RHO& rho, PI& pi) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    V pre2 = V{0.5} * descriptors::invCs2<V,DESCRIPTOR>() * omega0 / rho; // strain rate tensor prefactor
    pre2 *= pre2;
    V gamma{};
    if constexpr (DESCRIPTOR::template provides<descriptors::SHEAR_RATE_MAGNITUDE>()) {
      gamma = cell.template getField<descriptors::SHEAR_RATE_MAGNITUDE>();
    }
    else {
      if constexpr (DESCRIPTOR::template provides<descriptors::FORCE>()) {
        // Cannot be done in just one line, it gives error - I don't know why. Davide Dapelo
        const auto force = cell.template getField<descriptors::FORCE>();
        gamma = util::sqrt(V{2}*pre2*lbm<DESCRIPTOR>::computePiNeqNormSqr(cell, force));
      }
      else {
        gamma = util::sqrt(V{2}*pre2*lbm<DESCRIPTOR>::computePiNeqNormSqr(cell));
      }
    }

    V nuNew = MODEL::computeViscosity(parameters, gamma);
    V newOmega = V{1} / (nuNew*descriptors::invCs2<V,DESCRIPTOR>() + V{0.5});
    return newOmega;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = false;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
    {
      V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeAllMomenta(cell, rho, u, pi);
      const V oldOmega = cell.template getField<descriptors::OMEGA>();
      const V newOmega = computeOmega(cell, parameters, oldOmega, rho, pi);
      cell.template setField<descriptors::OMEGA>(newOmega);
      parameters.template set<descriptors::OMEGA>(newOmega);
      return CollisionO().apply(cell, parameters);
    }
  };
};

/// Viscosity model Carreau-Yasuda
struct CarreauYasuda
{
  // MODEL::parameters containing only used parameters
  using parameters = meta::list<N,MU_INF,MU_ZERO,LAMBDA,A>;

  static std::string getName()
  {
    return "visco::Carreau-Yasuda";
  }

  // Compute the dynamic viscosity from the shear rate
  template <typename T, typename PARAMETERS>
  static T computeViscosity(PARAMETERS& parameters, T gamma) any_platform
  {
    T n = parameters.template get<N>();
    T mu_inf = parameters.template get<MU_INF>();
    T mu_0 = parameters.template get<MU_ZERO>();
    T lambda = parameters.template get<LAMBDA>();
    T a = parameters.template get<A>();
    return mu_inf + (mu_0 - mu_inf) * util::pow( 1. + util::pow(lambda * gamma, a), (n-1.)/(a)); // Carreau-Yasuda model
  }
};

/// Viscosity model Casson
struct Casson
{
  // MODEL::parameters containing only used parameters
  using parameters = meta::list<K_ZERO,K_ONE>;

  static std::string getName()
  {
    return "visco::Casson";
  }

  // Compute the dynamic viscosity from the shear rate
  template <typename T, typename PARAMETERS>
  static T computeViscosity(PARAMETERS& parameters, T gamma) any_platform
  {
    T k0 = parameters.template get<K_ZERO>();
    T k1 = parameters.template get<K_ONE>();
    if (std::numeric_limits<T>::epsilon()*T(1e1) > gamma) {// to avoid numerical instability
      gamma = std::numeric_limits<T>::epsilon()*T(1e1);
    }
    return (1./(gamma)) * util::pow( k0 + k1 * util::sqrt(gamma) , 2.); // Casson model
  }
};

}


/// BGK collision using Carrau-Yasuda viscosity model
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using CarreauYasudaBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  visco::OmegaFromCell<collision::BGK,visco::CarreauYasuda>
>;

/// BGK collision using Carrau-Yasuda viscosity model with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using CarreauYasudaForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  visco::OmegaFromCell<collision::BGK,visco::CarreauYasuda>,
  forcing::Guo<momenta::ForcedWithStress>
>;

/// BGK collision using Casson viscosity model
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using CassonBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  visco::OmegaFromCell<collision::BGK,visco::Casson>
>;

/// BGK collision using Casson viscosity model with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using CassonForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  visco::OmegaFromCell<collision::BGK,visco::Casson>,
  forcing::Guo<momenta::ForcedWithStress>
>;

// *INDENT-ON*

}

#endif
