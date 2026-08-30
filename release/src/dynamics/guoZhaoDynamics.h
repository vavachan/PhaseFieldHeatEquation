/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Specific dynamics classes for Guo and Zhao (2002) porous model, with
 * which a Cell object can be instantiated -- header file.
 */
#ifndef LB_GUOZHAO_DYNAMICS_H
#define LB_GUOZHAO_DYNAMICS_H

#include <type_traits>
#include <functional>

#include "momenta/interface.h"
#include "interface.h"
#include "core/util.h"

namespace olb {

namespace guoZhao {


template <typename DESCRIPTOR>
struct guoZhao_equilibrium {
  /// Computation of equilibrium distribution, second order in u
  template <typename RHO, typename U, typename USQR, typename EPSILON, typename V=RHO>
  static V secondOrder(int iPop, const RHO& rho, const U& u, const USQR& uSqr, const EPSILON& epsilon) any_platform
  {
    V c_u{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho
           * descriptors::t<V,DESCRIPTOR>(iPop)
           * ( V{1}
               + descriptors::invCs2<V,DESCRIPTOR>() * c_u
               + descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * V{0.5} * c_u * c_u / epsilon
               - descriptors::invCs2<V,DESCRIPTOR>() * V{0.5} * uSqr / epsilon)
           - descriptors::t<V,DESCRIPTOR>(iPop);
  }
};

template <typename DESCRIPTOR>
struct guoZhao_lbm {
  /// Add a force term after BGK collision
  template <typename CELL, typename RHO, typename U, typename OMEGA, typename EPSILON, typename K, typename NU, typename FORCE, typename V=typename CELL::value_t>
  static void addExternalForce(CELL& cell, const RHO& rho, const U& u, const OMEGA& omega, const EPSILON& epsilon, const K& k, const NU& nu, const FORCE& force) any_platform
  {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V c_u = V();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      c_u *= descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()/epsilon;
      V forceTerm = V();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   (epsilon*(V)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<V,DESCRIPTOR>()/epsilon
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<V,DESCRIPTOR>(iPop);
      forceTerm *= V(1) - omega/V(2);
      forceTerm *= rho;
      cell[iPop] += forceTerm;
    }
  }

  /// Updates the force terms with the reaction from the porous matrix, before applying it
  template <typename CELL, typename U, typename EPSILON, typename K, typename NU, typename FORCE, typename BODY_F, typename V=typename CELL::value_t>
  static void updateExternalForce(CELL& cell, const U& u, const EPSILON& epsilon, const K& k, const NU& nu, FORCE& force, BODY_F& bodyF) any_platform
  {
    const V uMag = util::sqrt( util::normSqr<V,DESCRIPTOR::d>(u) );
    const V Fe = 0;//1.75/util::sqrt(150.*util::pow(epsilon,3));

    // Linear Darcy term, nonlinear Forchheimer term and body force
    for (int iDim=0; iDim <DESCRIPTOR::d; iDim++) {
      force[iDim] = -u[iDim]*epsilon*nu/k - epsilon*Fe/util::sqrt(k)*uMag*u[iDim] + bodyF[iDim]*epsilon;
    }
  }
};

struct GuoZhaoSecondOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "GuoZhaoSecondOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      const V epsilon  = cell.template getField<descriptors::EPSILON>();
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = guoZhao_equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr, epsilon);
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

template <template <typename> typename Forced = momenta::Forced>
struct GuoZhaoForcing {
  static std::string getName() {
    return "GuoEpsilonForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V epsilon  = cell.template getField<descriptors::EPSILON>();
      const V k        = cell.template getField<descriptors::K>();
      const V nu       = cell.template getField<descriptors::NU>();
      auto force       = cell.template getField<descriptors::FORCE>();
      auto bodyF       = cell.template getField<descriptors::BODY_FORCE>();
      guoZhao_lbm<DESCRIPTOR>::updateExternalForce(cell, u,epsilon, k, nu, force, bodyF);
      guoZhao_lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, epsilon, k, nu, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using GuoZhaoBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  guoZhao::GuoZhaoSecondOrder,
  collision::BGK,
  guoZhao::GuoZhaoForcing<momenta::GuoZhaoForced>
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyGuoZhaoBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  guoZhao::GuoZhaoSecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>,
  guoZhao::GuoZhaoForcing<momenta::GuoZhaoForcedWithStress>
>;


}

#endif
