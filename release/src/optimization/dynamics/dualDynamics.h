/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
 *                2024 Julius Jessberger, Shota Ito, Adrian Kummerlaender
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
 * The description of optimization algorthims -- header file.
 */
/** \file
 * A collection of dynamics classes for dual LB methods
 * (e.g. dual BGK) with which a Cell object can be
 * instantiated -- header file.
 */

#ifndef DUAL_DYNAMICS_H
#define DUAL_DYNAMICS_H

#include "utilities/vectorHelpers.h"

// All OpenLB code is contained in this namespace.
namespace olb {

namespace collision {

struct DualPorousBGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "DualPorousBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      // Revert for backwards-in-time propagation
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }

      // Preparation
      const auto pop_f =   cell.template getField<opti::F>();
      const auto dJdF  =   cell.template getFieldPointer<opti::DJDF>();
      const V d      =  cell.template getField<descriptors::POROSITY>();
      const V omega = parameters.template get<descriptors::OMEGA>();

      // Forward density and velocity
      V rho_f { };
      V u_f[DESCRIPTOR::d] { };
      //this->computeRhoU( pop_f.data(), rho_f, u_f);
      rho_f = V();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u_f[iD] = V();
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rho_f += pop_f[iPop];
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          u_f[iD] += pop_f[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
      }
      rho_f += (V)1;
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u_f[iD] /= rho_f;
      }

      // Adjoint equilibrium
      V pheq[DESCRIPTOR::q] { };
      for (int i=0; i < DESCRIPTOR::q; ++i) {
        pheq[i] = V{0};
        for (int j=0; j < DESCRIPTOR::q; ++j) {
          V feq_j = equilibrium<DESCRIPTOR>::secondOrder(j, rho_f, u_f)
              + descriptors::t<V,DESCRIPTOR>(j);
          V dot_ij = V();
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            dot_ij += ( descriptors::c<DESCRIPTOR>(j,iD) - d*u_f[iD] )
              * ( descriptors::c<DESCRIPTOR>(i,iD) - u_f[iD] );
          }
          pheq[i] += cell[j]*feq_j*( V{1} + descriptors::invCs2<V,DESCRIPTOR>()*d*dot_ij );
        }
        pheq[i] /= rho_f;
      }

      // Collision
      for (int i=0; i < DESCRIPTOR::q; ++i) {
        cell[i] = cell[i] - omega*( cell[i] - pheq[i] ) + dJdF[i];
      }

      // Statistics
      V rho_phi, u_phi[DESCRIPTOR::d];
      MomentaF().computeRhoU( cell, rho_phi, u_phi );
      V uSqr_phi = util::normSqr<V,DESCRIPTOR::d>( u_phi );

      // Undo revert
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }
      return {V(1) + cell[0], uSqr_phi};
    };
  };
};

struct DualForcedBGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "DualForcedBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }

      // Preparation
      auto dJdF = cell.template getFieldPointer<opti::DJDF>();
      auto f = cell.template getField<opti::F>();
      auto force = cell.template getFieldPointer<descriptors::FORCE>();
      const V omega = parameters.template get<descriptors::OMEGA>();

      //dualLbDynamicsHelpers<T,DESCRIPTOR>::computeRhoU(f.data(), rho_f, u_f);
      V rho { };
      V u[DESCRIPTOR::d] { };
      //this->computeRhoU( pop_f.data(), rho_f, u_f);
      rho = V();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] = V();
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rho += f[iPop];
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          u[iD] += f[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
      }
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] /= rho;
      }

      V eq_phi[DESCRIPTOR::q];
      V force_phi[DESCRIPTOR::q];

      // Force
      V F1_phi[DESCRIPTOR::q][DESCRIPTOR::q];
      V F2_phi[DESCRIPTOR::q][DESCRIPTOR::q];
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        V f_c = V();
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          f_c += force[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
        for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
          V sum = V();
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            sum += (f_c*descriptors::c<DESCRIPTOR>(iPop,iD)-force[iD]/(V)descriptors::invCs2<V,DESCRIPTOR>())*descriptors::c<DESCRIPTOR>(jPop,iD);
          }
          F1_phi[iPop][jPop] = descriptors::t<V,DESCRIPTOR>(iPop)*descriptors::invCs2<V,DESCRIPTOR>()*f_c;
          F2_phi[iPop][jPop] = descriptors::t<V,DESCRIPTOR>(iPop)*descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()*sum;
        }
      }

      // Collision preperation
      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        eq_phi[jPop] = V();
        force_phi[jPop] = V();
        V eq = V();
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          V eq_tmp = olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u) + descriptors::t<V,DESCRIPTOR>(iPop);
          V sum = V();
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            sum += (u[iD] - descriptors::c<DESCRIPTOR>(jPop,iD))*(u[iD] - descriptors::c<DESCRIPTOR>(iPop,iD));
          }
          eq = (descriptors::invCs2<V,DESCRIPTOR>()*sum + V(1))/rho*eq_tmp;
          eq_phi[jPop] += cell[iPop]*eq;
          force_phi[jPop] += cell[iPop]*(F1_phi[iPop][jPop] + (1.-.5*omega)*F2_phi[iPop][jPop]);
        }
      }

      // Collision
      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        cell[jPop] = cell[jPop] - omega*(cell[jPop]-eq_phi[jPop]) + force_phi[jPop] + dJdF[jPop];
      }

      // Incrementing statistic values for convergence
      V phi2 = V();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        phi2 += cell[iPop]*cell[iPop];
      }
      //cell.revert();
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }

      return {V(1) + cell[0], phi2};
    };
  };
};

}

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DualPorousBGKDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::DualPorousBGK
>;

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DualForcedBGKDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::DualForcedBGK
>;

} // namespace olb

#endif
