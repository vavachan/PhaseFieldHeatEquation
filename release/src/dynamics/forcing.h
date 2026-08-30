/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_FORCING_H
#define DYNAMICS_FORCING_H

#include "lbm.h"
#include "descriptor/fields.h"
#include "core/latticeStatistics.h"

namespace olb {

/// Dynamics combination rules for various forcing schemes
namespace forcing {

/// Dynamics combination rule implementing the forcing scheme by Guo et al.
template <template <typename> typename Forced = momenta::Forced>
struct Guo {
  static std::string getName() {
    return "GuoForcing";
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
      //cell.template setField<descriptors::VELOCITY>(u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Guo et al.
//template <template <typename> typename Sourced = momenta::Sourced>
struct AdeGuo {
  static std::string getName() {
    return "AdeGuoForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      const auto u = cell.template getField<descriptors::VELOCITY>();
      const V rho = MomentaF().computeRho(cell);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto source = cell.template getField<descriptors::SOURCE>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V sourceTerm{};
      sourceTerm = source;
      sourceTerm *= descriptors::t<V,DESCRIPTOR>(iPop);
      sourceTerm *= (V{1} - omega * V{0.5});
      cell[iPop] += sourceTerm;
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Guo et al. with barycentric velocity
template <template <typename> typename Forced = momenta::Identity>
struct MCGuo {
  static std::string getName() {
    return "MultiComponentGuoForcing";
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
      const auto force = cell.template getField<descriptors::FORCE>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

template <template <typename> typename Forced = momenta::Forced>
struct Liang {
  static std::string getName() {
    return "LiangForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    constexpr static bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>
      // hacky workaround until vectorizability is deduced by introspection
      && !std::is_same_v<MOMENTA, momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V u[DESCRIPTOR::d];
      MomentaF().computeU(cell, u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      const auto rho = cell.template getField<descriptors::RHO>();
      const auto gradRho = cell.template getField<descriptors::NABLARHO>();
      lbm<DESCRIPTOR>::addLiangForce(cell, rho, gradRho, u, omega, force);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

template <template <typename> typename Forced = momenta::Forced>
struct LiangTRT {
  static std::string getName() {
    return "LiangTRTForcing";
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

    constexpr static bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>
      // hacky workaround until vectorizability is deduced by introspection
      && !std::is_same_v<MOMENTA, momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V u[DESCRIPTOR::d];
      MomentaF().computeU(cell, u);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      const auto rho = cell.template getField<descriptors::RHO>();
      const auto gradRho = cell.template getField<descriptors::NABLARHO>();

      V fTerm[DESCRIPTOR::q], fTermPlus[DESCRIPTOR::q], fTermMinus[DESCRIPTOR::q];

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        V c_u{};
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        }
        fTerm[iPop] = 0.;
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          fTerm[iPop] += descriptors::c<DESCRIPTOR>(iPop,iD) * force[iD] * rho +c_u * descriptors::c<DESCRIPTOR>(iPop,iD) * gradRho[iD];
        }
        fTerm[iPop] *= descriptors::t<V,DESCRIPTOR>(iPop)*descriptors::invCs2<V,DESCRIPTOR>();
      }

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fTermPlus[iPop]  = V{0.5} * (fTerm[iPop] + fTerm[descriptors::opposite<DESCRIPTOR>(iPop)]);
        fTermMinus[iPop] = V{0.5} * (fTerm[iPop] - fTerm[descriptors::opposite<DESCRIPTOR>(iPop)]);
      }

      V tau = V{1} / omega;
      const V tau2 = parameters.template get<collision::ITRT::TAU_MINUS>();
      const auto nablaRho = cell.template getField<descriptors::NABLARHO>();
      V ratio = util::norm<DESCRIPTOR::d>(nablaRho) / parameters.template get<collision::ITRT::MAXNABLARHO>();
      const V omega2 = V{1} / (tau + ratio * (tau2 - tau));

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += ( V{1} - omega * V{0.5} ) * fTermPlus[iPop]
                    + ( V{1} - omega2 * V{0.5} ) * fTermMinus[iPop];
      }

      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

struct AllenCahn {
  static std::string getName() {
    return "AllenCahnForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V phi = MomentaF().computeRho(cell);
      CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      const auto u = cell.template getField<descriptors::VELOCITY>();
      const auto source = cell.template getField<descriptors::SOURCE>();
      lbm<DESCRIPTOR>::addAllenCahnForce(cell, omega, force);
      lbm<DESCRIPTOR>::addAllenCahnSource(cell, omega, source);
      return {phi, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};


struct WellBalancedCahnHilliard {
  static std::string getName() {
    return "WellBalancedCahnHilliardForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V phi = MomentaF().computeRho(cell);
      CollisionO().apply(cell, parameters);
      const auto source = cell.template getField<descriptors::SOURCE>();
      const auto u = cell.template getField<descriptors::VELOCITY>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        V c_c{};
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          c_c += descriptors::c<DESCRIPTOR>(iPop,iD)*descriptors::c<DESCRIPTOR>(iPop,iD);
        }
        V sourceTerm{};
        sourceTerm += -V{2} + 0.5*descriptors::invCs2<V,DESCRIPTOR>()*c_c;
        sourceTerm *= descriptors::t<V,DESCRIPTOR>(iPop)*source;
        cell[iPop] += sourceTerm;
      }
      return {phi, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Kupershtokh et al.
struct Kupershtokh {
  static std::string getName() {
    return "KupershtokhForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V u[DESCRIPTOR::d];
      MomentaF().computeU(cell, u);
      const auto statistic = CollisionO().apply(cell, parameters);
      const auto force = cell.template getField<descriptors::FORCE>();
      V uPlusDeltaU[DESCRIPTOR::d];
      for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
        uPlusDeltaU[iVel] = u[iVel] + force[iVel];
      }
      V fEq[DESCRIPTOR::q] { };
      V fEq2[DESCRIPTOR::q] { };
      EquilibriumF().compute(cell, statistic.rho, u, fEq);
      EquilibriumF().compute(cell, statistic.rho, uPlusDeltaU, fEq2);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += fEq2[iPop] - fEq[iPop];
      }
      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by A.J. Wagner, Phys. Rev. E (2006)
struct Wagner {
  static std::string getName() {
    return "WagnerForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      const V omega = cell.template getField<descriptors::OMEGA>();
      const auto force = cell.template getField<descriptors::FORCE>();
      const V d2_rho = cell.template getField<descriptors::SCALAR>();
      Vector<V,DESCRIPTOR::d * DESCRIPTOR::d> psi_w{};

      // Forcing scheme
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        V c_u{};
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        }
        c_u *= descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
        V forceTerm{};
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          forceTerm +=
            (   (descriptors::c<DESCRIPTOR>(iPop,iD) - u[iD]) * descriptors::invCs2<V,DESCRIPTOR>()
                + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
            ) * force[iD];
        }

        V delta_ij = 0;
        for (int n = 0; n < DESCRIPTOR::d; ++n) {
          for (int m = 0; m < DESCRIPTOR::d; ++m) {
            int ind = n*DESCRIPTOR::d + m;

            delta_ij = (n == m) ? 1 : 0;

            psi_w[ind] = (V(1.)-omega/V(4.))*force[n]*force[m] + omega/V(12.)*( d2_rho )/rho*delta_ij;

            // Computing the source term
            forceTerm += V(9.) * ( descriptors::c<DESCRIPTOR>(iPop)[n]
                * descriptors::c<DESCRIPTOR>(iPop)[m] - V(1./3.)*delta_ij ) / V(2.)
                * psi_w[ind];
          }
        }

        forceTerm *= descriptors::t<V,DESCRIPTOR>(iPop);
        forceTerm *= rho;
        cell[iPop] += forceTerm;
      }

      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Shan and Chen
class ShanChen {
private:
  template <typename EQUILIBRIUM>
  struct VelocityShiftedEquilibrium {
    static std::string getName() {
      return "ShanChen<" + EQUILIBRIUM::getName() + ">";
    }

    template <typename DESCRIPTOR, typename MOMENTA>
    struct type {
      using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
      using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

      template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
      CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
        V uSqr{};
        for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
          uSqr += u[iVel] * u[iVel];
        }
        EquilibriumF().compute(cell, rho, u, fEq);
        return {rho, uSqr};
      };

      template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
      CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
        V rho, u[DESCRIPTOR::d];
        MomentaF().computeRhoU(cell, rho, u);
        const auto force = cell.template getFieldPointer<descriptors::FORCE>();
        for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
          u[iVel] += force[iVel] / parameters.template get<descriptors::OMEGA>();
        }
        return compute(cell, rho, u, fEq);
      };
    };
  };

public:
  static std::string getName() {
    return "ShanChenForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename VelocityShiftedEquilibrium<EQUILIBRIUM>::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = combined_momenta<DESCRIPTOR,MOMENTA>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,VelocityShiftedEquilibrium<EQUILIBRIUM>>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d] { };
      MomentaF().computeRhoU(cell, rho, u);
      CollisionO().apply(cell, parameters);
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

struct LinearVelocity {
  static std::string getName() {
    return "LinearVelocityForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
      MomentaF().computeAllMomenta(cell, rho, u, pi);
      auto force = cell.template getFieldPointer<descriptors::FORCE>();
      constexpr int nDim = DESCRIPTOR::d;
      V forceSave[nDim];
      // adds a+Bu to force, where
      //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
      //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
      auto v = cell.template getFieldPointer<descriptors::V12>();
      for (int iDim=0; iDim<nDim; ++iDim) {
        forceSave[iDim] = force[iDim];
        force[iDim] += v[iDim];
        for (int jDim=0; jDim<nDim; ++jDim) {
          force[iDim] += v[jDim + iDim*nDim + nDim]*u[jDim];
        }
      }
      for (int iVel=0; iVel<nDim; ++iVel) {
        u[iVel] += force[iVel] / V{2.};
      }

      auto statistics = CollisionO().apply(cell, parameters);
      V newOmega = parameters.template get<descriptors::OMEGA>();
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, newOmega, force);
      // Writing back to froce fector
      for (int iVel=0; iVel<nDim; ++iVel) {
        force[iVel] = forceSave[iVel];
      }
      return statistics;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Homogenized LBM modelling moving porous media
struct HLBM {
  static std::string getName() {
    return "HLBM<Kupershtokh>";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template wrap_momentum<
    momenta::MovingPorousMomentum
  >::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      MomentaF().computeRhoU(cell, rho, u);

      const V porosity = cell.template getField<descriptors::POROSITY>();
      auto statistic = CollisionO().apply(cell, parameters);

      if (porosity < 1) {
        // This is Kuperstokh forcing
        Vector<V,DESCRIPTOR::d> uPlus = u;
        uPlus += (V{1} - porosity)
               * (cell.template getField<descriptors::VELOCITY>() - u);

        V fEq[DESCRIPTOR::q] { };
        V fEq2[DESCRIPTOR::q] { };
        EquilibriumF().compute(cell, statistic.rho, u, fEq);
        EquilibriumF().compute(cell, statistic.rho, uPlus, fEq2);
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] += fEq2[iPop] - fEq[iPop];
        }
        return {-1,-1};
      } else {
        return statistic;
      }
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

struct ForcedHLBM {
  static std::string getName() {
    return "ForcedHLBM<Kupershtokh>";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<
    typename MOMENTA::template wrap_momentum<
      momenta::MovingPorousMomentum
    >
  >::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      Vector<V,DESCRIPTOR::d> u{};
      MomentaF().computeU(cell, u);

      auto statistic = CollisionO().apply(cell, parameters);

      const V porosity = cell.template getField<descriptors::POROSITY>();
      const auto force = cell.template getField<descriptors::FORCE>();

      auto fHLBM = (1 - porosity) * (cell.template getField<descriptors::VELOCITY>() - u);

      // This is Kuperstokh forcing
      Vector<V,DESCRIPTOR::d> uPlus = u + force + fHLBM;

      V fEq[DESCRIPTOR::q] { };
      V fEq2[DESCRIPTOR::q] { };
      EquilibriumF().compute(cell, statistic.rho, u, fEq);
      EquilibriumF().compute(cell, statistic.rho, uPlus, fEq2);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += fEq2[iPop] - fEq[iPop];
      }

      if (porosity < 1) {
        return {-1,-1};
      } else {
        V rho{};
        MomentaF().computeRhoU(cell, rho, u);
        return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
      }
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

}

#endif
