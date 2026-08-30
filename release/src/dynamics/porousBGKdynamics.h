/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender, Nando Suntoyo
 *
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

#ifndef POROUS_BGK_DYNAMICS_H
#define POROUS_BGK_DYNAMICS_H

#include "interface.h"
#include "collision.h"
#include "dynamics/collisionLES.h"

namespace olb {

/// Porous BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::BGK
>;


namespace particles {

template <typename DESCRIPTOR, typename CELL,
          typename V = typename CELL::value_t>
inline void resetParticleRelatedFields(CELL& cell) noexcept
{
  if constexpr (DESCRIPTOR::template provides<descriptors::POROSITY>()) {
    cell.template setField<descriptors::POROSITY>(
        descriptors::POROSITY::template getInitialValue<V,DESCRIPTOR>() );
  }
  if constexpr (DESCRIPTOR::template provides<
                    descriptors::VELOCITY_DENOMINATOR>()) {
    cell.template setField<descriptors::VELOCITY_DENOMINATOR>(
        descriptors::VELOCITY_DENOMINATOR::template getInitialValue<V,DESCRIPTOR>() );
  }
  if constexpr (DESCRIPTOR::template provides<
                    descriptors::VELOCITY_NUMERATOR>()) {
    cell.template setField<descriptors::VELOCITY_NUMERATOR>(
        descriptors::VELOCITY_NUMERATOR::template getInitialValue<V,DESCRIPTOR>() );
  }
  if constexpr (DESCRIPTOR::template provides<descriptors::VELOCITY_SOLID>()) {
    cell.template setField<descriptors::VELOCITY_SOLID>(
        descriptors::VELOCITY_SOLID::template getInitialValue<V,DESCRIPTOR>() );
  }
}

template <typename DESCRIPTOR, typename CELL,
          typename V = typename CELL::value_t>
inline void resetParticleContactRelatedFields(CELL& cell) noexcept
{
  if constexpr (DESCRIPTOR::template provides<
                    descriptors::CONTACT_DETECTION>()) {
    cell.template setField<descriptors::CONTACT_DETECTION>(
        descriptors::CONTACT_DETECTION::template getInitialValue<V,DESCRIPTOR>() );
  }
}

template <typename DESCRIPTOR, typename CELL,
          typename V = typename CELL::value_t>
inline void resetAllParticleRelatedFields(CELL& cell) noexcept
{
  resetParticleRelatedFields<DESCRIPTOR,CELL,V>(cell);
  resetParticleContactRelatedFields<DESCRIPTOR,CELL,V>(cell);
}

} // namespace particles

namespace collision {

/* Implementation of the BGK collision for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field,
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
template <typename COLLISION, bool isStatic=false>
struct PorousParticle {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static constexpr bool is_vectorizable = false;

  static std::string getName() {
    return "PorousParticle<" + COLLISION::getName() + ">" ;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    constexpr static bool is_vectorizable = false;

    template <typename CELL, typename pVELOCITY, typename V=typename CELL::value_t>
    void calculate (CELL& cell, pVELOCITY& pVelocity) {
      if constexpr (isStatic) {
        for (int i=0; i<DESCRIPTOR::d; i++)  {
          pVelocity[i] -= (1.-(cell.template getField<descriptors::POROSITY>())) * pVelocity[i];
        }
      } else {
        for (int i=0; i<DESCRIPTOR::d; i++)  {
          pVelocity[i] += (1.-cell.template getField<descriptors::POROSITY>())
                        * (cell.template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(i)
                        / cell.template getField<descriptors::VELOCITY_DENOMINATOR>() - pVelocity[i]);
        }
      }
    }

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);

      auto statistic = CollisionO().apply(cell, parameters);

      V velDenominator = cell.template getFieldComponent<descriptors::VELOCITY_DENOMINATOR>(0);

      // use Kuperstokh forcing by default
      V uPlus[DESCRIPTOR::d]{ };
      V diff[DESCRIPTOR::q]{ };
      V fEqPlus[DESCRIPTOR::q]{ };
      V fEq[DESCRIPTOR::q]{ };

      if (velDenominator > std::numeric_limits<V>::epsilon()) {
        for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
          uPlus[iDim] = u[iDim];
        }
        calculate(cell, uPlus);
        if constexpr (!isStatic) {
          particles::resetParticleRelatedFields<DESCRIPTOR,CELL,V>(cell);
        }

        EquilibriumF().compute(cell, rho, uPlus, fEqPlus);
        EquilibriumF().compute(cell, rho, u, fEq);
        for (int tmp_iPop=0; tmp_iPop < DESCRIPTOR::q; tmp_iPop++) {
          diff[tmp_iPop] += fEqPlus[tmp_iPop] - fEq[tmp_iPop];
          cell[tmp_iPop] += diff[tmp_iPop];
        }
      }

      particles::resetParticleContactRelatedFields<DESCRIPTOR,CELL,V>(cell);

      return statistic;
    }
  };
};

template <typename COLLISION>
struct MovingPorosity {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "MovingPorosity<" + COLLISION::getName() + ">" ;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      Vector<V,DESCRIPTOR::d> u{};
      MomentaF().computeU(cell, u);

      auto statistic = CollisionO().apply(cell, parameters);

      // use Kuperstokh forcing
      Vector<V,DESCRIPTOR::d> uPlus = u;
      uPlus += (V{1} - cell.template getField<descriptors::POROSITY>())
             * (cell.template getField<descriptors::VELOCITY>() - uPlus);

      // TODO replace by single EquilibriumF()::compute
      V fEq[DESCRIPTOR::q] { };
      V fEq2[DESCRIPTOR::q] { };
      EquilibriumF().compute(cell, statistic.rho, u, fEq);
      EquilibriumF().compute(cell, statistic.rho, uPlus, fEq2);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += fEq2[iPop] - fEq[iPop];
      }

      return statistic;
    }
  };
};

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template <typename COLLISION>
struct PSM {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static constexpr bool is_vectorizable = false;

  static std::string getName() {
    return "PSM<" + COLLISION::getName() + ">" ;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      const V epsilon = 1. - cell.template getField<descriptors::POROSITY>();
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V paramA = V{1.} / omega - V{0.5};
      MomentaF().computeRhoU(cell, rho, u);
      // velocity at the boundary
      auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();
      if (epsilon < 1e-5) {
        return CollisionO().apply(cell, parameters);
      }
      else {
        // speed up paramB
        V paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
        // speed up paramC
        V paramC = (1. - paramB);
        V omega_s[DESCRIPTOR::q];
        V cell_tmp[DESCRIPTOR::q];

        V fEq[DESCRIPTOR::q] { };
        V fEq_s[DESCRIPTOR::q] { };
        EquilibriumF().compute(cell, rho, u, fEq);
        EquilibriumF().compute(cell, rho, u_s, fEq_s);
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell_tmp[iPop] = cell[iPop];
          //// To be reimplemented sensibly in case one needs it
          //if constexpr (MODE == 2) {
            omega_s[iPop] = (fEq_s[iPop] - cell[iPop]) + (V{1} - omega) * (cell[iPop] - fEq[iPop]);
          //} else if constexpr (MODE == 3) {
          //  omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - fEq_s[descriptors::opposite<DESCRIPTOR>(iPop)])
          //                  - (cell[iPop] - fEq_s[iPop]);
          //}
        }
        CollisionO().apply(cell, parameters);
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
          cell[iPop] += paramB * omega_s[iPop];
        }
        for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
          u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
        }
      }

      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    }
  };
};

// Implementation of the BGK collision step for subgridscale particles
// extended collision step, computes local drag in a given direction
template <typename COLLISION>
struct SubgridParticle {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "SubgridParticle<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V porosity = cell.template getField<descriptors::POROSITY>();
      auto extVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();
      //  if (porosity[0] != 0) {
      //    cout << "extVelocity: " << extVelocity[0] << " " <<  extVelocity[1] << " " <<  extVelocity[2] << " " << std::endl;
      //    cout << "porosity: " << porosity[0] << std::endl;
      //  }
      for (int i=0; i<DESCRIPTOR::d; i++)  {
        u[i] *= (1.-porosity);
        u[i] += extVelocity[i];
      }

      V uSqr = CollisionO().apply(cell, parameters);

      //statistics.incrementStats(rho, uSqr);
      cell.template setField<descriptors::POROSITY>(0);
      cell.template setField<descriptors::VELOCITY_NUMERATOR>(0);
      cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0);

      return {rho, uSqr};
    }
  };
};

// Implementation of the BGK collision for Zeta-Field formulation (Geng2019)
struct DBBParticleBGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static constexpr bool is_vectorizable = false;

  static std::string getName() {
    return "DBBParticleBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d],eta[DESCRIPTOR::q],uPlus[DESCRIPTOR::d],tmp_cell[(DESCRIPTOR::q+1)/2];
      MomentaF().computeRhoU(cell, rho, u);
      const V omega = parameters.template get<descriptors::OMEGA>();
      V tmpMomentumLoss[DESCRIPTOR::d]{ };
      auto velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
      auto zeta           = cell.template getFieldPointer<descriptors::ZETA>();
      auto velDenominator = cell.template getField<descriptors::VELOCITY_DENOMINATOR>();

      if (velDenominator > 1) {
        rho /= velDenominator;
      }
      for (int tmp_iPop=1; 2*tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
        eta[tmp_iPop]=6.*descriptors::t<V,DESCRIPTOR>(tmp_iPop)*rho*(descriptors::c<DESCRIPTOR>(tmp_iPop,0)*(velNumerator[0])+descriptors::c<DESCRIPTOR>(tmp_iPop,1)*(velNumerator[1]));
        tmp_cell[tmp_iPop]=(zeta[tmp_iPop]*(-cell[tmp_iPop]+cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]+eta[tmp_iPop]));
        cell[tmp_iPop]+=tmp_cell[tmp_iPop]/(1.+2.*(zeta[tmp_iPop]));
        cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]-=tmp_cell[tmp_iPop]/(1.+2.*(zeta[tmp_iPop]));
        zeta[tmp_iPop] = 0.;
        zeta[descriptors::opposite<DESCRIPTOR>(tmp_iPop)] = 0.;
      }

      cell.template setField<descriptors::POROSITY>(1.);
      cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0);

      MomentaF().computeRhoU(cell, rho, uPlus);
      V uPlusSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, uPlus, omega);

      V diff[DESCRIPTOR::q] = {};
      V fEqPlus[DESCRIPTOR::q] = {};
      V fEq[DESCRIPTOR::q] = {};
      EquilibriumF().compute(cell, rho, uPlus, fEqPlus);
      EquilibriumF().compute(cell, rho, u, fEq);
      for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
        diff[tmp_iPop] += fEqPlus[tmp_iPop] - fEq[tmp_iPop];

        for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
          tmpMomentumLoss[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
        }
      }

      for (int i_dim=0; i_dim<DESCRIPTOR::d; i_dim++) {
        velNumerator[i_dim] = tmpMomentumLoss[i_dim];
      }

      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    }
  };
};

// Implementation of the HBGK collision step for a porosity model enabling
// drag computation for many particles
// including the Krause turbulence modell
template <typename COLLISION>
struct KrauseH {
  using parameters = meta::list<descriptors::OMEGA, LES::SMAGORINSKY>;

  static std::string getName(){
    return "KrauseH<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      V newOmega[DESCRIPTOR::d];
      V _fieldTmp[4];
      _fieldTmp[0] = V{1};
      _fieldTmp[1] = V{};
      _fieldTmp[2] = V{};
      _fieldTmp[3] = V{};

      const V smagoConst = parameters.template get<collision::LES::SMAGORINSKY>();
      V preFactor = smagoConst*smagoConst
                    * descriptors::invCs2<V,DESCRIPTOR>()*descriptors::invCs2<V,DESCRIPTOR>()
                    * 2*util::sqrt(2);

      MomentaF().computeRhoU(cell, rho, u);

      // compute newOmega
      V omega = parameters.template get<descriptors::OMEGA>();
      V uSqr = u[0]*u[0];
      for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
        uSqr += u[iDim]*u[iDim];
      }
      /// Molecular realaxation time
      V tau_mol = 1./omega;
      V fEq[DESCRIPTOR::q] { };
      EquilibriumF().compute(cell, rho, u, fEq);
      for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
        V fNeq = util::fabs(cell[iPop] - fEq[iPop]);
        /// Turbulent realaxation time
        V tau_turb = 0.5*(util::sqrt(tau_mol*tau_mol+(preFactor*fNeq))-tau_mol);
        /// Effective realaxation time
        V tau_eff = tau_mol + tau_turb;
        newOmega[iPop] = 1./tau_eff;
      }
      parameters.template set<descriptors::OMEGA>(newOmega);


      V vel_denom = cell.template getField<descriptors::VELOCITY_DENOMINATOR>();
      if (vel_denom > std::numeric_limits<V>::epsilon()) {
        V porosity = cell.template getField<descriptors::POROSITY>(); // prod(1-smoothInd)
        auto vel_num = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
        porosity = 1.-porosity; // 1-prod(1-smoothInd)
        for (int i=0; i<DESCRIPTOR::d; i++)  {
          u[i] += porosity * (vel_num[i] / vel_denom - u[i]);
        }
      }
      uSqr = CollisionO().apply(cell, parameters);

      cell.template setField<descriptors::POROSITY>(_fieldTmp[0]);
      cell.template setField<descriptors::VELOCITY_NUMERATOR>({_fieldTmp[1], _fieldTmp[2]});
      cell.template setField<descriptors::VELOCITY_DENOMINATOR>(_fieldTmp[3]);

      return {rho, uSqr};
    }
  };
};

/// Implementation of the BGK collision step for a small particles enabling
/// two way coupling
template <typename COLLISION>
struct SmallParticle {
  using parameters = meta::list<descriptors::OMEGA, LES::SMAGORINSKY>;

  static std::string getName(){
    return "SmallParticle<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V porosity = cell.template getField<descriptors::POROSITY>();
      auto localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

      //cout << porosity[0]  << " " <<   localVelocity[0]<< " " <<   localVelocity[1]<< " " <<   localVelocity[2]<<std::endl;
      for (int i=0; i<DESCRIPTOR::d; i++)  {
        u[i] *= porosity;
        u[i] += localVelocity[i];
      }
      V uSqr = CollisionO().apply(cell, parameters);

      return {rho, uSqr};
    }
  };
};

} // namespace collision

namespace dynamics {

struct ExposePorousParticleMomenta {
  static std::string getName() {
    return "ExposePorousParticleMomenta";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::PorousParticle<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}


/// Porous particle BGK collision for moving particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,false>,
  dynamics::ExposePorousParticleMomenta
>;

/// Porous particle BGK collision for static particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using StaticPorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,true>,
  dynamics::ExposePorousParticleMomenta
>;

/// BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyPorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::PorousParticle<collision::BGK>>
>;

// BGK collision step for subgridscale particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SubgridParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SubgridParticle<collision::BGK>
>;

// BGK collision for Zeta-Field formulation (Geng2019)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DBBParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::DBBParticleBGK
>;

/// HBGK collision step for a porosity model enabling drag computation for many particles
/// including the Krause turbulence modell
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using KrauseHBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::KrauseH<collision::BGK>
>;

/// BGK collision step for a small particles enabling two way coupling
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmallParticleBGKdymaics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmallParticle<collision::BGK>
>;

/// Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PSMBGKdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PSM<collision::BGK>
>;

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR>
struct ForcedPSMBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,momenta::ForcedPSMBulkTuple> {

  using MomentaF = typename momenta::ForcedPSMBulkTuple::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::SecondOrder::template type<DESCRIPTOR,momenta::ForcedPSMBulkTuple>;

  using parameters = meta::list<descriptors::OMEGA>;

  constexpr static bool is_vectorizable = false;

  template <typename NEW_T>
  using exchange_value_type = ForcedPSMBGKdynamics<NEW_T,DESCRIPTOR>;

  std::type_index id() override {
    return typeid(ForcedPSMBGKdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ForcedPSMBGKdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform
  {
    V omega = parameters.template get<descriptors::OMEGA>();
    V rho, u[DESCRIPTOR::d], uSqr;
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);
    auto force = cell.template getField<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * V{0.5};
    }
    auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();

    V epsilon = V{1} - cell.template getField<descriptors::POROSITY>();
    if (epsilon < 1e-5) {
      uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
    }
    else {
      V paramA = V{1} / omega - V{0.5};
      V paramB = (epsilon * paramA) / ((V{1} - epsilon) + paramA);
      V paramC = (V{1} - paramB);

      V omega_s[DESCRIPTOR::q];
      V cell_tmp[DESCRIPTOR::q];

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell_tmp[iPop] = cell[iPop];
        omega_s[iPop] = (equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u_s) - cell[iPop])
                      + (V{1} - omega) * (cell[iPop] - equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u));
      }

      uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
        cell[iPop] += paramB * omega_s[iPop];
      }
      for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
        u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
      }
      uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    }
    return {rho, uSqr};
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "ForcedPSMBGKdynamics<" + MomentaF().getName() + ">";
  };

};

}

#endif
