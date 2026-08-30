/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_DYNAMICS_H
#define DYNAMICS_DYNAMICS_H

#include "interface.h"
#include "modifiers.h"

#include "core/util.h"
#include "core/latticeStatistics.h"
#include "descriptor/descriptor.h"

#include "momenta/interface.h"
#include "momenta/aliases.h"

#include "collision.h"
#include "equilibrium.h"
#include "forcing.h"

#include "collisionModifiers.h"

// TODO: Check this later
#include "collisionLES.h"

namespace olb {

/// Dynamics for "dead cells" doing nothing
template <typename T, typename DESCRIPTOR>
using NoDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  // Return 1 for the 0th moment, 0 for all others
  momenta::Tuple<
    momenta::OneDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::None,
  collision::None
>;

/// Dynamics for "dead cells" doing nothing. Variant with density=0
template <typename T, typename DESCRIPTOR>
using NoDynamicsWithZero = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::None,
  equilibria::None,
  collision::None
>;

/// Dynamics for "dead cells" with fixed density
template <typename T, typename DESCRIPTOR>
using NoDynamicsWithFixedDensity = dynamics::Tuple<
  T, DESCRIPTOR,
  // Return 1 for the 0th moment, 0 for all others
  momenta::Tuple<
    momenta::FixedDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::None,
  collision::None
>;

/// Common BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using BGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  std::conditional_t<   (DESCRIPTOR::d == 3 && DESCRIPTOR::q == 7)
                     || (DESCRIPTOR::d == 2 && DESCRIPTOR::q == 5),
                        equilibria::FirstOrder,
                        equilibria::SecondOrder>,
  collision::BGK
>;

//BGK collision with third order equilibrium function
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ThirdOrderBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ThirdOrder,
  collision::BGK
>;

// HRR collision https://hal.science/hal-02114308
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ThirdOrderHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ThirdOrder,
  collision::HRR
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ThirdOrderRLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ThirdOrder,
  collision::RLBThirdOrder
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedThirdOrderRLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ThirdOrder,
  collision::RLBThirdOrder,
  forcing::GuoThirdOrder<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedThirdOrderHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ThirdOrder,
  collision::HRR,
  forcing::GuoThirdOrder<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR>
using ThirdOrderHRLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
  momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
  momenta::BulkStress,
  momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::RLBThirdOrder
>;

template <typename T, typename DESCRIPTOR>
using ForcedThirdOrderHRLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
  momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
  momenta::BulkStress,
  momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::RLBThirdOrder,
  forcing::GuoThirdOrder<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR>
using ThirdOrderHHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
  momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
  momenta::BulkStress,
  momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::HRR
>;

template <typename T, typename DESCRIPTOR>
using ForcedThirdOrderHHRRdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
  momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
  momenta::BulkStress,
  momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::HRR,
  forcing::GuoThirdOrder<momenta::Forced>
>;

/// Pressure-corrected BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ConstRhoBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ConstRhoBGK
>;

/// BGK collision step with external force (Guo)
/**
 * Guo Z, Zheng C, Shi B (2002) Discrete lattice effects on the
 * forcing term in the lattice Boltzmann method. Phys Rev E 65(4)
 * DOI: 10.1103/PhysRevE.65.046308
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::Guo<momenta::Forced>
>;

// BGK collision step with external force (Guo) for multiple component lattices
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::ExternalVelocityTuple>
using MultiComponentForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::MCGuo<momenta::Identity>
>;

/// BGK collision step with external force (Wagner)
/**
 * Wagner, A (2006) Thermodynamic consistency of liquid-gas
 * lattice Boltzmann simulations. Phys Rev E 74, 056703
 * DOI: 10.1103/PhysRevE.74.056703
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedWagnerBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>,
  forcing::Wagner
>;

/// BGK collision step with external force (Kupershtokh)
/**
 * Kupershtokh A, Medvedev D, Karpov D (2009) On equations
 * of state in a lattice Boltzmann method. Comput Math Appl 58(5)
 * DOI: 10.1016/j.camwa.2009.02.024
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedKupershtokhBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::Kupershtokh
>;

/// BGK collision step with external force (Shan and Chen)
/**
 * Shan X, Chen H (1993) Lattice Boltzmann model for simulating
 * flows with multiple phases and components. Phys Rev E. 47 (3)
 * DOI: 10.1103/PhysRevE.47.1815
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedShanChenBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::ShanChen
>;

/// Incompressible BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using IncBGKdynamics = dynamics::Tuple<T,DESCRIPTOR,MOMENTA,equilibria::Incompressible,collision::BGK>;

/// Incompressible BGK collision step with external force
/**
 * Using Guo forcing on incompressible BGK
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedIncBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Incompressible,
  collision::BGK,
  forcing::Guo<momenta::Forced>
>;

/// Incompressible BGK collision step with relaxation frequency 1 / TAU_EFF and external force
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalTauForcedIncBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Incompressible,
  collision::OmegaFromCellTauEff<collision::BGK>,
  forcing::Guo<momenta::Forced>
>;

/// Multi-phase incompressible BGK collision step with relaxation frequency 1 / TAU_EFF and external force
/**
 * Liang, H., Xu, J., Chen, J., Wang, H., Chai, Z., & Shi, B. (2018).
 * Phase-field-based lattice Boltzmann modeling of large-density-ratio two-phase flows.
 * Physical Review E, 97(3), 033309.
 * DOI: 10.1103/PhysRevE.97.033309
 **/
template <typename T, typename DESCRIPTOR>
using MultiPhaseIncompressbileBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>,
  equilibria::MPIncompressible,
  collision::OmegaFromCellTauEff<collision::BGK>,
  forcing::Liang<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR>
using MultiPhaseIncompressbileTRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>,
  equilibria::MPIncompressible,
  collision::OmegaFromCellTauEff<collision::TRT>,
  forcing::Liang<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR>
using MultiPhaseIncompressbileInterfaceTRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>,
  equilibria::MPIncompressible,
  collision::OmegaFromCellTauEff<collision::ITRT>,
  forcing::LiangTRT<momenta::Forced>
>;

template <typename T, typename DESCRIPTOR>
using MultiPhaseIncompressbileSmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>,
  equilibria::MPIncompressible,
  collision::OmegaFromCellTauEff<collision::IncompressibleSmagorinskyEffectiveOmega<collision::BGK>>,
  forcing::Liang<momenta::ForcedWithIncompressibleStress>
>;

template <typename T, typename DESCRIPTOR>
using MultiPhaseIncompressbileSmagorinskyTRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>,
  equilibria::MPIncompressible,
  collision::OmegaFromCellTauEff<collision::IncompressibleSmagorinskyEffectiveOmega<collision::TRT>>,
  forcing::LiangTRT<momenta::ForcedWithIncompressibleStress>
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::ExternalVelocityTuple>
using AllenCahnBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::BGK,
  forcing::AllenCahn
>;

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::ExternalVelocityTuple>
using WellBalancedCahnHilliardBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::CahnHilliardZerothOrder,
  collision::BGK,
  forcing::WellBalancedCahnHilliard
>;


/// Regularized BGK collision step
/**
 * This model is substantially more stable than plain BGK, and has roughly
 * the same efficiency. However, it cuts out the modes at higher Knudsen
 * numbers and can not be used in the regime of rarefied gases.
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using RLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::RLB
>;

/// Regularized BGK collision followed by any other Dynamics
/**
 * Commonly used to model boundaries consisting of a specific boundary
 * momenta and standard dynamics::Tuple-declared dynamics.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA>
class CombinedRLBdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename CORRECTED_DYNAMICS::EquilibriumF;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CombinedRLBdynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA
  >;

  template <typename M>
  using exchange_momenta = CombinedRLBdynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedRLBdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<CombinedRLBdynamics>>();
  }

  template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n];
    MomentaF().computeAllMomenta(cell,rho,u,pi);
    V fEq[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, rho, u, fEq);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop] + equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
    }

    return CORRECTED_DYNAMICS().collide(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "CombinedRLBdynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

/// TRT collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using TRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::TRT
>;

/// TRT collision step with external force
/**
 * Using Guo forcing on incompressible BGK
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedTRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::TRT,
  forcing::Guo<momenta::Forced>
>;

/// Bounce Back boundary dynamics
/**
 * This is a very popular way to implement no-slip boundary conditions,
 * due to the simplicity and due to the independence of the boundary's
 * orientation.
 */
template <typename T, typename DESCRIPTOR>
using BounceBack = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::FixedDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::SecondOrder,
  collision::Revert
>;

/// Bounce Back boundary dynamics with bulk density
/**
 * Different from BounceBack these dynamics compute the
 * 0th moment from the cell's populations instead of
 * applying a fixed value.
 **/
template <typename T, typename DESCRIPTOR>
using BounceBackBulkDensity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::SecondOrder,
  collision::Revert
>;

template <typename T, typename DESCRIPTOR>
using BounceBackIncompressible = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::FixedDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::MPIncompressible,
  collision::Revert
>;

template <typename T, typename DESCRIPTOR>
using BounceBackBulkDensityADE = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::FirstOrder,
  collision::Revert
>;

template <typename T, typename DESCRIPTOR>
using BounceBackBulkDensityWellBalanced = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::CahnHilliardZerothOrder,
  collision::Revert
>;

/// Bounce Back boundary dynamics with Nguyen-Ladd velocity correction
/**
 * This is a very popular way to implement no-slip boundary conditions,
 * due to the simplicity and due to the independence of the boundary's
 * orientation.
 */
template <typename T, typename DESCRIPTOR>
using BounceBackVelocity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::ExternalVelocityFixedDensityTuple,
  equilibria::SecondOrder,
  collision::NguyenLaddCorrection<collision::Revert>
>;

/**
 * Corresponds to macro Robin boundary, micro Fresnel surface
 * Motivated by Hiorth et al. 2008; doi 10.1002/fld.1822
 */
template <typename T, typename DESCRIPTOR>
using PartialBounceBack = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<momenta::FixedDensity,momenta::ZeroMomentum,momenta::ZeroStress,momenta::DefineSeparately>,
  equilibria::FirstOrder,
  collision::PartialBounceBack
>;

/// Poisson dynamics
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::PoissonTuple>
using PoissonDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ZerothOrder,
  collision::Poisson
>;

/// P1 dynamics
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::P1Tuple>
using P1Dynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::P1,
  collision::P1
>;

/// Models a density sink by enforcing a zero distribution on the cell
template <typename T, typename DESCRIPTOR>
struct ZeroDistributionDynamics final : public dynamics::CustomCollision<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >
> {
  using EquilibriumF = typename equilibria::None::template type<DESCRIPTOR,momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >>;

  using parameters = meta::list<>;

  template <typename NEW_T>
  using exchange_value_type = ZeroDistributionDynamics<NEW_T,DESCRIPTOR>;

  template <typename M>
  using exchange_momenta = ZeroDistributionDynamics<T,DESCRIPTOR>;

  std::type_index id() override {
    return typeid(ZeroDistributionDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ZeroDistributionDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop);
    }
    return {-1, -1};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = 0;
    }
  };

  std::string getName() const override {
    return "ZeroDistributionDynamics";
  };
};

template <typename T, typename DESCRIPTOR, typename MOMENTA>
using ChopardDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Chopard,
  collision::BGK
>;

/// First order equilibrium boundary dynamics
/**
 * Applies the first order equilibrium distribution on every
 * time step using fixed density and velocity moments.
 */
template <typename T, typename DESCRIPTOR>
using EquilibriumBoundaryFirstOrder = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::EquilibriumBoundaryTuple,
  equilibria::FirstOrder,
  collision::FixedEquilibrium
>;

/// Second order equilibrium boundary dynamics
/**
 * Applies the second order equilibrium distribution on every
 * time step using fixed density and velocity moments.
 */
template <typename T, typename DESCRIPTOR>
using EquilibriumBoundarySecondOrder = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::EquilibriumBoundaryTuple,
  equilibria::SecondOrder,
  collision::FixedEquilibrium
>;

/// VANS BGK collision step with external force
/**
 * WIP: computeVANSRhoU doesn't exist (and should not)
 * Forcing to be realized using standard combinations rules (Likely Guo in this case)
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
struct ForcedVANSBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  using parameters = meta::list<descriptors::OMEGA>;

  template <typename NEW_T>
  using exchange_value_type = ForcedVANSBGKdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  std::type_index id() override {
    return typeid(ForcedVANSBGKdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ForcedVANSBGKdynamics>>();
  }

  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    T rho;
    T porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * T{0.5};
    }
  }

  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    T porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * T{0.5};
    }
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) {
    V rho, u[DESCRIPTOR::d];
    V porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * V{0.5};
      u[iVel] /= porosity;
    }
    rho *= porosity;
    const auto omega = parameters.template get<descriptors::OMEGA>();
    V uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbm<DESCRIPTOR>::addExternalForce(cell, u, omega, rho);
    return {rho / porosity, uSqr * porosity * porosity};
  }

};




/* This is essentially a no dynamics
 */
template < typename T, typename DESCRIPTOR, typename MOMENTA = momenta::BulkTuple, typename EQUILIBRIUM=equilibria::ZerothOrder>
struct DiffusionDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA>{

	using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

	using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

	using parameters = meta::list <descriptors::OMEGA>;

	template <typename M>
		using exchange_momenta = DiffusionDynamics<T,DESCRIPTOR,M>;
	template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
		CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
			V u[DESCRIPTOR::d];
			auto theta = MomentaF().computeRho(cell);
			MomentaF().computeU(cell,u);
			/* We don't do anything here */ 
			//particles::resetParticleRelatedFields<DESCRIPTOR,CELL,V>(cell);
			return {theta,util::normSqr<V,DESCRIPTOR::d>(u)};
		};

	//T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const any_platform override {
	//	return descriptors::t<T,DESCRIPTOR>(iPop) * rho - descriptors::t<T,DESCRIPTOR>(iPop);
	//};

	void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
		EquilibriumF().compute(cell, rho, u, fEq);
	};
	std::type_index id() override{
		return typeid(DiffusionDynamics);
	};
	AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
		return block.template getData<OperatorParameters<DiffusionDynamics>>();
	}
	std::string getName() const override {
		return "DiffusionDynamics";
	};
};

}

#endif

#include "legacy/dynamics.h"
