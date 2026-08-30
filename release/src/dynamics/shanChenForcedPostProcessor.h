/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Czelusniak,
 *                     Maximilian Schecher, Adrian Kummerlaender
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

#ifndef SHAN_CHEN_FORCED_POST_PROCESSOR_H
#define SHAN_CHEN_FORCED_POST_PROCESSOR_H

#define THIRD_COMPONENT
//#define FOURTH_COMPONENT

#include "core/blockStructure.h"
#include "core/postProcessing.h"

namespace olb {

namespace multicomponent_velocity {

struct Baricentric {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<>;

  int getPriority() const {
    return 0;
  }

  template <typename V, typename DESCRIPTOR, typename CELL, typename PARAMETERS>
  static void compute(CELL& cellA, CELL& cellB, PARAMETERS& params) any_platform {

    V rhoA = cellA.template getFieldComponent<descriptors::STATISTIC>(0);
    V rhoB = cellB.template getFieldComponent<descriptors::STATISTIC>(0);

    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    auto forceA = cellA.template getField<descriptors::FORCE>();
    auto forceB = cellB.template getField<descriptors::FORCE>();

    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoA + rhoB;

    Vector<V,DESCRIPTOR::d> uTot = ( uA + V(0.5)*rhoA*forceA
                                   + uB + V(0.5)*rhoB*forceB ) / rhoTot;

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);

    return;
  };
};

struct ShanChen {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct OMEGA_A : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_B : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<OMEGA_A,OMEGA_B>;

  int getPriority() const {
    return 0;
  }

  template <typename V, typename DESCRIPTOR, typename CELL, typename PARAMETERS>
  static void compute(CELL& cellA, CELL& cellB, PARAMETERS& params) any_platform {

    V rhoA = cellA.template getFieldComponent<descriptors::STATISTIC>(0);
    V rhoB = cellB.template getFieldComponent<descriptors::STATISTIC>(0);

    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    V omegaA = params.template get<OMEGA_A>();
    V omegaB = params.template get<OMEGA_B>();
    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoA*omegaA + rhoB*omegaB;

    Vector<V,DESCRIPTOR::d> uTot = (uA*omegaA + uB*omegaB) / rhoTot;

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);

    return;
  };
};

}

/**
* Multiphysics class for coupling between different lattices.
*/

struct RhoStatistics  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    auto statistic = cell.template getField<descriptors::STATISTIC>();
    statistic[0] = cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(statistic);
  }
};

template <typename POTENTIAL, unsigned N_COMPONENTS>
struct RhoPsiStatistics  {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct TEMPERATURE : public descriptors::FIELD_BASE<1> { };
  struct G           : public descriptors::FIELD_BASE<1> { };
  struct K           : public descriptors::FIELD_BASE<1> { };
  struct A           : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct B           : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct MM          : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct TCRIT       : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct DEVI        : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct ALPHA       : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
  struct GI          : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
  struct GII         : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };

  using parameters = meta::list<TEMPERATURE,G,K,A,B,MM,TCRIT,DEVI,ALPHA,GI,GII>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;

    auto t = parameters.template get<TEMPERATURE>();
    auto g = parameters.template get<G>();
    auto k = parameters.template get<K>();
    auto a = parameters.template get<A>();
    auto b = parameters.template get<B>();
    auto M = parameters.template get<MM>();
    auto T_c = parameters.template get<TCRIT>();
    auto m = parameters.template get<DEVI>();
    auto alpha = parameters.template get<ALPHA>();
    auto g_I = parameters.template get<GI>();
    auto g_II = parameters.template get<GII>();
    Vector<V,N_COMPONENTS> rhoField{};

    auto& cell1 = cells.template get<names::Component1>();
    rhoField[0] = cell1.computeRho();
    cell1.template setField<descriptors::STATISTIC>(rhoField[0]);

    auto& cell2 = cells.template get<names::Component2>();
    rhoField[1] = cell2.computeRho();
    cell2.template setField<descriptors::STATISTIC>(rhoField[1]);

#ifdef THIRD_COMPONENT
    auto& cell3 = cells.template get<names::Component3>();
    rhoField[2] = cell3.computeRho();
    cell3.template setField<descriptors::STATISTIC>(rhoField[2]);
#endif
#ifdef FOURTH_COMPONENT
    auto& cell4 = cells.template get<names::Component4>();
    rhoField[3] = cell4.computeRho();
    cell4.template setField<descriptors::STATISTIC>(rhoField[3]);
#endif

    V rhoTot = 0;
    for (unsigned n = 0; n < N_COMPONENTS; ++n) {
      rhoTot += rhoField[n];
    }
    V p = POTENTIAL().computeP(rhoField, t, a, b, T_c, m, alpha, g_I, g_II, M);
    cell1.template setField<descriptors::SCALAR>(p);
    V phi = 6.*(k*p - rhoTot/3.)/g;
    cell1.template setField<descriptors::PSI>(phi);
  }
};

// =========================================================================//
// ===========Shan Chen coupling without wall interaction===================//
// =========================================================================//
template <typename POTENTIAL, typename VELOCITY>
struct PseudopotentialForcedCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct G : public descriptors::FIELD_BASE<1> { };

  using parameters = typename meta::merge<
    meta::merge< typename POTENTIAL::parameters,
                 typename VELOCITY::parameters
               >,
    meta::list<G>
    >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    auto& cellA = cells.template get<names::A>();
    auto& cellB = cells.template get<names::B>();

    V rhoA = cellA.template getFieldComponent<descriptors::STATISTIC>(0);
    V rhoB = cellB.template getFieldComponent<descriptors::STATISTIC>(0);

    VELOCITY().template compute<V, DESCRIPTOR>(cellA, cellB, parameters);

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoBlockContribution{};
    Vector<V,DESCRIPTOR::d> rhoPartnerContribution{};
    V psiA = POTENTIAL().compute(rhoA, parameters);
    V psiB = POTENTIAL().compute(rhoB, parameters);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V rhoA_i = nextCellA.template getFieldComponent<descriptors::STATISTIC>(0);
      V rhoB_i = nextCellB.template getFieldComponent<descriptors::STATISTIC>(0);
      V psiA_i = POTENTIAL().compute(rhoA_i, parameters);
      V psiB_i = POTENTIAL().compute(rhoB_i, parameters);
      rhoBlockContribution   += psiB * psiA_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      rhoPartnerContribution += psiA * psiB_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    auto externalBlockForce   = cellA.template getField<descriptors::EXTERNAL_FORCE>();
    auto externalPartnerForce = cellB.template getField<descriptors::EXTERNAL_FORCE>();

    auto g = parameters.template get<G>();

    cellA.template setField<descriptors::FORCE>(externalBlockForce
        - g*rhoPartnerContribution/rhoA);
    cellB.template setField<descriptors::FORCE>(externalPartnerForce
        - g*rhoBlockContribution/rhoB);
  }

};

template <typename POTENTIAL>
struct ShanChenForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct RHO0 : public descriptors::FIELD_BASE<2> { };
  struct G    : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_A : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_B : public descriptors::FIELD_BASE<1> { };

  using parameters = typename meta::merge<
    typename POTENTIAL::parameters,
    meta::list<RHO0,G,OMEGA_A,OMEGA_B>
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    auto& cellA = cells.template get<names::A>();
    auto& cellB = cells.template get<names::B>();

    auto rho0 = parameters.template get<RHO0>();

    Vector<V,2> rhoField{
      cellA.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[0],
      cellB.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[1]
    };
    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    V omegaA = parameters.template get<OMEGA_A>();
    V omegaB = parameters.template get<OMEGA_B>();
    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoField[0]*omegaA
             + rhoField[1]*omegaB;

    Vector<V,DESCRIPTOR::d> uTot = (uA*rho0[0]*omegaA + uB*rho0[1]*omegaB) / rhoTot;

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoBlockContribution{};
    Vector<V,DESCRIPTOR::d> rhoPartnerContribution{};
    V psi1 = POTENTIAL().compute(rhoField[0], parameters);
    V psi2 = POTENTIAL().compute(rhoField[1], parameters);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V rhoA = POTENTIAL().compute(nextCellA.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[0], parameters);
      V rhoB = POTENTIAL().compute(nextCellB.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[1], parameters);
      rhoBlockContribution   += psi2 * rhoA * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      rhoPartnerContribution += psi1 * rhoB * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    auto externalBlockForce   = cellA.template getField<descriptors::EXTERNAL_FORCE>();
    auto externalPartnerForce = cellB.template getField<descriptors::EXTERNAL_FORCE>();

    auto g = parameters.template get<G>();

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);
    cellA.template setField<descriptors::FORCE>(externalBlockForce
        - g*rhoPartnerContribution/rhoField[0]);
    cellB.template setField<descriptors::FORCE>(externalPartnerForce
        - g*rhoBlockContribution/rhoField[1]);
  }

};

template <typename POTENTIAL>
struct ShanChenForcedHlbmCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct RHO0 : public descriptors::FIELD_BASE<2> { };
  struct G    : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_A : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_B : public descriptors::FIELD_BASE<1> { };

  using parameters = typename POTENTIAL::parameters::template decompose_into<
    meta::list<RHO0,G,OMEGA_A,OMEGA_B>::template include
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    auto& cellA = cells.template get<names::A>();
    auto& cellB = cells.template get<names::B>();

    auto rho0 = parameters.template get<RHO0>();

    Vector<V,2> rhoField{
      cellA.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[0],
      cellB.template getFieldComponent<descriptors::STATISTIC>(0) * rho0[1]
    };
    Vector<V,DESCRIPTOR::d> uA{};
    lbm<DESCRIPTOR>::computeJ(cellA, uA);
    Vector<V,DESCRIPTOR::d> uB{};
    lbm<DESCRIPTOR>::computeJ(cellB, uB);

    V omegaA = parameters.template get<OMEGA_A>();
    V omegaB = parameters.template get<OMEGA_B>();
    // Computation of the common velocity, shared among the two populations
    V rhoTot = rhoField[0]*omegaA
             + rhoField[1]*omegaB;

    const V p = cellA.template getField<descriptors::POROSITY>();
    Vector<V,DESCRIPTOR::d> uCell{};
    cellA.computeU(uCell.data());
    Vector<V,DESCRIPTOR::d> uTot = p     * ((uA*rho0[0]*omegaA + uB*rho0[1]*omegaB) / rhoTot)
                                 + (1-p) * (uCell);

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoBlockContribution{};
    Vector<V,DESCRIPTOR::d> rhoPartnerContribution{};
    V psi1 = POTENTIAL().compute(rhoField[0], parameters);
    V psi2 = POTENTIAL().compute(rhoField[1], parameters);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V rhoA = POTENTIAL().compute(nextCellA.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[0], parameters);
      V rhoB = POTENTIAL().compute(nextCellB.template getFieldComponent<descriptors::STATISTIC>(0)*rho0[1], parameters);
      rhoBlockContribution   += psi2 * rhoA * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      rhoPartnerContribution += psi1 * rhoB * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    Vector<V,DESCRIPTOR::d> externalBlockForce   = cellA.template getField<descriptors::EXTERNAL_FORCE>();
    Vector<V,DESCRIPTOR::d> externalPartnerForce = cellB.template getField<descriptors::EXTERNAL_FORCE>();

    auto g = parameters.template get<G>();

    cellA.template setField<descriptors::VELOCITY>(uTot);
    cellB.template setField<descriptors::VELOCITY>(uTot);
    cellA.template setField<descriptors::FORCE>(externalBlockForce
        - g*rhoPartnerContribution/rhoField[0]);
    cellB.template setField<descriptors::FORCE>(externalPartnerForce
        - g*rhoBlockContribution/rhoField[1]);
  }

};

// Pseudopotential Single-component-multi-phase force
template <typename POTENTIAL>
struct PseudopotentialForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = typename POTENTIAL::parameters::template include<>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Component1>::descriptor_t;

    const auto g = parameters.template get<typename POTENTIAL::G>();
    const auto sig = parameters.template get<typename POTENTIAL::KAPPAP>();

    V rho = 0;

    auto& cell1 = cells.template get<names::Component1>();
    rho = cell1.template getFieldComponent<descriptors::STATISTIC>(0);

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> totalForce{};
    V p = POTENTIAL().computeP(rho, parameters);
    V psi = util::sqrt(util::abs(6.*(p - rho/3.)));
    const int symmNumer = DESCRIPTOR::d * ( DESCRIPTOR::d + 1 ) / 2;
    Vector<V,DESCRIPTOR::d> M1{};
    Vector<V,symmNumer> M2{};

    V cs2 = 1./descriptors::invCs2<V,DESCRIPTOR>();
    V cs4 = 1./descriptors::invCs2<V,DESCRIPTOR>()/descriptors::invCs2<V,DESCRIPTOR>();

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      V rho_i;
      auto nextCell1 = cell1.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rho_i = nextCell1.template getFieldComponent<descriptors::STATISTIC>(0);

      V p_i = POTENTIAL().computeP(rho_i, parameters);
      V psi_i = util::sqrt(util::abs(6.*(p_i - rho_i/3.)));

      for (int iD_1 = 0; iD_1 < DESCRIPTOR::d; ++iD_1) {

        M1[iD_1] += descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::c<DESCRIPTOR>(iPop,iD_1) * psi_i;

        for (int iD_2 = iD_1; iD_2 < DESCRIPTOR::d; ++iD_2) {

          const int iD_Tensor = iD_1 * DESCRIPTOR::d - ( iD_1 - 1 ) * iD_1 / 2  + iD_2 - iD_1;

          M2[iD_Tensor] += descriptors::t<V,DESCRIPTOR>(iPop)
                          * ( descriptors::c<DESCRIPTOR>(iPop,iD_1) * descriptors::c<DESCRIPTOR>(iPop,iD_2)
                              - cs2 * ( iD_1 == iD_2 )
                            ) * psi_i;

        }

      }

    }

    for (int iD_1 = 0; iD_1 < DESCRIPTOR::d; ++iD_1) {

      totalForce[iD_1] = - g * psi * M1[iD_1];

      for (int iD_2 = 0; iD_2 < DESCRIPTOR::d; ++iD_2) {

        int i = util::min( iD_1, iD_2 );
        int j = util::max( iD_1, iD_2 );

        const int iD_Tensor = i * DESCRIPTOR::d - ( i - 1 ) * i / 2  + j - i;
        const int iD_Symm = iD_2 * DESCRIPTOR::d - ( iD_2 - 1 ) * iD_2 / 2;

        totalForce[iD_1] += ( sig*cs2 - cs2 - 1. ) * M1[iD_2] / cs2 * M2[iD_Tensor] / cs4 / 6.
                          - ( sig - 1. ) * cs2 * M1[iD_1] / cs2 * M2[iD_Symm] / cs4 / 6.;

      }

    }

    // Saving total force and pressure from EOS
    auto externalForce1 = cell1.template getField<descriptors::EXTERNAL_FORCE>();
    cell1.template setField<descriptors::FORCE>(externalForce1 + totalForce/rho);
    cell1.template setField<descriptors::SCALAR>(p);

  }

};

// =========================================================================//
// ==================Single-component-multi-phase coupling==================//
// =========================================================================//
// Free-Energy Chemical Potential
/*
For now, this Free-Energy will be placed here, move it to a separeted
folder in the future.
This free-energy is based on the work "A.J. Wagner, Phys. Rev. E (2006)".
*/
template<typename EQUATION_OF_STATE>
struct ChemicalPotentialPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = typename EQUATION_OF_STATE::parameters::template include<>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Component1>::descriptor_t;

    const auto kappa = parameters.template get<typename EQUATION_OF_STATE::KAPPA>();

    V rho;
    auto& cell = cells.template get<names::Component1>();
    rho = cell.template getFieldComponent<descriptors::STATISTIC>(0);

    // Computing chemical potential given by EOS
    V mu = EQUATION_OF_STATE().computeMU( rho, parameters );
    // Computing the density derivatives
    V d2_rho = 0.;
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      V rho_i;
      auto nextCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      rho_i = nextCell.template getFieldComponent<descriptors::STATISTIC>(0);

      d2_rho += V(6.) * descriptors::t<V,DESCRIPTOR>(iPop) * ( rho_i - rho );

    }

    mu += - kappa*d2_rho;

    cell.template setField<descriptors::CHEM_POTENTIAL>(mu);
    cell.template setField<descriptors::SCALAR>(d2_rho);

  }

};

// Free-energy single component force
/*
For now, this Free-Energy will be placed here, move it to a separeted
folder in the future.
This free-energy is based on the work "A.J. Wagner, Phys. Rev. E (2006)".
*/
struct FreeEnergyForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  //using parameters = meta::list<KAPPA,TCRIT>;
  using parameters = meta::list<>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Component1>::descriptor_t;

    auto& cell = cells.template get<names::Component1>();

    V rho = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> dmu{}, drho{};

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {

      auto nextCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto mu_i = nextCell.template getField<descriptors::CHEM_POTENTIAL>();
      auto rho_i = nextCell.template getFieldComponent<descriptors::STATISTIC>(0);

      for(int iD = 0; iD < DESCRIPTOR::d; ++iD){
        dmu[iD] += V(3.) * descriptors::t<V,DESCRIPTOR>(iPop)
                * descriptors::c<DESCRIPTOR>(iPop)[iD] * mu_i;
        drho[iD] += V(3.) * descriptors::t<V,DESCRIPTOR>(iPop)
                 * descriptors::c<DESCRIPTOR>(iPop)[iD] * rho_i;
      }

    }

    Vector<V,DESCRIPTOR::d> interactionForce{};

    for(int iD = 0; iD < DESCRIPTOR::d; ++iD){
        interactionForce[iD] = - rho * dmu[iD] + drho[iD]/V(3.);
    }

    // Computation of the common velocity, shared among the two populations, consisting
    // of u and the momentum difference due to interaction potential plus external force
    cell.template setField<descriptors::FORCE>(interactionForce/rho);

  }

};


// =========================================================================//
// ==================Multi-component-multi-phase coupling===================//
// =========================================================================//

/** Multi-Component-Multi-Phase Shan-Chen force with thermodynamic equation
 * of state based on
 *
 * 1. Czelusniak L E, Mapelli V P, Guzella M S, Cabezas-Gómez L, Wagner A J (2020)
 * Force approach for the pseudopotential lattice Boltzmann method. Physical
 * Review E, 102(3), 033307.
 * DOI: 10.1103/PhysRevE.102.033307
 * 2. Peng C, Ayala L F, Ayala O M (2021) A thermodynamically consistent pseudo-
 * potential lattice Boltzmann model for multi-component, multiphase, partially
 * miscible mixtures. Journal of Computational Physics, 429, 110018.
 * DOI: 10.1016/j.jcp.2020.110018
 **/
template <unsigned N_COMPONENTS>
struct MCMPForcedPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct CHI         : public descriptors::FIELD_BASE<N_COMPONENTS> { };
  struct G           : public descriptors::FIELD_BASE<1> { };
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct EPSILON     : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<CHI,G,SIGMA,EPSILON>;

  int getPriority() const {
    return 0;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Component1>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Component1>::descriptor_t;

    auto chi = parameters.template get<CHI>();
    auto g = parameters.template get<G>();
    auto sig = parameters.template get<SIGMA>();
    auto eps = parameters.template get<EPSILON>();

    Vector<V,10> rhoField{};
    V rhoTot = 0;
    Vector<V,DESCRIPTOR::d> u_sum{};
    auto& cell1 = cells.template get<names::Component1>();
    rhoField[0] = cell1.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u1{};
    lbm<DESCRIPTOR>::computeJ(cell1, u1);
    rhoTot += rhoField[0];
    u_sum += u1;

    auto& cell2 = cells.template get<names::Component2>();
    rhoField[1] = cell2.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u2{};
    lbm<DESCRIPTOR>::computeJ(cell2, u2);
    rhoTot += rhoField[1];
    u_sum += u2;

#ifdef THIRD_COMPONENT
    auto& cell3 = cells.template get<names::Component3>();
    rhoField[2] = cell3.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u3{};
    lbm<DESCRIPTOR>::computeJ(cell3, u3);
    rhoTot += rhoField[2];
    u_sum += u3;
#endif
#ifdef FOURTH_COMPONENT
    auto& cell4 = cells.template get<names::Component4>();
    rhoField[3] = cell4.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> u4{};
    lbm<DESCRIPTOR>::computeJ(cell4, u4);
    rhoTot += rhoField[3];
    u_sum += u4;
#endif

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> totalForce{};
    Vector<V,DESCRIPTOR::d> A_ij{};
    V phi = cell1.template getField<descriptors::PSI>();
    V psi = util::sqrt(util::abs(phi));

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      auto nextCell1 = cell1.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      V phi_i = nextCell1.template getField<descriptors::PSI>();
      V psi_i = 0;
      if (phi_i >= 0){
        psi_i = util::sqrt(util::abs(phi_i));
      }
      else {
        psi_i = -1*util::sqrt(util::abs(phi_i));
      }
      totalForce += -g * psi * psi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop);
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        auto nextCell1 = cell1.neighbor(descriptors::c<DESCRIPTOR>(jPop));
        V phi_j = nextCell1.template getField<descriptors::PSI>();
        V psi_j = 0;
        if (phi_j >= 0){
          psi_j = util::sqrt(util::abs(phi_j));
        }
        else {
          psi_j = -1*util::sqrt(util::abs(phi_j));
        }
        Vector<V,DESCRIPTOR::d> c1{};
        Vector<V,DESCRIPTOR::d> c2{};
        for (int m = 0; m < DESCRIPTOR::d; ++m) {
          Vector<V,DESCRIPTOR::d> C1{};
          V C2 = 0.;
          for (int n = 0; n < DESCRIPTOR::d; ++n) {
            C1[n] += descriptors::c<DESCRIPTOR>(jPop)[m] * descriptors::c<DESCRIPTOR>(jPop)[n];
            C2 += descriptors::c<DESCRIPTOR>(jPop)[n] * descriptors::c<DESCRIPTOR>(jPop)[n] - 1./descriptors::invCs2<V,DESCRIPTOR>();
            if (n==m){
              C1[n] -= 1./descriptors::invCs2<V,DESCRIPTOR>();
            }
          }
          for (int n = 0; n < DESCRIPTOR::d; ++n) {
            c1[m] += descriptors::c<DESCRIPTOR>(iPop)[n] * C1[n];
          }
          c2[m] += descriptors::c<DESCRIPTOR>(iPop)[m] * C2;
        }
        A_ij = 0.5*g*descriptors::t<V,DESCRIPTOR>(iPop)*descriptors::t<V,DESCRIPTOR>(jPop)*descriptors::invCs2<V,DESCRIPTOR>()*
                 ((3./2.*eps-(sig-1))*c1 + (sig-1)*c2);
        totalForce += A_ij * psi_i * psi_j;
      }

    }
    // Computation of the common velocity, shared among the two populations, consisting
    // of u and the momentum difference due to interaction potential plus external force
    auto externalForce1 = cell1.template getField<descriptors::EXTERNAL_FORCE>();
    Vector<V,DESCRIPTOR::d> uTot = (u_sum + totalForce/2.) / rhoTot + externalForce1/2.;
    cell1.template setField<descriptors::FORCE>(externalForce1 + chi[0]*totalForce/rhoField[0]);
    cell1.template setField<descriptors::VELOCITY>(uTot);

    auto externalForce2 = cell2.template getField<descriptors::EXTERNAL_FORCE>();
    cell2.template setField<descriptors::FORCE>(externalForce2 + chi[1]*totalForce/rhoField[1]);
    cell2.template setField<descriptors::VELOCITY>(uTot);

#ifdef THIRD_COMPONENT
    auto externalForce3 = cell3.template getField<descriptors::EXTERNAL_FORCE>();
    cell3.template setField<descriptors::FORCE>(externalForce3 + chi[2]*totalForce/rhoField[2]);
    cell3.template setField<descriptors::VELOCITY>(uTot);
#endif
#ifdef FOURTH_COMPONENT
    auto externalForce4 = cell4.template getField<descriptors::EXTERNAL_FORCE>();
    cell4.template setField<descriptors::FORCE>(externalForce4 + chi[3]*totalForce/rhoField[3]);
    cell4.template setField<descriptors::VELOCITY>(uTot);
#endif

  }

};

}

#endif
