/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_H

#include "core/operator.h"

namespace olb {

/// TotalEnthalpyPhaseChange between a Navier-Stokes and an Advection-Diffusion lattice
struct TotalEnthalpyPhaseChangeCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct T_S             : public descriptors::FIELD_BASE<1> { };
  struct T_L             : public descriptors::FIELD_BASE<1> { };
  struct CP_S            : public descriptors::FIELD_BASE<1> { };
  struct CP_L            : public descriptors::FIELD_BASE<1> { };
  struct L               : public descriptors::FIELD_BASE<1> { };
  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T_COLD          : public descriptors::FIELD_BASE<1> { };
  struct DELTA_T         : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<T_S,T_L,CP_S,CP_L,L,FORCE_PREFACTOR,T_COLD,DELTA_T>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAD = cells.template get<names::Temperature>();

    auto cp_s = parameters.template get<CP_S>();
    auto cp_l = parameters.template get<CP_L>();
    auto T_s = parameters.template get<T_S>();
    auto T_l = parameters.template get<T_L>();
    auto l = parameters.template get<L>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    auto T0 = parameters.template get<T_COLD>();
    //auto deltaT = parameters.template get<DELTA_T>();

    const V H_s  = cp_s * T_s;
    const V H_l  = cp_l * T_l + l;
    V enthalpy = cellAD.computeRho();
    V temperature, liquid_fraction;

    if (enthalpy <= H_s) {
      temperature = T_s - (H_s - enthalpy) / cp_s;
      liquid_fraction = 0.;
    }
    else if (enthalpy >= H_l) {
      temperature = T_l + (enthalpy - H_l) / cp_l;
      liquid_fraction = 1.;
    }
    else {
      temperature = (H_l - enthalpy) / (H_l - H_s) * T_s + (enthalpy - H_s) / (H_l - H_s) * T_l;
      liquid_fraction = (enthalpy - H_s) / l;
    }
    //std::cout << "temperature: " << temperature << std::endl;
    //std::cout << "liquid fraction: " << liquid_fraction << std::endl;
    cellNS.template setField<descriptors::POROSITY>(liquid_fraction);
    cellAD.template setField<descriptors::TEMPERATURE>(temperature);

    // computation of the bousinessq force
    auto force = cellNS.template getFieldPointer<descriptors::FORCE>();
    V temperatureDifference = temperature - T0;
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNS.computeU(u);
    cellAD.template setField<descriptors::VELOCITY>(u);
  }
};



/// Coupling between a Navier-Stokes and an Advection-Diffusion lattice
struct NavierStokesAdvectionDiffusionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,T0>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }
};

/// Velocity coupling between Navier-Stokes and an Advection-Diffusion lattice
struct NavierStokesAdvectionDiffusionVelocityCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Concentration0>();

    // Velocity coupling
    V u[DESCRIPTOR::d] { };
    cellNSE.computeU(u);
    cellADE.template setField<descriptors::VELOCITY>(u);
  }
};

/// AD coupling with Boussinesq bouancy for Smagorinsky-LES
struct SmagorinskyBoussinesqCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct PR_TURB : public descriptors::FIELD_BASE<1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,SMAGORINSKY_PREFACTOR,PR_TURB,T0,OMEGA_NSE,OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // computation of the bousinessq force
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    V temperatureDifference = cellADE.computeRho() - parameters.template get<T0>();

    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }

    // Velocity coupling
    auto u = cellADE.template getField<descriptors::VELOCITY>();
    // tau coupling
    auto tauNS = cellNSE.template getFieldPointer<descriptors::TAU_EFF>();
    auto tauAD = cellADE.template getFieldPointer<descriptors::TAU_EFF>();

    V rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cellNSE.computeAllMomenta(rho, u.data(), pi);
    cellADE.template setField<descriptors::VELOCITY>(u);
    V PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    V PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    V tau_mol_NS = V{1} / parameters.template get<OMEGA_NSE>();
    V tau_mol_AD = V{1} / parameters.template get<OMEGA_ADE>();
    /// Turbulent realaxation time
    V smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    V tau_turb_NS = V{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    /// Effective realaxation time
    tauNS[0] = tau_mol_NS+tau_turb_NS;

    V prTurb = parameters.template get<PR_TURB>();
    V tauTurbADPrefactor = descriptors::invCs2<V,DESCRIPTOR_ADE>() / descriptors::invCs2<V,DESCRIPTOR>() / prTurb;
    V tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
    tauAD[0] = tau_mol_AD+tau_turb_AD;
  }

};


/// Reaction coupling for three species homogeneous bulk reaction
template<typename T>
struct ReactionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct LATTICE_REACTION_COEFF : public descriptors::FIELD_BASE<1> { };
  struct STOCH_COEFF : public descriptors::FIELD_BASE<3> { };
  struct REACTION_ORDER : public descriptors::FIELD_BASE<3> { };

  using parameters = meta::list<LATTICE_REACTION_COEFF, STOCH_COEFF, REACTION_ORDER>;

  template<typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    auto stochCoeff = parameters.template get<STOCH_COEFF>();
    auto reactionOrder = parameters.template get<REACTION_ORDER>();
    T source = parameters.template get<LATTICE_REACTION_COEFF>();
    {
      T conc = cells.template get<names::Concentration0>().computeRho();
      source *= util::pow(conc, reactionOrder[0]);
    }
    {
      T conc = cells.template get<names::Concentration1>().computeRho();
      source *= util::pow(conc, reactionOrder[1]);
    }
    {
      T conc = cells.template get<names::Concentration2>().computeRho();
      source *= util::pow(conc, reactionOrder[2]);
    }
    {
      cells.template get<names::Concentration0>().template setField<descriptors::SOURCE>(stochCoeff[0]*source);
      cells.template get<names::Concentration1>().template setField<descriptors::SOURCE>(stochCoeff[1]*source);
      cells.template get<names::Concentration2>().template setField<descriptors::SOURCE>(stochCoeff[2]*source);
    }
  }
};


/// Reaction Coupling for the In-Bulk appraoch of lognitudinalMixing3d example
template<typename T>
struct LongitudinalMixingReactionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct REACTION_CONSTANT : public descriptors::FIELD_BASE<1> { };
  struct EQUILIBRIUM : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<REACTION_CONSTANT, EQUILIBRIUM>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
        T forwardConstant = parameters.template get<REACTION_CONSTANT>();
        T Ceq = parameters.template get<EQUILIBRIUM>();

        auto allow_source = cells.template get<names::Concentration0>().template getField<descriptors::GAMMA>(); //allow source only at surface
        T concC = cells.template get<names::Concentration0>().computeRho(); //current concentration
        T source = -forwardConstant*(concC - Ceq); //difference between current Conc and equilibrium concentration determines source
        cells.template get<names::Concentration0>().template setField<descriptors::SOURCE>(source*allow_source);
  }
};



/// LES-ADE coupling for multiple reactions
template<typename T, int numComp>
struct LESReactionCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct LATTICE_REACTION_COEFF : public descriptors::FIELD_BASE<1> { };
  struct STOCH_COEFF : public descriptors::FIELD_BASE<numComp> { };
  struct REACTION_ORDER : public descriptors::FIELD_BASE<numComp> { };
  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct SCHMIDT : descriptors::FIELD_BASE<numComp> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGAS_ADE : public descriptors::FIELD_BASE<numComp> { };

  using parameters = meta::list<LATTICE_REACTION_COEFF, STOCH_COEFF, REACTION_ORDER, SMAGORINSKY_PREFACTOR, SCHMIDT, OMEGA_NSE, OMEGAS_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Concentration0>::descriptor_t;
    // Velocity coupling
    auto u = cells.template get<names::Concentration0>().template getField<descriptors::VELOCITY>();
    T rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u.data(), pi);
    cells.template get<names::Concentration0>().template setField<descriptors::VELOCITY>(u);
    cells.template get<names::Concentration1>().template setField<descriptors::VELOCITY>(u);
    cells.template get<names::Concentration2>().template setField<descriptors::VELOCITY>(u);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    T tau_mol_NS = T{1} / parameters.template get<OMEGA_NSE>();
    auto omegasAD = parameters.template get<OMEGAS_ADE>();
    T tau_mol_AD0 = T{1} / omegasAD[0];
    T tau_mol_AD1 = T{1} / omegasAD[1];
    T tau_mol_AD2 = T{1} / omegasAD[2];
    /// Turbulent realaxation time
    T smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    T tau_turb_NS = T{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    // Schmidt number stabilization
    auto Sc = parameters.template get<SCHMIDT>();
    T tauTurbADPrefactor0 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[0];
    T tauTurbADPrefactor1 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[1];
    T tauTurbADPrefactor2 = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / Sc[2];
    T tau_turb_AD0 = tau_turb_NS * tauTurbADPrefactor0;
    T tau_turb_AD1 = tau_turb_NS * tauTurbADPrefactor1;
    T tau_turb_AD2 = tau_turb_NS * tauTurbADPrefactor2;
    cells.template get<names::Concentration0>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD0 + tau_turb_AD0));
    cells.template get<names::Concentration1>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD1 + tau_turb_AD1));
    cells.template get<names::Concentration2>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD2 + tau_turb_AD2));

    auto stochCoeff = parameters.template get<STOCH_COEFF>();
    auto reactionOrder = parameters.template get<REACTION_ORDER>();
    T source = parameters.template get<LATTICE_REACTION_COEFF>();
    {
      T conc = cells.template get<names::Concentration0>().computeRho();
      source *= util::pow(conc, reactionOrder[0]);
    }
    {
      T conc = cells.template get<names::Concentration1>().computeRho();
      source *= util::pow(conc, reactionOrder[1]);
    }
    {
      T conc = cells.template get<names::Concentration2>().computeRho();
      source *= util::pow(conc, reactionOrder[2]);
    }
    {
      cells.template get<names::Concentration0>().template setField<descriptors::SOURCE>(stochCoeff[0]*source);
      cells.template get<names::Concentration1>().template setField<descriptors::SOURCE>(stochCoeff[1]*source);
      cells.template get<names::Concentration2>().template setField<descriptors::SOURCE>(stochCoeff[2]*source);
    }
  }

};

/// Porous ADE correction term
template<typename T>
struct PorousADECorrection {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct DIFFUSION : public descriptors::FIELD_BASE<1> { };
  //struct VELOCITY : public descriptors::FIELD_BASE<3> { };

  using parameters = meta::list<DIFFUSION>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    T diffusion = parameters.template get<DIFFUSION>();

    auto& cell = cells.template get<names::Concentration0>();

    auto vel = cell.template getField<descriptors::VELOCITY>();

    T porOld = cell.template getField<descriptors::SCALAR>();
    T porosity = cell.template getField<descriptors::POROSITY>();
    T porosityXP = cell.neighbor({1,0,0}).template getField<descriptors::POROSITY>();
    T porosityXM = cell.neighbor({-1,0,0}).template getField<descriptors::POROSITY>();
    T porosityYP = cell.neighbor({0,1,0}).template getField<descriptors::POROSITY>();
    T porosityYM = cell.neighbor({0,-1,0}).template getField<descriptors::POROSITY>();
    T porosityZP = cell.neighbor({0,0,1}).template getField<descriptors::POROSITY>();
    T porosityZM = cell.neighbor({0,0,-1}).template getField<descriptors::POROSITY>();
    T dxPor = 0.5 * (porosityXP - porosityXM);
    T dyPor = 0.5 * (porosityYP - porosityYM);
    T dzPor = 0.5 * (porosityZP - porosityZM);
    T conc = cell.computeRho();
    T concPlusX = cell.neighbor({1,0,0}).computeRho();
    T concMinusX = cell.neighbor({-1,0,0}).computeRho();
    T concPlusY = cell.neighbor({0,1,0}).computeRho();
    T concMinusY = cell.neighbor({0,-1,0}).computeRho();
    T concPlusZ = cell.neighbor({0,0,1}).computeRho();
    T concMinusZ = cell.neighbor({0,0,-1}).computeRho();
    T dxConc = 0.5 * (concPlusX - concMinusX);
    T dyConc = 0.5 * (concPlusY - concMinusY);
    T dzConc = 0.5 * (concPlusZ - concMinusZ);
    //T source = cell.template getField<descriptors::SOURCE>();
    T source = (-conc*(porosity - porOld)/porosity -conc*(vel[0]*dxPor + vel[1]*dyPor + vel[2]*dzPor)/porosity + diffusion*(dxConc*dxPor + dyConc*dyPor + dzConc*dzPor)/porosity);
    cell.template setField<descriptors::SOURCE>(source);
    cell.template setField<descriptors::SCALAR>(porosity);
  }

};



/// LES-ADE coupling with Schmidt number stabilization
template<typename T>
struct LESADECoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct SMAGORINSKY_PREFACTOR : public descriptors::FIELD_BASE<1> { };
  struct SCHMIDT : descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<SMAGORINSKY_PREFACTOR, SCHMIDT, OMEGA_NSE, OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Concentration0>::descriptor_t;
    // Velocity coupling
    auto u = cells.template get<names::Concentration0>().template getField<descriptors::VELOCITY>();
    T rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u.data(), pi);
    cells.template get<names::Concentration0>().template setField<descriptors::VELOCITY>(u);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
    /// Molecular realaxation time
    T tau_mol_NS = T{1} / parameters.template get<OMEGA_NSE>();
    T tau_mol_AD = T{1} / parameters.template get<OMEGA_ADE>();
    /// Turbulent realaxation time
    T smagoPrefactor = parameters.template get<SMAGORINSKY_PREFACTOR>();
    T tau_turb_NS = T{0.5}*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
    // Schmidt number stabilization
    T tauTurbADPrefactor = descriptors::invCs2<T,DESCRIPTOR_ADE>() / descriptors::invCs2<T,DESCRIPTOR>() / parameters.template get<SCHMIDT>();
    T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
    cells.template get<names::Concentration0>().template setField<descriptors::OMEGA>(T{1} / (tau_mol_AD + tau_turb_AD));
  }

};

/// VANS-ADE coupling
template<typename T>
struct VANSADECoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct PARTICLE_DIAMETER : public descriptors::FIELD_BASE<1> { };
  struct VISCOSITY : public descriptors::FIELD_BASE<1> { };
  struct CONV_VEL : public descriptors::FIELD_BASE<1> { };
  struct CONV_DENS : public descriptors::FIELD_BASE<1> { };
  struct PART_DENS : public descriptors::FIELD_BASE<1> { };
  struct DT : public descriptors::FIELD_BASE<1> { };
  struct CONV_FORCE : public descriptors::FIELD_BASE<1> { };
  struct CONV_MASS : public descriptors::FIELD_BASE<1> { };
  struct EARTH_ACC : public descriptors::FIELD_BASE<3> { };

  using parameters = meta::list<PARTICLE_DIAMETER,VISCOSITY,CONV_VEL,CONV_DENS,PART_DENS,DT,CONV_FORCE,CONV_MASS,EARTH_ACC>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    T particleDiam = parameters.template get<PARTICLE_DIAMETER>();
    T visc = parameters.template get<VISCOSITY>();
    T convVel = parameters.template get<CONV_VEL>();
    T convDens = parameters.template get<CONV_DENS>();
    T convForce = parameters.template get<CONV_FORCE>();
    T convMass = parameters.template get<CONV_MASS>();
    T partDens = parameters.template get<PART_DENS>();
    T dt = parameters.template get<DT>();
    auto earthAcc = parameters.template get<EARTH_ACC>();
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Concentration0>();

    auto u_p = cellADE.template getField<descriptors::VELOCITY>();
    auto u_pXp = cellADE.neighbor({1,0,0}).template getField<descriptors::VELOCITY2>();
    auto u_pXm = cellADE.neighbor({-1,0,0}).template getField<descriptors::VELOCITY2>();
    auto u_pYp = cellADE.neighbor({0,1,0}).template getField<descriptors::VELOCITY2>();
    auto u_pYm = cellADE.neighbor({0,-1,0}).template getField<descriptors::VELOCITY2>();
    auto u_pZp = cellADE.neighbor({0,0,1}).template getField<descriptors::VELOCITY2>();
    auto u_pZm = cellADE.neighbor({0,0,-1}).template getField<descriptors::VELOCITY2>();
    auto u_f = cellADE.template getField<descriptors::VELOCITY>();
    T por = (1. - cellADE.computeRho());
    T rho { };
    cellNSE.computeRhoU(rho, u_f.data());

    T relU[3] = {0.};
    relU[0] = (u_f[0] - u_p[0]) * convVel;
    relU[1] = (u_f[1] - u_p[1]) * convVel;
    relU[2] = (u_f[2] - u_p[2]) * convVel;
    T norm = util::sqrt(relU[0]*relU[0] + relU[1]*relU[1] + relU[2]*relU[2]);
    T Re_p = norm * particleDiam / visc;
    T C_D = 0.;
    if(Re_p != 0.){
    if(Re_p <= 1000.){
      C_D = 24.*(1. + 0.15*util::pow(Re_p,0.687))/Re_p;
    }else{
      C_D = 0.44;
    }}
    T F_D[3] = {0.};
    //F_D[0] = C_D * 3./4. * rho * convDens * norm * relU[0] / particleDiam / partDens;
    //F_D[1] = C_D * 3./4. * rho * convDens * norm * relU[1] / particleDiam / partDens;
    //F_D[2] = C_D * 3./4. * rho * convDens * norm * relU[2] / particleDiam / partDens;
    //dragCoeff = (9.*converter_.getPhysViscosity()*converter_.getPhysDensity()*converter_.getConversionFactorTime()) / (2.*pRho_*pRadius_*pRadius_);
    F_D[0] = 9.*visc*convDens/partDens/particleDiam*relU[0];
    F_D[1] = 9.*visc*convDens/partDens/particleDiam*relU[1];
    F_D[2] = 9.*visc*convDens/partDens/particleDiam*relU[2];

    T mass = partDens * 3.14/4. * particleDiam * particleDiam;
    if(cellADE.computeRho() > 0.001){
    u_p[0] += dt * (F_D[0]/*mass*/ + earthAcc[0]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pXp - u_pXm);
    u_p[1] += dt * (F_D[1]/*mass*/ + earthAcc[1]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pYp - u_pYm);
    u_p[2] += dt * (F_D[2]/*mass*/ + earthAcc[2]*(partDens - convDens)/partDens) / convVel - u_p * 0.5 * (u_pZp - u_pZm);
    cellADE.template setField<descriptors::VELOCITY>(u_p);
    }

    T porPlusX = (1. - cellADE.neighbor({1,0,0}).computeRho());
    T porMinusX = (1. - cellADE.neighbor({-1,0,0}).computeRho());
    T porPlusY = (1. - cellADE.neighbor({0,1,0}).computeRho());
    T porMinusY = (1. - cellADE.neighbor({0,-1,0}).computeRho());
    T porPlusZ = (1. - cellADE.neighbor({0,0,1}).computeRho());
    T porMinusZ = (1. - cellADE.neighbor({0,0,-1}).computeRho());

    T porDX2 = porPlusX -2.*por + porMinusX;
    T porDY2 = porPlusY -2.*por + porMinusY;
    T porDZ2 = porPlusZ -2.*por + porMinusZ;

    int choice = 0;
    if( porDX2 != 0.) choice += 1;
    if( porDY2 != 0.) choice += 1;
    if( porDZ2 != 0.) choice += 1;

    T coeff = (choice == 1) *1./4. + (choice == 2) * 1./6. + (choice == 3) * 5./36.;
    T pressCorr = por + coeff * (porDX2 + porDY2 + porDZ2);
    if(pressCorr < 0.001){ pressCorr = 1.; }
    cellNSE.template setField<descriptors::PRESSCORR>(pressCorr);
    T press = rho / pressCorr / descriptors::invCs2<T,DESCRIPTOR>();
    T pressCorrForce[3] = {0.};
    pressCorrForce[0] = press * 0.5 * (porPlusX - porMinusX);
    pressCorrForce[1] = press * 0.5 * (porPlusY - porMinusY);
    pressCorrForce[2] = press * 0.5 * (porPlusZ - porMinusZ);

    T force[3] = {0.};
    if(cellADE.computeRho() > 0.001){
    force[0] = pressCorrForce[0] / rho - (1. - por) * por * F_D[0] / convForce * convMass * partDens / rho + por * earthAcc[0] / convVel * dt;
    force[1] = pressCorrForce[1] / rho - (1. - por) * por * F_D[1] / convForce * convMass * partDens / rho + por * earthAcc[1] / convVel * dt;
    force[2] = pressCorrForce[2] / rho - (1. - por) * por * F_D[2] / convForce * convMass * partDens / rho + por * earthAcc[2] / convVel * dt;
    }
    cellNSE.template setField<descriptors::FORCE>(force);
  }

};

/// granular flow
template<typename T>
struct GranularCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FRICTION_ANGLE : public descriptors::FIELD_BASE<1> { };
  struct RHO0 : public descriptors::FIELD_BASE<1> { };
  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<FRICTION_ANGLE,RHO0,FORCE_PREFACTOR>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u, pi);
    // Stress tensor
    T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
    }
    T PiNeqNorm = util::sqrt(PiNeqNormSqr);

    if(PiNeqNorm != T{0.}){
      auto omega = cells.template get<names::NavierStokes>().template getFieldPointer<descriptors::OMEGA>();
      PiNeqNorm /= (-rho/omega[0]/descriptors::invCs2<T,DESCRIPTOR>());
      T angle = parameters.template get<FRICTION_ANGLE>();
      T ps = ( T{0.5}*omega[0] - T{1} ) * T{1./3.} * ( pi[0] + pi[3] + pi[5] );
      T granularTau = ps/rho * util::sin(T{3.14}*angle/T{180}) / PiNeqNorm * descriptors::invCs2<T,DESCRIPTOR>() + T{0.5};
      omega[0] = T{1} / granularTau;
    }

    // computation of the bousinessq force
    auto force = cells.template get<names::NavierStokes>().template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    T densityDifference = (rho - parameters.template get<RHO0>()) / rho;
    for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
      force[iD] = forcePrefactor[iD] * densityDifference;
    }
  }

};

/// Poisson-Nernst-Planck coupling
template<typename T>
struct PNPCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct DX : public descriptors::FIELD_BASE<1> { };
  struct NPVELCOEFF : public descriptors::FIELD_BASE<1> { };
  struct POISSONCOEFF : public descriptors::FIELD_BASE<1> { };
  struct OMEGA : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<DX,NPVELCOEFF,POISSONCOEFF,OMEGA>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    T dX = parameters.template get<DX>();
    T velCoeff = parameters.template get<NPVELCOEFF>();
    T poissonCoeff = parameters.template get<POISSONCOEFF>();
    T omega = parameters.template get<OMEGA>();
    using DESCRIPTOR = typename CELLS::template value_t<names::Concentration0>::descriptor_t;

    auto& cellNP = cells.template get<names::Concentration0>();
    auto& cellP = cells.template get<names::Concentration1>();
    auto& cellNP2 = cells.template get<names::Concentration2>();

    T dxPsi = 0.;
    T dyPsi = 0.;
    T dzPsi = 0.;

    for(int iPop = 0; iPop < DESCRIPTOR::q; iPop++){
      dxPsi += descriptors::c<DESCRIPTOR>(iPop, 0) * cellP[iPop];
      dyPsi += descriptors::c<DESCRIPTOR>(iPop, 1) * cellP[iPop];
      dzPsi += descriptors::c<DESCRIPTOR>(iPop, 2) * cellP[iPop];
    }

    dxPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);
    dyPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);
    dzPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);

    T vel[3] ={0.};
    vel[0] -= velCoeff * dxPsi;
    vel[1] -= velCoeff * dyPsi;
    vel[2] -= velCoeff * dzPsi;
    cellNP.template setField<descriptors::VELOCITY>(vel);

    T vel2[3] ={0.};
    vel2[0] += velCoeff * dxPsi;
    vel2[1] += velCoeff * dyPsi;
    vel2[2] += velCoeff * dzPsi;
    cellNP2.template setField<descriptors::VELOCITY>(vel2);

    T concentration = cellNP.computeRho();
    if(util::abs(concentration) > 1.){ concentration = 0.; }
    T concentration2 = cellNP2.computeRho();
    if(util::abs(concentration2) > 1.){ concentration2 = 0.; }
    T poissonSource = poissonCoeff * (concentration - concentration2);
    cellP.template setField<descriptors::SOURCE>(poissonSource);
  }
};

/// Naver-Stokes-Poisson-Nernst-Planck coupling
template<typename T>
struct NSPNPCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct DX : public descriptors::FIELD_BASE<1> { };
  struct NPVELCOEFF : public descriptors::FIELD_BASE<1> { };
  struct POISSONCOEFF : public descriptors::FIELD_BASE<1> { };
  struct FORCECOEFF : public descriptors::FIELD_BASE<1> { };
  struct DTADE : public descriptors::FIELD_BASE<1> { };
  struct DTNSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<DX,NPVELCOEFF,POISSONCOEFF,FORCECOEFF,DTADE,DTNSE,OMEGA>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    T dX = parameters.template get<DX>();
    T velCoeff = parameters.template get<NPVELCOEFF>();
    T poissonCoeff = parameters.template get<POISSONCOEFF>();
    T forceCoeff = parameters.template get<FORCECOEFF>();
    T dtADE = parameters.template get<DTADE>();
    T dtNSE = parameters.template get<DTNSE>();
    T omega = parameters.template get<OMEGA>();

    using DESCRIPTOR = typename CELLS::template value_t<names::Concentration0>::descriptor_t;
    auto& cellNP = cells.template get<names::Concentration0>();
    auto& cellNP2 = cells.template get<names::Concentration1>();
    auto& cellP = cells.template get<names::Temperature>();
    auto& cellNSE = cells.template get<names::NavierStokes>();

    T dxPsi = 0.;
    T dyPsi = 0.;
    T dzPsi = 0.;

    for(int iPop = 0; iPop < DESCRIPTOR::q; iPop++){
      dxPsi += descriptors::c<DESCRIPTOR>(iPop, 0) * cellP[iPop];
      dyPsi += descriptors::c<DESCRIPTOR>(iPop, 1) * cellP[iPop];
      dzPsi += descriptors::c<DESCRIPTOR>(iPop, 2) * cellP[iPop];
    }

    dxPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);
    dyPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);
    dzPsi *= (-1.*omega*descriptors::invCs2<T,DESCRIPTOR>()/dX);

    auto vel = cellNP.template getField<descriptors::VELOCITY>();
    cellNSE.computeU(vel.data());
    vel[0] *= dtADE / dtNSE;
    vel[1] *= dtADE / dtNSE;
    vel[2] *= dtADE / dtNSE;

    vel[0] -= velCoeff * dxPsi;
    vel[1] -= velCoeff * dyPsi;
    vel[2] -= velCoeff * dzPsi;
    cellNP.template setField<descriptors::VELOCITY>(vel);

    auto vel2 = cellNP2.template getField<descriptors::VELOCITY>();
    cellNSE.computeU(vel2.data());
    vel2[0] *= dtADE / dtNSE;
    vel2[1] *= dtADE / dtNSE;
    vel2[2] *= dtADE / dtNSE;

    vel2[0] += velCoeff * dxPsi;
    vel2[1] += velCoeff * dyPsi;
    vel2[2] += velCoeff * dzPsi;
    cellNP2.template setField<descriptors::VELOCITY>(vel2);

    T concentration = cellNP.computeRho();
    if(util::abs(concentration) > 1.e5){ concentration = 0.; }
    T concentration2 = cellNP2.computeRho();
    if(util::abs(concentration2) > 1.e5){ concentration2 = 0.; }
    T poissonSource = poissonCoeff * (concentration - concentration2);
    cellP.template setField<descriptors::SOURCE>(poissonSource);

    auto force = cellNSE.template getField<descriptors::FORCE>();
    force = {0.,0.,0.};
    force[0] = forceCoeff * (concentration - concentration2);
    //force[1] = -forceCoeff * (concentration - concentration2) * dyPsi;
    //force[2] = -forceCoeff * (concentration - concentration2) * dzPsi;
    cellNSE.template setField<descriptors::FORCE>(force);
  }
};

/*
  LATTICE_UdV_SUM from superLatticeFieldReduction0 is  ∫UxdV
*/
template<typename V>
struct TurbulentChannelForce {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct LATTICE_U_SUM : public descriptors::FIELD_BASE<0,1> { };
  struct CHAR_LATTICE_U : public descriptors::FIELD_BASE<1> { };
  struct CHAR_LATTICE_L : public descriptors::FIELD_BASE<1> { };
  struct LATTICE_UTAU : public descriptors::FIELD_BASE<1> { };
  struct LATTICE_CHANNEL_VOLUME : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<LATTICE_U_SUM,CHAR_LATTICE_U,CHAR_LATTICE_L,LATTICE_UTAU,LATTICE_CHANNEL_VOLUME>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    //every value is the lattice(non-dimensional) unit
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    auto& cell = cells.template get<names::NavierStokes>();
    const V charlat_u = parameters.template get<CHAR_LATTICE_U>();
    const V charlat_l = parameters.template get<CHAR_LATTICE_L>();
    auto U_SUM = parameters.template get<LATTICE_U_SUM>();
    const V volume = parameters.template get<LATTICE_CHANNEL_VOLUME>();
    const V uave      = U_SUM[0] / volume;// uave = <Ux> = 1/V ∫UxdV
    const V utau   = parameters.template get<LATTICE_UTAU>();
    V tmpForce[DESCRIPTOR::d]  = {};
    tmpForce[0] = (utau * utau / charlat_l + (charlat_u - uave) * charlat_u / charlat_l);
    tmpForce[1] = (V) 0.0;
    tmpForce[2] = (V) 0.0;
    cell.template setField<descriptors::FORCE>(tmpForce);
  }
};

template<typename V>
struct PhaseFieldHeatEquationCouplingFirstPart{
	static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

	struct LAMBDA: public descriptors::FIELD_BASE<1> { };
	struct TAU : public descriptors::FIELD_BASE<1> { };
	struct DELTA_T : public descriptors::FIELD_BASE<1> { };
	struct OMEGAPFE : public descriptors::FIELD_BASE<1> { };

	using parameters = meta::list<LAMBDA, TAU, DELTA_T, OMEGAPFE>;

	template <typename CELLS, typename PARAMETERS>
		void apply(CELLS& cells, PARAMETERS& parameters) any_platform
		{
			//using V = typename CELLS::template value_t<names::A>::value_t;
			using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

			auto& cellPFE = cells.template get<names::A>();
			auto& cellHE = cells.template get<names::B>();


			auto porosity = cellPFE. template getField<descriptors::POROSITY>();
			V phi = cellPFE.computeRho();
			V theta = cellHE.computeRho(); // the current theta and phi 
			V oneMinusPhiSq = 1.- phi*phi;

			//V Q =  (phi - (porosity*parameters.template get<LAMBDA>()*(theta-1))*(1-phi*phi)-1*(1-porosity)*(1-phi*phi))*(1-phi*phi);  
			V Q =  (phi - (porosity*parameters.template get<LAMBDA>()*(theta-1))*(oneMinusPhiSq))*(oneMinusPhiSq)-5*(1-porosity)*2.*(1+phi);  

			Q = Q*parameters.template get<DELTA_T>()/parameters.template get<TAU>();	

			V omegaPFE = parameters.template get<OMEGAPFE>();
			cellPFE.template setField<descriptors::PREV_PHI>(phi);

			for(int iPop =0; iPop <DESCRIPTOR::q; iPop++)
			{
				cellPFE[iPop] = cellPFE[iPop] - omegaPFE*(cellPFE[iPop]-descriptors::t<V,DESCRIPTOR>(iPop)*(phi-1.))+descriptors::t<V,DESCRIPTOR>(iPop)*Q;
			}


		}
};

template<typename V>
struct PhaseFieldHeatEquationCouplingSecondPart{
	static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

	struct thermalConductivityWater : public descriptors::FIELD_BASE<1> { };
	struct thermalConductivityIce   : public descriptors::FIELD_BASE<1> { };
	struct thermalConductivitySoil  : public descriptors::FIELD_BASE<1> { };
	struct heatCapacityWater        : public descriptors::FIELD_BASE<1> { };
	struct heatCapacityIce          : public descriptors::FIELD_BASE<1> { };
	struct density 			: public descriptors::FIELD_BASE<1> { };
	struct etaPhi 			: public descriptors::FIELD_BASE<1> { };
	struct deltax 			: public descriptors::FIELD_BASE<1> { };
	struct deltat			: public descriptors::FIELD_BASE<1> { };

	using parameters = meta::list<
		thermalConductivityWater,
		thermalConductivityIce,
		thermalConductivitySoil,
		heatCapacityWater,
		heatCapacityIce,
		density,
		etaPhi,
		deltax,
		deltat
		>;

	template <typename CELLS, typename PARAMETERS>
		void apply(CELLS& cells, PARAMETERS& parameters) any_platform
		{
			//using V = typename CELLS::template value_t<names::A>::value_t;
			using DESCRIPTOR_PFM = typename CELLS::template value_t<names::A>::descriptor_t;
			using DESCRIPTOR_HE = typename CELLS::template value_t<names::B>::descriptor_t;

			auto& cellPFE = cells.template get<names::A>();
			auto& cellHE = cells.template get<names::B>();

			V phi = cellPFE.computeRho();
			V theta = cellHE.computeRho(); // the current theta and phi 

			//V gradPhi[3] = {V{0}, V{0}, V{0}};
			//V gradTheta[3] = {V{0}, V{0}, V{0}};

			//descriptors::invCs2<T,DESCRIPTOR>(); 1/e^2 

			V prevphi = cellPFE.template getFieldPointer <descriptors::PREV_PHI>()[0];
			auto porosity = cellPFE. template getField<descriptors::POROSITY>();

			V kW = parameters.template get<thermalConductivityWater>();
			V kI = parameters.template get<thermalConductivityIce>();
			V kS = parameters.template get<thermalConductivitySoil>();

			V cW = parameters.template get<heatCapacityWater>();
			V cI = parameters.template get<heatCapacityIce>();
			V rho = parameters.template get<density>();
			V eta_phi = parameters.template get<etaPhi>();
			V dx = parameters.template get<deltax>();
			V dt = parameters.template get<deltat>();

			V C_phi = ((1+phi)/2.*cI+(1-phi)/2.*cW);
			V K_phi = (1-porosity)*kS + porosity*((1+phi)/2.*kI+(1-phi)/2.*kW);
			V effectiveDiffusivity = K_phi/(rho*C_phi);
			V etaTheta = descriptors::invCs2<V,DESCRIPTOR_HE>()*effectiveDiffusivity*(dt/(dx*dx))+1./2.;
			V omegaEff = 1./etaTheta;
			//for(int iPop =0; iPop < TDESCRIPTOR::q; iPop++)
			//{
			//	V w_i = descriptors::t<V,TDESCRIPTOR>(iPop);
			//	V cx = descriptors::c<DESCRIPTOR>(iPop, 0);
			//	V cy = descriptors::c<DESCRIPTOR>(iPop, 1);
			//	V cz = descriptors::c<DESCRIPTOR>(iPop, 2);

			//	V eqPhi   = w_i*phi;
			//	V eqTheta = w_i*theta;

			//	V neqPhi = cellPFE[iPop] - eqPhi;
			//	V neqTheta = cellHE[iPop] - eqTheta;

			//	gradPhi[0] += cx * neqPhi;
			//	gradPhi[1] += cy * neqPhi;
			//	gradPhi[2] += cz * neqPhi;

			//	gradTheta[0] += cx * neqTheta;
			//	gradTheta[1] += cy * neqTheta;
			//	gradTheta[2] += cz * neqTheta;
			//}

			//V scalePhi   = V{-1.} / (eta_phi *dx)*descriptors::invCs2<T,DESCRIPTOR>();
			//V scaleTheta = V{-1.} / (etaTheta*dx)*descriptors::invCs2<T,DESCRIPTOR>();

			//for (int d = 0; d < 3; d++) {
			//	gradPhi[d] *= scalePhi;
			//	gradTheta[d] *= scaleTheta;
			//}

			//cellPFE.template setField<descriptors::GRADIENT_PHI>({gradPhi[0],gradPhi[1],gradPhi[2]});
			//cellPFE.template setField<descriptors::GRADIENT_THETA>({gradTheta[0],gradTheta[1],gradTheta[2]});


			//V gradPhi_dot_gradTheta = 0.;
			//for (int d = 0; d < 3; d++) {
			//	gradPhi_dot_gradTheta = gradPhi_dot_gradTheta +  gradPhi[d]*gradTheta[d];
			//}
			for(int iPop =0; iPop < DESCRIPTOR_HE::q; iPop++)
			{
				cellHE[iPop] = cellHE[iPop] 
					- omegaEff*(cellHE[iPop]-descriptors::t<V,DESCRIPTOR_HE>(iPop)*(theta-V{1.}))
					+ descriptors::t<V,DESCRIPTOR_HE>(iPop)*V{0.5}*(phi-prevphi);
				//+ descriptors::t<V,TDESCRIPTOR>(iPop)*V{0.5}*effectiveDiffusivity*(rho*cI-rho*cW)/(rho*C_phi)*gradPhi_dot_gradTheta*dt;
			}
		}
	};
}
#endif
