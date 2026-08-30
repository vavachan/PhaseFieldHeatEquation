/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *                2022 Nando Suntoyo, Adrian Kummerlaender
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_H
#define ADVECTION_DIFFUSION_DYNAMICS_H

#include "descriptor/descriptor.h"
#include "dynamics/dynamics.h"
#include "core/unitConverter.h"
#include "collisionMRT.h"

namespace olb {

namespace TotalEnthalpy {
  struct T_S               : public descriptors::FIELD_BASE<1> { };
  struct T_L               : public descriptors::FIELD_BASE<1> { };
  struct CP_S              : public descriptors::FIELD_BASE<1> { };
  struct CP_L              : public descriptors::FIELD_BASE<1> { };
  struct LAMBDA_S          : public descriptors::FIELD_BASE<1> { };
  struct LAMBDA_L          : public descriptors::FIELD_BASE<1> { };
  struct L                 : public descriptors::FIELD_BASE<1> { };
}

namespace Light {
  struct ABSORPTION               : public descriptors::FIELD_BASE<1> { };
  struct SCATTERING               : public descriptors::FIELD_BASE<1> { };
  struct ANISOMATRIX              : public descriptors::FIELD_BASE_CUSTOM_SIZE {
    template <unsigned D, unsigned Q, unsigned...>
    static constexpr unsigned size() {
      return Q*Q;
    }

    template <typename DESCRIPTOR>
    static constexpr unsigned size() {
      return DESCRIPTOR::q * DESCRIPTOR::q;
    }

    template <typename T, typename DESCRIPTOR>
    static constexpr auto getInitialValue() {
      return Vector<value_type<T>, DESCRIPTOR::template size<ANISOMATRIX>()>{};
    }

    template <typename T, typename DESCRIPTOR, typename FIELD>
    static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
      return true;
    }
  };
}


struct AdvectionDiffusionExternalVelocityCollision {
  static std::string getName() {
    return "AdvectionDifffusionExternalVelocityCollision";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename COLLISION::template type<
    DESCRIPTOR,
    momenta::Tuple<
      typename MOMENTA::density,
      momenta::FixedVelocityMomentum,
      typename MOMENTA::stress,
      typename MOMENTA::definition
    >,
    EQUILIBRIUM
  >;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
using AdvectionDiffusionRLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::AdvectionDiffusionRLB,
  AdvectionDiffusionExternalVelocityCollision
>;

template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
class CombinedAdvectionDiffusionRLBdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CombinedAdvectionDiffusionRLBdynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA
  >;

  template <typename M>
  using exchange_momenta = CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedAdvectionDiffusionRLBdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<CombinedAdvectionDiffusionRLBdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V jNeq[DESCRIPTOR::d] { };

    const V rho = MomentaF().computeRho(cell);
    const auto u = cell.template getField<descriptors::VELOCITY>();
    MomentaF().computeJ(cell, jNeq);

    for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
      jNeq[iD] -= u[iD] * rho;
    }

    V fEq[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, rho, u, fEq);
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop] + equilibrium<DESCRIPTOR>::template fromJneqToFneq<V>(iPop, jNeq);
    }

    return typename DYNAMICS::CollisionO().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "CombinedAdvectionDiffusionRLBdynamics<" + MomentaF().getName() + ">";
  };

};

// ========= the BGK advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
using AdvectionDiffusionBGKdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::BGK,
  AdvectionDiffusionExternalVelocityCollision
>;


// ========= the TRT advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
using AdvectionDiffusionTRTdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::TRT,
  AdvectionDiffusionExternalVelocityCollision
>;


// ======= BGK advection diffusion dynamics with source term  ======//
// following Seta, T. (2013). Implicit temperature-correction-based
// immersed-boundary thermal lattice Boltzmann method for the simulation
// of natural convection. Physical Review E, 87(6), 063304.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct SourcedAdvectionDiffusionBGKdynamics final : public dynamics::CustomCollision<
  T,DESCRIPTOR,
  momenta::Tuple<
    momenta::SourcedDensity<typename MOMENTA::density>,
    typename MOMENTA::momentum,
    typename MOMENTA::stress,
    typename MOMENTA::definition
  >
> {
  using MomentaF = typename momenta::Tuple<
    momenta::SourcedDensity<typename MOMENTA::density>,
    typename MOMENTA::momentum,
    typename MOMENTA::stress,
    typename MOMENTA::definition
  >::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,MOMENTA>;

  using parameters = meta::list<descriptors::OMEGA>;

  template <typename NEW_T>
  using exchange_value_type = SourcedAdvectionDiffusionBGKdynamics<
    NEW_T, DESCRIPTOR,
    MOMENTA
  >;

  template<typename M>
  using exchange_momenta = SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(SourcedAdvectionDiffusionBGKdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<SourcedAdvectionDiffusionBGKdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    const auto u = cell.template getField<descriptors::VELOCITY>();
    const V temperature = MomentaF().computeRho(cell);
    const V omega = parameters.template get<descriptors::OMEGA>();

    const V uSqr = lbm<DESCRIPTOR>::adeBgkCollision(cell, temperature, u, omega);
    const V sourceMod = cell.template getField<descriptors::SOURCE>() * (V{1} - V{0.5} * omega);

    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      cell[iPop] += sourceMod * descriptors::t<T,DESCRIPTOR>(iPop);
    }

    return {temperature, uSqr};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "SourcedAdvectionDiffusionBGKdynamics<" + MomentaF().getName() + ">";
  };
};


// ======= BGK advection diffusion dynamics with source term  ======//
// following Seta, T. (2013). Implicit temperature-correction-based
// immersed-boundary thermal lattice Boltzmann method for the simulation
// of natural convection. Physical Review E, 87(6), 063304.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct SourcedLimitedAdvectionDiffusionBGKdynamics final : public dynamics::CustomCollision<
  T,DESCRIPTOR,
  momenta::Tuple<
    momenta::SourcedDensity<typename MOMENTA::density>,
    typename MOMENTA::momentum,
    typename MOMENTA::stress,
    typename MOMENTA::definition
  >
> {
  using MomentaF = typename momenta::Tuple<
    momenta::SourcedDensity<typename MOMENTA::density>,
    typename MOMENTA::momentum,
    typename MOMENTA::stress,
    typename MOMENTA::definition
  >::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,MOMENTA>;

  using parameters = meta::list<descriptors::OMEGA>;

  constexpr static bool is_vectorizable = false;

  template <typename NEW_T>
  using exchange_value_type = SourcedLimitedAdvectionDiffusionBGKdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = SourcedLimitedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(SourcedLimitedAdvectionDiffusionBGKdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<SourcedLimitedAdvectionDiffusionBGKdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    const auto u = cell.template getField<descriptors::VELOCITY>();
    V temperature = MomentaF().computeRho(cell);
    if (temperature < V{1.e-8}) temperature = V{1.e-8};
    const V omega = cell.template getField<descriptors::OMEGA>();

    const V uSqr = lbm<DESCRIPTOR>::adeBgkCollision(cell, temperature, u, omega);
    const V sourceMod = cell.template getField<descriptors::SOURCE>() * (V{1} - V{0.5} * omega);

    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      cell[iPop] += sourceMod * descriptors::t<T,DESCRIPTOR>(iPop);
    }

    return {temperature, uSqr};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "SourcedLimitedAdvectionDiffusionBGKdynamics<" + MomentaF().getName() + ">";
  };
};

// ======= BGK advection diffusion dynamics for solid-liquid phase change  ======//
// following Huang, R. (2015). Phase interface effects in the total
// enthalpy-based lattice Boltzmann model for solid–liquid phase change.
// Journal of Computational Physics, 294, 345-362.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct TotalEnthalpyAdvectionDiffusionBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  static constexpr bool is_vectorizable = false;

  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,MOMENTA>;

  using parameters = meta::list<
    descriptors::OMEGA,
    TotalEnthalpy::T_S,
    TotalEnthalpy::T_L,
    TotalEnthalpy::CP_S,
    TotalEnthalpy::CP_L,
    TotalEnthalpy::LAMBDA_S,
    TotalEnthalpy::LAMBDA_L,
    TotalEnthalpy::L
  >;

  template <typename NEW_T>
  using exchange_value_type = TotalEnthalpyAdvectionDiffusionBGKdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = TotalEnthalpyAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(TotalEnthalpyAdvectionDiffusionBGKdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<TotalEnthalpyAdvectionDiffusionBGKdynamics>>();
  }

  template<typename V, typename PARAMETERS, typename ENTHALPY>
  V computeTemperature(const PARAMETERS& parameters, const ENTHALPY& enthalpy) const any_platform
  {
    using namespace TotalEnthalpy;

    const V cp_s = parameters.template get<CP_S>();
    const V cp_l = parameters.template get<CP_L>();
    const V T_s  = parameters.template get<T_S>();
    const V T_l  = parameters.template get<T_L>();
    const V l    = parameters.template get<L>();
    const V H_s  = cp_s * T_s;
    const V H_l  = cp_l * T_l + l;
    V temperature{};

    if (enthalpy <= H_s) {
      temperature = T_s - (H_s - enthalpy) / cp_s;
    }
    else if (enthalpy >= H_l) {
      temperature = T_l + (enthalpy - H_l) / cp_l;
    }
    else {
      temperature = (H_l - enthalpy) / (H_l - H_s) * T_s + (enthalpy - H_s) / (H_l - H_s) * T_l;
    }
    return temperature;
  }

  template<typename V, typename PARAMETERS, typename ENTHALPY>
  V computeLiquidFraction(const PARAMETERS& parameters, const ENTHALPY& enthalpy) const any_platform
  {
    using namespace TotalEnthalpy;

    const V cp_s = parameters.template get<CP_S>();
    const V cp_l = parameters.template get<CP_L>();
    const V T_s  = parameters.template get<T_S>();
    const V T_l  = parameters.template get<T_L>();
    const V l    = parameters.template get<L>();
    const V H_s  = cp_s * T_s;
    const V H_l  = cp_l * T_l + l;
    V liquid_fraction{};

    if (enthalpy <= H_s) {
      liquid_fraction = 0.;
    }
    else if (enthalpy >= H_l) {
      liquid_fraction = 1.;
    }
    else {
      liquid_fraction = (enthalpy - H_s) / l;
    }
    return liquid_fraction;
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    using namespace TotalEnthalpy;

    const V lambda_s = parameters.template get<LAMBDA_S>();
    const V lambda_l = parameters.template get<LAMBDA_L>();
    const V cp_s     = parameters.template get<CP_S>();
    const V cp_l     = parameters.template get<CP_L>();
    const V cp_ref   = V{2} * cp_s * cp_l / (cp_s + cp_l);

    const V enthalpy = MomentaF().computeRho( cell );
    const V temperature = computeTemperature<V>( parameters, enthalpy );
    const V liquid_fraction = computeLiquidFraction<V>( parameters, enthalpy );
    const V lambda = (V{1} - liquid_fraction) * lambda_s + liquid_fraction * lambda_l;
    const V cp = (V{1} - liquid_fraction) * cp_s + liquid_fraction * cp_l;
    const V omega = V{1} / ( lambda / cp_ref * descriptors::invCs2<T,DESCRIPTOR>() + V{0.5} );

    const auto u = cell.template getFieldPointer<descriptors::VELOCITY>();

    const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

    const V f_eq = enthalpy - cp_ref * temperature
                   + cp * temperature * descriptors::t<T,DESCRIPTOR>(0) * ( cp_ref / cp
                       - descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * uSqr )
                   - descriptors::t<T,DESCRIPTOR>(0);
    cell[0] *= V{1} - omega;
    cell[0] += omega * f_eq;
    for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      V c_u{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      const V f_eq = cp * temperature * descriptors::t<T,DESCRIPTOR>(iPop) * ( cp_ref / cp + descriptors::invCs2<T,DESCRIPTOR>() * c_u
                     + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * c_u *c_u
                     - descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * uSqr )
                     - descriptors::t<T,DESCRIPTOR>(iPop);
      cell[iPop] *= V{1} - omega;
      cell[iPop] += omega * f_eq;
    }
    return {enthalpy, uSqr};
  };

  std::string getName() const override {
    return "TotalEnthalpyAdvectionDiffusionBGKdynamics<" + MomentaF().getName() + ">";
  };
};

// ======= TRT advection diffusion dynamics for solid-liquid phase change  ======//
// following Huang, R. (2015). Phase interface effects in the total
// enthalpy-based lattice Boltzmann model for solid–liquid phase change.
// Journal of Computational Physics, 294, 345-362.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct TotalEnthalpyAdvectionDiffusionTRTdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  static constexpr bool is_vectorizable = false;

  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,MOMENTA>;

  using parameters = meta::list<
    descriptors::OMEGA,
    collision::TRT::MAGIC,
    TotalEnthalpy::T_S,
    TotalEnthalpy::T_L,
    TotalEnthalpy::CP_S,
    TotalEnthalpy::CP_L,
    TotalEnthalpy::LAMBDA_S,
    TotalEnthalpy::LAMBDA_L,
    TotalEnthalpy::L
  >;

  template <typename NEW_T>
  using exchange_value_type = TotalEnthalpyAdvectionDiffusionTRTdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = TotalEnthalpyAdvectionDiffusionTRTdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(TotalEnthalpyAdvectionDiffusionTRTdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<TotalEnthalpyAdvectionDiffusionTRTdynamics>>();
  }

  template<typename V, typename PARAMETERS, typename ENTHALPY>
  V computeTemperature(const PARAMETERS& parameters, const ENTHALPY& enthalpy) const any_platform
  {
    using namespace TotalEnthalpy;

    const V cp_s = parameters.template get<CP_S>();
    const V cp_l = parameters.template get<CP_L>();
    const V T_s  = parameters.template get<T_S>();
    const V T_l  = parameters.template get<T_L>();
    const V l    = parameters.template get<L>();
    const V H_s  = cp_s * T_s;
    const V H_l  = cp_l * T_l + l;
    V temperature{};

    if (enthalpy <= H_s) {
      temperature = T_s - (H_s - enthalpy) / cp_s;
    }
    else if (enthalpy >= H_l) {
      temperature = T_l + (enthalpy - H_l) / cp_l;
    }
    else {
      temperature = (H_l - enthalpy) / (H_l - H_s) * T_s + (enthalpy - H_s) / (H_l - H_s) * T_l;
    }
    return temperature;
  }

  template<typename V, typename PARAMETERS, typename ENTHALPY>
  V computeLiquidFraction(const PARAMETERS& parameters, const ENTHALPY& enthalpy) const any_platform
  {
    using namespace TotalEnthalpy;

    const V cp_s = parameters.template get<CP_S>();
    const V cp_l = parameters.template get<CP_L>();
    const V T_s  = parameters.template get<T_S>();
    const V T_l  = parameters.template get<T_L>();
    const V l    = parameters.template get<L>();
    const V H_s  = cp_s * T_s;
    const V H_l  = cp_l * T_l + l;
    V liquid_fraction{};

    if (enthalpy <= H_s) {
      liquid_fraction = 0.;
    }
    else if (enthalpy >= H_l) {
      liquid_fraction = 1.;
    }
    else {
      liquid_fraction = (enthalpy - H_s) / l;
    }
    return liquid_fraction;
  }

  template<typename V, typename PARAMETERS, typename RHO, typename U>
  V computeEquilibrium(int iPop, const PARAMETERS& parameters, RHO& rho, U& u) const any_platform
  {
    using namespace TotalEnthalpy;

    const V temperature = computeTemperature<V>(parameters, rho);
    const V liquid_fraction = computeLiquidFraction<V>(parameters, rho);

    const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    const V cp_s = parameters.template get<CP_S>();
    const V cp_l = parameters.template get<CP_L>();
    const V cp = (V{1} - liquid_fraction) * cp_s + liquid_fraction * cp_l;
    const V cp_ref = V{2} * cp_s * cp_l / (cp_s + cp_l);

    V c_u{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }

    V f_eq{};
    if (iPop == 0) {
      f_eq = rho - cp_ref * temperature
             + cp * temperature * descriptors::t<T,DESCRIPTOR>(0) * ( cp_ref / cp
                 - descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * uSqr )
             - descriptors::t<T,DESCRIPTOR>(0);
    }
    else {
      f_eq = cp * temperature * descriptors::t<T,DESCRIPTOR>(iPop) * ( cp_ref / cp + descriptors::invCs2<T,DESCRIPTOR>() * c_u
             + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * c_u *c_u
             - descriptors::invCs2<T,DESCRIPTOR>() * V{0.5} * uSqr )
             - descriptors::t<T,DESCRIPTOR>(iPop);
    }

    return f_eq;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    using namespace TotalEnthalpy;

    const V lambda_s = parameters.template get<LAMBDA_S>();
    const V lambda_l = parameters.template get<LAMBDA_L>();
    const V cp_s     = parameters.template get<CP_S>();
    const V cp_l     = parameters.template get<CP_L>();
    const V cp_ref   = V{2} * cp_s * cp_l / (cp_s + cp_l);

    const V enthalpy = MomentaF().computeRho( cell );
    const V liquid_fraction = computeLiquidFraction<V>( parameters, enthalpy );
    const V lambda = (V{1} - liquid_fraction) * lambda_s + liquid_fraction * lambda_l;
    const V omega = V{1} / ( lambda / cp_ref * descriptors::invCs2<T,DESCRIPTOR>() + V{0.5} );
    const V magic_parameter = parameters.template get<collision::TRT::MAGIC>();
    const V omega2 = V{1} / (magic_parameter/(V{1}/omega-V{0.5})+V{0.5});

    const auto u = cell.template getField<descriptors::VELOCITY>();

    const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);

    V fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q] { };
    V fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q] { };

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fPlus[iPop] = 0.5 * ( cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
      fMinus[iPop] = 0.5 * ( cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
      fEq[iPop] = computeEquilibrium<V>(iPop, parameters, enthalpy, u);
    }

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
      fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
    }

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] -= omega2 * (fPlus[iPop] - fEqPlus[iPop]) + omega * (fMinus[iPop] - fEqMinus[iPop]);
    }

    return {enthalpy, uSqr};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "TotalEnthalpyAdvectionDiffusionTRTdynamics<" + MomentaF().getName() + ">";
  };
};

namespace descriptors {

struct INTERFACE_THICKNESS : public descriptors::FIELD_BASE<1> { };

}

// ======= BGK advection diffusion dynamics for phase field equation  ======//
// following Fakhari, Abbas, et al. (2017). Improved locality of the phase-field
// lattice-Boltzmann model for immiscible fluids at high density ratios.
// Physical Review E 96.5, 053301.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct PhaseFieldAdvectionDiffusionBGKdynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,momenta::BulkTuple>;

  using parameters = meta::list<descriptors::OMEGA,descriptors::INTERFACE_THICKNESS>;

  template <typename NEW_T>
  using exchange_value_type = PhaseFieldAdvectionDiffusionBGKdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = PhaseFieldAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(PhaseFieldAdvectionDiffusionBGKdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldAdvectionDiffusionBGKdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) {
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V interface_thickness = parameters.template get<descriptors::INTERFACE_THICKNESS>();

    const V phi = MomentaF().computeRho( cell );
    const auto u = cell.template getField<descriptors::VELOCITY>();
    const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    const auto mobility = (V{1} / omega - V{0.5}) / descriptors::invCs2<V,DESCRIPTOR>();

    const auto n = cell.template getFieldPointer<descriptors::INTERPHASE_NORMAL>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V c_u{};
      V c_n{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        c_n += descriptors::c<DESCRIPTOR>(iPop,iD)*n[iD];
      }
      V f_eq = equilibrium<DESCRIPTOR>::firstOrder(iPop, phi, u);
      f_eq += descriptors::t<V,DESCRIPTOR>(iPop) * mobility * descriptors::invCs2<V,DESCRIPTOR>()
            * (V{4} * phi * (V{1} - phi) / interface_thickness) * c_n;

      cell[iPop] *= V{1} - omega;
      cell[iPop] += omega * f_eq;
    }

    return {phi, uSqr};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "PhaseFieldAdvectionDiffusionBGKdynamics<" + MomentaF().getName() + ">";
  };
};

// ========= the BGK advection diffusion Stokes drag dynamics  ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct ParticleAdvectionDiffusionBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
public:
  using MomentaF = typename momenta::Tuple<
    typename MOMENTA::density,
    momenta::FixedVelocityMomentum,
    typename MOMENTA::stress,
    typename MOMENTA::definition
  >::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::FirstOrder::template type<DESCRIPTOR,MOMENTA>;

  using parameters = meta::list<descriptors::OMEGA,descriptors::LATTICE_TIME>;

  template <typename NEW_T>
  using exchange_value_type = ParticleAdvectionDiffusionBGKdynamics<NEW_T,DESCRIPTOR,MOMENTA>;

  template <typename M>
  using exchange_momenta = ParticleAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(ParticleAdvectionDiffusionBGKdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ParticleAdvectionDiffusionBGKdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    const V    omega  = parameters.template get<descriptors::OMEGA>();
    const auto time = parameters.template get<descriptors::LATTICE_TIME>();
    const auto u = (time % 2 == 0) ? cell.template getField<descriptors::VELOCITY>()
                                   : cell.template getField<descriptors::VELOCITY2>();
    V rho = MomentaF().computeRho(cell);
    V uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    return {rho, uSqr};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "ParticleAdvectionDiffusionBGKdynamics<" + MomentaF().getName() + ">";
  };

};

//==================================================================//
//============= BGK Model for Advection diffusion anisotropic ===//
//==================================================================//
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct RTLBMdynamicsMcHardy final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA>
{
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::ZerothOrder::template type<DESCRIPTOR,MOMENTA>;
  using parameters = meta::list<Light::ABSORPTION,Light::SCATTERING,Light::ANISOMATRIX>;

  template <typename NEW_T>
  using exchange_value_type = RTLBMdynamicsMcHardy<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = RTLBMdynamicsMcHardy<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(RTLBMdynamicsMcHardy);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<RTLBMdynamicsMcHardy>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    const V absorption  = parameters.template get<Light::ABSORPTION>();
    const V scattering  = parameters.template get<Light::SCATTERING>();
    auto anisoMatrix  = parameters.template get<Light::ANISOMATRIX>();
    V irradiance = MomentaF().computeRho(cell);
    Vector<V, DESCRIPTOR::q> feq;
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      for ( int jPop = 0; jPop < DESCRIPTOR::q; ++jPop ) {
        V fFull = cell[jPop] + descriptors::t<V,DESCRIPTOR>(jPop);
        feq[iPop] += fFull * anisoMatrix[jPop+iPop];
      }
      feq[iPop] *= descriptors::t<V,DESCRIPTOR>(iPop);
    }

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      V normC = util::norm<DESCRIPTOR::d>(descriptors::c<DESCRIPTOR>(iPop));
      V fFull = cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop);
      cell[iPop] = fFull
                 - normC * (absorption + scattering) * ( fFull- feq[iPop] )
                 - absorption * normC * fFull
                 - descriptors::t<V,DESCRIPTOR>(iPop);
    }

    return {irradiance, V(0)};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T irradiance, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, irradiance, u, fEq);
  };

  std::string getName() const override {
    return "RTLBMdynamicsMcHardy<" + MomentaF().getName() + ">";
  };
};


//=======================================================================================//
//======= BGK Model for Advection diffusion anisotropic with Runge Kutta scheme =========//
//=======================================================================================//

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
struct RTLBMdynamicsMcHardyRK final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA>
{
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::ZerothOrder::template type<DESCRIPTOR,MOMENTA>;
  using parameters = meta::list<Light::ABSORPTION,Light::SCATTERING,Light::ANISOMATRIX>;

  template <typename NEW_T>
  using exchange_value_type = RTLBMdynamicsMcHardyRK<NEW_T,DESCRIPTOR,MOMENTA>;

  template<typename M>
  using exchange_momenta = RTLBMdynamicsMcHardyRK<T,DESCRIPTOR,M>;

  std::type_index id() override {
    return typeid(RTLBMdynamicsMcHardyRK);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<RTLBMdynamicsMcHardyRK>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    const V absorption  = parameters.template get<Light::ABSORPTION>();
    const V scattering  = parameters.template get<Light::SCATTERING>();
    auto anisoMatrix  = parameters.template get<Light::ANISOMATRIX>();
    V irradiance = MomentaF().computeRho(cell);
    Vector<V, DESCRIPTOR::q> feq;
    Vector<V, DESCRIPTOR::q> f_pre_collision;

    // calculate f_pre_collision
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      f_pre_collision[iPop] = cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // preparation for RK
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Runge-Kutta 4
    computeEquilibriumAniso(cell, feq, anisoMatrix);
    auto k1 = doCollision(cell, feq, absorption, scattering);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = f_pre_collision[iPop] + 0.5 * k1[iPop];
    }

    computeEquilibriumAniso(cell, feq, anisoMatrix);
    auto k2 = doCollision(cell, feq, absorption, scattering);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = f_pre_collision[iPop] + 0.5 * k2[iPop];
    }

    computeEquilibriumAniso(cell, feq, anisoMatrix);
    auto k3 = doCollision(cell, feq, absorption, scattering);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = f_pre_collision[iPop] + k3[iPop];
    }

    computeEquilibriumAniso(cell, feq, anisoMatrix);
    auto k4 = doCollision(cell, feq, absorption, scattering);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = f_pre_collision[iPop] + (1.0 / 6.0) * (k1[iPop] + 2 * k2[iPop] + 2 * k3[iPop] + k4[iPop])
                 - descriptors::t<V,DESCRIPTOR>(iPop);
    }

    return {irradiance, V(0)};
  };

  template <typename CELL, typename MATRIX, typename V=typename CELL::value_t>
  void computeEquilibriumAniso(CELL& cell, Vector<V, DESCRIPTOR::q>& feq, MATRIX& anisoMatrix) const {
//    feq.fill(V());
    for (int i = 0; i < DESCRIPTOR::q; ++i) {
      feq[i] = V();
    }
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        feq[iPop] += cell[jPop] * anisoMatrix[jPop+iPop];
      }
      feq[iPop] *= descriptors::t<V, DESCRIPTOR>(iPop);
    }
  };

  template <typename CELL, typename VECTOR, typename V=typename CELL::value_t>
  VECTOR doCollision(CELL& cell, VECTOR& feq, V absorption, V scattering) const {
    VECTOR k;
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      V normC = util::norm<DESCRIPTOR::d>(descriptors::c<DESCRIPTOR>(iPop));
      V fFull = cell[iPop];
      k[iPop] = -normC * (absorption + scattering) * fFull
                + normC * scattering * feq[iPop];
    }
    return k;
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T irradiance, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, irradiance, u, fEq);
  };

  std::string getName() const override {
    return "RTLBMdynamicsMcHardyRK<" + MomentaF().getName() + ">";
  };
};




// ========= the MRT advection diffusion dynamics ========//
/// This approach is based on the multi-distribution LBM model.
/// The coupling is done using the Boussinesq approximation
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
using AdvectionDiffusionMRTdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::MRT,
  AdvectionDiffusionExternalVelocityCollision
>;

template <typename T, typename DESCRIPTOR>
using NoCollideDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::BulkTuple,
  equilibria::None,
  collision::None
>;
} // namespace olb

#endif
