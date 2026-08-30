/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-2024 Adrian Kummerlaender
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

#ifndef DYNAMICS_MODIFIERS_H
#define DYNAMICS_MODIFIERS_H

namespace olb {

namespace dynamics {

/// Set PARAMETER of DYNAMICS from CELL (for CustomCollision-based DYNAMICS)
/**
 * Allows for e.g. overriding DYNAMICS-wide relaxation frequency by per-cell values:
 *
 * \code
 * template <typename T, typename DESCRIPTOR>
 * using PerCellOmegaSourcedAdvectionDiffusionBGKdynamics = dynamics::ParameterFromCell<
 *   descriptors::OMEGA,
 *   SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>
 * >;
 * \endcode
 **/
template <typename PARAMETER, typename DYNAMICS>
struct ParameterFromCell final : public CustomCollision<
  typename DYNAMICS::value_t,
  typename DYNAMICS::descriptor_t,
  typename DYNAMICS::MomentaF::abstract
> {
  using value_t = typename DYNAMICS::value_t;
  using descriptor_t = typename DYNAMICS::descriptor_t;

  using MomentaF = typename DYNAMICS::MomentaF;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = ParameterFromCell<PARAMETER,typename DYNAMICS::template exchange_value_type<NEW_T>>;

  template<typename M>
  using exchange_momenta = ParameterFromCell<PARAMETER,typename DYNAMICS::template exchange_momenta<M>>;

  std::type_index id() override {
    return typeid(ParameterFromCell);
  };

  AbstractParameters<value_t,descriptor_t>&
  getParameters(BlockLattice<value_t,descriptor_t>& block) override {
    return block.template getData<OperatorParameters<ParameterFromCell>>();
  }

  template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    parameters.template set<PARAMETER>(
      cell.template getField<PARAMETER>());
    return DYNAMICS().collide(cell, parameters);
  };

  void computeEquilibrium(ConstCell<value_t,descriptor_t>& cell, value_t rho, const value_t u[descriptor_t::d], value_t fEq[descriptor_t::q]) const override {
    DYNAMICS().computeEquilibrium(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "ParameterFromCell<" + DYNAMICS().getName() + ">";
  };

};

template <concepts::IntrospectableDynamics DYNAMICS>
struct StoreStatisticInField final : public CustomCollision<
  typename DYNAMICS::value_t,
  typename DYNAMICS::descriptor_t,
  typename DYNAMICS::MomentaF::abstract
> {
  using value_t = typename DYNAMICS::value_t;
  using descriptor_t = typename DYNAMICS::descriptor_t;

  using MomentaF = typename DYNAMICS::MomentaF;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = StoreStatisticInField<typename DYNAMICS::template exchange_value_type<NEW_T>>;

  template<typename M>
  using exchange_momenta = StoreStatisticInField<typename DYNAMICS::template exchange_momenta<M>>;

  std::type_index id() override {
    return typeid(StoreStatisticInField);
  };

  AbstractParameters<value_t,descriptor_t>&
  getParameters(BlockLattice<value_t,descriptor_t>& block) override {
    return block.template getData<OperatorParameters<StoreStatisticInField>>();
  }

  template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    auto statistic = DYNAMICS().collide(cell, parameters);
    cell.template setField<descriptors::STATISTIC>({statistic.rho, statistic.uSqr});
    return statistic;
  };

  void computeEquilibrium(ConstCell<value_t,descriptor_t>& cell, value_t rho, const value_t u[descriptor_t::d], value_t fEq[descriptor_t::q]) const override {
    DYNAMICS().computeEquilibrium(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "StoreStatisticInField<" + DYNAMICS().getName() + ">";
  };

};

}

}

#endif
