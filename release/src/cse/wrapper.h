/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender, Shota Ito
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

#ifndef CSE_WRAPPER_H
#define CSE_WRAPPER_H

#include "dynamics/interface.h"

namespace olb {

namespace dynamics {

/// Marker for the default non-optimized specialization of olb::CSE<>
struct not_cse_optimized_tag { };

/// To be specialized for automatically generated CSE-optimized DYNAMICS::collide
/**
 * Falls back to unoptimized collision by default
 **/
template <typename DYNAMICS>
struct CSE : public not_cse_optimized_tag {
  template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    return DYNAMICS().collide(cell, parameters);
  };
};

/// Exposes whether a auto-generated CSE specialization is available
template <typename DYNAMICS>
using is_cse_optimized = std::integral_constant<
  bool,
  !std::is_base_of_v<not_cse_optimized_tag, CSE<DYNAMICS>>
>;

/// Evaluates to true iff a auto-generated CSE specialization is available
template <typename DYNAMICS>
constexpr static bool is_cse_optimized_v = is_cse_optimized<DYNAMICS>::value;

}

/// Wrapper for auto-generated CSE-optimized DYNAMICS
template <concepts::IntrospectableDynamics DYNAMICS>
class CSE : public Dynamics<typename DYNAMICS::value_t,
                            typename DYNAMICS::descriptor_t> {
private:
  using T = typename DYNAMICS::value_t;
  using DESCRIPTOR = typename DYNAMICS::descriptor_t;

public:
  using value_t = typename DYNAMICS::value_t;
  using descriptor_t = typename DYNAMICS::descriptor_t;

  using MomentaF = typename DYNAMICS::MomentaF;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CSE<typename DYNAMICS::template exchange_value_type<NEW_T>>;

  template <typename NEW_MOMENTA>
  using exchange_momenta = CSE<typename DYNAMICS::template exchange_momenta<NEW_MOMENTA>>;

  constexpr static bool is_vectorizable = dynamics::is_vectorizable_v<DYNAMICS>;
  constexpr static bool is_optimized    = dynamics::is_cse_optimized_v<DYNAMICS>;

  std::type_index id() override {
    return typeid(CSE);
  }

  std::string getName() const override {
    return "dynamics::cse<" + DYNAMICS().getName() + ">";
  };

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    return dynamics::CSE<DYNAMICS>().collide(cell, parameters);
  };

  /// Return true iff FIELD is a parameter
  template <typename FIELD>
  constexpr bool hasParameter() const {
    return DYNAMICS().template hasParameter<FIELD>();
  };

  /// Interim workaround for accessing dynamics parameters in legacy post processors
  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return DYNAMICS().getParameters(block);
  }

  void initialize(Cell<T,DESCRIPTOR>& cell) override {
    DYNAMICS().initialize(cell);
  };
  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    DYNAMICS().computeEquilibrium(cell, rho, u, fEq);
  };
  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override {
    return DYNAMICS().computeRho(cell);
  };
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    DYNAMICS().computeU(cell, u);
  };
  void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const override {
    DYNAMICS().computeJ(cell, j);
  };
  void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                     T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    DYNAMICS().computeStress(cell, rho, u, pi);
  };
  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    DYNAMICS().computeRhoU(cell, rho, u);
  };
  void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                         T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    DYNAMICS().computeAllMomenta(cell, rho, u, pi);
  };
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override {
    DYNAMICS().defineRho(cell, rho);
  };
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override {
    DYNAMICS().defineU(cell, u);
  };
  void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override {
    DYNAMICS().defineRhoU(cell, rho, u);
  };
  void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                        T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) override {
    DYNAMICS().defineAllMomenta(cell, rho, u, pi);
  };
  void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    DYNAMICS().inverseShiftRhoU(cell, rho, u);
  }
};


namespace operators {

/// Marker for the default non-optimized specialization of olb::CSE<>
struct not_cse_optimized_tag { };

/// To be specialized for automatically generated CSE-optimized OPERATOR::apply
/**
 * Falls back to unoptimized operator by default
 **/
template <typename OPERATOR, typename DESCRIPTOR>
struct CSE_O : public not_cse_optimized_tag { };

template <typename OPERATOR, typename DESCRIPTOR>
requires (OPERATOR::scope == OperatorScope::PerCell)
struct CSE_O<OPERATOR,DESCRIPTOR> : public not_cse_optimized_tag {
  template <concepts::Cell CELL>
  void apply(CELL& cell) any_platform {
    OPERATOR().apply(cell);
  };
};

template <typename OPERATOR, typename DESCRIPTOR>
requires (OPERATOR::scope == OperatorScope::PerCellWithParameters)
struct CSE_O<OPERATOR,DESCRIPTOR> : public not_cse_optimized_tag {
  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    OPERATOR().apply(cell, parameters);
  };
};

/// Exposes whether a auto-generated CSE specialization is available
template <typename OPERATOR, typename DESCRIPTOR>
using is_cse_optimized = std::integral_constant<
  bool,
  (!std::is_base_of_v<not_cse_optimized_tag, CSE_O<OPERATOR,DESCRIPTOR>> &&
   concepts::CellOperator<OPERATOR>)
>;

/// Evaluates to true iff a auto-generated CSE specialization is available
template <typename OPERATOR, typename DESCRIPTOR>
constexpr static bool is_cse_optimized_v = is_cse_optimized<OPERATOR,DESCRIPTOR>::value;

}

/// Wrapper for auto-generated CSE-optimized OPERATOR
template <typename OPERATOR, typename DESCRIPTOR>
struct CSE_O { };

template <typename OPERATOR, typename DESCRIPTOR>
requires (OPERATOR::scope == OperatorScope::PerCell)
struct CSE_O<OPERATOR,DESCRIPTOR> {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return OPERATOR().getPriority();
  }

  template <concepts::Cell CELL>
  void apply(CELL& cell) any_platform {
    operators::CSE_O<OPERATOR,DESCRIPTOR>().apply(cell);
  };
};

template <typename OPERATOR, typename DESCRIPTOR>
requires (OPERATOR::scope == OperatorScope::PerCellWithParameters)
struct CSE_O<OPERATOR,DESCRIPTOR> {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return OPERATOR().getPriority();
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    operators::CSE_O<OPERATOR,DESCRIPTOR>().apply(cell, parameters);
  };
};

}

#endif
