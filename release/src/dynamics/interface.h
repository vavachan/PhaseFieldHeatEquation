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

#ifndef DYNAMICS_INTERFACE_H
#define DYNAMICS_INTERFACE_H

#include "lbm.h"

#include "core/concepts.h"

#include "momenta/interface.h"
#include "momenta/aliases.h"

namespace olb {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR> class ConstCell;
template <typename T, typename DESCRIPTOR> class Cell;
template <typename T, typename DESCRIPTOR> class BlockLattice;

/// Return value of any collision
/**
 * Used for lattice statistics if enabled
 **/
template <typename T>
struct CellStatistic {
  T rho;
  T uSqr;

  operator bool() const any_platform {
    return rho != T{-1} && uSqr != T{-1};
  }
};

/// Interface for per-cell dynamics
template <typename T, typename DESCRIPTOR>
struct Dynamics {
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  virtual ~Dynamics() any_platform { }

  /// Expose unique type-identifier for RTTI
  virtual std::type_index id() = 0;
  /// Return human-readable name
  virtual std::string getName() const {
    return "Dynamics";
  };

  /// Initialize dynamics-specific data for cell
  virtual void initialize(Cell<T,DESCRIPTOR>& cell) { };

  /// Parameters access for legacy post processors
  virtual AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) = 0;

  /// Perform purely-local collision step on Cell interface (legacy, to be deprecated)
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) {
    throw std::bad_function_call();
  };

  /// Compute particle density
  virtual T computeRho(ConstCell<T,DESCRIPTOR>& cell) const = 0;
  /// Compute fluid velocity
  virtual void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const = 0;
  /// Compute fluid momentum
  virtual void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const = 0;
  /// Compute stress tensor
  virtual void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                             T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const = 0;
  /// Compute fluid velocity and particle density
  virtual void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const = 0;
  /// Compute all momenta  up to second order
  virtual void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                                 T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const = 0;

  /// Set particle density
  virtual void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) = 0;
  /// Set fluid velocity
  virtual void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) = 0;
  /// Define fluid velocity and particle density
  virtual void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) = 0;
  /// Define all momenta up to second order
  virtual void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                                T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) = 0;

  /// Return iPop equilibrium for given first and second momenta
  virtual void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell,
                                  T rho,
                                  const T u[DESCRIPTOR::d],
                                  T fEq[DESCRIPTOR::q]) const = 0;
  /// Initialize to equilibrium distribution
  void iniEquilibrium(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) {
    T rhoShift { rho };
    T uShift[DESCRIPTOR::d] { };
    util::copyN(uShift, u, DESCRIPTOR::d);
    inverseShiftRhoU(cell, rhoShift, uShift);
    T fEq[DESCRIPTOR::q] { };
    computeEquilibrium(cell, rhoShift, uShift, fEq);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop];
    }
  };
  /// Initialize cell to equilibrium and non-equilibrum part
  void iniRegularized(Cell<T,DESCRIPTOR>& cell,
                      T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR>::n])  {
    iniEquilibrium(cell, rho, u);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
    }
  };
  /// Calculate population momenta s.t. the physical momenta are reproduced by the computeRhoU
  // This is relevant for correct initialization for some forcing/ sourcing schemes
  virtual void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const { };

};


namespace concepts {

/// Equilibrium functor
template <typename EQUILIBRIUM>
concept EquilibriumF = requires(EQUILIBRIUM eqF,
                                placeholder::Cell cell,
                                placeholder::Parameters parameters,
                                placeholder::Cell::value_t* fEq) {
  { eqF.compute(cell, parameters, fEq) } -> std::same_as<CellStatistic<placeholder::Cell::value_t>>;
};

/// Equilibrium element of dynamics::Tuple
template <typename EQUILIBRIUM>
concept EquilibriumElement = requires(EQUILIBRIUM eq) {
  // declares a list of parameters
  requires std::is_base_of_v<meta::list_base, typename EQUILIBRIUM::parameters>;
  // has human-readable name
  { EQUILIBRIUM::getName() } -> std::same_as<std::string>;
  // is concretizable into a equilibrium functor given a descriptor and a momenta tuple
  requires EquilibriumF<typename EQUILIBRIUM::template type<placeholder::Descriptor, momenta::BulkTuple>>;
};

namespace placeholder {

struct Equilibrium {
  using parameters = meta::list<>;

  static std::string getName() {
    return "Placeholder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      return {0, 0};
    };
  };
};

}

/// Collision operator
template <typename COLLISION>
concept CollisionO = requires(COLLISION op,
                              placeholder::Cell cell,
                              placeholder::Parameters parameters) {
  { op.apply(cell, parameters) } -> std::same_as<CellStatistic<placeholder::Cell::value_t>>;
};

/// Collision element of dynamics::Tuple
template <typename COLLISION>
concept CollisionElement = requires(COLLISION op) {
  // declares a list of parameters
  requires std::is_base_of_v<meta::list_base, typename COLLISION::parameters>;
  // has human-readable name
  { COLLISION::getName() } -> std::same_as<std::string>;
  // is concretizable into a collision operator given a descriptor, equilibrium and momenta tuple
  requires CollisionO<typename COLLISION::template type<
    placeholder::Descriptor,
    momenta::BulkTuple,
    placeholder::Equilibrium
  >>;
};

namespace placeholder {

struct Collision {
  using parameters = meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "Placeholder";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      return {0, 0};
    };
  };
};

}

/// Combinator rule for dynamics::Tuple
template <typename RULE>
concept CombinationRule = requires() {
  // has human-readable name
  { RULE::getName() } -> std::same_as<std::string>;

  requires EquilibriumF<typename RULE::template combined_equilibrium<
    placeholder::Descriptor,
    momenta::BulkTuple, // TODO replace
    placeholder::Equilibrium
  >>;
  requires CollisionO<typename RULE::template combined_collision<
    placeholder::Descriptor,
    momenta::BulkTuple, // TODO replace
    placeholder::Equilibrium,
    placeholder::Collision
  >>;
  requires std::is_base_of_v<meta::list_base, typename RULE::template combined_parameters<
    placeholder::Descriptor,
    momenta::BulkTuple, // TODO replace
    placeholder::Equilibrium,
    placeholder::Collision
  >>;
};

template <typename DYNAMICS>
concept IntrospectableDynamics = requires {
  typename DYNAMICS::template exchange_value_type<Expr>;
};

}

namespace dynamics {

/// Default combination rule used by dynamics::Tuple
/**
 * Use MOMENTA, EQUILIBRIUM and COLLISION elements as they are given
 **/
struct DefaultCombination {
  static std::string getName() {
    return "Default";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// DYNAMICS is not explicitly marked as unvectorizable
template <typename DYNAMICS, typename = void>
struct is_vectorizable : std::true_type { };

/// DYNAMICS is explicitly marked as unvectorizable
/**
 * Required as a workaround to enable non-SIMD fallback in ConcreteBlockCollisionO
 **/
template <typename DYNAMICS>
struct is_vectorizable<
  DYNAMICS,
  std::enable_if_t<!DYNAMICS::is_vectorizable>
> : std::false_type { };

/// Evaluates to true iff DYNAMICS is not explicitly marked as unvectorizable
template <typename DYNAMICS>
static constexpr bool is_vectorizable_v = is_vectorizable<DYNAMICS>::value;

/// Dynamics constructed as a tuple of momenta, equilibrium and collision
/**
 * Optionally also a combination rule thereof (e.g. a forcing scheme)
 **/
template <
  concepts::BaseType           T,
  concepts::LatticeDescriptor  DESCRIPTOR,
  typename MOMENTA,
  concepts::EquilibriumElement EQUILIBRIUM,
  concepts::CollisionElement   COLLISION,
  concepts::CombinationRule    COMBINATION_RULE = DefaultCombination
>
struct Tuple final : public Dynamics<T,DESCRIPTOR> {
  using MomentaF     = typename COMBINATION_RULE::template combined_momenta<DESCRIPTOR,MOMENTA>;
  using EquilibriumF = typename COMBINATION_RULE::template combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
  using CollisionO   = typename COMBINATION_RULE::template combined_collision<DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;

  using parameters = typename COMBINATION_RULE::template combined_parameters<DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;

  constexpr static bool is_vectorizable = is_vectorizable_v<CollisionO>;

  template <typename NEW_T>
  using exchange_value_type = Tuple<
    NEW_T, DESCRIPTOR,
    MOMENTA, EQUILIBRIUM, COLLISION,
    COMBINATION_RULE
  >;

  template <typename NEW_MOMENTA>
  using exchange_momenta = Tuple<
    T, DESCRIPTOR,
    NEW_MOMENTA, EQUILIBRIUM, COLLISION,
    COMBINATION_RULE
  >;

  template <concepts::CombinationRule NEW_RULE>
  using exchange_combination_rule = Tuple<
    T, DESCRIPTOR,
    MOMENTA, EQUILIBRIUM, COLLISION,
    NEW_RULE
  >;

  template <template<typename> typename WRAPPER>
  using wrap_collision = Tuple<
    T, DESCRIPTOR,
    MOMENTA, EQUILIBRIUM, WRAPPER<COLLISION>,
    COMBINATION_RULE
  >;

  std::type_index id() override {
    return typeid(Tuple);
  }

  std::string getName() const override {
    return "dynamics::Tuple<"
      + MomentaF().getName() + ","
      + EQUILIBRIUM::getName() + ","
      + COLLISION::getName() + ","
      + COMBINATION_RULE::getName() +
    ">";
  };

  /// Apply purely-local collision step to a generic CELL
  template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    // Copy parameters to enable changes in composed dynamics
    auto params = parameters.template copyAs<V>();
    return CollisionO().apply(cell, params);
  };

  /// Return true iff FIELD is a parameter
  template <typename FIELD>
  constexpr bool hasParameter() const {
    return parameters::template contains<FIELD>();
  };

  /// Interim workaround for accessing dynamics parameters in legacy post processors
  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<Tuple>>();
  }

  /// Initialize MOMENTA-specific data for cell
  void initialize(Cell<T,DESCRIPTOR>& cell) override {
    MomentaF().initialize(cell);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override {
    return MomentaF().computeRho(cell);
  };
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeU(cell, u);
  };
  void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const override {
    MomentaF().computeJ(cell, j);
  };
  void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                     T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeStress(cell, rho, u, pi);
  };
  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeRhoU(cell, rho, u);
  };
  void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                         T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeAllMomenta(cell, rho, u, pi);
  };

  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override {
    MomentaF().defineRho(cell, rho);
  };
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineU(cell, u);
  };
  void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineRhoU(cell, rho, u);
  };
  void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                        T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) override {
    MomentaF().defineAllMomenta(cell, rho, u, pi);
  };
  void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};

template <typename T, typename DESCRIPTOR, typename MOMENTA>
struct CustomCollision : public Dynamics<T,DESCRIPTOR> {
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  void initialize(Cell<T,DESCRIPTOR>& cell) override {
    MomentaF().initialize(cell);
  };

  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override {
    return MomentaF().computeRho(cell);
  };
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeU(cell, u);
  };
  void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const override {
    MomentaF().computeJ(cell, j);
  };
  void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                     T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeStress(cell, rho, u, pi);
  };
  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeRhoU(cell, rho, u);
  };
  void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                         T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeAllMomenta(cell, rho, u, pi);
  };

  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override {
    MomentaF().defineRho(cell, rho);
  };
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineU(cell, u);
  };
  void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineRhoU(cell, rho, u);
  };
  void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                        T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) override {
    MomentaF().defineAllMomenta(cell, rho, u, pi);
  };
  void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};

/// DYNAMICS is not explicitly marked as requiring parameters outside DYNAMICS::collide
template <typename DYNAMICS, typename = void>
struct has_parametrized_momenta : std::false_type { };

/// DYNAMICS is explicitly marked as requiring parameters outside DYNAMICS::collide
/**
 * Required to support ForcedPSMBGKdynamics
 **/
template <typename DYNAMICS>
struct has_parametrized_momenta<
  DYNAMICS,
  std::enable_if_t<DYNAMICS::has_parametrized_momenta>
> : std::true_type { };

/// Evaluates to true iff DYNAMICS doesn't require parameters outside DYNAMICS::collide
template <typename DYNAMICS>
static constexpr bool has_parametrized_momenta_v = has_parametrized_momenta<DYNAMICS>::value;


/// Evaluates combination rule OUTER on result of combination rule INNER (prototype)
template <typename OUTER, typename INNER>
class RuleComposition {
private:
  template <typename MOMENTA>
  struct inner_momenta {
    template <typename DESCRIPTOR>
    using type = typename INNER::template combined_momenta<DESCRIPTOR,MOMENTA>;

    template <template<typename> typename F>
    using wrap_momentum = inner_momenta<typename MOMENTA::template wrap_momentum<F>>;
  };

  template <typename EQUILIBRIUM>
  struct inner_equilibrium {
    template <typename DESCRIPTOR, typename MOMENTA>
    using type = typename INNER::template combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
  };

  template <typename COLLISION>
  struct inner_collision {
    using parameters = typename COLLISION::parameters;

    template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
    using type = typename INNER::template combined_collision<DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;
  };

public:
  static std::string getName() {
    return "RuleComposition<" + OUTER().getName() + "," + INNER().getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename OUTER::template combined_momenta<DESCRIPTOR, inner_momenta<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename OUTER::template combined_equilibrium<
    DESCRIPTOR,
    inner_momenta<MOMENTA>,
    inner_equilibrium<EQUILIBRIUM>
  >;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename OUTER::template combined_collision<
    DESCRIPTOR,
    inner_momenta<MOMENTA>,
    inner_equilibrium<EQUILIBRIUM>,
    inner_collision<COLLISION>
  >;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename INNER::template combined_parameters<
    DESCRIPTOR, MOMENTA, EQUILIBRIUM, COLLISION
  >::template decompose_into<
    OUTER::template combined_parameters<
      DESCRIPTOR,
      inner_momenta<MOMENTA>,
      inner_equilibrium<EQUILIBRIUM>,
      inner_collision<COLLISION>
    >::template include
   >;
};

}

}

#endif
