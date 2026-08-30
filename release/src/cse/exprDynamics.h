/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Shota Ito
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

#ifndef CSE_EXPR_DYNAMICS_H
#define CSE_EXPR_DYNAMICS_H

#include "cse/symbolGenerator.h"

namespace olb {

namespace cse {

struct ExprDensity {
  template <typename TYPE, typename CELL, typename R>
  void compute(CELL& cell, R& rho) any_platform
  {
    using V = typename CELL::value_t;
    auto scalar = cell.template getFieldPointer<cse::SCALAR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = scalar->next();
    rho = V(cse::symbol(name));

    std::string body = prefix + ".computeRho";
    std::string call = cse::assignment(cse::symbol(name),
                                       cse::function(cse::symbol(body),cse::symbol(body+name)));

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(call);
  }

  template <typename TYPE, typename CELL, typename R>
  void define(CELL& cell, const R& rho) any_platform
  {
    auto scalar = cell.template getFieldPointer<cse::SCALAR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = scalar->next();
    std::string decl = cse::assignment(cse::symbol(name),
                                       cse::function(cse::symbol("DECLARE_WITH_ARGUMENT:V "+name),rho.describe()));

    std::string body = prefix + ".defineRho";
    std::string call = cse::assignment(cse::symbol(body),
                                       cse::function(cse::symbol(body),cse::symbol(name)));

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(decl+","+call);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename RHO>
  void inverseShift(CELL& cell, RHO& rho) any_platform {};

  static std::string getName(){
    return "ExprDensity";
  }
};

struct ExprMomentum {
  // compute the momentum
  template <typename TYPE, typename CELL, typename J, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, J& j) any_platform
  {
    auto vector = cell.template getFieldPointer<cse::VECTOR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = vector->next();
    std::string body = prefix + ".computeJ";
    std::string call = cse::assignment(cse::symbol(name),
                                       cse::function(cse::symbol(body),cse::symbol(body+name)));
    std::string elem;
    for (auto D=0; D < DESCRIPTOR::d; ++D) {
      elem += ",";
      elem += cse::assignment(cse::symbol(name+"["+std::to_string(D)+"]"),
                              cse::function(cse::symbol("EXPOSE_ARRAY_ELEMENT:"+name+"["+std::to_string(D)+"]"),cse::symbol(name)));
      j[D] = V(cse::symbol(name+"["+std::to_string(D)+"]"));
    }

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(call+elem);
  }

  // compute the velocity
  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void computeU(CELL& cell, U& u) any_platform
  {
    auto vector = cell.template getFieldPointer<cse::VECTOR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = vector->next();
    std::string body = prefix + ".computeU";
    std::string call = cse::assignment(cse::symbol(name),
                                       cse::function(cse::symbol(body),cse::symbol(body+name)));
    std::string elem;
    for (auto D=0; D < DESCRIPTOR::d; ++D) {
      elem += ",";
      elem += cse::assignment(cse::symbol(name+"["+std::to_string(D)+"]"),
                              cse::function(cse::symbol("EXPOSE_ARRAY_ELEMENT:"+name+"["+std::to_string(D)+"]"),cse::symbol(name)));
      u[D] = V(cse::symbol(name+"["+std::to_string(D)+"]"));
    }

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(call+elem);
  }

  // define the velocity
  template <typename TYPE, typename CELL, typename U>
  void define(CELL& cell, const U& u) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto vector = cell.template getFieldPointer<cse::VECTOR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = vector->next();
    std::string decl = cse::assignment(cse::symbol(name),
                                       cse::symbol("DECLARE:V "+name+ " [DESCRIPTOR::d]"));
    std::string elem;
    for (auto D=0; D < DESCRIPTOR::d; ++D) {
      elem += ",";
      elem += cse::assignment(cse::symbol(name+"["+std::to_string(D)+"]"),
                              cse::function(cse::symbol("FILL_ARRAY_ELEMENT:"+name+"["+std::to_string(D)+"]"),cse::symbol(name)+","+u[D].describe()));
    }
    elem += ",";
    std::string body = prefix + ".defineU";
    std::string call = cse::assignment(cse::symbol(body),
                                       cse::function(cse::symbol(body),cse::symbol(name)));

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(decl+elem+call);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void initialize(CELL& cell) any_platform {};

  template <typename TYPE, typename CELL, typename U>
  void inverseShift(CELL& cell, U& u) any_platform {};

  static std::string getName() {
    return "ExprMomentum";
  }
};

struct ExprStress {
  // Currently, it is not possible during extraction to detect either computeStress or computeAllMomenta is calling
  // this method, which is why computeStress is replaced by computeAllMomenta during extraction.
  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void compute(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    auto tensor = cell.template getFieldPointer<cse::TENSOR_SYMBOL_GENERATOR>()[0];
    std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

    std::string name = tensor->next();
    std::string rho_name = cse::strip_symbol(rho.describe());
    std::string u_name = cse::extract_c_array_name(cse::strip_symbol(u[0].describe()));
    std::string collected_arguments_as_symbols = cse::symbol(rho_name) + "," +
                                                 cse::symbol(u_name);

    std::string body = prefix + ".computeStress";
    std::string call = cse::assignment(cse::symbol(name),
                                       cse::function(cse::symbol(body),collected_arguments_as_symbols));
    std::string elem;
    for (auto D=0; D < util::TensorVal<DESCRIPTOR>::n; ++D) {
      elem += ",";
      elem += cse::assignment(cse::symbol(name+"["+std::to_string(D)+"]"),
                              cse::function(cse::symbol("EXPOSE_ARRAY_ELEMENT:"+name+"["+std::to_string(D)+"]"),cse::symbol(name)));
      pi[D] = V(cse::symbol(name+"["+std::to_string(D)+"]"));
    }

    cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(call+elem);
  }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void define(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    // Not defined in the cell interface yet, so should not be called independently
    throw std::bad_function_call();
  }

  static std::string getName() {
    return "ExprStress";
  }
};

struct ExprEquilibrium {
  using parameters = meta::list<>;

  static std::string getName() {
    return "ExprEquilibrium";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      auto scalar = cell.template getFieldPointer<cse::SCALAR_SYMBOL_GENERATOR>()[0];
      auto vector = cell.template getFieldPointer<cse::VECTOR_SYMBOL_GENERATOR>()[0];
      auto population = cell.template getFieldPointer<cse::POPULATION_SYMBOL_GENERATOR>()[0];
      std::string prefix = cell.template getField<cse::NEIGHBOR_PREFIX>().describe();

      std::string name = population->next();
      std::string collected_arguments_as_symbols;
      std::string pre_decl;

      collected_arguments_as_symbols += cse::symbol(prefix) + ",";
      if (scalar->contains_symbol(cse::strip_symbol(rho.describe()))) {
        collected_arguments_as_symbols += cse::symbol(*(scalar->find_symbol(cse::strip_symbol(rho.describe())))) + ",";
      } else {
        std::string new_name = scalar->next();
        pre_decl += cse::assignment(cse::symbol(new_name),
                                    cse::function(cse::symbol("DECLARE_WITH_ARGUMENT:V "+new_name),rho.describe()));
        pre_decl += ",";
        collected_arguments_as_symbols += cse::symbol(new_name) + ",";
      }
      if (vector->contains_symbol(cse::extract_c_array_name(cse::strip_symbol(u[0].describe())))) {
        std::string array_name = cse::extract_c_array_name(cse::strip_symbol(u[0].describe()));
        collected_arguments_as_symbols += cse::symbol(*(vector->find_symbol(array_name))) + ",";
      } else {
        std::string new_name = vector->next();
        pre_decl += cse::assignment(cse::symbol(new_name),
                                    cse::symbol("DECLARE:V "+new_name+ " [DESCRIPTOR::d]"));
        pre_decl += ",";
        for (auto D=0; D < DESCRIPTOR::d; ++D) {
          pre_decl += cse::assignment(cse::symbol(name+"["+std::to_string(D)+"]"),
                                  cse::function(cse::symbol("FILL_ARRAY_ELEMENT:"+name+"["+std::to_string(D)+"]"),cse::symbol(name)+","+u[D].describe()));
          pre_decl += ",";
        }
        collected_arguments_as_symbols += cse::symbol(new_name);
      }

      std::string body = prefix + ".getDynamics().computeEquilibrium";
      std::string call = cse::assignment(cse::symbol(name),
                                         cse::function(cse::symbol(body),collected_arguments_as_symbols));

      std::string elem;
      for (auto iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        elem += ",";
        elem += cse::assignment(cse::symbol(name+"["+std::to_string(iPop)+"]"),
                                cse::function(cse::symbol("EXPOSE_ARRAY_ELEMENT:"+name+"["+std::to_string(iPop)+"]"),cse::symbol(name)));
        fEq[iPop] = V(cse::symbol(name + "[" + std::to_string(iPop) + "]"));
      }

      cell.template getFieldPointer<cse::FUNCTION_CALLS>()[0]->emplace_back(pre_decl+call+elem);
      return {V(0), V(0)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      // This overload is currently not used by operators
      throw std::bad_function_call();
      return {0, 0};
    };
  };
};

/* Dummy dynamics struct for overloading methods during expression
 * extraction.
 *
 * Currently, the extraction is limited to one call of each method
 * per cell. Calls on neighbors is working as intended.
 */
template <typename DESCRIPTOR>
using ExprDynamics = dynamics::Tuple<
  Expr, DESCRIPTOR,
  momenta::Tuple<
    cse::ExprDensity,
    cse::ExprMomentum,
    cse::ExprStress,
    momenta::DefineSeparately
  >,
  cse::ExprEquilibrium,
  collision::ParameterFromCell<descriptors::OMEGA, collision::None>
>;

}

}

#endif
