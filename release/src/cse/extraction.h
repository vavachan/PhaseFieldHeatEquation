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

#ifndef CSE_EXTRACTION_H
#define CSE_EXTRACTION_H

#include <type_traits>
#include "core/expr.h"
#include "cse/symbolGenerator.h"
#include "cse/exprDynamics.h"

namespace olb {

namespace cse {

/// Extract expression tree of the collision operator for DYNAMICS
template <typename DYNAMICS>
std::stringstream extractExpressionTree() {
  std::stringstream tree;
  using T = Expr;
  using DESCRIPTOR = DYNAMICS::descriptor_t;
  using PARAMETERS = DYNAMICS::parameters;

  // Instatiate a 1x1 BlockLattice with Expr type
  Vector<int,DESCRIPTOR::d> size(1);
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> exprLattice(size, 0);
  // This is required as invalid operators for Expr are used in CellStatistics
  exprLattice.setStatisticsEnabled(false);

  // Set dynamics to the cells
  exprLattice.setDynamics(0, meta::id<dynamics::StoreStatisticInField<DYNAMICS>>{});

  // Initialize expressions for the populations
  Vector<int,DESCRIPTOR::d> pos(0);
  auto cell = exprLattice.get(pos);

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = T(symbol("cell[" + std::to_string(iPop) + "]"));
  }

  std::vector<T> params;

  // Initialize expressions for the accessed fields
  for (auto field : introspection::getFieldsAccessedByDynamics<T,DESCRIPTOR,DYNAMICS>()) {
    if (auto dim = field.dimension()) {
      std::vector<T> initExpr;
      for (unsigned iD=0; iD < dim; ++iD) {
        initExpr.emplace_back(T(symbol("cell.template getFieldComponent<" + remove_array_from_name(field.name()) + ">(" + std::to_string(iD) + ")")));
        params.emplace_back(T(symbol("cell.template getFieldComponent<" + remove_array_from_name(field.name()) + ">(" + std::to_string(iD) + ")")));
      }
      field.setPlaceholderExpression(cell, initExpr);
    }
  }

  // Retrieve parameters from Dynamics and initialize unique expressions
  PARAMETERS::for_each([&](auto field) {
    using FIELD = typename decltype(field)::type;
    FieldD<T,DESCRIPTOR,FIELD> f;
    if (f.d == 1) {
      T expr = T(symbol("parameters.template get<"+std::string(fields::name<FIELD>())+">()"));
      exprLattice.template setParameter<FIELD>(expr);
      params.push_back(expr);
    } else {
      for (unsigned d=0; d < f.d; ++d) {
        f[d] = T(symbol("parameters.template get<"+std::string(fields::name<FIELD>())+">()["+std::to_string(d)+"]"));
        params.push_back(f[d]);
      }
      exprLattice.template setParameter<FIELD>(f);
    }
  });

  // Execute one colliding step to obtain post-collision expressions
  exprLattice.collide();

  // Obtain fields which were updated and overwritten by the operator
  std::vector<std::pair<std::string,std::vector<T>>> updatedFields;
  for (auto field : introspection::getFieldsAccessedByDynamics<T,DESCRIPTOR,DYNAMICS>()) {
    if (auto dim = field.dimension()) {
      auto postExpr = field.getPlaceholderExpression(cell);
      bool updated = false;
      // check if fields are updated
      for (unsigned iD=0; iD < dim; ++iD) {
        updated = !(postExpr[iD].describe() == "Symbol(\"cell.template getFieldComponent<"
                                                + remove_array_from_name(field.name())
                                                + ">(" + std::to_string(iD) +  ")\")");
        if (updated) {
          break;
        }
      }
      // for any updated field, add the field type and expressions to the list
      if (updated) {
        updatedFields.emplace_back(std::make_pair(remove_array_from_name(field.name()), std::move(postExpr)));
      }
    }
  }

  tree << "# Dynamics info: \n";
  tree << "collisionO = \"" << fields::name<DYNAMICS>() << "\"\n";

  // Retrieve expressions for the post-collision populations
  tree << "# Expressions of the post-collision populations\n";
  tree << "cell_assignments = [\n";
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    tree << "  " << assignment(symbol("cell["+std::to_string(iPop)+"]"),cell[iPop].describe()) << ",\n";
  }
  tree << "]\n";

  // Retrieve parameter expressions
  tree << "# Expressions of the used parameters\n";
  tree << "params_symbols = [\n";
  for (std::size_t i=0; i < params.size(); ++i) {
    tree << "  " << params[i].describe() << ",\n";
  }
  tree << "]\n";

  // Retrieve expressions for the cell statistics
  tree << "# Expressions of the return values\n";
  tree << "return_assignments = [\n";
  tree << "  " << assignment(symbol("rho"),cell.template getFieldComponent<descriptors::STATISTIC>(0).describe()) << ",\n";
  tree << "  " << assignment(symbol("uSqr"),cell.template getFieldComponent<descriptors::STATISTIC>(1).describe()) << ",\n";
  tree << "]\n";

  // Retrieve updated field expressions
  tree << "# Expressions of the updated fields\n";
  tree << "fields_assignments = [\n";
  for (std::size_t i=0; i < updatedFields.size(); ++i) {
    for (std::size_t iD=0; iD < updatedFields[i].second.size(); ++iD) {
      tree << "  " << assignment(symbol("cell.template getFieldPointer<" + updatedFields[i].first + ">()[" + std::to_string(iD) + "]"),
                                 updatedFields[i].second[iD].describe()) << ",\n";
    }
  }
  tree << "]\n";
  return tree;
}

/// Extract expression tree of the post-processor OPERATOR
template<typename OPERATOR, typename DESCRIPTOR>
std::stringstream extractExpressionTree() {
  std::stringstream tree;
  using T = Expr;
  const int neighbor_radius = 3; // can be adjusted if necessary

  // Instatiate a 5x5 BlockLattice with Expr type
  Vector<int,DESCRIPTOR::d> size(neighbor_radius*2+1);
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> exprLattice(size, 0);
  // This is required as invalid operators for Expr are used in CellStatistics
  exprLattice.setStatisticsEnabled(false);

  // Set post-processor to the center cell
  Vector<int,DESCRIPTOR::d> center(neighbor_radius);
  exprLattice.addPostProcessor(typeid(stage::CSE), center, meta::id<OPERATOR>{});

  // First initialize neighbor-related prefix expressions in fields
  auto cell = exprLattice.get(center);
  auto variable_scalar = std::make_shared<SymbolGenerator>("v_S");
  auto variable_vector = std::make_shared<SymbolGenerator>("v_V");
  auto variable_tensor = std::make_shared<SymbolGenerator>("v_T");
  auto variable_population = std::make_shared<SymbolGenerator>("v_P");
  auto calls = std::make_shared<std::vector<Expr>>();

  exprLattice.forCoreSpatialLocations([&](Vector<int,DESCRIPTOR::d> loc) {
    // Set dummy dynamics to extract method calls based on dynamics
    auto cellId = exprLattice.getCellId(loc);
    exprLattice.setDynamics(cellId, DynamicsPromise(meta::id<ExprDynamics<DESCRIPTOR>>{}));

    // Provide access to symbol generators and function call list
    auto loc_cell = exprLattice.get(loc);
    loc_cell.template getFieldPointer<SCALAR_SYMBOL_GENERATOR>()[0] = variable_scalar.get();
    loc_cell.template getFieldPointer<VECTOR_SYMBOL_GENERATOR>()[0] = variable_vector.get();
    loc_cell.template getFieldPointer<TENSOR_SYMBOL_GENERATOR>()[0] = variable_tensor.get();
    loc_cell.template getFieldPointer<POPULATION_SYMBOL_GENERATOR>()[0] = variable_population.get();
    loc_cell.template getFieldPointer<FUNCTION_CALLS>()[0] = calls.get();

    // Fill function call prefix for calling functions from neighbors
    std::string arg("cell");
    auto dist = loc - center;
    if (norm(dist) > 0) {
      arg +=".neighbor({";
      for (std::size_t d=0; d < dist.size(); ++d) {
        if (d < dist.size() - 1) {
          arg += std::to_string(dist[d]) + ",";
        } else {
          arg += std::to_string(dist[d]);
        }
      }
      arg += "})";
    }
    cell.neighbor(dist).template setField<cse::NEIGHBOR_PREFIX>(T(arg));
  });

  // Initialize expressions for all cells
  exprLattice.forCoreSpatialLocations([&](Vector<int,DESCRIPTOR::d> loc) {
    auto dist = loc - center;
    std::string prefix = cell.neighbor(dist).template getField<cse::NEIGHBOR_PREFIX>().describe();
    auto loc_cell = cell.neighbor(dist);

    // Initialize expressions for the populations
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell.neighbor(dist)[iPop] = T(symbol(prefix + "[" + std::to_string(iPop) + "]"));
    }

    // Initialize expressions for the accessed fields
    for (auto field : introspection::getFieldsAccessedByOperator<T,DESCRIPTOR,OPERATOR>()) {
      if (auto dim = field.dimension()) {
        std::vector<T> initExpr;
        for (unsigned iD=0; iD < dim; ++iD) {
          initExpr.emplace_back(T(symbol( prefix + ".template getFieldComponent<"
                                                 + remove_array_from_name(field.name())
                                                 + ">(" + std::to_string(iD) + ")")));
        }
        field.setPlaceholderExpression(loc_cell, initExpr);
      }
    }
  });

  // Retrieve parameters from Operator and initialize unique expressions
  std::vector<T> params;
  if constexpr (OPERATOR::scope == OperatorScope::PerCellWithParameters) {
    using PARAMETERS = OPERATOR::parameters;
    PARAMETERS::for_each([&](auto field) {
      using FIELD = typename decltype(field)::type;
      FieldD<T,DESCRIPTOR,FIELD> f;
      if (f.d == 1) {
        T expr = T(symbol("parameters.template get<"+std::string(fields::name<FIELD>())+">()"));
        exprLattice.template setParameter<FIELD>(expr);
        params.push_back(expr);
      } else {
        for (unsigned d=0; d < f.d; ++d) {
          f[d] = T(symbol("parameters.template get<"+std::string(fields::name<FIELD>())+">()["+std::to_string(d)+"]"));
          params.push_back(f[d]);
        }
        exprLattice.template setParameter<FIELD>(f);
      }
    });
  }
  // Set placeholder for omega maunally as this is often done in the app
  T exprOmega = T(symbol("cell.getDynamics().getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN())"));
  exprLattice.template setParameter<descriptors::OMEGA>(exprOmega);
  params.push_back(exprOmega);

  // Execute one colliding step to obtain post-collision expressions
  exprLattice.postProcess(typeid(stage::CSE));

  // Obtain fields which were updated and overwritten by the operator
  // Currently, this is only checked for the cell on which the operator is applied
  std::vector<std::pair<std::string,std::vector<T>>> updatedFields;
  for (auto field : introspection::getFieldsAccessedByOperator<T,DESCRIPTOR,OPERATOR>()) {
    if (auto dim = field.dimension()) {
      auto postExpr = field.getPlaceholderExpression(cell);
      bool updated = false;
      // check if fields are updated
      for (unsigned iD=0; iD < dim; ++iD) {
        updated = !(postExpr[iD].describe() == "Symbol(\"cell.template getFieldComponent<"
                                                + remove_array_from_name(field.name())
                                                + ">(" + std::to_string(iD) +  ")\")");
        if (updated) {
          break;
        }
      }
      // for any updated field, add the field type and expressions to the list
      if (updated) {
        updatedFields.emplace_back(std::make_pair(remove_array_from_name(field.name()), std::move(postExpr)));
      }
    }
  }

  tree << "# Operator info: \n";
  tree << "operatorO = \"" << fields::name<OPERATOR>() << "\"\n";
  tree << "# Desriptor info: \n";
  tree << "descriptor = \"" << fields::name<DESCRIPTOR>() << "\"\n";
  tree << "# isOperatorWithParameter: \n";
  tree << "isOperatorWithParameter = " << ((OPERATOR::scope == OperatorScope::PerCellWithParameters) ?
                                               "true" : "false") << "\n";

  // Retrieve expressions for the post-collision populations
  tree << "# Expressions of the post-collision populations\n";
  tree << "cell_assignments = [\n";
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    tree << "  " << assignment(symbol("cell["+std::to_string(iPop)+"]"),cell[iPop].describe()) << ",\n";
  }
  tree << "]\n";

  // Retrieve parameter expressions
  tree << "# Expressions of the used parameters\n";
  tree << "params_symbols = [\n";
  for (std::size_t i=0; i < params.size(); ++i) {
    tree << "  " << params[i].describe() << ",\n";
  }
  tree << "]\n";

  // Retrieve updated field expressions
  tree << "# Expressions of the updated fields\n";
  tree << "fields_assignments = [\n";
  for (std::size_t i=0; i < updatedFields.size(); ++i) {
    for (std::size_t iD=0; iD < updatedFields[i].second.size(); ++iD) {
      tree << "  " << assignment(symbol("cell.template getFieldPointer<" + updatedFields[i].first + ">()[" + std::to_string(iD) + "]"),
                                 updatedFields[i].second[iD].describe()) << ",\n";
    }
  }
  tree << "]\n";

  // Print dynamics based cell methods if called
  tree << "# Symbols of the called dynamics based functions\n";
  tree << "function_calls = [\n";
  for (auto call : *calls) {
    tree << "  " << call.describe() << ",\n";
  }
  tree << "]\n";
  return tree;
}

}

}

#endif
