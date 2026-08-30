/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef OLB_CORE_INTROSPECTION_H
#define OLB_CORE_INTROSPECTION_H

#include "concepts.h"

#include <set>
#include <optional>
#include <regex>

#include "olbInit.h"

#include "cse/symbolGenerator.h"

namespace olb {

template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;

namespace introspection {

template <typename DYNAMICS>
[[gnu::noinline]]
std::optional<std::size_t> getArithmeticOperationCount() {
  if constexpr (concepts::IntrospectableDynamics<DYNAMICS>) {
    try {
      using DESCRIPTOR = typename DYNAMICS::descriptor_t;
      ConcreteBlockLattice<Expr,DESCRIPTOR,Platform::CPU_SISD> exprLattice(1, 0);
      exprLattice.setStatisticsEnabled(false);
      exprLattice.setIntrospectability(false);
      exprLattice.setDynamics(0, meta::id<typename DYNAMICS::template exchange_value_type<Expr>>{});
      Expr::reset();
      exprLattice.collide();
      return Expr::count();
    }
    catch (std::domain_error& error) {
      return std::nullopt;
    }
  } else {
    return std::nullopt;
  }
}

template <typename OPERATOR, typename DESCRIPTOR>
[[gnu::noinline]]
std::optional<std::size_t> getArithmeticOperationCount() {
  try {
    const int r = 3;
    Vector<int,DESCRIPTOR::d> extend(r*2+1);
    Vector<int,DESCRIPTOR::d> center(r);
    ConcreteBlockLattice<Expr,DESCRIPTOR,Platform::CPU_SISD> exprLattice(extend, 0);
    exprLattice.setStatisticsEnabled(false);
    exprLattice.setIntrospectability(false);
    if (OPERATOR::scope == OperatorScope::PerBlock) {
      exprLattice.addPostProcessor(typeid(stage::CSE),
                                   meta::id<OPERATOR>{});
    } else {
      exprLattice.addPostProcessor(typeid(stage::CSE),
                                   center,
                                   meta::id<OPERATOR>{});
    }
    Expr::reset();
    exprLattice.postProcess(typeid(stage::CSE));
    return Expr::count();
  }
  catch (std::domain_error& error) {
    return std::nullopt;
  }
}

template <typename DYNAMICS>
[[gnu::noinline]]
std::optional<std::size_t> getComplexity() {
  if constexpr (concepts::IntrospectableDynamics<DYNAMICS>) {
    try {
      using DESCRIPTOR = typename DYNAMICS::descriptor_t;
      ConcreteBlockLattice<Expr,DESCRIPTOR,Platform::CPU_SISD> exprLattice(1, 0);
      exprLattice.setStatisticsEnabled(false);
      exprLattice.setIntrospectability(false);
      exprLattice.setDynamics(0, meta::id<typename DYNAMICS::template exchange_value_type<Expr>>{});
      exprLattice.collide();
      std::size_t size = 0;
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        size += exprLattice.get(0)[iPop].size();
      }
      return size;
    }
    catch (std::domain_error& error) {
      return std::nullopt;
    }
  } else {
    return std::nullopt;
  }
}

template <concepts::IntrospectableDynamics DYNAMICS>
[[gnu::noinline]]
bool isOptimizable() {
  try {
    using DESCRIPTOR = typename DYNAMICS::descriptor_t;
    ConcreteBlockLattice<Expr,DESCRIPTOR,Platform::CPU_SISD> exprLattice(1, 0);
    exprLattice.setStatisticsEnabled(false);
    exprLattice.setIntrospectability(false);
    exprLattice.setDynamics(0, meta::id<typename DYNAMICS::template exchange_value_type<Expr>>{});
    exprLattice.collide();
  }
  catch (std::domain_error& error) {
    return false;
  }
  return true;
}

template <typename OPERATOR, typename DESCRIPTOR>
[[gnu::noinline]]
bool isOptimizable() {
  try {
    const int r = 3;
    Vector<int,DESCRIPTOR::d> extend(r*2+1);
    Vector<int,DESCRIPTOR::d> center(r);
    ConcreteBlockLattice<Expr,DESCRIPTOR,Platform::CPU_SISD> exprLattice(extend, 0);
    exprLattice.setStatisticsEnabled(false);
    exprLattice.setIntrospectability(false);
    if (OPERATOR::scope == OperatorScope::PerBlock) {
      exprLattice.addPostProcessor(typeid(stage::CSE),
                                   meta::id<OPERATOR>{});
    } else {
      exprLattice.addPostProcessor(typeid(stage::CSE),
                                   center,
                                   meta::id<OPERATOR>{});
    }
    exprLattice.postProcess(typeid(stage::CSE));
  }
  catch (std::domain_error& error) {
    return false;
  }
  return true;
}

template <typename DYNAMICS>
requires (!concepts::IntrospectableDynamics<DYNAMICS>)
bool isOptimizable() {
  return false;
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
[[gnu::noinline]]
std::set<FieldTypePromise<T,DESCRIPTOR>> getFieldsAccessedByDynamics() {
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> lattice(1, 0);
  lattice.setStatisticsEnabled(false);
  lattice.setIntrospectability(false);

  std::set<FieldTypePromise<T,DESCRIPTOR>> preCollide;
  lattice.getData().forEach([&](auto& anyFieldType) {
    preCollide.emplace(anyFieldType.getPromise());
  });

  lattice.defineDynamics(0, meta::id<DYNAMICS>());
  lattice.collide();

  std::set<FieldTypePromise<T,DESCRIPTOR>> postCollide;
  lattice.getData().forEach([&](auto& anyFieldType) {
    postCollide.emplace(anyFieldType.getPromise());
  });

  std::set<FieldTypePromise<T,DESCRIPTOR>> accessed;
  std::set_difference(postCollide.begin(), postCollide.end(),
                      preCollide.begin(), preCollide.end(),
                      std::inserter(accessed, accessed.end()));
  return accessed;
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
[[gnu::noinline]]
std::set<FieldTypePromise<T,DESCRIPTOR>> getFieldsAccessedByCollision() {
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> lattice(1, 0);
  lattice.setStatisticsEnabled(false);
  lattice.setIntrospectability(false);

  lattice.defineDynamics(0, meta::id<DYNAMICS>());

  std::set<FieldTypePromise<T,DESCRIPTOR>> preCollide;
  lattice.getData().forEach([&](auto& anyFieldType) {
    preCollide.emplace(anyFieldType.getPromise());
  });

  lattice.collide();

  std::set<FieldTypePromise<T,DESCRIPTOR>> postCollide;
  lattice.getData().forEach([&](auto& anyFieldType) {
    postCollide.emplace(anyFieldType.getPromise());
  });

  std::set<FieldTypePromise<T,DESCRIPTOR>> accessed;
  std::set_difference(postCollide.begin(), postCollide.end(),
                      preCollide.begin(), preCollide.end(),
                      std::inserter(accessed, accessed.end()));
  return accessed;
}

template <typename T, typename DESCRIPTOR, typename OPERATOR>
[[gnu::noinline]]
std::set<FieldTypePromise<T,DESCRIPTOR>> getFieldsAccessedByOperator() {
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> lattice(1, 3);
  lattice.setStatisticsEnabled(false);
  lattice.setIntrospectability(false);

  std::set<FieldTypePromise<T,DESCRIPTOR>> prePostProcess;
  lattice.getData().forEach([&](auto& anyFieldType) {
    prePostProcess.emplace(anyFieldType.getPromise());
  });

  if (OPERATOR::scope == OperatorScope::PerBlock) {
    lattice.addPostProcessor(typeid(stage::CSE),
                             meta::id<OPERATOR>{});
  } else {
    lattice.addPostProcessor(typeid(stage::CSE),
                             0,
                             meta::id<OPERATOR>{});
  }
  lattice.postProcess(typeid(stage::CSE));

  std::set<FieldTypePromise<T,DESCRIPTOR>> postPostProcess;
  lattice.getData().forEach([&](auto& anyFieldType) {
    postPostProcess.emplace(anyFieldType.getPromise());
  });

  std::set<FieldTypePromise<T,DESCRIPTOR>> accessed;
  std::set_difference(postPostProcess.begin(), postPostProcess.end(),
                      prePostProcess.begin(), prePostProcess.end(),
                      std::inserter(accessed, accessed.end()));
  return accessed;
}

/// Returns estimate of number of scalar reads and writes for DYNAMICS::collide
template <typename DYNAMICS>
[[gnu::noinline]]
std::optional<std::size_t> getMemoryBandwidthOfCollision() {
  using T = Expr;
  using DESCRIPTOR = DYNAMICS::descriptor_t;

  try {
    Vector<int,DESCRIPTOR::d> size(1);
    ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> exprLattice(size, 0);
    exprLattice.setIntrospectability(false);
    exprLattice.setStatisticsEnabled(false);

    exprLattice.setDynamics(0, meta::id<typename DYNAMICS::template exchange_value_type<T>>{});

    std::size_t bandwidth = 0;

    Vector<int,DESCRIPTOR::d> pos(0);
    auto cell = exprLattice.get(pos);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = "cell[" + std::to_string(iPop) + "]";
    }
    // Assume that all populations are read
    bandwidth += DESCRIPTOR::q;

    auto accessedFields = getFieldsAccessedByCollision<typename DYNAMICS::value_t,DESCRIPTOR,DYNAMICS>();

    // Initialize expressions for the accessed fields
    for (auto field : accessedFields) {
      if (auto dim = field.dimension()) {
        // Assume that all touched fields are read
        bandwidth += *dim;
        std::vector<T> initExpr;
        for (unsigned iD=0; iD < dim; ++iD) {
          initExpr.emplace_back("cell[" + cse::remove_array_from_name(field.name()) + "][" + std::to_string(iD) + "]");
        }
        field.setPlaceholderExpression(cell, initExpr);
      }
    }

    // Execute one colliding step to obtain post-collision expressions
    exprLattice.collide();

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      if (!cell[iPop].isSymbol("cell[" + std::to_string(iPop) + "]")) {
        // Count number of population updates
        bandwidth += 1;
      }
    }

    // Obtain fields which were updated and overwritten by the operator
    for (auto field : accessedFields) {
      if (auto dim = field.dimension()) {
        auto postExpr = field.getPlaceholderExpression(cell);
        for (unsigned iD=0; iD < dim; ++iD) {
          if (!postExpr[iD].isSymbol("cell["
                                     + cse::remove_array_from_name(field.name())
                                     + "][" + std::to_string(iD) +  "]")) {
            // Count number of field component updates
            bandwidth += 1;
          }
        }
      }
    }

    return bandwidth;
  }
  catch (std::domain_error& error) {
    return std::nullopt;
  }
}

/// Return reduced and reasonable newline-separated version of given full dynamics name
[[gnu::noinline]]
std::string getSimplifiedDynamicsName(std::string name) {
  std::regex removeBaseTypeAndDescriptorRe("(float|double|Expr),\\s?descriptors::D[23]Q[0-9]+<([^<>]*|<([^<>]*|<[^<>]*>)*>)*>,?\\s?");
  std::regex lineBreakRe("((dynamics|collision|forcing|equilibria)::[a-zA-Z]+|momenta::Tuple<)");

  std::string stage0;
  std::string stage1;
  std::regex_replace(std::back_inserter(stage0),
                     name.begin(), name.end(), removeBaseTypeAndDescriptorRe,
                     "");
  std::regex_replace(std::back_inserter(stage1),
                     stage0.begin(), stage0.end(), lineBreakRe,
                     "\n$&");
  return stage1;
}


}

}

#endif
