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

#ifndef BLOCK_POST_PROCESSOR_MAP_H
#define BLOCK_POST_PROCESSOR_MAP_H

#include <functional>
#include <typeindex>

#include "operator.h"
#include "postProcessing.h"
#include "cse/wrapper.h"
#include "introspection.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename OPERATOR>
struct ConcretizableBlockO {

using value_t = T;

using base_t = AbstractBlockO;

template <Platform PLATFORM>
using type = ConcreteBlockO<T,DESCRIPTOR,PLATFORM,OPERATOR,OPERATOR::scope>;

};

/// Factory for instances of a specific OPERATOR type
/**
 * Factory callable for ConcreteBlockO<OPERATOR> is constructed at
 * PostProcessorPromise construction time. Recipients accepting such
 * _promised post processors_ are not obligated to actually _realize_
 * this promise.
 *
 * Analogously to DynamicsPromise, this structure is needed to bridge the gap between
 * high level virtual interfaces and efficient platform-specific implementations with
 * full type knowledge.
 **/
template <typename T, typename DESCRIPTOR>
class PostProcessorPromise {
protected:
  const std::type_index _id;
  const int _priority;
  const OperatorScope _scope;
  const std::string _name;

  std::function<AbstractBlockO*(Platform)> _constructor;
  /// Returns the set of all accessed fields
  std::function<std::set<FieldTypePromise<T,DESCRIPTOR>>()> _accessedFields;

  const bool _isOptimizationAvailable;
  std::function<std::optional<std::size_t>()> _inspectArithmeticOperationCount;
  std::function<std::optional<bool>()> _inspectOptimizability;

public:
  template <typename OPERATOR>
  PostProcessorPromise(meta::id<OPERATOR> id = meta::id<OPERATOR>{}):
    _id(typeid(OPERATOR)),
    _priority(OPERATOR().getPriority()),
    _scope(OPERATOR::scope),
    _name{fields::name<OPERATOR>()},
    _constructor([](Platform platform) -> AbstractBlockO* {
      if constexpr (operators::is_cse_optimized<OPERATOR,DESCRIPTOR>::value) {
        return constructUsingConcretePlatform<ConcretizableBlockO<T,DESCRIPTOR,CSE_O<OPERATOR,DESCRIPTOR>>>(platform);
      } else {
        return constructUsingConcretePlatform<ConcretizableBlockO<T,DESCRIPTOR,OPERATOR>>(platform);
      }
    }),
    _accessedFields([]() -> std::set<FieldTypePromise<T,DESCRIPTOR>> {
      if constexpr (!std::is_same_v<T,Expr>) {
        return introspection::getFieldsAccessedByOperator<T,DESCRIPTOR,OPERATOR>();
      } else {
        throw std::domain_error("Can not introspect the introspection");
      }
    }),
    _isOptimizationAvailable{operators::is_cse_optimized<OPERATOR,DESCRIPTOR>::value},
    _inspectArithmeticOperationCount([]() -> std::optional<std::size_t> {
#ifdef FEATURE_INSPECT_POST_PROCESSORS
      if constexpr (!std::is_same_v<T,Expr>) {
        if constexpr (operators::is_cse_optimized<OPERATOR,DESCRIPTOR>::value) {
          return introspection::getArithmeticOperationCount<CSE_O<OPERATOR,DESCRIPTOR>>();
        } else {
          return introspection::getArithmeticOperationCount<OPERATOR,DESCRIPTOR>();
        }
      } else {
        throw std::domain_error("Can not introspect the introspection - we, sadly, are not in LISP land.");
      }
#endif // FEATURE_INSPECT_POST_PROCESSORS
      return std::nullopt;
    }),
    _inspectOptimizability([]() -> std::optional<bool> {
#ifdef FEATURE_INSPECT_POST_PROCESSORS
      if constexpr (!std::is_same_v<T,Expr>) {
        return introspection::isOptimizable<OPERATOR,DESCRIPTOR>();
      } else {
        return false;
      }
#else
      return std::nullopt;
#endif // FEATURE_INSPECT_POST_PROCESSORS
    })
  { }

  /// Returns type index of the promised OPERATOR
  std::type_index id() const {
    return _id;
  }

  std::string name() const {
    return _name;
  }

  int priority() const {
    return _priority;
  }

  OperatorScope scope() const {
    return _scope;
  }

  std::set<FieldTypePromise<T,DESCRIPTOR>> accessedFields() const {
    return _accessedFields();
  }

  bool hasOptimizedVersion() const {
    return _isOptimizationAvailable;
  }

  std::optional<std::size_t> getArithmeticOperationCount() const {
    return _inspectArithmeticOperationCount();
  }

  std::optional<bool> isOptimizable() const {
    return _inspectOptimizability();
  }

  bool operator<(const PostProcessorPromise<T,DESCRIPTOR>& rhs) const {
    return id() < rhs.id();
  }

  template <Platform PLATFORM>
  BlockO<T,DESCRIPTOR,PLATFORM>* realize() {
    return static_cast<BlockO<T,DESCRIPTOR,PLATFORM>*>(
      _constructor(PLATFORM));
  };

};

template <typename PP>
PostProcessorPromise(meta::id<PP>) -> PostProcessorPromise<typename PP::value_t,
                                                           typename PP::descriptor_t>;

/// Map of post processors of a single priority and stage
/**
 * Maintained in ConcreteBlockLattice
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class BlockPostProcessorMap {
private:
  /// Lattice reference to be passed to BlockO::setup and used for apply
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _lattice;

  /// Type indexed map of BlockO instances
  /**
   * Only add new BlockO instances to this map via resolve
   **/
  std::map<std::type_index,
           std::tuple<PostProcessorPromise<T,DESCRIPTOR>,
                      std::unique_ptr<BlockO<T,DESCRIPTOR,PLATFORM>>>> _map;

  /// Resolve promised post processor into concrete BlockO instance
  BlockO<T,DESCRIPTOR,PLATFORM>& resolve(PostProcessorPromise<T,DESCRIPTOR>&& promise)
  {
    auto iter = _map.find(promise.id());
    if (iter == _map.end()) {
      iter = _map.emplace(std::piecewise_construct,
                          std::forward_as_tuple(promise.id()),
                          std::forward_as_tuple(promise, promise.template realize<PLATFORM>())).first;
      auto& [_, op] = iter->second;
      op->setup(_lattice);
      if (_lattice.isIntrospectable()) {
        auto accessedFields = promise.accessedFields();
        for (FieldTypePromise<T,DESCRIPTOR> field : accessedFields) {
          field.ensureAvailabilityIn(_lattice);
        }
      }
      return *op;
    } else {
      auto& [_, op] = iter->second;
      return *op;
    }
  }

public:
  BlockPostProcessorMap(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>* lattice):
    _lattice{*lattice}
  { }

  /// Schedule post processor for application at iCell
  void add(std::size_t iCell, PostProcessorPromise<T,DESCRIPTOR>&& promise)
  {
    resolve(std::forward<decltype(promise)>(promise)).set(iCell, true);
  }

  /// Add post processor to map, do nothing if it already exists
  void add(PostProcessorPromise<T,DESCRIPTOR>&& promise)
  {
    resolve(std::forward<decltype(promise)>(promise));
  }

  /// Returns true if map contains post processor
  bool contains(PostProcessorPromise<T,DESCRIPTOR>&& promise) const
  {
    return _map.find(promise.id()) != _map.end();
  }

  /// Returns number of cells assigned to promised post processor
  /**
   * Doesn't allocate, intended for introspection
   **/
  std::size_t getWeight(const PostProcessorPromise<T,DESCRIPTOR>& promise) const {
    auto iter = _map.find(promise.id());
    if (iter != _map.end()) {
      const auto& [_, op] = std::get<1>(*iter);
      return op->weight();
    } else {
      return 0;
    }
  }

  /// Returns set of all post processors maintained this map
  std::set<PostProcessorPromise<T,DESCRIPTOR>> getAll() const
  {
    std::set<PostProcessorPromise<T,DESCRIPTOR>> operators;
    for (auto& [id, value] : _map) {
      const auto& [promise, _] = value;
      operators.emplace(promise);
    }
    return operators;
  }

  /// Apply all managed post processors to lattice
  /**
   * All post processors within a single BlockPostProcessorMap should be expected
   * to be executed in parallel as far as feasible. E.g. post processors on CPU
   * blocks are parallelized using OpenMP if enabled while per-cell operators
   * on GPU blocks are parallelized both on the level of kernels and by running
   * multiple kernels asynchronously in separate streams.
   *
   * In general, interdependent post processors should be separated into different
   * priorities (not doing so also "works" in many cases but leads to non-deterministic
   * data relationships).
   **/
  void apply()
  {
    for (auto& [_, value] : _map) {
      auto& [promise, postProcessor] = value;
      postProcessor->apply(_lattice);
    }

    #ifdef PLATFORM_GPU_CUDA
    gpu::cuda::device::synchronize();
    #endif
  }

};


}

#endif
