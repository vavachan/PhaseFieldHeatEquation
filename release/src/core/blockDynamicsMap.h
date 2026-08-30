/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef BLOCK_DYNAMICS_MAP_H
#define BLOCK_DYNAMICS_MAP_H

#include <functional>
#include <typeindex>
#include <algorithm>
#include <optional>
#include <unordered_set>

#include "dynamics/dynamics.h"
#include "cse/wrapper.h"

#include "platform/cpu/cell.h"

#include "introspection.h"

namespace olb {

template<typename T, typename DESCRIPTOR> struct AbstractCollisionO;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> struct BlockCollisionO;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;
template<typename T,                      Platform PLATFORM> class ConcreteBlockMask;
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS> class ConcreteBlockCollisionO;

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
struct ConcretizableBlockCollisionO {

using value_t = T;

using base_t = AbstractCollisionO<T,DESCRIPTOR>;

template <Platform PLATFORM>
using type = ConcreteBlockCollisionO<T,DESCRIPTOR,PLATFORM,DYNAMICS>;

};

/// Factory for instances of a specific Dynamics type
/**
 * Factory callables for DYNAMICS and ConcreteBlockCollisionO<DYNAMICS> are constructed
 * at DynamicsPromise construction time. Recipients accepting such _promised dynamics_
 * are not obligated to actually _realize_ this promise.
 *
 * This structure is needed to bridge the gap between high level virtual interfaces and
 * efficient platform-specific implementations with full type knowledge.
 * i.e. DynamicsPromise can carry full dynamics type information through the virtual
 * layer surrounding concrete platform-specialized block lattices (needed as virtual
 * template methods are not supported by C++).
 **/
template <typename T, typename DESCRIPTOR>
class DynamicsPromise {
protected:
  const std::type_index _id;

  std::function<Dynamics<T,DESCRIPTOR>*()>                   _makeDynamics;
  std::function<AbstractCollisionO<T,DESCRIPTOR>*(Platform)> _makeOperator;

  const std::string _name;

  /// Returns the set of all accessed fields
  std::function<std::set<FieldTypePromise<T,DESCRIPTOR>>()> _accessedFields;

  const bool _isOptimizationAvailable;
  std::function<std::optional<std::size_t>()> _inspectArithmeticOperationCount;
  std::function<std::optional<std::size_t>()> _inspectMemoryBandwidth;
  std::function<std::optional<bool>()> _inspectOptimizability;
  std::function<std::optional<std::size_t>()> _inspectComplexity;

public:
  template <typename DYNAMICS>
  DynamicsPromise(meta::id<DYNAMICS> id = meta::id<DYNAMICS>{}):
    _id(typeid(DYNAMICS)),
    _makeDynamics([]() -> Dynamics<T,DESCRIPTOR>* {
      return new DYNAMICS();
    }),
    _makeOperator([](Platform platform) -> AbstractCollisionO<T,DESCRIPTOR>* {
      if constexpr (concepts::IntrospectableDynamics<DYNAMICS>) {
        if constexpr (dynamics::is_cse_optimized<DYNAMICS>::value) {
         return constructUsingConcretePlatform<ConcretizableBlockCollisionO<T,DESCRIPTOR,CSE<DYNAMICS>>>(platform);
        } else {
          return constructUsingConcretePlatform<ConcretizableBlockCollisionO<T,DESCRIPTOR,DYNAMICS>>(platform);
        }
      } else {
        return constructUsingConcretePlatform<ConcretizableBlockCollisionO<T,DESCRIPTOR,DYNAMICS>>(platform);
      }
    }),
    _name{fields::name<DYNAMICS>()},
    _accessedFields([]() -> std::set<FieldTypePromise<T,DESCRIPTOR>> {
      if constexpr (!std::is_same_v<T,Expr>) {
        if constexpr (std::is_same_v<NoDynamics<T,DESCRIPTOR>, DYNAMICS>) {
          return {};
        } else {
          return introspection::getFieldsAccessedByDynamics<T,DESCRIPTOR,DYNAMICS>();
        }
      } else {
        throw std::domain_error("Can not introspect the introspection");
      }
    }),
    _isOptimizationAvailable{dynamics::is_cse_optimized<DYNAMICS>::value},
    _inspectArithmeticOperationCount([]() -> std::optional<std::size_t> {
#ifdef FEATURE_INSPECT_DYNAMICS
      if constexpr (!std::is_same_v<T,Expr>) {
        if constexpr (dynamics::is_cse_optimized<DYNAMICS>::value) {
          return introspection::getArithmeticOperationCount<CSE<DYNAMICS>>();
        } else {
          return introspection::getArithmeticOperationCount<DYNAMICS>();
        }
      } else {
        throw std::domain_error("Can not introspect the introspection - we, sadly, are not in LISP land.");
      }
#endif // FEATURE_INSPECT_DYNAMICS
      return std::nullopt;
    }),
    _inspectMemoryBandwidth([]() -> std::optional<std::size_t> {
#ifdef FEATURE_INSPECT_DYNAMICS
      if constexpr (!std::is_same_v<T,Expr>) {
        if constexpr (dynamics::is_cse_optimized<DYNAMICS>::value) {
          return introspection::getMemoryBandwidthOfCollision<CSE<DYNAMICS>>();
        } else {
          return introspection::getMemoryBandwidthOfCollision<DYNAMICS>();
        }
      } else {
        throw std::domain_error("Can not introspect the introspection");
      }
#endif // FEATURE_INSPECT_DYNAMICS
      return std::nullopt;
    }),
    _inspectOptimizability([]() -> std::optional<bool> {
#ifdef FEATURE_INSPECT_DYNAMICS
      if constexpr (!std::is_same_v<T,Expr>) {
        return introspection::isOptimizable<DYNAMICS>();
      } else {
        return false;
      }
#else
      return std::nullopt;
#endif // FEATURE_INSPECT_DYNAMICS
    }),
    _inspectComplexity([]() -> std::optional<std::size_t> {
#ifdef FEATURE_INSPECT_DYNAMICS
      if constexpr (!std::is_same_v<T,Expr>) {
        if constexpr (dynamics::is_cse_optimized<DYNAMICS>::value) {
          return introspection::getComplexity<CSE<DYNAMICS>>();
        } else {
          return introspection::getComplexity<DYNAMICS>();
        }
      } else {
        throw std::domain_error("Can not introspect the introspection - we, sadly, are not in LISP land.");
      }
#endif // FEATURE_INSPECT_DYNAMICS
      return std::nullopt;
    })
  { }

  /// Returns type index of the promised DYNAMICS
  std::type_index id() const {
    return _id;
  }

  /// Returns name of promised DYNAMICS
  std::string name() const {
    return _name;
  }

  /// Returns set of accessed fields
  std::set<FieldTypePromise<T,DESCRIPTOR>> accessedFields() const {
    return _accessedFields();
  }

  /// Returns true if a CSE version is available
  bool hasOptimizedVersion() const {
    return _isOptimizationAvailable;
  }

  /// Returns whether a CSE version can be generated (if this can be determined)
  std::optional<bool> isOptimizable() const {
    return _inspectOptimizability();
  }

  /// Returns FLOPs per collision (if this can be determined)
  std::optional<std::size_t> getArithmeticOperationCount() const {
    return _inspectArithmeticOperationCount();
  }

  /// Returns memory bandwidth per collision (if this can be determined)
  std::optional<std::size_t> getMemoryBandwidth() const {
    if (auto bandwidth = _inspectMemoryBandwidth()) {
      return *bandwidth * sizeof(T);
    } else {
      return std::nullopt;
    }
  }

  /// Returns size of expanded expression tree (only for internal CSE generation use)
  std::optional<std::size_t> getComplexity() const {
    return _inspectComplexity();
  }

  bool operator<(const DynamicsPromise<T,DESCRIPTOR>& rhs) const {
    return id() < rhs.id();
  }

  /// Returns new instance of the promised DYNAMICS
  Dynamics<T,DESCRIPTOR>* realize() {
    return _makeDynamics();
  };

  /// Returns new instance of the collision operator for promised DYNAMICS
  template <Platform PLATFORM>
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>* realize() {
    AbstractCollisionO<T,DESCRIPTOR>* op = _makeOperator(PLATFORM);
    return static_cast<BlockCollisionO<T,DESCRIPTOR,PLATFORM>*>(op);
  };

};

template <typename DYNAMICS>
DynamicsPromise(meta::id<DYNAMICS>) -> DynamicsPromise<typename DYNAMICS::value_t,
                                                       typename DYNAMICS::descriptor_t>;


/// Mask describing the subdomain on which to apply the collision step
/**
 * By default this is the entire non-overlap area of the block lattice
 **/
struct CollisionSubdomainMask {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteBlockMask<T,PLATFORM>;
};

/// Map between cell indices and concrete dynamics
/**
 * Central class for managing and applying dynamics of / to a block lattice.
 *
 * Maintains both a map of cell indices to dynamics and
 * of dynamics to a mask of all assigned cell indices.
 *
 * Automatically determines the dominant dynamics and applies
 * the collision step to all cells via the dominant dynamics'
 * ConcreteBlockCollisionO.
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class BlockDynamicsMap {
private:
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _lattice;

  std::map<std::type_index,
           std::tuple<DynamicsPromise<T,DESCRIPTOR>,
                      std::unique_ptr<BlockCollisionO<T,DESCRIPTOR,PLATFORM>>>> _map;

  std::unique_ptr<Dynamics<T,DESCRIPTOR>*[]>                 _dynamicsOfCells;
  std::unique_ptr<BlockCollisionO<T,DESCRIPTOR,PLATFORM>*[]> _operatorOfCells;

  /// Subdomain on which to apply collisions
  ConcreteBlockMask<T,PLATFORM>& _coreMask;
  /// Pointer to collision operator with highest cell fraction
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>* _dominantCollisionO;

  /// Returns (Dynamics, CollisionO, Mask) tuple for promise
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>& resolve(DynamicsPromise<T,DESCRIPTOR>&& promise)
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
  /// Constructor for a BlockDynamicsMap
  BlockDynamicsMap(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice):
    _lattice(lattice),
    _dynamicsOfCells(new Dynamics<T,DESCRIPTOR>*                [_lattice.getNcells()] { nullptr }),
    _operatorOfCells(new BlockCollisionO<T,DESCRIPTOR,PLATFORM>*[_lattice.getNcells()] { nullptr }),
    _coreMask(lattice.template getData<CollisionSubdomainMask>()),
    _dominantCollisionO(nullptr)
  {
    _lattice.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> lattice) {
      _coreMask.set(_lattice.getCellId(lattice), true);
    });
  }

  /// Returns a pointer to dynamics matching the promise
  /**
   * Promise is only realized if dynamics were not previously constructed
   **/
  Dynamics<T,DESCRIPTOR>* get(DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    return resolve(std::forward<decltype(promise)>(promise)).getDynamics();
  }

  const Dynamics<T,DESCRIPTOR>* get(std::size_t iCell) const
  {
    return _dynamicsOfCells[iCell];
  }

  Dynamics<T,DESCRIPTOR>* get(std::size_t iCell)
  {
    return _dynamicsOfCells[iCell];
  }

  /// Assigns promised dynamics to cell index and updates block masks
  void set(std::size_t iCell, DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    auto& collisionO = resolve(std::forward<decltype(promise)>(promise));
    if (_operatorOfCells[iCell] != nullptr) {
      _operatorOfCells[iCell]->set(iCell, false, !_coreMask[iCell]);
    }
    collisionO.set(iCell, true, !_coreMask[iCell]);
    _operatorOfCells[iCell] = &collisionO;
    if (Dynamics<T,DESCRIPTOR>* dynamics = collisionO.getDynamics()) {
      _dynamicsOfCells[iCell] = dynamics;
    }
    _dominantCollisionO = nullptr;
  }

  /// Executes local collision step for entire non-overlap area of lattice
  /**
   * When using the default CollisionDispatchStrategy::Dominant, the most frequently
   * assigned dynamics resp. collision operator is applied and responsible for also
   * colliding cells excluded by its masks using their respective dynamics.
   *
   * Correspondingly, the alternative CollisionDispatchStrategy::Individual applies
   * all dynamics separately using a list-based approach.
   **/
  void collide(CollisionDispatchStrategy strategy)
  {
    switch (strategy) {
    case CollisionDispatchStrategy::Dominant:
      if (!_dominantCollisionO) {
        _dominantCollisionO = std::get<1>(std::max_element(
          _map.begin(),
          _map.end(),
          [](const auto& lhs, const auto& rhs) -> bool {
            return std::get<1>(lhs.second)->weight()
                 < std::get<1>(rhs.second)->weight();
          })->second).get();
      }
      _dominantCollisionO->apply(_lattice, _coreMask, strategy);
      break;

    case CollisionDispatchStrategy::Individual:
      for (auto& [id, value] : _map) {
        auto& [promise, collisionO] = value;
        if (collisionO->weight() > 0) {
          collisionO->apply(_lattice, _coreMask, strategy);
        }
      }
      break;

    default:
      throw std::runtime_error("Invalid collision dispatch strategy");
      break;
    }
  }

  /// Returns number of cells assigned to promised dynamics
  /**
   * Doesn't allocate, intended for introspection
   **/
  std::size_t getWeight(const DynamicsPromise<T,DESCRIPTOR>& promise) const {
    auto iter = _map.find(promise.id());
    if (iter != _map.end()) {
      const auto& [_, collisionO] = std::get<1>(*iter);
      return collisionO->weight();
    } else {
      return 0;
    }
  }

  /// Returns set of all dynamics maintained this map
  std::set<DynamicsPromise<T,DESCRIPTOR>> getAll() const {
    std::set<DynamicsPromise<T,DESCRIPTOR>> dynamics;
    for (auto& [id, value] : _map) {
      const auto& [promise, collisionO] = value;
      dynamics.emplace(promise);
    }
    return dynamics;
  }

};


}

#endif
