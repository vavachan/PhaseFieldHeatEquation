/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
 *
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

#ifndef REFINEMENT_OPERATOR_PROMISE_H
#define REFINEMENT_OPERATOR_PROMISE_H

#include "operator.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class BlockRefinementOperatorPromise {
private:
  std::function<void(BlockRefinementContextD<T,DESCRIPTOR>&)> _executor;
  std::function<std::set<FieldTypePromise<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>>>()> _requiredData;
  bool _requiresContextNeigborAccess;
  std::type_index _operatorId;

public:
  template <typename OPERATOR>
  BlockRefinementOperatorPromise(meta::id<OPERATOR>):
    _executor([](BlockRefinementContextD<T,DESCRIPTOR>& context) {
      switch (context.getPlatform()) {
#ifdef PLATFORM_CPU_SISD
      case Platform::CPU_SISD:
        ConcreteBlockRefinementO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR>().apply(
          context.template asConcrete<Platform::CPU_SISD>());
        break;
#endif
#ifdef PLATFORM_CPU_SIMD
      case Platform::CPU_SIMD:
        ConcreteBlockRefinementO<T,DESCRIPTOR,Platform::CPU_SIMD,OPERATOR>().apply(
          context.template asConcrete<Platform::CPU_SIMD>());
        break;
#endif
#ifdef PLATFORM_GPU_CUDA
      case Platform::GPU_CUDA:
        ConcreteBlockRefinementO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR>().apply(
          context.template asConcrete<Platform::GPU_CUDA>());
        break;
#endif
      default:
        throw std::runtime_error("Refinement executor not implemented for PLATFORM");
      }
    }),
    _requiredData([]() -> std::set<FieldTypePromise<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>>> {
      std::set<FieldTypePromise<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>>> fields;
      OPERATOR::data::for_each([&fields](auto id) {
        fields.emplace(meta::id<Array<decltype(id.get())>>{});
      });
      return fields;
    }),
    _requiresContextNeigborAccess{
      OPERATOR::data::template contains<fields::refinement::CONTEXT_NEIGHBORS>()
    },
    _operatorId{typeid(OPERATOR)}
  { }

  std::set<FieldTypePromise<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>>> requiredData() const {
    return _requiredData();
  }

  bool requiresContextNeighborAccess() const {
    return _requiresContextNeigborAccess;
  }

  std::type_index id() const {
    return _operatorId;
  }

  void apply(BlockRefinementContextD<T,DESCRIPTOR>& context) {
    for (auto promise : requiredData()) {
      promise.ensureAvailabilityIn(context.getData());
    }
    _executor(context);
  }

};

}

#endif
