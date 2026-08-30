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

#ifndef REFINEMENT_OPERATOR_H
#define REFINEMENT_OPERATOR_H

#include "operatorScope.h"
#include "fields.h"

namespace olb {

template <typename T, typename DESCRIPTOR> class BlockD;

/// Context for the execution of block refinement operators
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockRefinementContextD;

/// Executor of a concrete block refinement operator in a given context
template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename COUPLER>
struct ConcreteBlockRefinementO;

/// Untyped promise to execute some block refinement operator
template <typename T, typename DESCRIPTOR>
class BlockRefinementOperatorPromise;

/// Interface for the contextual data of a block refinement
template <typename T, typename DESCRIPTOR>
struct BlockRefinementContextD {
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;

  using Data = olb::BlockD<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>>;

  virtual Platform getPlatform() const = 0;
  virtual void setProcessingContext(ProcessingContext context) = 0;

  /// Return reference to contextual data
  virtual Data& getData() = 0;

  /// Connect the coarse and fine lattice at position latticeR
  virtual void add(LatticeR<DESCRIPTOR::d> latticeR) = 0;
  /// Return index where contextual data of connection point latticeR is stored
  virtual std::optional<std::size_t> getDataIndex(LatticeR<DESCRIPTOR::d> latticeR) const = 0;

  /// Apply promised algorithm step to the previously added connection points
  /**
   * Scope of operator must match scope of context
   **/
  virtual void apply(BlockRefinementOperatorPromise<T,DESCRIPTOR>&& promise) = 0;

  template <Platform PLATFORM>
  auto& asConcrete() {
    return static_cast<ConcreteBlockRefinementContextD<T,DESCRIPTOR,PLATFORM>&>(*this);
  }

};

}

#endif
