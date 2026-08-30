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

#ifndef REFINEMENT_SUPER_LATTICE_REFINEMENT_H
#define REFINEMENT_SUPER_LATTICE_REFINEMENT_H

#include "operatorScope.h"
#include "operatorPromise.h"

#include "loadBalancer.h"

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockRefinementContextD final : public BlockRefinementContextD<T,DESCRIPTOR> {
private:
  /// Context data for coupling cells only
  ConcreteBlockD<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>,PLATFORM>& _data;
  /// Map between fine lattice cell index and context data cell index
  std::map<CellID, CellID> _dataIndex;
  /// True if data index needs to be recomputed
  bool _dataIndexChanged;

  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _cLattice;
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _fLattice;

  void updateContextNeighborIndices() {
    using namespace fields::refinement;
    const auto& cellIdFine = _data.template getField<CELL_ID_FINE>();
    for (std::size_t i=0; i < cellIdFine.getSize(); ++i) {
      auto latticeR = _fLattice.getLatticeR(cellIdFine[0][i]);
      FieldD<T,DESCRIPTOR,CONTEXT_NEIGHBORS> neighbors(std::numeric_limits<CellID>::max());
      for (unsigned iN=0; iN < CONTEXT_NEIGHBORS::count<DESCRIPTOR>(); ++iN) {
        if (auto index = getDataIndex(latticeR + CONTEXT_NEIGHBORS::c<DESCRIPTOR>(iN))) {
          neighbors[iN] = *index;
        }
      }
      _data.template getField<CONTEXT_NEIGHBORS>().set(i, neighbors);
    }
  }

public:
  ConcreteBlockRefinementContextD(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& cLattice,
                                  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& fLattice,
                                  ConcreteBlockD<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>,PLATFORM>& data)
  : _data{data}
  , _dataIndex{}
  , _dataIndexChanged{true}
  , _cLattice{cLattice}
  , _fLattice{fLattice}
  {
    auto size = _data.getCore();
    size[0] = 0;
    for (unsigned iD=1; iD < DESCRIPTOR::d; ++iD) {
      size[iD] = 1;
    }
    _data.resize(size);
  }

  Platform getPlatform() const override {
    return PLATFORM;
  }

  void setProcessingContext(ProcessingContext context) override {
    _data.setProcessingContext(context);
  }

  auto& getFineLattice() {
    return _fLattice;
  }

  auto& getCoarseLattice() {
    return _cLattice;
  }

  auto& getConcreteData() {
    return _data;
  }

  typename BlockRefinementContextD<T,DESCRIPTOR>::Data& getData() override {
    return _data;
  }

  void add(LatticeR<DESCRIPTOR::d> latticeR) override {
    CellID iCell = _fLattice.getCellId(latticeR);
    if (_dataIndex.find(iCell) == _dataIndex.end()) {
      _dataIndexChanged = true;
      _dataIndex[iCell] = _data.getNcells();

      auto size = _data.getCore();
      size[0] += 1;
      _data.resize(size);

      _data.template getField<fields::refinement::CELL_ID_FINE>().set(_dataIndex[iCell], iCell);

      bool coIncident = true;
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        coIncident &= latticeR[iD] % 2 == 0;
      }
      if (coIncident) {
        _data.template getField<fields::refinement::CELL_ID_COARSE>().set(
          _dataIndex[iCell], _cLattice.getCellId(latticeR / 2));
      } else {
        _data.template getField<fields::refinement::CELL_ID_COARSE>().set(_dataIndex[iCell], 0);
      }
    }
  }

  std::optional<std::size_t> getDataIndex(LatticeR<DESCRIPTOR::d> latticeR) const override {
    auto iter = _dataIndex.find(_fLattice.getCellId(latticeR));
    if (iter != _dataIndex.end()) {
      return std::get<1>(*iter);
    } else {
      return std::nullopt;
    }
  }

  void apply(BlockRefinementOperatorPromise<T,DESCRIPTOR>&& promise) override {
    if (_dataIndexChanged && promise.requiresContextNeighborAccess()) {
      updateContextNeighborIndices();
      _data.template getField<fields::refinement::CONTEXT_NEIGHBORS>()
           .setProcessingContext(ProcessingContext::Simulation);
      _dataIndexChanged = false;
    }
    promise.apply(*this);
  }

};

template <typename T, typename DESCRIPTOR>
class SuperLatticeRefinement {
private:
  SuperD<T,refinement::DATA_DESCRIPTOR<DESCRIPTOR>> _data;

  SuperLattice<T,DESCRIPTOR>& _cLattice;
  SuperLattice<T,DESCRIPTOR>& _fLattice;

  const RefinedLoadBalancer<T,DESCRIPTOR::d>& _balancer;

  std::vector<std::unique_ptr<BlockRefinementContextD<T,DESCRIPTOR>>> _context;

  template <Platform PLATFORM>
  auto constructConcreteBlockContext(int fineC) {
    const int coarseC = _balancer.cloc(fineC);
    return std::make_unique<ConcreteBlockRefinementContextD<T,DESCRIPTOR,PLATFORM>>(
        _cLattice.getBlock(coarseC).template asConcrete<PLATFORM>()
      , _fLattice.getBlock(fineC).template asConcrete<PLATFORM>()
      , _data.getBlock(fineC).template asConcrete<PLATFORM>());
  }

  std::map<std::type_index,
           std::unique_ptr<SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>>> _cCommunicator;
  std::map<std::type_index,
           std::unique_ptr<SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>>> _fCommunicator;

public:
  SuperLatticeRefinement(SuperLattice<T,DESCRIPTOR>& cLattice,
                         SuperLattice<T,DESCRIPTOR>& fLattice,
                         const LoadBalancer<T>& balancer)
  : _data(balancer)
  , _cLattice(cLattice)
  , _fLattice(fLattice)
  , _balancer(dynamic_cast<const RefinedLoadBalancer<T,DESCRIPTOR::d>&>(balancer))
  {
    for (int iC = 0; iC < _balancer.size(); ++iC) {
      callUsingConcretePlatform(_balancer.platform(iC),
                                [&](auto platform) {
        _context.emplace_back(constructConcreteBlockContext<platform.value>(iC));
      });
    }
  }

  BlockRefinementContextD<T,DESCRIPTOR>& getBlock(int iC) {
    return *_context[iC];
  }

  /// Apply coupling operation on all blocks
  void apply(BlockRefinementOperatorPromise<T,DESCRIPTOR>&& promise) {
    {
      auto iter = _cCommunicator.find(promise.id());
      if (iter != _cCommunicator.end()) {
        std::get<1>(*iter)->communicate();
      }
    }
    {
      auto iter = _fCommunicator.find(promise.id());
      if (iter != _fCommunicator.end()) {
        std::get<1>(*iter)->communicate();
      }
    }

    #ifdef PARALLEL_MODE_OMP
    #pragma omp taskloop
    #endif
    for (int iC = 0; iC < _balancer.size(); ++iC) {
      _context[iC]->apply(std::forward<BlockRefinementOperatorPromise<T,DESCRIPTOR>&&>(promise));
    }
  }

  void setProcessingContext(ProcessingContext context) {
    for (int iC = 0; iC < _balancer.size(); ++iC) {
      _context[iC]->setProcessingContext(context);
    }
  }

  SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>& getCoarseCommunicator(
    BlockRefinementOperatorPromise<T,DESCRIPTOR>&& promise) {
    auto iter = _cCommunicator.find(promise.id());
    if (iter == _cCommunicator.end()) {
      iter = std::get<0>(_cCommunicator.emplace(promise.id(),
                                                std::make_unique<SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>>(_cLattice)));
    }
    return *std::get<1>(*iter);
  }

  SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>& getFineCommunicator(
    BlockRefinementOperatorPromise<T,DESCRIPTOR>&& promise) {
    auto iter = _fCommunicator.find(promise.id());
    if (iter == _fCommunicator.end()) {
      iter = std::get<0>(_fCommunicator.emplace(promise.id(),
                                                std::make_unique<SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>>(_fLattice)));
    }
    return *std::get<1>(*iter);
  }

};

}

#endif
