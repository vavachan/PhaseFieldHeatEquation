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

#ifndef SISD_OPERATOR_H
#define SISD_OPERATOR_H

#include "mask.h"

#include "core/meta.h"
#include "core/operator.h"

#include "core/platform/cpu/cell.h"
#include "core/latticeStatistics.h"

#include "dynamics/dynamics.h"

namespace olb {

namespace cpu {

/// Implementations of scalar CPU specifics
namespace sisd {

/// Implementation of cpu::Dynamics for concrete DYNAMICS on SISD blocks
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* _parameters;

public:
  ConcreteDynamics(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* parameters):
    _parameters{parameters} {
  }

  CellStatistic<T> collide(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell) override {
    return DYNAMICS().collide(cell, *_parameters);
  }

  T computeRho(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell) override {
    return typename DYNAMICS::MomentaF().computeRho(cell);
  }
  void computeU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T* u) override {
    typename DYNAMICS::MomentaF().computeU(cell, u);
  }
  void computeJ(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T* j) override {
    typename DYNAMICS::MomentaF().computeJ(cell, j);
  }
  void computeRhoU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u) override {
    typename DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  }
  void computeStress(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  }
  void computeAllMomenta(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
  }

  T getOmegaOrFallback(T fallback) override {
    if constexpr (DYNAMICS::parameters::template contains<descriptors::OMEGA>()) {
      return _parameters->template get<descriptors::OMEGA>();
    } else {
      return fallback;
    }
    __builtin_unreachable();
  }

  void computeEquilibrium(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T rho, T* u, T* fEq) override {
    typename DYNAMICS::EquilibriumF().compute(cell, rho, u, fEq);
  };

  void defineRho(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho) override {
    typename DYNAMICS::MomentaF().defineRho(cell, rho);
  }

  void defineU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T* u) override {
    typename DYNAMICS::MomentaF().defineU(cell, u);
  }

  void defineRhoU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u) override {
    typename DYNAMICS::MomentaF().defineRhoU(cell, rho, u);
  }

  void defineAllMomenta(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().defineAllMomenta(cell, rho, u, pi);
  }

  void inverseShiftRhoU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u) override {
    typename DYNAMICS::MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};

}

}

/// Application of the collision step on a concrete SISD block
/**
 * ConcreteBlockCollisionO::apply allows for applying DYNAMICS
 * on a given cell index or an entire block. The latter option
 * accepts a ConcreteBlockMask instance describing the subset
 * of cells for which DYNAMICS is to be called directly
 * instead of via roundtrip.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::CPU_SISD,DYNAMICS> final
  : public BlockCollisionO<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  std::unique_ptr<DYNAMICS> _dynamics;
  std::unique_ptr<cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD>> _concreteDynamics;

  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* _parameters;
  ConcreteBlockMask<T,Platform::CPU_SISD>* _mask;

  cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD>** _dynamicsOfCells;

  std::vector<CellID> _cells;
  bool _modified;

  /// Apply DYNAMICS using its mask and fall back to dynamic dispatch for others
  /**
   * Loop excludes overlap areas of block as collisions are never applied there.
   **/
  void applyDominant(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block,
                     ConcreteBlockMask<T,Platform::CPU_SISD>&               subdomain)
  {
    auto& parameters = *_parameters;
    auto& mask = *_mask;
    typename LatticeStatistics<T>::Aggregatable statistics{};
    #ifdef PARALLEL_MODE_OMP
    #pragma omp declare reduction(+ : typename LatticeStatistics<T>::Aggregatable : omp_out += omp_in) initializer (omp_priv={})
    #endif

    if constexpr (DESCRIPTOR::d == 3) {
      if (block.statisticsEnabled()) {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(dynamic,1) reduction(+ : statistics)
        #endif
        for (int iX=0; iX < block.getNx(); ++iX) {
          for (int iY=0; iY < block.getNy(); ++iY) {
            std::size_t iCell = block.getCellId(iX,iY,0);
            for (int iZ=0; iZ < block.getNz(); ++iZ) {
              cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
              if (auto cellStatistic = mask[iCell] ? DYNAMICS().collide(cell, parameters)
                                                   : _dynamicsOfCells[iCell]->collide(cell)) {
                statistics.increment(cellStatistic.rho, cellStatistic.uSqr);
              }
              iCell += 1;
            }
          }
        }
      } else {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for (int iX=0; iX < block.getNx(); ++iX) {
          for (int iY=0; iY < block.getNy(); ++iY) {
            std::size_t iCell = block.getCellId(iX,iY,0);
            for (int iZ=0; iZ < block.getNz(); ++iZ) {
              cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
              if (mask[iCell]) [[likely]] {
                DYNAMICS().collide(cell, parameters);
              } else {
                _dynamicsOfCells[iCell]->collide(cell);
              }
              iCell += 1;
            }
          }
        }
      }
    } else {
      if (block.statisticsEnabled()) {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(dynamic,1) reduction(+ : statistics)
        #endif
        for (int iX=0; iX < block.getNx(); ++iX) {
          std::size_t iCell = block.getCellId(iX,0);
          for (int iY=0; iY < block.getNy(); ++iY) {
            cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
            if (auto cellStatistic = mask[iCell] ? DYNAMICS().collide(cell, parameters)
                                                 : _dynamicsOfCells[iCell]->collide(cell)) {
              statistics.increment(cellStatistic.rho, cellStatistic.uSqr);
            }
            iCell += 1;
          }
        }
      } else {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(dynamic,1)
        #endif
        for (int iX=0; iX < block.getNx(); ++iX) {
          std::size_t iCell = block.getCellId(iX,0);
          for (int iY=0; iY < block.getNy(); ++iY) {
            cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
            if (mask[iCell]) [[likely]] {
              DYNAMICS().collide(cell, parameters);
            } else {
              _dynamicsOfCells[iCell]->collide(cell);
            }
            iCell += 1;
          }
        }
      }
    }

    if (block.statisticsEnabled()) {
      block.getStatistics().incrementStats(statistics);
    }
  }

  /// Apply only DYNAMICS, do not apply others
  void applyIndividual(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block,
                       ConcreteBlockMask<T,Platform::CPU_SISD>&               subdomain)
  {
    // Update cell list from mask
    if (_modified) {
      _cells.clear();
      for (CellID iCell=0; iCell  < block.getNcells(); ++iCell) {
        if (_mask->operator[](iCell)) {
          _cells.push_back(iCell);
        }
      }
      _modified = false;
    }

    auto& parameters = *_parameters;
    typename LatticeStatistics<T>::Aggregatable statistics{};

    #ifdef PARALLEL_MODE_OMP
    #pragma omp declare reduction(+ : typename LatticeStatistics<T>::Aggregatable : omp_out += omp_in) initializer (omp_priv={})
    #endif

    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static) reduction(+ : statistics)
    #endif
    for (std::size_t i=0; i < _cells.size(); ++i) {
      std::size_t iCell = _cells[i];
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
      if (auto cellStatistic = DYNAMICS().collide(cell, parameters)) {
        statistics.increment(cellStatistic.rho, cellStatistic.uSqr);
      }
    }

    block.getStatistics().incrementStats(statistics);
  }

public:
  ConcreteBlockCollisionO():
    _dynamics(new DYNAMICS()),
    _parameters(nullptr),
    _mask(nullptr),
    _cells(0),
    _modified(true)
  { }

  std::type_index id() const override
  {
    return typeid(DYNAMICS);
  }

  std::size_t weight() const override
  {
    return _mask->weight();
  }

  void set(CellID iCell, bool state, bool overlap) override
  {
    /// Only unmask cells that actually do something
    //if constexpr (!std::is_same_v<DYNAMICS,NoDynamics<T,DESCRIPTOR>>) {
      if (!overlap) {
        _mask->set(iCell, state);
      }
    //}
    if (state) {
      _dynamicsOfCells[iCell] = _concreteDynamics.get();
    }
    _modified = true;
  }

  Dynamics<T,DESCRIPTOR>* getDynamics() override
  {
    return _dynamics.get();
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    _parameters = &block.template getData<OperatorParameters<DYNAMICS>>().parameters;
    _mask = &block.template getData<DynamicsMask<DYNAMICS>>();
    if constexpr (dynamics::has_parametrized_momenta_v<DYNAMICS>) {
      _dynamics->setMomentaParameters(_parameters);
    }

    _concreteDynamics.reset(new cpu::sisd::ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(_parameters));
    // Fetch pointer to concretized dynamic-dispatch field
    _dynamicsOfCells = block.template getField<cpu::DYNAMICS<T,DESCRIPTOR,Platform::CPU_SISD>>()[0].data();
  }

  /// Apply collision on subdomain of block
  /**
   * This assumes that `subdomain` is the core mask of BlockDynamicsMap.
   **/
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block,
             ConcreteBlockMask<T,Platform::CPU_SISD>&               subdomain,
             CollisionDispatchStrategy                              strategy) override
  {
    switch (strategy) {
    case CollisionDispatchStrategy::Dominant:
      return applyDominant(block, subdomain);
    case CollisionDispatchStrategy::Individual:
      return applyIndividual(block, subdomain);
    default:
      throw std::runtime_error("Invalid collision dispatch strategy");
    }
  }

};


/// Application of a cell-wise OPERATOR on a concrete scalar CPU block
template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerCell> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  std::vector<CellID> _cells;
  bool _modified;

public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  std::size_t weight() const override
  {
    return _cells.size();
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _cells.emplace_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  { }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    if (_modified) {
      std::sort(_cells.begin(), _cells.end());
      _cells.erase(std::unique(_cells.begin(), _cells.end()), _cells.end());
      _modified = false;
    }
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, 0);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) firstprivate(cell)
      #endif
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell);
      }
    }
  }

};


template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerCellWithParameters> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  std::vector<CellID> _cells;
  bool _modified;

  ParametersOfOperatorD<T,DESCRIPTOR,OPERATOR>* _parameters;

public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  std::size_t weight() const override
  {
    return _cells.size();
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _cells.emplace_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    _parameters = &block.template getData<OperatorParameters<OPERATOR>>().parameters;
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    if (_modified) {
      std::sort(_cells.begin(), _cells.end());
      _cells.erase(std::unique(_cells.begin(), _cells.end()), _cells.end());
      _modified = false;
    }
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, 0);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) firstprivate(cell)
      #endif
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell, *_parameters);
      }
    }
  }

};


/// Application of a block-wise OPERATOR on a concrete scalar CPU block
/**
 * e.g. StatisticsPostProcessor
 **/
template <typename T, typename DESCRIPTOR, concepts::BlockOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerBlock> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  std::size_t weight() const override
  {
    return 0;
  }

  void set(CellID iCell, bool state) override
  {
    throw std::logic_error("BlockO::set not supported for OperatorScope::PerBlock");
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    OPERATOR().setup(block);
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    OPERATOR().apply(block);
  }

};

}

#endif
