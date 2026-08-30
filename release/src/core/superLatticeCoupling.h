/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_COUPLING_H
#define SUPER_LATTICE_COUPLING_H

#include "operator.h"
#include "superLattice.h"

#include "utilities/typeMap.h"
#include "solver/names.h"

namespace olb {


template <typename COUPLER, typename COUPLEES, Platform PLATFORM>
class ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,OperatorScope::PerCell>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              PLATFORM>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  typename AbstractCouplingO<COUPLEES>::ParametersD _parameters;

  std::unique_ptr<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                    PLATFORM>> _mask;

  void execute(typename AbstractCouplingO<COUPLEES>::LatticeR latticeR)
  {
    CellID iCell = _lattices.template get<0>()->getCellId(latticeR);
    auto cells = _lattices.exchange_values([&](auto name) -> auto {
      return cpu::Cell{*_lattices.get(name), iCell};
    });
    COUPLER().apply(cells);
  }

public:
  template <typename LATTICES>
  ConcreteBlockCouplingO(LATTICES&& lattices):
    _lattices{lattices}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void set(CellID iCell, bool state) override
  {
    if (!_mask) {
      _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                                 PLATFORM>>(
        _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
      );
    }
    _mask->set(iCell, state);
  }

  void execute() override
  {
    using loc = typename AbstractCouplingO<COUPLEES>::LatticeR::value_t;
    auto* lattice = _lattices.template get<0>();
    if (_mask) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) collapse(1)
      #endif
      for (loc iX=0; iX < lattice->getNx(); ++iX) {
        for (loc iY=0; iY < lattice->getNy(); ++iY) {
          if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
            for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
              if (_mask->operator[](lattice->getCellId(iX,iY,iZ))) {
                execute({iX,iY,iZ});
              }
            }
          } else {
            if (_mask->operator[](lattice->getCellId(iX,iY))) {
              execute({iX,iY});
            }
          }
        }
      }
    } else {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) collapse(1)
      #endif
      for (loc iX=0; iX < lattice->getNx(); ++iX) {
        for (loc iY=0; iY < lattice->getNy(); ++iY) {
          if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
            for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
              execute({iX,iY,iZ});
            }
          } else {
            execute({iX,iY});
          }
        }
      }
    }
  }

};

template <typename COUPLER, typename COUPLEES, Platform PLATFORM>
class ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,OperatorScope::PerCellWithParameters>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              PLATFORM>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  typename COUPLER::parameters::template decompose_into<
    AbstractCouplingO<COUPLEES>::ParametersD::template include_fields
  > _parameters;

  std::unique_ptr<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                    PLATFORM>> _mask;

  void execute(typename AbstractCouplingO<COUPLEES>::LatticeR latticeR)
  {
    CellID iCell = _lattices.template get<0>()->getCellId(latticeR);
    auto cells = _lattices.exchange_values([&](auto name) -> auto {
      return cpu::Cell{*_lattices.get(name), iCell};
    });
    COUPLER().apply(cells, _parameters);
  }

public:
  template <typename LATTICES>
  ConcreteBlockCouplingO(LATTICES&& lattices):
    _lattices{lattices}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void set(CellID iCell, bool state) override
  {
    if (!_mask) {
      _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                                 PLATFORM>>(
        _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
      );
    }
    _mask->set(iCell, state);
  }

  void execute() override
  {
    using loc = typename AbstractCouplingO<COUPLEES>::LatticeR::value_t;
    auto* lattice = _lattices.template get<0>();
    if (_mask) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) collapse(1)
      #endif
      for (loc iX=0; iX < lattice->getNx(); ++iX) {
        for (loc iY=0; iY < lattice->getNy(); ++iY) {
          if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
            for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
              if (_mask->operator[](lattice->getCellId(iX,iY,iZ))) {
                execute({iX,iY,iZ});
              }
            }
          } else {
            if (_mask->operator[](lattice->getCellId(iX,iY))) {
              execute({iX,iY});
            }
          }
        }
      }
    } else {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static) collapse(1)
      #endif
      for (loc iX=0; iX < lattice->getNx(); ++iX) {
        for (loc iY=0; iY < lattice->getNy(); ++iY) {
          if constexpr (AbstractCouplingO<COUPLEES>::descriptor_t::d == 3) {
            for (loc iZ=0; iZ < lattice->getNz(); ++iZ) {
              execute({iX,iY,iZ});
            }
          } else {
            execute({iX,iY});
          }
        }
      }
    }
  }

};


/// Coupling operator COUPLER on named COUPLEES
template <typename COUPLER, typename COUPLEES>
class SuperLatticeCoupling {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = SuperLattice<typename VALUED_DESCRIPTOR::value_t,
                                      typename VALUED_DESCRIPTOR::descriptor_t>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  std::vector<std::unique_ptr<AbstractCouplingO<COUPLEES>>> _block;

  template <Platform PLATFORM>
  auto constructConcreteBlockCoupling(int iC)
  {
    auto block = _lattices.exchange_values([&](auto name) -> auto {
      using NAME = typename decltype(name)::type;
      using T = typename COUPLEES::template value<NAME>::value_t;
      using DESCRIPTOR = typename COUPLEES::template value<NAME>::descriptor_t;
      return dynamic_cast<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>*>(
        &_lattices.get(name)->getBlock(iC));
    });
    return std::make_unique<ConcreteBlockCouplingO<COUPLEES,PLATFORM,COUPLER,COUPLER::scope>>(block);
  }

public:
  template <typename... MAP>
  SuperLatticeCoupling(COUPLER,
                       MAP&&... args) {
    auto map = std::make_tuple(&args...);
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      constexpr unsigned idx = COUPLEES::keys_t::template index<name_t>();
      _lattices.template set<name_t>(std::get<2*idx+1>(map));
    });

    auto& load = _lattices.template get<0>()->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      Platform reference = _lattices.template get<0>()->getBlock(iC).getPlatform();
      _lattices.for_each([&](auto name, auto lattice) {
        if (lattice->getBlock(iC).getPlatform() != reference) {
          throw std::runtime_error("Platforms of coupled block lattices must match");
        }
      });
      callUsingConcretePlatform(reference,
                                [&](auto platform) {
        _block.emplace_back(constructConcreteBlockCoupling<platform.value>(iC));
      });
    }
  }

  /// Execute coupling operation on all blocks
  void execute()
  {
    auto& load = _lattices.template get<0>()->getLoadBalancer();
    #ifdef PARALLEL_MODE_OMP
    #pragma omp taskloop
    #endif
    for (int iC = 0; iC < load.size(); ++iC) {
      _block[iC]->execute();
    }
  }

  /// Set coupling parameter FIELD
  template <typename FIELD>
  void setParameter(typename AbstractCouplingO<COUPLEES>::template FieldD<FIELD>&& field)
  {
    auto& load = _lattices.template get<0>()->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      _block[iC]->getParameters().template set<FIELD>(std::forward<decltype(field)>(field));
    }
  }

  /// Get coupling parameter FIELD for block 0
  template <typename FIELD>
  auto getParameter() const -> typename AbstractCouplingO<COUPLEES>::template FieldD<FIELD>
  {
    return _block[0]->getParameters().template get<FIELD>();
  }

  /// Return block-level abstract coupling operator
  /**
   * e.g. to set block-wise mask or parameter:
   *
   * superCoupling.getBlock(iC).set(latticeR, false);
   **/
  AbstractCouplingO<COUPLEES>& getBlock(int iC)
  {
    return *_block[iC];
  }

  /// Restrict coupling operation to indicated locations
  /**
   * By default the entire (non-overlap) lattice area is coupled.
   **/
  void restrictTo(FunctorPtr<SuperIndicatorF<typename AbstractCouplingO<COUPLEES>::value_t,
                                             AbstractCouplingO<COUPLEES>::descriptor_t::d>>&& indicator)
  {
    using DESCRIPTOR = typename COUPLEES::values_t::template get<0>::descriptor_t;
    auto& load = _lattices.template get<0>()->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      const auto& blockL = _lattices.template get<0>()->getBlock(iC);
      auto& blockI = indicator->getBlockIndicatorF(iC);
      blockL.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
        _block[iC]->set(blockL.getCellId(latticeR), blockI(latticeR));
      });
    }
  }

};

template <typename COUPLER, typename... MAP>
SuperLatticeCoupling(COUPLER, MAP&&...) -> SuperLatticeCoupling<
  COUPLER,
  typename meta::map<MAP...>::template map_values<descriptors::extract_valued_descriptor_t>
>;

template <typename COUPLER, typename... MAP>
auto constructSharedCoupling(COUPLER, MAP&&... map) {
  return std::make_shared<SuperLatticeCoupling<
    COUPLER,
    typename meta::map<MAP...>::template map_values<descriptors::extract_valued_descriptor_t>
  >>(COUPLER{}, std::forward<decltype(map)>(map)...);
}

}

#endif
