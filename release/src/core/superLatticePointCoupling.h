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

#ifndef SUPER_PARTICLE_LATTICE_COUPLING_H
#define SUPER_PARTICLE_LATTICE_COUPLING_H

#include "operator.h"
#include "superLattice.h"

#include "utilities/typeMap.h"
#include "solver/names.h"

namespace olb {

template <typename COUPLER, typename COUPLEES, Platform PLATFORM>
class ConcreteBlockPointCouplingO<COUPLEES,PLATFORM,COUPLER,OperatorScope::PerCellWithParameters>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = std::conditional_t<
    VALUED_DESCRIPTOR::descriptor_t::template provides<descriptors::POPULATION>(),
    ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                         typename VALUED_DESCRIPTOR::descriptor_t,
                         PLATFORM>,
    ConcreteBlockD<typename VALUED_DESCRIPTOR::value_t,
                          typename VALUED_DESCRIPTOR::descriptor_t,
                          PLATFORM>
  >*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  Cuboid<typename AbstractCouplingO<COUPLEES>::value_t,
         AbstractCouplingO<COUPLEES>::descriptor_t::d>& _cuboid;

  typename AbstractCouplingO<COUPLEES>::ParametersD::template include<
    typename meta::merge<
      typename COUPLER::parameters,
      meta::list<
        fields::converter::PHYS_DELTA_X,
        fields::BLOCK_LOWER
      >
    >
  > _parameters;

public:
  template <typename LATTICES>
  ConcreteBlockPointCouplingO(LATTICES&& lattices,
                              Cuboid<typename AbstractCouplingO<COUPLEES>::value_t,
                              AbstractCouplingO<COUPLEES>::descriptor_t::d>& cuboid):
    _lattices{lattices},
    _cuboid{cuboid},
    _parameters{}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void execute() override
  {
    using V = typename AbstractCouplingO<COUPLEES>::value_t;
    using DESCRIPTOR = typename AbstractCouplingO<COUPLEES>::descriptor_t;
    auto* points = _lattices.get(meta::id<names::Points>{});
    const auto originR = _parameters.template get<fields::BLOCK_LOWER>();
    const auto deltaX = _parameters.template get<fields::converter::PHYS_DELTA_X>();
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static)
    #endif
    for (CellID iP=0; iP < points->getNcells(); ++iP) {
      const auto physR = points->template getField<fields::PHYS_R>().get(iP);
      const LatticeR<DESCRIPTOR::d> latticeR = util::floor((physR - originR) / deltaX + V{0.5});
      const CellID iCell = _lattices.template get<0>()->getCellId(latticeR);
      auto cells = _lattices.exchange_values([&](auto name) -> auto {
        if constexpr (std::is_same_v<typename decltype(name)::type, names::Points>) {
          return cpu::Row{*_lattices.get(name), iP};
        } else {
          return cpu::Cell{*_lattices.get(name), iCell};
        }
      });
      COUPLER().apply(cells, _parameters);
    }
  }

  void set(CellID iCell, bool state) override {
    throw std::bad_function_call();
  }

};


/// Coupling operator COUPLER on named COUPLEES
template <typename COUPLER, typename COUPLEES>
class SuperLatticePointCoupling {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = std::conditional_t<
    VALUED_DESCRIPTOR::descriptor_t::template provides<descriptors::POPULATION>(),
    SuperLattice<typename VALUED_DESCRIPTOR::value_t,
                 typename VALUED_DESCRIPTOR::descriptor_t>,
    SuperD<typename VALUED_DESCRIPTOR::value_t,
                  typename VALUED_DESCRIPTOR::descriptor_t>
  >*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  std::vector<std::unique_ptr<AbstractCouplingO<COUPLEES>>> _block;

  template <Platform PLATFORM>
  auto constructConcreteBlockPointCoupling(int iC)
  {
    auto block = _lattices.exchange_values([&](auto name) -> auto {
      using NAME = typename decltype(name)::type;
      using T = typename COUPLEES::template value<NAME>::value_t;
      using DESCRIPTOR = typename COUPLEES::template value<NAME>::descriptor_t;
      if constexpr (DESCRIPTOR::template provides<descriptors::POPULATION>()) {
        return dynamic_cast<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>*>(
          &_lattices.get(name)->getBlock(iC));
      } else {
        return dynamic_cast<ConcreteBlockD<T,DESCRIPTOR,PLATFORM>*>(
          &_lattices.get(name)->getBlock(iC));
      }
    });
    return std::make_unique<ConcreteBlockPointCouplingO<COUPLEES,PLATFORM,COUPLER,COUPLER::scope>>(
      block,
      _lattices.template get<0>()->getCuboidDecomposition().get(
        _lattices.template get<0>()->getLoadBalancer().glob(iC)));
  }

public:
  template <typename... MAP>
  SuperLatticePointCoupling(COUPLER, MAP&&... args)
  {
    auto map = std::make_tuple(&args...);
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      constexpr unsigned idx = COUPLEES::keys_t::template index<name_t>();
      _lattices.template set<name_t>(std::get<2*idx+1>(map));
    });

    auto& load = _lattices.template get<0>()->getLoadBalancer();
    auto& cGeometry = _lattices.template get<0>()->getCuboidDecomposition();

    for (int iC = 0; iC < load.size(); ++iC) {
      Platform reference = _lattices.template get<0>()->getBlock(iC).getPlatform();

      _lattices.for_each([&](auto name, auto lattice) {
        if (lattice->getBlock(iC).getPlatform() != reference) {
          throw std::runtime_error("Platforms of coupled block lattices must match");
        }
      });

      callUsingConcretePlatform(reference,
                                [&](auto platform) {
        _block.emplace_back(constructConcreteBlockPointCoupling<platform.value>(iC));
      });
    }

    for (int iC = 0; iC < load.size(); ++iC) {
      const auto origin = cGeometry.get(load.glob(iC)).getOrigin();
      const auto deltaX = cGeometry.get(load.glob(iC)).getDeltaR();
      _block[iC]->getParameters().template set<fields::BLOCK_LOWER>(origin);
      _block[iC]->getParameters().template set<fields::converter::PHYS_DELTA_X>(deltaX);
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

};

/// SuperPointLatticeCoupling
template <typename COUPLER, typename... MAP>
SuperLatticePointCoupling(COUPLER, MAP&&...)
  -> SuperLatticePointCoupling<
       COUPLER,
       typename meta::map<MAP...>::template map_values<descriptors::extract_valued_descriptor_t>
>;

}

#endif
