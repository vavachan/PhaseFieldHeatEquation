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

#ifndef GPU_CUDA_OPERATOR_H
#define GPU_CUDA_OPERATOR_H

#include "core/operator.h"

#include "dynamics/dynamics.h"

#include "device.h"
#include "mask.h"
#include "context.h"
#include "dynamics.h"

namespace olb {

namespace gpu::cuda {

template <typename T, typename DESCRIPTOR, typename... DYNAMICS>
std::function<void(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>&)> getFusedCollisionO();

}


/// Application of the collision step on a concrete CUDA block
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS> final
  : public BlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  std::unique_ptr<DYNAMICS> _dynamics;

  ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,typename DYNAMICS::parameters>* _parameters;
  ConcreteBlockMask<T,Platform::GPU_CUDA>* _mask;

  gpu::cuda::Dynamics<T,DESCRIPTOR>** _dynamicsOfCells;

  gpu::cuda::device::unique_ptr<gpu::cuda::Dynamics<T,DESCRIPTOR>> _deviceDynamics;

  gpu::cuda::Column<CellID> _cells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

  /// Apply DYNAMICS using its mask and fall back to dynamic dispatch for others
  void applyDominant(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
                     ConcreteBlockMask<T,Platform::GPU_CUDA>&               subdomain);
  /// Apply only DYNAMICS, do not apply others
  void applyIndividual(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
                       ConcreteBlockMask<T,Platform::GPU_CUDA>&               subdomain);

public:
  ConcreteBlockCollisionO();

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
    if constexpr (!std::is_same_v<DYNAMICS,NoDynamics<T,DESCRIPTOR>>) {
      if (!overlap) {
        _mask->set(iCell, state);
      }
    }
    if (state) {
      _dynamicsOfCells[iCell] = _deviceDynamics.get();
    }
    _modified = true;
  }

  Dynamics<T,DESCRIPTOR>* getDynamics() override
  {
    return _dynamics.get();
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override;

  /// Apply collision to subdomain of block using strategy
  /**
   * The subdomain argument is currently assumed to be the core mask of BlockDynamicsMap
   **/
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
             ConcreteBlockMask<T,Platform::GPU_CUDA>&               subdomain,
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


/// Application of a cell-wise OPERATOR on a concrete CUDA block
template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCell> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  gpu::cuda::Column<CellID> _cells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

public:
  ConcreteBlockO();

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
      _cells.push_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    _modified = false;
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override;

};


/// Application of a parametrized cell-wise OPERATOR on a concrete CUDA block
template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCellWithParameters> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  gpu::cuda::Column<CellID> _cells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

  ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,typename OPERATOR::parameters>* _parameters;

public:
  ConcreteBlockO();

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
      _cells.push_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    _modified = false;
    _parameters = &block.template getData<OperatorParameters<OPERATOR>>();
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override;

};


/// Application of a block-wise OPERATOR on a concrete CUDA block
/**
 * e.g. StatisticsPostProcessor
 **/
template <typename T, typename DESCRIPTOR, concepts::BlockOperator OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerBlock> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
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

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    OPERATOR().setup(block);
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    OPERATOR().apply(block);
  }

};


/// Application of a block-wise COUPLER on concrete CUDA COUPLEES
template <typename COUPLER, typename COUPLEES>
class ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCell>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              Platform::GPU_CUDA>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  typename AbstractCouplingO<COUPLEES>::ParametersD _parameters;

  std::unique_ptr<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                    Platform::GPU_CUDA>> _mask;

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

  void set(CellID iCell, bool state) override;

  void execute() override;

};

/// Application of a block-wise COUPLER on concrete CUDA COUPLEES with parameters
template <typename COUPLER, typename COUPLEES>
class ConcreteBlockCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                                              typename VALUED_DESCRIPTOR::descriptor_t,
                                              Platform::GPU_CUDA>*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  typename COUPLER::parameters::template decompose_into<
    AbstractCouplingO<COUPLEES>::ParametersD::template include_fields
  > _parameters;

  std::unique_ptr<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                    Platform::GPU_CUDA>> _mask;

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

  void set(CellID iCell, bool state) override;

  void execute() override;

};

/// Application of a block-wise particle COUPLER on concrete CUDA COUPLEES with parameters
template <typename COUPLER, typename COUPLEES>
class ConcreteBlockPointCouplingO<COUPLEES,Platform::GPU_CUDA,COUPLER,OperatorScope::PerCellWithParameters>
  final : public AbstractCouplingO<COUPLEES> {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = std::conditional_t<
    VALUED_DESCRIPTOR::descriptor_t::template provides<descriptors::POPULATION>(),
    ConcreteBlockLattice<typename VALUED_DESCRIPTOR::value_t,
                         typename VALUED_DESCRIPTOR::descriptor_t,
                         Platform::GPU_CUDA>,
    ConcreteBlockD<typename VALUED_DESCRIPTOR::value_t,
                          typename VALUED_DESCRIPTOR::descriptor_t,
                          Platform::GPU_CUDA>
  >*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  Cuboid<
    typename AbstractCouplingO<COUPLEES>::value_t,
    AbstractCouplingO<COUPLEES>::descriptor_t::d
  >& _cuboid;

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
    _cuboid{cuboid}
  { }

  std::type_index id() const override {
    return typeid(COUPLER);
  }

  typename AbstractCouplingO<COUPLEES>::AbstractParameters& getParameters() override {
    return _parameters;
  }

  void set(CellID iCell, bool state) override;

  void execute() override;

};

}

#include "communicator.h"
#include "statistics.h"

#endif
