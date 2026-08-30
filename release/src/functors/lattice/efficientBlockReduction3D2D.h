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

#ifndef EFFICIENT_BLOCK_REDUCTION_3D2D_H
#define EFFICIENT_BLOCK_REDUCTION_3D2D_H

#include "core/blockData.h"
#include "core/vector.h"
#include "blockBaseF2D.h"
#include "superBaseF2D.h"
#include "superBaseF3D.h"

#include "utilities/hyperplane3D.h"
#include "utilities/hyperplaneLattice3D.h"
#include "utilities/blockDataSyncMode.h"
#include "utilities/blockDataReductionMode.h"

#include "core/superLatticePointCoupling.h"

#include <tuple>
#include <map>

namespace olb {

struct VelocityNormF {
  using parameters = meta::list<fields::converter::PHYS_VELOCITY>;

  using result_field = descriptors::SCALAR;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    Vector<V,DESCRIPTOR::d> u{};
    cell.computeU(u.data());
    u *= parameters.template get<fields::converter::PHYS_VELOCITY>();
    return norm(u);
  }
};

template <typename FUNCTOR>
struct BlockReduction3D2DO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = FUNCTOR::parameters;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    auto point = cells.template get<names::Points>();
    auto cell = cells.template get<names::Lattice>();
    auto result = FUNCTOR().compute(cell, parameters);
    point.template setField<typename FUNCTOR::result_field>(result);
  }
};

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
class EfficientBlockReduction3D2D final : public HyperplaneLattice3D<T>, public BlockDataF2D<T,T> {
private:
  SuperLattice<T,DESCRIPTOR>& _sLattice;

  /// Data fields to hold the reduced data
  std::unique_ptr<BlockData<2,T,T>> _blockDataMemory;
  /// Plane points scheduled for storage in _blockData
  /// i.e. Plane points whose physical location intersects the mother cuboid
  ///      and is nearest to a rank-local cuboid
  std::vector<std::tuple<int,int,int>> _rankLocalSubplane;
  std::map<int,int> _rankLocalSubplaneSize;
  /// Synchronization mode, see BlockDataSyncMode enum for further information.
  /// This value only matters when PARALLEL_MODE_MPI is defined.
  const BlockDataSyncMode _syncMode;

  SuperD<T,descriptors::D3<fields::PHYS_R,typename FUNCTOR::result_field>> _sliceD;
  SuperLatticePointCoupling<
    BlockReduction3D2DO<FUNCTOR>,
    meta::map<names::Lattice, descriptors::VALUED_DESCRIPTOR<T,DESCRIPTOR>,
              names::Points, descriptors::VALUED_DESCRIPTOR<T,descriptors::D3<fields::PHYS_R,typename FUNCTOR::result_field>>>
  > _couplingO;

  void updateBlockAnalytical(BlockData<2,T,T>& block);

public:
  EfficientBlockReduction3D2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                              const UnitConverter<T,DESCRIPTOR>& converter,
                              const HyperplaneLattice3D<T>& lattice,
                              BlockDataSyncMode syncMode=BlockDataSyncMode::ReduceAndBcast);

  /// Initialize rank-local list of plane points to be stored in _blockData
  void initialize();
  /// Updates and writes the data to _blockData using _rankLocalSubplane
  void update();
  /// Overload of virtual function from class BlockF2D
  BlockStructureD<2>& getBlockStructure() override;
  /// \return reference to the rank local list of discrete plane points, cuboid ids
  const std::vector<std::tuple<int,int,int>>& getRankLocalSubplane() const;

};

}

#endif
