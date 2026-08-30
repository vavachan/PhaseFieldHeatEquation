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

#ifndef EFFICIENT_BLOCK_REDUCTION_3D2D_HH
#define EFFICIENT_BLOCK_REDUCTION_3D2D_HH

#include "efficientBlockReduction3D2D.h"

#include <limits>
#include "utilities/omath.h"

#include "utilities/vectorHelpers.h"
#include "communication/mpiManager.h"

namespace olb {


template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::updateBlockAnalytical(BlockData<2,T,T>& block)
{
  LoadBalancer<T>& load = _sLattice.getLoadBalancer();
  std::map<int,int> rankLocalSubplaneIndex;
  for (auto [iX, iY, iC] : _rankLocalSubplane) {
    auto& sliceBlock = _sliceD.getBlock(load.loc(iC));
    auto iPoint = rankLocalSubplaneIndex[iC]++;
    auto point = sliceBlock.get(iPoint);

    auto result = point.template getField<typename FUNCTOR::result_field>();
    for (int iD=0; iD < result.d; ++iD) {
      block.get({iX, iY}, iD) = result[iD];
    }
  }
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::EfficientBlockReduction3D2D(
  SuperLattice<T,DESCRIPTOR>&        sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  const HyperplaneLattice3D<T>&      lattice,
  BlockDataSyncMode                  syncMode)
  : HyperplaneLattice3D<T>(lattice),
    BlockDataF2D<T,T>(lattice.getNx(), lattice.getNy(),
                      descriptors::D3<>::template size<typename FUNCTOR::result_field>()),
    _sLattice(sLattice),
    _syncMode(syncMode),
    _sliceD(sLattice.getLoadBalancer()),
    _couplingO(BlockReduction3D2DO<FUNCTOR>{},
               names::Lattice{}, sLattice,
               names::Points{}, _sliceD)
{
  this->getName() = "planeReduction()";

  // intialize list of relevant rank local points making up the reduced plane
  initialize();

  _couplingO.template setParameter<fields::converter::PHYS_VELOCITY>(converter.getConversionFactorVelocity());
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::initialize()
{
  const auto& geometry = _sLattice.getCuboidDecomposition();
  LoadBalancer<T>&           load     = _sLattice.getLoadBalancer();

  for (int iC=0; iC < geometry.size(); ++iC) {
    _rankLocalSubplaneSize[iC] = 0;
  }

  for ( int iX = 0; iX < this->getNx(); ++iX ) {
    for ( int iY = 0; iY < this->getNy(); ++iY ) {
      const Vector<T,3> physR = this->getPhysR(iX, iY);

      // Schedule plane point for storage if its physical position intersects the
      // mother cuboid and the cuboid of the nearest lattice position is local to
      // the current rank:
      if (auto iC = geometry.getC(physR)) {
        if (load.isLocal(*iC)) {
          _rankLocalSubplane.emplace_back(iX, iY, *iC);
          _rankLocalSubplaneSize[*iC]++;
        }
      }
    }
  }

  std::map<int,int> rankLocalSubplaneIndex;
  for (auto [iC, size] : _rankLocalSubplaneSize) {
    if (load.isLocal(iC)) {
      auto& block = _sliceD.getBlock(load.loc(iC));
      block.resize({size,1,1});
      rankLocalSubplaneIndex[iC] = 0;
      std::cout << "Local subplane " << iC << " on rank " << singleton::mpi().getRank() << " is " << size << std::endl;
    }
  }

  for (auto [iX, iY, iC] : _rankLocalSubplane) {
    auto& block = _sliceD.getBlock(load.loc(iC));
    auto iPoint = rankLocalSubplaneIndex[iC]++;
    auto point = block.get(iPoint);
    const Vector<T,3> physR = this->getPhysR(iX, iY);
    point.template setField<fields::PHYS_R>(physR);
  }

  _sliceD.setProcessingContext(ProcessingContext::Simulation);
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::update()
{
  _couplingO.execute();
  //_sliceD.template setProcessingContext<Array<typename FUNCTOR::result_field>>(ProcessingContext::Evaluation);
  _sliceD.setProcessingContext(ProcessingContext::Evaluation);

#ifdef PARALLEL_MODE_MPI
  std::unique_ptr<BlockData<2,T,T>> localBlockData(
    new BlockData<2,T,T>({{this->getNx(), this->getNy()}, 0},
                         descriptors::D3<>::template size<typename FUNCTOR::result_field>()));

  updateBlockAnalytical(*localBlockData);

  switch ( _syncMode ) {
  case BlockDataSyncMode::ReduceAndBcast:
    singleton::mpi().reduce(*localBlockData, this->getBlockData(), MPI_SUM);
    singleton::mpi().bCast(this->getBlockData());
    break;
  case BlockDataSyncMode::ReduceOnly:
    singleton::mpi().reduce(*localBlockData, this->getBlockData(), MPI_SUM);
    break;
  case BlockDataSyncMode::None:
    if (this->_owning) {
      delete this->_blockData;
    }
    this->_blockData = localBlockData.release();
    this->_owning = true;
    break;
  }
#else
  updateBlockAnalytical(this->getBlockData());
#endif
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
BlockStructureD<2>& EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::getBlockStructure()
{
  return *this->_blockData;
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
const std::vector<std::tuple<int,int,int>>& EfficientBlockReduction3D2D<T,DESCRIPTOR,FUNCTOR>::getRankLocalSubplane() const
{
  return this->_rankLocalSubplane;
}


} // end namespace olb

#endif
