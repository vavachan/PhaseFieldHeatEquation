/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CORE_SUPER_LATTICE_POINT_EXTRACTION_HH
#define CORE_SUPER_LATTICE_POINT_EXTRACTION_HH

#include "superLatticePointExtraction.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::updateExtractionResultD(
  FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,typename FUNCTOR::result_field>& resultD)
{
  LoadBalancer<T>& load = _sLattice.getLoadBalancer();
  std::map<int,int> rankLocalPointsIndex;
  for (const auto& [iC, iExtracted, physR] : _rankLocalPoints) {
    auto& block = _sampleD.getBlock(load.loc(iC));
    auto iPoint = rankLocalPointsIndex[iC]++;
    auto point = block.get(iPoint);

    resultD.setField(iExtracted, point.template getField<typename FUNCTOR::result_field>());
  }
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::SuperLatticePointExtraction(
  SuperLattice<T,DESCRIPTOR>&                sLattice,
  const std::vector<Vector<T,DESCRIPTOR::d>>& points)
  : _sLattice(sLattice),
    _sampleD(sLattice.getLoadBalancer()),
    _couplingO(PointExtractionO<FUNCTOR>{},
               names::Lattice{}, sLattice,
               names::Points{}, _sampleD),
    _extractionPointD(points.size()),
    _extractionResultD(points.size())
{
  const auto& geometry  = _sLattice.getCuboidDecomposition();
  LoadBalancer<T>& load = _sLattice.getLoadBalancer();

  for (int iC=0; iC < geometry.size(); ++iC) {
    _rankLocalPointsSize[iC] = 0;
  }

  for (std::size_t iG=0; iG < points.size(); ++iG) {
    const auto& physR = points[iG];
    _extractionPointD.setField(iG, physR);
    if (auto iC = geometry.getC(physR)) {
      if (load.isLocal(*iC)) {
        _rankLocalPoints.emplace_back(*iC, iG, physR);
        _rankLocalPointsSize[*iC]++;
      }
    }
  }

  std::map<int,int> rankLocalPointsIndex;
  for (auto [iC, size] : _rankLocalPointsSize) {
    if (load.isLocal(iC)) {
      auto& block = _sampleD.getBlock(load.loc(iC));
      block.resize({size,1,1});
      rankLocalPointsIndex[iC] = 0;
    }
  }

  for (auto& [iC, iG, physR] : _rankLocalPoints) {
    auto& block = _sampleD.getBlock(load.loc(iC));
    auto iPoint = rankLocalPointsIndex[iC]++;
    auto point = block.get(iPoint);
    point.template setField<fields::PHYS_R>(physR);
  }

  _sampleD.setProcessingContext(ProcessingContext::Simulation);
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::update()
{
  _couplingO.execute();
  _sampleD.setProcessingContext(ProcessingContext::Evaluation);

#ifdef PARALLEL_MODE_MPI
  FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,typename FUNCTOR::result_field> localResultD(
    _extractionResultD.getSize());

  updateExtractionResultD(localResultD);

  for (unsigned iD=0; iD < DESCRIPTOR::template size<typename FUNCTOR::result_field>(); ++iD) {
    singleton::mpi().reduce(localResultD[iD].data(),
                            _extractionResultD[iD].data(),
                            _extractionResultD[iD].size(),
                            MPI_SUM);
  }
#else
  updateExtractionResultD(_extractionResultD);
#endif
}

}

#endif
