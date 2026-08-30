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

#ifndef REFINEMENT_PLATFORM_CPU_OPERATOR_H
#define REFINEMENT_PLATFORM_CPU_OPERATOR_H

#include "refinement/operator.h"
#include "refinement/coarseCell.h"

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename OPERATOR>
requires (isPlatformCPU(PLATFORM))
struct ConcreteBlockRefinementO<T,DESCRIPTOR,PLATFORM,OPERATOR> {
  void apply(ConcreteBlockRefinementContextD<T,DESCRIPTOR,PLATFORM>& data) {
    const auto& cellIdCoarse = data.getConcreteData().template getField<fields::refinement::CELL_ID_COARSE>();
    const auto& cellIdFine = data.getConcreteData().template getField<fields::refinement::CELL_ID_FINE>();
    auto parameters = data.getConcreteData().template getData<OperatorParameters<OPERATOR>>().parameters;

           if constexpr (OPERATOR::scope == refinement::OperatorScope::PerCoarseCell) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for (CellID i=0; i < cellIdFine.getSize(); ++i) {
        if (cellIdCoarse[0][i] != 0) {
          cpu::Cell cCell{data.getCoarseLattice(), cellIdCoarse[0][i]};
          cpu::Cell fCell{data.getFineLattice(), cellIdFine[0][i]};
          cpu::Row  cData{data.getConcreteData(), i};
          OPERATOR().apply(cCell, fCell, cData, parameters);
        }
      }
    } else if constexpr (OPERATOR::scope == refinement::OperatorScope::PerFineCell) {
      if constexpr (OPERATOR::data::template contains<fields::refinement::CONTEXT_NEIGHBORS>()) {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(static)
        #endif
        for (CellID i=0; i < cellIdFine.getSize(); ++i) {
          cpu::Cell fCell{data.getFineLattice(), cellIdFine[0][i]};
          refinement::CoarseCell cCell{cpu::Cell{data.getCoarseLattice(), cellIdCoarse[0][i]},
                                       data.getFineLattice().getLatticeR(cellIdFine[0][i])};
          refinement::ContextDataWithNeighbors cData{cpu::Row{data.getConcreteData(), i}};
          OPERATOR().apply(cCell, fCell, cData, parameters);
        }
      } else {
        #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel for schedule(static)
        #endif
        for (CellID i=0; i < cellIdFine.getSize(); ++i) {
          cpu::Cell fCell{data.getFineLattice(), cellIdFine[0][i]};
          refinement::CoarseCell cCell{cpu::Cell{data.getCoarseLattice(), cellIdCoarse[0][i]},
                                       data.getFineLattice().getLatticeR(cellIdFine[0][i])};
          refinement::ContextData cData{cpu::Row{data.getConcreteData(), i}};
          OPERATOR().apply(cCell, fCell, cData, parameters);
        }
      }
    }
  }
};

}

#endif
