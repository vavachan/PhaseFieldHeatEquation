/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Yuji (Sam) Shimojima, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_FIELD_REDUCTION_O_H
#define SUPER_LATTICE_FIELD_REDUCTION_O_H

#include "core/superD.h"

#include "blockLatticeFieldReductionO.h"

namespace olb {

namespace stage::reduction {

struct ReductionField {};

}

template <typename T, typename DESCRIPTOR, typename FIELD, typename REDUCTION_OP,
          typename CONDITION = reduction::ConditionTrue>
class SuperLatticeFieldReductionO final {
public:
  template <unsigned D>
  struct REDUCED_FIELD : public descriptors::SPATIAL_DESCRIPTOR<D, FIELD> {};

private:
  SuperLattice<T, DESCRIPTOR>&                             _sLattice;
  std::unique_ptr<SuperD<T, REDUCED_FIELD<DESCRIPTOR::d>>> _superReducedFieldD;
  static constexpr unsigned                                _fieldDimension = FieldD<T, DESCRIPTOR, FIELD>::d;
  FieldD<T, DESCRIPTOR, FIELD>                             _reductedField {};

  bool _rankDoesReduction;

  void _gatherField(const FieldD<T, DESCRIPTOR, FIELD>& f)
  {
    for (unsigned iD = 0; iD < f.d; iD++) {
      _reductedField[iD] = REDUCTION_OP {}(_reductedField[iD], f[iD]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;

  void _mpiGatherField(MPI_Op op)
  {
    std::vector<T> localFields(_fieldDimension, T {});
    for (std::size_t iD = 0; iD < _fieldDimension; ++iD) {
      localFields[iD] = _reductedField[iD];
    }

    std::vector<T> globalFields(_fieldDimension, T {});
    singleton::mpi().allreduce(localFields.data(), globalFields.data(), globalFields.size(), op, _mpiCommunicator);

    for (std::size_t iD = 0; iD < _fieldDimension; ++iD) {
      _reductedField[iD] = globalFields[iD];
    }
  }

#endif // PARALLEL_MODE_MPI

public:
  ~SuperLatticeFieldReductionO() {};

  SuperLatticeFieldReductionO(SuperLattice<T, DESCRIPTOR>& sLattice)
      : _sLattice(sLattice)
      , _superReducedFieldD(new SuperD<T, REDUCED_FIELD<DESCRIPTOR::d>>(sLattice.getLoadBalancer()))
      , _rankDoesReduction {singleton::mpi().isMainProcessor()}

  {
    OstreamManager clout(std::cout, "SuperLatticeFieldReductionO");
    auto&          load = _sLattice.getLoadBalancer();

    {
      for (int iC = 0; iC < load.size(); ++iC) {
        auto& block = _sLattice.getBlock(iC);
        block.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
          block.get(loc).template setField<field::reduction::TAG_CORE>(
              (int)1); //this tag exracts core region. The main purpose is for GPU.
        });
        block.setProcessingContext(ProcessingContext::Simulation);
      }
    }

    clout << "Set up operators and communication" << std::endl;

    sLattice.template addPostProcessor<stage::reduction::ReductionField>(
        meta::id<BlockLatticeFieldReductionO<FIELD, REDUCTION_OP, CONDITION>> {});
    {
      auto& c = _sLattice.getCommunicator(stage::reduction::ReductionField {});

      c.template requestField<FIELD>();
      c.requestOverlap(1);
      c.exchangeRequests();
    }

#ifdef PARALLEL_MODE_MPI
    for (int iC = 0; iC < load.size(); ++iC) {

      _rankDoesReduction |= sLattice.getBlock(iC).hasPostProcessor(
          typeid(stage::reduction::ReductionField),
          meta::id<BlockLatticeFieldReductionO<FIELD, REDUCTION_OP, CONDITION>> {});
    }
    MPI_Comm_split(MPI_COMM_WORLD, _rankDoesReduction ? 0 : MPI_UNDEFINED, singleton::mpi().getRank(),
                   &_mpiCommunicator);
#endif

    clout << "Set operator parameters" << std::endl;
    for (int iC = 0; iC < load.size(); ++iC) {
      auto& block         = _sLattice.getBlock(iC);
      auto& elementsBlock = _superReducedFieldD->getBlock(iC);
      elementsBlock.resize(1);
      block.template setParameter<fields::array_of<FIELD>>(elementsBlock.template getField<FIELD>());
    }

    _sLattice.getCommunicator(stage::reduction::ReductionField {}).communicate();
    _superReducedFieldD->setProcessingContext(ProcessingContext::Simulation);
  }

  bool rankDoesReduction() const { return _rankDoesReduction; }

  FieldD<T, DESCRIPTOR, FIELD> compute()
  {
    OstreamManager clout(std::cout, "SuperLatticeFieldReductionO");
    _sLattice.executePostProcessors(stage::reduction::ReductionField {});
    for (unsigned iD = 0; iD < _fieldDimension; iD++) {
      _reductedField[iD] = REDUCTION_OP{}.reset(_reductedField[iD]);
    }

    for (int iC = 0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block      = _superReducedFieldD->getBlock(iC);
      auto& blockField = block.template getField<FIELD>();
      block.setProcessingContext(ProcessingContext::Evaluation);

      FieldD<T, DESCRIPTOR, FIELD> f([&](std::size_t iD) {
        return blockField[iD][0];
      });
      _gatherField(f);
    }

#ifdef PARALLEL_MODE_MPI
    _mpiGatherField(REDUCTION_OP{}.forMPI());
#endif
    return _reductedField;
  }
};

} // namespace olb

#endif
