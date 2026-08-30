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

#ifndef FSI_SUPER_POROUS_ELEMENT_REDUCTION_O_H
#define FSI_SUPER_POROUS_ELEMENT_REDUCTION_O_H

#include "fields.h"
#include "operators.h"
#include "integral.h"

#include "core/superD.h"

namespace olb {

namespace stage {

namespace fsi {

struct CollectFromFluid { };

}

}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
class SuperPorousElementReductionO final {
private:
  /// Fluid lattice
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  /// FSI region indicator
  FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>> _indicatorF;
  /// Reduced element data obtained from the fluid lattice
  std::unique_ptr<SuperD<T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>> _superReducedElementsD;

  /// Reduced per-element fields
  std::map<unsigned, ParametersD<T,DESCRIPTOR,FIELDS...>> _fields;

  std::size_t _nElements;

  bool _rankDoesFSI;

  #ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
  #endif

public:
  SuperPorousElementReductionO(SuperLattice<T,DESCRIPTOR>& sLattice,
                               FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicatorF):
    _sLattice(sLattice),
    _indicatorF{std::move(indicatorF)},
    _superReducedElementsD(new SuperD<
      T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>(sLattice.getLoadBalancer())),
    _rankDoesFSI{singleton::mpi().isMainProcessor()}
  {
    OstreamManager clout(std::cout, "SuperPorousElementReductionO");
    auto& load = _sLattice.getLoadBalancer();

    // Schedule boundary force collection for FSI domain + overlap for block migration
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = _indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
            if (block.isPadding(latticeR, 2)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     meta::id<GrowPaddingLayerO>{});
            }
          }
        });
      }
    }
    sLattice.template addPostProcessor<stage::fsi::CollectFromFluid>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<IntegratePorousElementFieldsO<FIELDS...>>{});

    {
      auto& c = sLattice.getCommunicator(stage::fsi::CollectFromFluid{});
      c.template requestField<descriptors::POPULATION>();
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    #ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC < load.size(); ++iC) {
      _rankDoesFSI |= sLattice.getBlock(iC).hasPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                                             meta::id<IntegratePorousElementFieldsO<FIELDS...>>{});
    }
    MPI_Comm_split(MPI_COMM_WORLD, _rankDoesFSI ? 0 : MPI_UNDEFINED, singleton::mpi().getRank(), &_mpiCommunicator);
    #endif

    {
      auto& c = sLattice.getCommunicator(stage::Evaluation{});
      c.template requestField<descriptors::POROSITY>();
      c.requestOverlap(1, std::forward<decltype(indicatorF)>(indicatorF));
      c.exchangeRequests();
    }

    {
      auto& c = sLattice.getCommunicator(stage::Full{});
      c.template requestField<descriptors::POROSITY>();
      c.template requestField<fields::fsi::ELEMENT_TAG>();
      c.exchangeRequests();
    }

    for (int iC=0; iC < load.size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      (block.template setParameter<fields::array_of<FIELDS>>(
        elementsBlock.template getField<FIELDS>()), ...);
    }
  }

  /// Resize reduction target buffer must be >= total global number of elements
  void resize(std::size_t nElements) {
    _nElements = nElements;
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      elementsBlock.template getField<fields::fsi::ELEMENT_TAG>().resize(nElements);
      (elementsBlock.template getField<FIELDS>().resize(nElements), ...);

      auto& block = _sLattice.getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      (block.template setParameter<fields::array_of<FIELDS>>(
        elementsBlock.template getField<FIELDS>()), ...);
    }
    _superReducedElementsD->setProcessingContext(ProcessingContext::Simulation);
  }

  void addCollectionO(PostProcessorPromise<T,DESCRIPTOR>&& collectionO) {
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = _indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
            if (block.isInsideCore(latticeR)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     std::forward<decltype(collectionO)>(collectionO));
            }
          }
        });
      }
    }

    _sLattice.template addPostProcessor<stage::Evaluation>(
      std::forward<decltype(_indicatorF)>(_indicatorF),
      std::forward<decltype(collectionO)>(collectionO));
  }

  bool rankDoesFSI() const
  {
    return _rankDoesFSI;
  }

  /// Globally integrate all element surface forces
  void apply()
  {
    OstreamManager clout(std::cout, "SuperPorousBoundaryFieldReductionO");

    _sLattice.executePostProcessors(stage::fsi::CollectFromFluid{});

    if (!_rankDoesFSI) {
      return;
    }

    _fields.clear();

    /// Gather local element fields
    int maxElement{};
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      auto& blockTags = block.template getField<fields::fsi::ELEMENT_TAG>();
      blockTags.setProcessingContext(ProcessingContext::Evaluation);

      ((block.template getField<FIELDS>().setProcessingContext(ProcessingContext::Evaluation)), ...);

      auto set = [&]<typename FIELD>(meta::id<FIELD>, std::size_t iElement, int tag) {
        auto& fieldArray = block.template getField<FIELD>();
        FieldD<T,DESCRIPTOR,FIELD> f([&](std::size_t iD) {
          return fieldArray[iD][iElement];
        });
        _fields[tag].template set<FIELD>(
          _fields[tag].template get<FIELD>() + f);
      };

      auto& blockLattice = _sLattice.getBlock(iC);
      const std::size_t nElements = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO<FIELDS...>>>()
                                                .template get<fields::fsi::REDUCED_ELEMENTS_COUNT>();
      for (std::size_t iElement=0; iElement < nElements; ++iElement) {
        const int tag = blockTags[0][iElement];
        if (tag > maxElement) {
          maxElement = tag;
        }
        if (tag > 0) {
          ((set(meta::id<FIELDS>{}, iElement, tag)), ...);
        }
      }
    }

    #ifdef PARALLEL_MODE_MPI
    int nTotalElements{};
    singleton::mpi().allreduce(&maxElement, &nTotalElements, 1, MPI_MAX, _mpiCommunicator);
    if (nTotalElements >= 0) {
      _nElements = nTotalElements;
    }

    /// Aggregate global element fields
    auto communicate = [&]<typename FIELD>(meta::id<FIELD>) {
      const unsigned fieldSize = DESCRIPTOR::template size<FIELD>();
      std::vector<T> localValues(fieldSize*(_nElements+1), T{});
      for (auto [tag, fields] : _fields) {
        FieldD<T,DESCRIPTOR,FIELD> value = fields.template get<FIELD>();
        for (std::size_t iD=0; iD < fieldSize; ++iD) {
          localValues[fieldSize*tag+iD] = value[iD];
        }
      }
      std::vector<T> globalValues(fieldSize*(_nElements+1), T{});
      singleton::mpi().allreduce(localValues.data(), globalValues.data(), globalValues.size(), MPI_SUM, _mpiCommunicator);
      for (std::size_t iElement=1; iElement <= _nElements; ++iElement) {
        _fields[iElement].template set<FIELD>(globalValues.data() + fieldSize*iElement);
      }
    };

    (communicate(meta::id<FIELDS>{}), ...);
    #endif
  }

  std::size_t getElementCount() const {
    return _nElements;
  }

  template <concepts::Field FIELD>
  auto getField(unsigned iElement) const {
    return _fields.at(iElement).template get<FIELD>();
  }

};

}

#endif
