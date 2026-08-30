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

#ifndef FSI_SUPER_POROUS_ELEMENT_EMBEDDING_O_H
#define FSI_SUPER_POROUS_ELEMENT_EMBEDDING_O_H

#include "fields.h"
#include "operators.h"
#include "elements/concept.h"

#include "core/data.h"

namespace olb {

namespace stage {

namespace fsi {

struct Initialize { };
struct DistributeToFluid { };

}

}

template <typename T, typename DESCRIPTOR>
class SuperPorousElementEmbeddingO final {
private:
  SuperLattice<T,DESCRIPTOR>& _sLattice;

  /// Element data
  std::unique_ptr<Data<T,descriptors::fsi::ELEMENTS<DESCRIPTOR::d>>> _elementsD;

  std::map<std::type_index, std::function<void(std::size_t)>> _exposeElementTypeParametersO;

  std::size_t getElementCount() const {
    return _elementsD->template get<Array<fields::fsi::ELEMENT_TAG>>()[0].size();
  }

  void setElementCount(std::size_t newSize) {
    _elementsD->resize(newSize);
    for (auto& [_, exposeO] : _exposeElementTypeParametersO) {
      exposeO(newSize);
    }
  }

public:
  SuperPorousElementEmbeddingO(SuperLattice<T, DESCRIPTOR>& sLattice):
    _sLattice(sLattice),
    _elementsD(makeSharedData<
      T,descriptors::fsi::ELEMENTS<DESCRIPTOR::d>>(sLattice.getLoadBalancer(), 0))
  {
    OstreamManager clout(std::cout, "SuperPorousElementEmbeddingO");

    //clout << "Expose physical cell locations" << std::endl;

    auto& load = _sLattice.getLoadBalancer();
    for (int iC=0; iC < load.size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      block.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();
      block.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        auto latticeR = loc.withPrefix(load.glob(iC));
        auto physR = _sLattice.getCuboidDecomposition().getPhysR(latticeR);
        block.get(loc).template setField<descriptors::LOCATION>(physR);
      });
    }
    _sLattice.template setProcessingContext<Array<descriptors::LOCATION>>(ProcessingContext::Simulation);

    _exposeElementTypeParametersO[typeid(SuperPorousElementEmbeddingO)] = [this](std::size_t nElements){
      descriptors::fsi::ELEMENTS<DESCRIPTOR::d>::fields_t::for_each([this,nElements](auto field) {
        using field_t = typename decltype(field)::type;
        auto& fieldArrayD = _elementsD->template get<Array<field_t>>();
        fieldArrayD.resize(nElements);
        _sLattice.template setParameter<fields::array_of<field_t>>(fieldArrayD);
      });
    };
  }

  template <concepts::PorosityElementF PorosityF>
  void registerElementType(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicatorF) {
    PorosityF::data::for_each([&](auto field) {
      using field_t = typename decltype(field)::type;
      if (!_elementsD->template provides<Array<field_t>>()) {
        _elementsD->template allocate<Array<field_t>>(getElementCount());
      }
    });

    _sLattice.template addPostProcessor<stage::fsi::Initialize>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<InitializePorosityO<PorosityF>>{});
    _sLattice.template addPostProcessor<stage::fsi::DistributeToFluid>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<UpdatePorosityO<PorosityF>>{});

    _exposeElementTypeParametersO[typeid(PorosityF)] = [this](std::size_t nElements){
      PorosityF::data::for_each([this,nElements](auto field) {
        using field_t = typename decltype(field)::type;
        auto& fieldArrayD = _elementsD->template get<Array<field_t>>();
        fieldArrayD.resize(nElements);
        _sLattice.template setParameter<fields::array_of<field_t>>(fieldArrayD);
      });
    };
  }

  template <concepts::Parameters PARAMETERS>
  void add(PARAMETERS& element) {
    std::size_t iElement = getElementCount();
    setElementCount(iElement+1);
    const std::size_t nElements = getElementCount();

    PARAMETERS::fields_t::for_each([&](auto field) {
      using field_t = typename decltype(field)::type;
      if (_elementsD->template provides<Array<field_t>>()) {
        _elementsD->template get<Array<field_t>>().set(
          iElement, element.template get<field_t>());
      }
    });

    if constexpr (PARAMETERS::fields_t::template contains<fields::fsi::ELEMENT_LOWER>()) {
      if constexpr (PARAMETERS::fields_t::template contains<fields::fsi::ELEMENT_UPPER>()) {
        _elementsD->template get<Array<fields::fsi::ELEMENT_UPPER>>().set(
          iElement, element.template get<fields::fsi::ELEMENT_UPPER>());
      } else {
        auto upper = element.template get<fields::fsi::ELEMENT_LOWER>()
                   + element.template get<fields::fsi::ELEMENT_REFERENCE_DELTA_X>()
                   * element.template get<fields::fsi::ELEMENT_REFERENCE_EXTENT>();
        _elementsD->template get<Array<fields::fsi::ELEMENT_UPPER>>().set(iElement, upper);
      }
    }

    if constexpr (PARAMETERS::fields_t::template contains<fields::fsi::ELEMENT_REFERENCE_DELTA_X>()) {
      _elementsD->template get<Array<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>().set(
        iElement, element.template get<fields::fsi::ELEMENT_REFERENCE_DELTA_X>());
      _elementsD->template get<Array<fields::fsi::ELEMENT_REFERENCE_POROSITY>>().set(
        iElement, element.template get<fields::fsi::ELEMENT_REFERENCE_POROSITY>());
      auto extent = element.template get<fields::fsi::ELEMENT_REFERENCE_EXTENT>();
      FieldD<T,DESCRIPTOR,fields::fsi::ELEMENT_REFERENCE_PROJECTION> projection{};
      if constexpr (DESCRIPTOR::d == 3) {
        projection = {extent[1]*extent[2], extent[2], 1};
      } else {
        projection = {extent[1], 1};
      }
      _elementsD->template get<Array<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>().set(iElement, projection);
    }

    _elementsD->setProcessingContext(ProcessingContext::Simulation);
    _sLattice.template setParameter<fields::fsi::ELEMENTS_COUNT>(nElements);
  }

  void initialize() {
    setElementCount(getElementCount());
    _elementsD->setProcessingContext(ProcessingContext::Simulation);
    _sLattice.executePostProcessors(stage::fsi::Initialize{});
  }

  template <typename FIELD>
  FieldD<T,DESCRIPTOR,FIELD> getField(std::size_t iElement) const {
    return _elementsD->template get<Array<FIELD>>().get(iElement);
  }

  template <typename FIELD>
  void setField(std::size_t iElement, FieldD<T,DESCRIPTOR,FIELD> data) {
    _elementsD->template get<Array<FIELD>>().set(iElement, data);
  }

  template <typename FIELD_TYPE>
  void setProcessingContext(ProcessingContext context) {
    _elementsD->template get<FIELD_TYPE>().setProcessingContext(context);
  }

  void apply() {
    _sLattice.executePostProcessors(stage::fsi::DistributeToFluid{});
  }

};

}

#endif
