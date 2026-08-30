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

#ifndef CORE_SUPER_LATTICE_POINT_EXTRACTION_H
#define CORE_SUPER_LATTICE_POINT_EXTRACTION_H

#include "vector.h"
#include "superLatticePointCoupling.h"

#include <tuple>
#include <map>

namespace olb {

template <typename FUNCTOR>
struct PointExtractionO {
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

/// Convenience class for extracting data associated to set of arbitrary physical points
template <typename T, typename DESCRIPTOR, typename FUNCTOR>
class SuperLatticePointExtraction final {
private:
  SuperLattice<T,DESCRIPTOR>& _sLattice;

  std::vector<std::tuple<int,std::size_t,Vector<T,DESCRIPTOR::d>>> _rankLocalPoints;
  std::map<int,int> _rankLocalPointsSize;

  SuperD<T,descriptors::D3<fields::PHYS_R,typename FUNCTOR::result_field>> _sampleD;
  SuperLatticePointCoupling<
    PointExtractionO<FUNCTOR>,
    meta::map<names::Lattice, descriptors::VALUED_DESCRIPTOR<T,DESCRIPTOR>,
              names::Points, descriptors::VALUED_DESCRIPTOR<T,descriptors::D3<fields::PHYS_R,typename FUNCTOR::result_field>>>
  > _couplingO;

  FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,fields::PHYS_R>                 _extractionPointD;
  FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,typename FUNCTOR::result_field> _extractionResultD;

  void updateExtractionResultD(
    FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,typename FUNCTOR::result_field>&);

public:
  SuperLatticePointExtraction(SuperLattice<T,DESCRIPTOR>& sLattice,
                              const std::vector<Vector<T,DESCRIPTOR::d>>& points);

  std::size_t getSize() const {
    return _extractionPointD.getSize();
  }

  Vector<T,DESCRIPTOR::d> getPhysR(std::size_t iExtracted) const {
    return _extractionPointD.getField(iExtracted);
  }

  auto getValue(std::size_t iExtracted) const {
    return _extractionResultD.getField(iExtracted);
  }

  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD>&& value) {
    _couplingO.template setParameter<FIELD>(std::move(value));
  }

  void update();

};

}

#endif
