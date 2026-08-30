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

#ifndef FSI_OPERATORS_H
#define FSI_OPERATORS_H

#include <ranges>

#include "fields.h"
#include "elements/concept.h"

#include "elements/couplingFaceF2D.h"

namespace olb {

/// Operator for discretizing FSI elements into HLBM porosities
template <concepts::PorosityElementF PorosityF>
struct InitializePorosityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::merge<
    typename PorosityF::parameters,
    meta::list<
      fields::array_of<fields::fsi::ELEMENT_TAG>,
      fields::fsi::ELEMENTS_COUNT
    >
  >;

  int getPriority() const {
    return 0;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    const auto physR = cell.template getField<descriptors::LOCATION>();
    const auto nElements = parameters.template get<fields::fsi::ELEMENTS_COUNT>();
    if (cell.template getField<descriptors::POROSITY>() == 1) {
      for (std::size_t iElement=0; iElement < nElements; ++iElement) {
        const auto porosity = PorosityF().compute(parameters, physR, iElement);
        const auto prevPorosity = cell.template getField<descriptors::POROSITY>();
        if (porosity < prevPorosity) {
          const auto newTag = PorosityF().tag(parameters, physR, iElement);
          cell.template setField<fields::fsi::ELEMENT_TAG>(newTag);
          cell.template setField<descriptors::POROSITY>(porosity);
          if constexpr (concepts::PorosityElementFcomputesU<PorosityF>) {
            const auto u = PorosityF().computeU(parameters, physR, iElement);
            cell.template setField<descriptors::VELOCITY>(u);
          }
          if constexpr (concepts::PorosityElementFcomputesY1<PorosityF>) {
            const auto y1 = PorosityF().computeY1(parameters, physR, iElement);
            cell.template setField<descriptors::Y1>(y1);
          }
        }
      }
      if (cell.template getField<descriptors::POROSITY>() == 1) {
        cell.template setField<fields::fsi::ELEMENT_TAG>(0);
        cell.template setField<descriptors::VELOCITY>(0);
      }
    }
  }
};

/// Operator for updating FSI element discretizations into HLBM porosities
template <concepts::PorosityElementF PorosityF>
struct UpdatePorosityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::merge<
    typename PorosityF::parameters,
    meta::list<fields::array_of<fields::fsi::ELEMENT_TAG>>
  >;

  int getPriority() const {
    return 0;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    if (tag > 0) {
      const int iElement = tag - 1;
      const auto physR = cell.template getField<descriptors::LOCATION>();

      const auto newTag      = PorosityF().tag(parameters, physR, iElement);
      const auto newPorosity = PorosityF().compute(parameters, physR, iElement);

      const bool isPorous = newPorosity < 1;
      cell.template setField<fields::fsi::ELEMENT_TAG>(isPorous * newTag);

      if constexpr (concepts::PorosityElementFcomputesU<PorosityF>) {
        const auto newU = PorosityF().computeU(parameters, physR, iElement);
        cell.template setField<descriptors::VELOCITY>(isPorous * newU);
      }
      if constexpr (concepts::PorosityElementFcomputesY1<PorosityF>) {
        const auto newY1 = PorosityF().computeY1(parameters, physR, iElement);
        cell.template setField<descriptors::Y1>(isPorous * newY1);
      }

      const bool isFrontier = !PorosityF().isInterior(parameters, physR, iElement);
      cell.template setField<descriptors::POROSITY>(isFrontier * newPorosity);
    }
  }
};

/// Operator for calculating per-cell momentum exchange torque in a HLBM context
struct CollectPorousBoundaryTorqueO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::converter::PHYS_LENGTH
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> force { };
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> torque { };
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && porosity < 1) {
      Vector<V,DESCRIPTOR::d> u{};
      cell.computeU(u.data());

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);

        const bool isNeighborFluid      = neighbor.template getField<descriptors::POROSITY>() == 1;

        const V pop1 = neighbor[iPop];
        const V pop2 = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        force -= isNeighborFluid * (pop1*(c-u) + pop2*(c+u));

        isAtSurface |= isNeighborFluid;
      }

      auto physR = cell.template getField<descriptors::LOCATION>();
      auto pivots = parameters.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
      const Vector<V,DESCRIPTOR::d> pivot{pivots, static_cast<std::size_t>(tag-1)};
      torque = crossProduct((physR - pivot) * parameters.template get<fields::converter::PHYS_LENGTH>(),
                            force);

      if (isAtSurface) {
        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        torque = 0;
      }
    }

    cell.template setField<fields::fsi::ELEMENT_TORQUE>(torque);
  }

};

/// Operator for calculating per-cell momentum exchange forces in a HLBM context
struct CollectPorousBoundaryForceAndTorqueO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::converter::PHYS_LENGTH
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> force { };
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> torque { };
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && porosity < 1) {
      Vector<V,DESCRIPTOR::d> u{};
      cell.computeU(u.data());

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);

        const bool isNeighborFluid      = neighbor.template getField<descriptors::POROSITY>() == 1;
        //const bool isIncreasingPorosity = neighbor.template getField<descriptors::POROSITY>() >  porosity;

        const V pop1 = neighbor[iPop];
        const V pop2 = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        force -= isNeighborFluid * (pop1*(c-u) + pop2*(c+u));

        isAtSurface |= isNeighborFluid;
      }

      auto physR = cell.template getField<descriptors::LOCATION>();
      auto pivots = parameters.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
      const Vector<V,DESCRIPTOR::d> pivot{pivots, static_cast<std::size_t>(tag-1)};
      torque = crossProduct((physR - pivot) * parameters.template get<fields::converter::PHYS_LENGTH>(),
                            force);

      if (isAtSurface) {
        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        force = 0;
        torque = 0;
      }
    }

    cell.template setField<fields::fsi::ELEMENT_FORCE>(force);
    cell.template setField<fields::fsi::ELEMENT_TORQUE>(torque);
  }

};

struct CollectPorousBoundaryForceO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::converter::PHYS_LENGTH
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> force { };
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && porosity < 1) {
      Vector<V,DESCRIPTOR::d> u{};
      cell.computeU(u.data());

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);

        const bool isNeighborFluid = neighbor.template getField<descriptors::POROSITY>() == 1;

        const V pop1 = neighbor[iPop];
        const V pop2 = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        force -= isNeighborFluid * (pop1*(c-u) + pop2*(c+u));

        isAtSurface |= isNeighborFluid;
      }

      if (isAtSurface) {
        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        force = 0;
      }
    }

    cell.template setField<fields::fsi::ELEMENT_FORCE>(force);
  }

};

/// Operator for calculating per-cell stress in a HLBM context
struct CollectPorousBoundaryStressO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<CouplingFaceF2D::NEIGHBORS>,
    fields::converter::PHYS_LENGTH,
    fields::converter::LATTICE_VISCOSITY
  >;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_STRESS> stress { };
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);

    bool isAtSurface = false;
    if (tag > 0 && cell.template getField<descriptors::POROSITY>() < 1) {
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);
        const bool isNeighborFluid = neighbor.template getField<descriptors::POROSITY>() == 1;
        isAtSurface |= isNeighborFluid;
      }

      if (isAtSurface) {
        V rho{};
        Vector<V,DESCRIPTOR::d> u{};
        Vector<V,util::TensorVal<DESCRIPTOR>::n> strainRate{};
        cell.computeAllMomenta(rho, u.data(), strainRate.data());

        const V viscosity = parameters.template get<fields::converter::LATTICE_VISCOSITY>();
        const V p = (rho - V{1}) / descriptors::invCs2<V,DESCRIPTOR>();

        const auto neighbors = parameters.template get<fields::array_of<CouplingFaceF2D::NEIGHBORS>>();
        const auto pivots = parameters.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();

        const int prevElement = neighbors[0][tag-1];
        const int nextElement = neighbors[1][tag-1];

        const Vector<V,DESCRIPTOR::d> prevR{pivots, static_cast<unsigned>(prevElement)};
        const Vector<V,DESCRIPTOR::d> nextR{pivots, static_cast<unsigned>(nextElement)};

        auto tangent = prevR-nextR;
        tangent /= norm(tangent);
        Vector normal{tangent[1], -tangent[0]};

        stress[0] = -p*normal[0];// + 2*rho*viscosity*(normal[0]*strainRate[0]+normal[1]*strainRate[1]);
        stress[1] = -p*normal[1];// + 2*rho*viscosity*(normal[0]*strainRate[1]+normal[1]*strainRate[2]);

        cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(tag);

        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      } else {
        stress = 0;
      }
    }

    cell.template setField<fields::fsi::ELEMENT_STRESS>(stress);
  }

};

/// Operator for allowing migration of padding layer
struct GrowPaddingLayerO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <concepts::Cell CELL>
  void apply(CELL& cell) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    cell.template setField<fields::fsi::REDUCED_ELEMENT_TAG>(0);
    const int tag = cell.template getField<fields::fsi::ELEMENT_TAG>();
    if (tag > 0 && cell.template getField<descriptors::POROSITY>() < 1) {
      bool isAtSurface = false;

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto c = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);
        const bool isNeighborFluid = neighbor.template getField<descriptors::POROSITY>() == 1;
        isAtSurface |= isNeighborFluid;
      }

      if (isAtSurface) {
        for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop);
          auto neighbor = cell.neighbor(c);
          if (neighbor.template getField<fields::fsi::ELEMENT_TAG>() == 0) {
            neighbor.template setField<fields::fsi::ELEMENT_TAG>(tag);
          }
        }
      }
    }
  }
};

}

#endif
