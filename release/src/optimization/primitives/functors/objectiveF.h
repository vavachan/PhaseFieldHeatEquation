/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef OBJECTIVE_F_H
#define OBJECTIVE_F_H

#include "../concept.h"
#include "core/matrixView.h"

namespace olb {

namespace functors {

// Computes "j = 0.5 * (phi - phi_ref)^2 / normalize", used for inverse problems
// FUNCTOR specifies how "phi" is computed in physical units and "phi_ref" is
// provided via fields (could be simulation or external data).
template <typename FUNCTOR>
struct L2DistanceF {
  using parameters = typename FUNCTOR::parameters::template include<descriptors::NORMALIZE>;

  using result_t = opti::J;

  struct Reference : public FUNCTOR::result_t { };

  // TODO: this is only required for the DualO which computes the derivative,
  // Thus make a normal operator concept and a differentiable operator concept, which additionally
  // requires differentiable_field_t to be exposed.
  using fields_t = typename FUNCTOR::fields_t::template include<result_t,
                                                                L2DistanceF<FUNCTOR>::Reference>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    V normalize = parameters.template get<descriptors::NORMALIZE>();

    auto phi = FUNCTOR{}.compute(cell, parameters);
    auto phi_ref = cell.template getField<Reference>();

    FieldD<V,DESCRIPTOR,result_t> j;
    for (int iDim=0; iDim < phi.getDim(); ++iDim) {
      j[0] += util::pow(phi[iDim] - phi_ref[iDim], 2);
    }
    j[0] *= 0.5 / normalize;

    return j;
  }
};

// Compute the derivative regarding the populations of any passed functor. This
// is used to compute automatically the derivatives of the objective functional.
template <typename FUNCTOR, typename DYNAMICS>
struct DualF {
  using parameters = typename FUNCTOR::parameters::template include<descriptors::DX>;

  using result_t = opti::DJDF;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using COMBINED_FIELDS = meta::merge<typename CELL::descriptor_t::fields_t, typename FUNCTOR::fields_t>;
    using DESCRIPTOR = typename COMBINED_FIELDS::template decompose_into<CELL::descriptor_t::template extend_by_fields>;
    using ADf = typename util::ADf<V,DESCRIPTOR::q>;
    using ADf_DYNAMICS = typename DYNAMICS::template exchange_value_type<ADf>;
    const V dx = parameters.template get<descriptors::DX>();

    FullCellD<ADf,DESCRIPTOR,ADf_DYNAMICS> adfCell;
    auto adfParams = parameters.template copyAs<ADf>();
    DESCRIPTOR::fields_t::for_each([&](auto id) {
      using FIELD = typename decltype(id)::type;
      adfCell.template setField<FIELD>(FieldD<ADf,DESCRIPTOR,FIELD>(cell.template getField<FIELD>()));
    });

    for (std::size_t iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      adfCell[iPop].setDiffVariable(iPop);
    }

    FieldD<ADf,DESCRIPTOR,typename FUNCTOR::result_t> j = FUNCTOR{}.compute(adfCell, adfParams);
    FieldD<V,DESCRIPTOR,result_t> djdf;
    for (std::size_t iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      djdf[iPop] = j[0].d(iPop);
      djdf[iPop] *= util::pow(dx, DESCRIPTOR::d); // scaling due to cell-wise contribution
    }
    return djdf;
  }
};

/// Generic functor to compute the jacobian for any evaluated functor in respect to
/// the variables stored in the field RESPECT_TO. Dynamics type is required to
/// provide a cell with full cell interface, e.g., to compute the derivative of a collision.
template <typename FUNCTOR, typename RESPECT_TO, typename DYNAMICS>
struct DerivativeF {
  using parameters = typename FUNCTOR::parameters::template include<descriptors::DX>;

  using result_t = descriptors::FIELD_MATRIX<typename FUNCTOR::result_t, RESPECT_TO>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;

    /// Collecting all accessed fields in order to include them all for correct derivative computation
    using COMBINED_FIELDS = meta::merge<typename CELL::descriptor_t::fields_t, typename FUNCTOR::fields_t>;
    using DESCRIPTOR = typename COMBINED_FIELDS::template decompose_into<CELL::descriptor_t::template extend_by_fields>;

    /// Forward AD type for derivative computation
    using ADf = typename util::ADf<V,DESCRIPTOR::template size<RESPECT_TO>()>;
    using ADf_DYNAMICS = typename DYNAMICS::template exchange_value_type<ADf>;
    const V dx = parameters.template get<descriptors::DX>();

    /// Full cell providing complete cell interface instantiated with forward AD type
    FullCellD<ADf,DESCRIPTOR,ADf_DYNAMICS> adfCell;
    auto adfParams = parameters.template copyAs<ADf>();
    DESCRIPTOR::fields_t::for_each([&](auto id) {
      using FIELD = typename decltype(id)::type;
      adfCell.template setField<FIELD>(FieldD<ADf,DESCRIPTOR,FIELD>(cell.template getField<FIELD>()));
    });

    for (std::size_t iDim=0; iDim<DESCRIPTOR::template size<RESPECT_TO>(); ++iDim) {
      adfCell.template getFieldComponent<RESPECT_TO>(iDim).setDiffVariable(iDim);
    }

    FieldD<ADf,DESCRIPTOR,typename FUNCTOR::result_t> y = FUNCTOR{}.compute(adfCell, adfParams);
    FieldD<V,DESCRIPTOR,typename DerivativeF::result_t> jacobian;
    /// Provide matrix-native view on the serial FieldD
    auto view = DerivativeF::result_t::template getMatrixView<V,DESCRIPTOR>(jacobian);
    for (unsigned row=0; row<view.rows; ++row) {
      for (unsigned col=0; col<view.cols; ++col) {
        view[row][col] = y[row].d(col);
        view[row][col] *= util::pow(dx, DESCRIPTOR::d); // scaling due to cell-wise contribution
      }
    }
    return jacobian;
  }
};

}

}

#endif
