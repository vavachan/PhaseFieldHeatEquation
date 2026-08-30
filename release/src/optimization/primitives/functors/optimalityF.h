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

#ifndef OPTIMALITY_F_H
#define OPTIMALITY_F_H

namespace olb {

namespace functors {

/// Generic functor to compute the optimality condition in adjoint based
/// simulations. The Jacobian of the collision operator can be computed by
/// the DerivativeF functor.
template <typename PRIMAL_DYNAMICS, typename CONTROLS>
struct OptimalityF {
  using parameters = typename PRIMAL_DYNAMICS::parameters;

  using result_t = opti::DJDALPHA<CONTROLS>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    /// Projecting optimization control to simulation variables
    const V dProjectionDalpha = cell.template getFieldComponent<opti::DPROJECTIONDALPHA<CONTROLS>>(0);
    auto phi = cell.template getField<descriptors::POPULATION>();
    auto dCDalpha = cell.template getField<opti::DCDALPHA<CONTROLS>>();

    /// Provides matrix-native view on the serial dCDalpha vector
    auto view = opti::DCDALPHA<CONTROLS>::template getTransposedMatrixView<V,DESCRIPTOR>(dCDalpha);
    FieldD<V,DESCRIPTOR,typename OptimalityF::result_t> djdalpha;
    for (unsigned row=0; row<view.rows; ++row) {
      for (unsigned col=0; col<view.cols; ++col) {
        djdalpha[row] -= (view[row][col]*phi[col]);
      }
      djdalpha[row] *= dProjectionDalpha;
    }
    return djdalpha;
  }

};

}

}

#endif
