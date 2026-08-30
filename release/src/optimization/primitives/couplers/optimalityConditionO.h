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

#ifndef OPTIMALITY_CONDITION_O_H
#define OPTIMALITY_CONDITION_O_H

#include <solver/names.h>
#include <core/superLatticeCoupling.h>

namespace olb {

namespace couplers {

// Operator which implements the optimality condition for a flow control
// problem where the porosity is controlled.
template <typename PRIMAL_DYNAMICS>
struct PorosityControlO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = typename PRIMAL_DYNAMICS::parameters;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) {
    using V = typename CELLS::template value_t<names::Primal>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Primal>::descriptor_t;
    auto& primal_cell = cells.template get<names::Primal>();
    auto& dual_cell = cells.template get<names::Dual>();
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V dProjectiondAlpha = dual_cell.template getField<opti::DPROJECTIONDALPHA>();

    V sensitivity{0};

    V rho_f{0}; V u_f[DESCRIPTOR::d]{0};
    rho_f = typename PRIMAL_DYNAMICS::MomentaF().computeRho(primal_cell);
    typename PRIMAL_DYNAMICS::MomentaF().computeU(primal_cell, u_f);

    const V d = dual_cell.template getField<descriptors::POROSITY>();
    for (std::size_t jPop=0; jPop<DESCRIPTOR::q; ++jPop) {
      const V phi_j = dual_cell[jPop];
      const V feq_j = equilibrium<DESCRIPTOR>::secondOrder(jPop, rho_f, u_f) + descriptors::t<V,DESCRIPTOR>(jPop);
      for (std::size_t iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
        sensitivity += phi_j*feq_j*(descriptors::c<DESCRIPTOR>(jPop, iDim) - d*u_f[iDim])*u_f[iDim]*dProjectiondAlpha;
      }
    }
    sensitivity *= -omega*descriptors::invCs2<V,DESCRIPTOR>();
    dual_cell.template setField<opti::SENSITIVITY>(sensitivity);
  }
};

}

}

#endif
