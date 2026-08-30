/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef SHAN_CHEN_FORCED_SINGLE_COMPONENT_COUPLING_H
#define SHAN_CHEN_FORCED_SINGLE_COMPONENT_COUPLING_H

#include "core/operator.h"

namespace olb {


template <typename T, typename DESCRIPTOR, typename POTENTIAL>
struct ShanChenForcedSingleComponentPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct G : public descriptors::FIELD_BASE<1> { };
  struct RHO0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA : public descriptors::FIELD_BASE<1> { };

  using parameters = typename POTENTIAL::parameters::template include<G,RHO0,OMEGA>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    //using DESCRIPTOR = typename CELL::descriptor_t;

    auto j = cell.template getField<descriptors::VELOCITY>();
    lbm<DESCRIPTOR>::computeJ(cell, j);
    cell.template setField<descriptors::VELOCITY>(j);

    const V g = parameters.template get<G>();
    const V rho0 = parameters.template get<RHO0>();
    const V omega = parameters.template get<OMEGA>();

    const V rho = cell.template getFieldComponent<descriptors::STATISTIC>(0) * rho0;
    // Computation of the common velocity, shared among the two populations
    const V rhoTot = rho * omega;

    Vector<V,DESCRIPTOR::d> uTot{};
    auto u = cell.template getField<descriptors::VELOCITY>();
    uTot = (u*rho0*omega) / rhoTot;

    // Computation of the interaction potential
    Vector<V,DESCRIPTOR::d> rhoContribution;
    const V psi = POTENTIAL().compute(rho, parameters);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V rho_ = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop))
                         .template getFieldComponent<descriptors::STATISTIC>(0)
                   * rho0;
      const V psi_ = POTENTIAL().compute(rho_, parameters);
      rhoContribution += psi * psi_
                       * descriptors::c<DESCRIPTOR>(iPop)
                       * descriptors::t<V,DESCRIPTOR>(iPop);
    }

    // Computation and storage of the final velocity, consisting
    //   of u and the momentum difference due to interaction
    //   potential plus external force
    auto force = cell.template getField<descriptors::EXTERNAL_FORCE>();
    cell.template setField<descriptors::VELOCITY>(uTot);
    cell.template setField<descriptors::FORCE>(force - g*rhoContribution/rho);
  }

};


}

#endif
