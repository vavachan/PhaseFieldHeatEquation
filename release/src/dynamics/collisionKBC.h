/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Louis Kronberg, Stephan Simonis
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

#ifndef KBC_COLLISION_H
#define KBC_COLLISION_H

#include "lbm.h"
#include "core/latticeStatistics.h"

namespace olb {

namespace collision {

/// Implementation of the KBC method. See 10.1103/PhysRevE.90.031302.
struct KBC {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "KBC";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {

    // not using a D3Q19 lattice is a compiler error, because currently this struct implements KBC only for D3Q19
    static_assert(DESCRIPTOR::d == 3 && DESCRIPTOR::q == 19, "KBC only implemented for D3Q19");

    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {

      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);

      const V omega = parameters.template get<descriptors::OMEGA>();

      V df[DESCRIPTOR::q];
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        df[iPop]   = cell[iPop] -  fEq[iPop];
        fEq[iPop] += descriptors::t<V,DESCRIPTOR>(iPop);
      }

      V M200 = 0 , M011 = 0;
      V M002 = 0, M110 = 0;
      V M020 = 0, M101 = 0;
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        V c0 = descriptors::c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, 0);
        V c1 = descriptors::c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, 1);
        V c2 = descriptors::c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, 2);

        M200 += c0 * c0 * df[iPop];
        M020 += c1 * c1 * df[iPop];
        M002 += c2 * c2 * df[iPop];

        M110 += c0 * c1 * df[iPop];
        M011 += c1 * c2 * df[iPop];
        M101 += c0 * c2 * df[iPop];
    }


      V Nxz = M200 - M002;
      V Nyz = M020 - M002;

      V Pxy = M110;
      V Pxz = M101;
      V Pyz = M011;


      V ds[DESCRIPTOR::q] = {
        0,
        Nxz/3 - Nyz/6, Nyz/3 - Nxz/6, - Nxz/6 - Nyz/6,
        Pxy/4, -Pxy/4, Pxz/4,
        -Pxz/4, Pyz/4, -Pyz/4,
        Nxz/3 - Nyz/6, Nyz/3 - Nxz/6, - Nxz/6 - Nyz/6,
        Pxy/4, -Pxy/4, Pxz/4,
        -Pxz/4, Pyz/4, -Pyz/4,
      };


      V dh[DESCRIPTOR::q];
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        dh[iPop] = df[iPop] - ds[iPop];
      }

      V beta = omega / 2;
      V dsdh(0);
      V dhdh(0);
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        dsdh += ds[iPop] * dh[iPop] / fEq[iPop];
        dhdh += dh[iPop] * dh[iPop] / fEq[iPop];
      }

      V gamma = dhdh < 1.0e-10 ? 2.0: 1.0/beta - (2 - 1./beta) * dsdh/dhdh;
      //V gamma = 2.0;

      cell.template setField<descriptors::GAMMA>(gamma);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
              cell[iPop] -=  beta * (2*ds[iPop] + gamma*dh[iPop]);
      }

      return statistic;
    };
  };
};

}
}
#endif
