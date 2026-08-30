
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

#ifndef FSI_DYNAMICS_H
#define FSI_DYNAMICS_H

#include "dynamics/forcing.h"

namespace olb {

namespace forcing {

namespace fsi {

/// Homogenized LBM modelling moving porous media and caching post-collision velocity for FSI
struct HLBM {
  static std::string getName() {
    return "HLBM<Kupershtokh>";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template wrap_momentum<
    momenta::MovingPorousMomentum
  >::template wrap_stress<
    momenta::MovingPorousStress
  >::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho{};
      Vector<V,DESCRIPTOR::d> u{};
      MomentaF().computeRhoU(cell, rho, u);

      const V porosity = cell.template getField<descriptors::POROSITY>();
      auto statistic = CollisionO().apply(cell, parameters);

      // This is Kuperstokh forcing
      Vector<V,DESCRIPTOR::d> uPlus = u;
      uPlus += (V{1} - porosity)
             * (cell.template getField<descriptors::VELOCITY>() - u);

      Vector<V,DESCRIPTOR::q> fEqPlus{};
      EquilibriumF().compute(cell, statistic.rho, uPlus, fEqPlus);

      Vector<V,DESCRIPTOR::q> fEq{};
      EquilibriumF().compute(cell, statistic.rho, u, fEq);

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += fEqPlus[iPop] - fEq[iPop];
      }

      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

}

}

#endif
