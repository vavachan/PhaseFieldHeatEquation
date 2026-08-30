/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Asher Zarth, Lukas Baron, Mathias J. Krause, Jonas Latt, Jan E. Marquardt
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

/** \file
 * BGK Dynamics for porous media -- header file.
 */
#ifndef POROUS_FORCED_BGK_DYNAMICS_H
#define POROUS_FORCED_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousForcedBGKdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::Guo<momenta::Forced>
>;

//////////////////////////////////////////////////////////////////////////////

/* Implementations of the BGK collision for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field,
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
/// Guo forced BGK collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticleGuoForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,false>,
  forcing::Guo<momenta::Forced>
>;

/// Guo forced BGK static collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using StaticPorousParticleGuoForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,true>,
  forcing::Guo<momenta::Forced>
>;

/// ShanChen forced BGK collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticleShanChenForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,false>,
  forcing::ShanChen
>;

/// ShanChen forced BGK static collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using StaticPorousParticleShanChenForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,true>,
  forcing::ShanChen
>;

namespace forcing {

template <bool isStatic = false>
struct PorousParticleKupershtokh {
  static std::string getName() {
    return "PorousParticleKupershtokhForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::Forced<momenta::PorousParticle<MOMENTA>>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,momenta::PorousParticle<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    constexpr static bool is_vectorizable = false;

    template <typename CELL, typename VELOCITY, typename V=typename CELL::value_t>
    void calculate(CELL& cell, VELOCITY& u) {
      if constexpr (isStatic) {
        for (int i=0; i<DESCRIPTOR::d; ++i)  {
          u[i] -= (V{1} - cell.template getField<descriptors::POROSITY>()) * u[i];
        }
      } else {
        for (int i=0; i<DESCRIPTOR::d; ++i)  {
          u[i] +=   (V{1} - cell.template getField<descriptors::POROSITY>())
                  * (  cell.template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(i)
                     / cell.template getField<descriptors::VELOCITY_DENOMINATOR>() - u[i]);
        }
      }
    }

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V u[DESCRIPTOR::d];
      V rho;
      MomentaF().computeRhoU(cell, rho, u);
      const auto velDenominator = cell.template getField<descriptors::VELOCITY_DENOMINATOR>();
      const auto statistic = CollisionO().apply(cell, parameters);
      const auto force = cell.template getField<descriptors::FORCE>();
      V uPlusDeltaU[DESCRIPTOR::d];
      for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
        uPlusDeltaU[iVel] = u[iVel] + force[iVel] / rho;
      }
      if (velDenominator > std::numeric_limits<V>::epsilon()) {
        calculate(cell, uPlusDeltaU);
        if constexpr (!isStatic) {
          // reset external field for next timestep
          cell.template setField<descriptors::POROSITY>(1.);
          cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0.);
          cell.template setField<descriptors::VELOCITY_NUMERATOR>({0.,0.,0.});
        }
      }
      V fEq[DESCRIPTOR::q] { };
      V fEq2[DESCRIPTOR::q] { };
      EquilibriumF().compute(cell, statistic.rho, u, fEq);
      EquilibriumF().compute(cell, statistic.rho, uPlusDeltaU, fEq2);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] += fEq2[iPop] - fEq[iPop];
      }

      // Reset contact helper if utilized
      if constexpr (DESCRIPTOR::template provides<descriptors::CONTACT_DETECTION>()) {
        cell.template setField<descriptors::CONTACT_DETECTION>(0);
      }

      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;

};

}

/// Kuperstokh forced BGK collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticleKupershtokhForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::PorousParticleKupershtokh<false>
>;

/// Kuperstokh forced BGK static collision for moving porous media (HLBM approach)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using StaticPorousParticleKupershtokhForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::PorousParticleKupershtokh<true>
>;

} // olb

#endif
