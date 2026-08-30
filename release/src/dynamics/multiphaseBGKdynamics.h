/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Luiz Eduardo Czelusniak
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
 * BGK Dynamics with adjusted omega -- header file.
 *
 */
#ifndef MULTIPHASE_BGK_DYNAMICS_H
#define MULTIPHASE_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

// *INDENT-OFF*

namespace multiphase {

struct OMEGA_VAPOR : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct OMEGA_LIQUID : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
 };

// The following is used for Herschel-Bulkley only
struct RHO_VAPOR : public descriptors::FIELD_BASE<1> { };
struct RHO_LIQUID : public descriptors::FIELD_BASE<1> { };

/// Compute and update cell-wise multiphase OMEGA
template <typename COLLISION>
struct OmegaFromCell {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, OMEGA_VAPOR, OMEGA_LIQUID, RHO_VAPOR, RHO_LIQUID
  >;

  static std::string getName()
  {
    return "multiphase::OmegaFromCell<" + COLLISION::getName() + ">";
  }

  template <typename CELL, typename PARAMETERS, typename RHO, typename V=typename CELL::value_t>
  static V computeOmega(CELL& cell, PARAMETERS& parameters, RHO& rho) any_platform
  {
    //using DESCRIPTOR = typename CELL::descriptor_t;
    V rho_v = parameters.template get<RHO_VAPOR>();
    V rho_l = parameters.template get<RHO_LIQUID>();
    V omega_v = parameters.template get<OMEGA_VAPOR>();
    V omega_l = parameters.template get<OMEGA_LIQUID>();
    V newOmega = ( rho - rho_l ) / ( rho_v - rho_l ) * omega_v
                + ( rho - rho_v ) / ( rho_l - rho_v ) * omega_l;
    return newOmega;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
    {
      const V rho = MomentaF().computeRho(cell);
      const V newOmega = computeOmega(cell, parameters, rho);
      cell.template setField<descriptors::OMEGA>(newOmega);
      parameters.template set<descriptors::OMEGA>(newOmega);
      return CollisionO().apply(cell, parameters);
    }
  };
};

}

/// BGK collision using Multiphase collision frequency
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using MultiphaseBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  multiphase::OmegaFromCell<collision::BGK>
>;

/// BGK collision using Multiphase collision frequency with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using MultiphaseForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  multiphase::OmegaFromCell<collision::BGK>,
  forcing::Guo<momenta::Forced>
>;

// *INDENT-ON*

}

#endif
