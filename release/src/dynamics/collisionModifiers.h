/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-2023 Adrian Kummerlaender
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

#ifndef DYNAMICS_COLLISION_MODIFIERS_H
#define DYNAMICS_COLLISION_MODIFIERS_H

#include "lbm.h"
#include "descriptor/fields.h"

namespace olb {

namespace collision {

/// Override COLLISION parameter PARAMETER with cell field PARAMETER
template <typename PARAMETER, typename COLLISION>
struct ParameterFromCell {
  using parameters = typename COLLISION::parameters::template include<PARAMETER>;

  static std::string getName() {
    return "ParameterFromCell<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      parameters.template set<PARAMETER>(
        cell.template getField<PARAMETER>());
      return CollisionO().apply(cell, parameters);
    }
  };
};

/// Override COLLISION parameter OMEGA with cell field OMEGA
template <typename COLLISION>
using OmegaFromCell = ParameterFromCell<descriptors::OMEGA, COLLISION>;

/// Override COLLISION parameter OMEGA with inverse of cell field TAU_EFF
template <typename COLLISION>
struct OmegaFromCellTauEff {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA
  >;

  static std::string getName() {
    return "OmegaFromCellTauEff<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      parameters.template set<descriptors::OMEGA>(
        V{1} / cell.template getField<descriptors::TAU_EFF>());
      return CollisionO().apply(cell, parameters);
    }
  };
};

/// Track time-averaged velocity of COLLISION into cell field AVERAGE_VELOCITY
template <typename COLLISION>
struct TrackAverageVelocity {
  using parameters = typename COLLISION::parameters::template include<descriptors::LATTICE_TIME>;

  static std::string getName() {
    return "TrackAverageVelocity<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      auto statistics = CollisionO().apply(cell, parameters);

      Vector<V,DESCRIPTOR::d> u;
      MomentaF().computeU(cell, u);

      auto iT = parameters.template get<descriptors::LATTICE_TIME>();
      auto uAvg = cell.template getField<descriptors::AVERAGE_VELOCITY>();
      cell.template setField<descriptors::AVERAGE_VELOCITY>((uAvg * (iT-1) + u) / iT);

      return statistics;
    }
  };
};

template <typename COLLISION>
struct StoreAndTrackAverageVelocity {
  using parameters = typename COLLISION::parameters::template include<descriptors::LATTICE_TIME>;

  static std::string getName() {
    return "StoreAndTrackAverageVelocity<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      auto statistics = CollisionO().apply(cell, parameters);

      Vector<V,DESCRIPTOR::d> u;
      MomentaF().computeU(cell, u);

      auto iT = parameters.template get<descriptors::LATTICE_TIME>();
      auto uAvg = cell.template getField<descriptors::AVERAGE_VELOCITY>();
      cell.template setField<descriptors::AVERAGE_VELOCITY>((uAvg * (iT-1) + u) / iT);
      cell.template setField<descriptors::VELOCITY2>(u);

      return statistics;
    }
  };
};

  /// Track time-averaged density of COLLISION into cell field AVERAGE_DENSITY
template <typename COLLISION>
struct TrackAverageDensity {
  using parameters = typename COLLISION::parameters::template include<descriptors::LATTICE_TIME>;

  static std::string getName() {
    return "TrackAverageDensity<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      auto statistics = CollisionO().apply(cell, parameters);

      auto d = MomentaF().computeRho(cell);
      auto iT = parameters.template get<descriptors::LATTICE_TIME>();
      auto dAvg = cell.template getField<descriptors::AVERAGE_DENSITY>();
      cell.template setField<descriptors::AVERAGE_DENSITY>((dAvg * (iT-1) + d) / iT);

      return statistics;
    }
  };
};

/// Save velocity of COLLISION into cell field VELOCITY
template <typename COLLISION>
struct SaveVelocity {
  using parameters = typename COLLISION::parameters;

  static std::string getName() {
    return "SaveVelocity<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    constexpr static bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      auto statistics = CollisionO().apply(cell, parameters);

      Vector<V,DESCRIPTOR::d> u;
      MomentaF().computeU(cell, u);
      cell.template setField<descriptors::VELOCITY>(u);

      return statistics;
    }
  };
};

  /// Track time-averaged TKE and velocity
template <typename COLLISION>
struct TrackAverageTKE {
  using parameters = typename COLLISION::parameters::template include<descriptors::LATTICE_TIME>;

  static std::string getName() {
    return "TrackAverageTKE<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    static constexpr bool is_vectorizable = dynamics::is_vectorizable_v<CollisionO>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      auto statistics = CollisionO().apply(cell, parameters);

      Vector<V,DESCRIPTOR::d> u;
      MomentaF().computeU(cell, u);

      auto iT = parameters.template get<descriptors::LATTICE_TIME>();
      auto uAvg = cell.template getField<descriptors::AVERAGE_VELOCITY>();
      cell.template setField<descriptors::AVERAGE_VELOCITY>((uAvg * (iT-1) + u) / iT);
      auto uDiff = u - uAvg;
      auto TKE = cell.template getField<descriptors::AVERAGE_TKE>();
      auto TKE_new = 0.5*(uDiff[0]*uDiff[0]+uDiff[1]*uDiff[1]+uDiff[2]*uDiff[2]);
      cell.template setField<descriptors::AVERAGE_TKE>((TKE * (iT-1) + TKE_new) / iT);
      return statistics;
    }
  };
};

}
}



#endif
