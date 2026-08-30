/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Andreas Schneider, Fedor Bukreev
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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef COLLISION_HRR_H
#define COLLISION_HRR_H

#include "descriptor/fields.h"
#include "core/latticeStatistics.h"
#include "functors/primitive/strainRateTensorFDM2D.h"
#include "functors/primitive/strainRateTensorFDM3D.h"

namespace olb {

template <typename DESCRIPTOR>
struct TensorIndices {
    static constexpr int xx = 0;
    static constexpr int xy = 1;
    static constexpr int xz = 2;
    static constexpr int yy = (DESCRIPTOR::d == 2) ? 2 : 3;
    static constexpr int yz = 4;
    static constexpr int zz = 5;
};

namespace collision {

struct HYBRID : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<HYBRID>()>(1);
  }
};
struct HYBRID_RHO : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<HYBRID_RHO>()>(1);
  }
};

struct RLBThirdOrder {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "RLBThirdOrder";
  }


  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      const V omega = parameters.template get<descriptors::OMEGA>();
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeStress(cell, pi);
      if constexpr (DESCRIPTOR::d == 2) {
        V aNeq3XXY = V(2) * u[0] * pi[TensorIndices<DESCRIPTOR>::xy] + u[1] * pi[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * pi[TensorIndices<DESCRIPTOR>::xy] + u[0] * pi[TensorIndices<DESCRIPTOR>::yy];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * pi[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * pi[TensorIndices<DESCRIPTOR>::yy])
                                                          + V(1.0) * reciSpeedOfSoundQc * hermite2XY * pi[TensorIndices<DESCRIPTOR>::xy]
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3XYY) * (aNeq3XXY + aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3XYY) * (aNeq3XXY - aNeq3XYY));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      else if constexpr (DESCRIPTOR::d == 3) {
        V aNeq3XXY = V(2) * u[0] * pi[TensorIndices<DESCRIPTOR>::xy] + u[1] * pi[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * pi[TensorIndices<DESCRIPTOR>::xy] + u[0] * pi[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3XXZ = V(2) * u[0] * pi[TensorIndices<DESCRIPTOR>::xz] + u[2] * pi[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XZZ = V(2) * u[2] * pi[TensorIndices<DESCRIPTOR>::xz] + u[0] * pi[TensorIndices<DESCRIPTOR>::zz];
        V aNeq3YYZ = V(2) * u[1] * pi[TensorIndices<DESCRIPTOR>::yz] + u[2] * pi[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3YZZ = V(2) * u[2] * pi[TensorIndices<DESCRIPTOR>::yz] + u[1] * pi[TensorIndices<DESCRIPTOR>::zz];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XZ = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2YZ = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2ZZ = descriptors::c<DESCRIPTOR>(iPop,2)*descriptors::c<DESCRIPTOR>(iPop,2) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3XXZ = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3XZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3YYZ = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3YZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,1);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * pi[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * pi[TensorIndices<DESCRIPTOR>::yy] + hermite2ZZ * pi[TensorIndices<DESCRIPTOR>::zz])
                                                          + V(1.0) * reciSpeedOfSoundQc * (hermite2XY * pi[TensorIndices<DESCRIPTOR>::xy] + hermite2XZ * pi[TensorIndices<DESCRIPTOR>::xz] + hermite2YZ * pi[TensorIndices<DESCRIPTOR>::yz])
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3YZZ) * (aNeq3XXY + aNeq3YZZ)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XZZ + hermite3XYY) * (aNeq3XZZ + aNeq3XYY)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3YYZ + hermite3XXZ) * (aNeq3YYZ + aNeq3XXZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3YZZ) * (aNeq3XXY - aNeq3YZZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XZZ - hermite3XYY) * (aNeq3XZZ - aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3YYZ - hermite3XXZ) * (aNeq3YYZ - aNeq3XXZ));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
};

struct HRR {
  using parameters = typename meta::list<descriptors::OMEGA, collision::HYBRID>;

  static std::string getName() {
    return "HRR";
  }


  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      using namespace olb::util::tensorIndices3D;
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V hybridPar = parameters.template get<collision::HYBRID>();
      const V complementHybridPar = V(1) - hybridPar;
      V kinVisc = (V(1)/omega - V(0.5))/descriptors::invCs2<V,DESCRIPTOR>();
      V invPreFactor = -V(1)/(omega * kinVisc * descriptors::invCs2<V,DESCRIPTOR>());
      auto strainRateTensorFD = cell.template getField<descriptors::TENSOR>();

      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeStress(cell, pi);

      V aHybrid[util::TensorVal<DESCRIPTOR>::n] {V(0)};
      for(int i = 0; i < util::TensorVal<DESCRIPTOR>::n; i++) {
        aHybrid[i] = hybridPar * pi[i] + complementHybridPar * invPreFactor * V(2) * kinVisc * rho * strainRateTensorFD[i];
      }

      if constexpr (DESCRIPTOR::d == 2) {
        V aNeq3XXY = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::yy];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * aHybrid[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * aHybrid[TensorIndices<DESCRIPTOR>::yy])
                                                          + V(1.0) * reciSpeedOfSoundQc * hermite2XY * aHybrid[TensorIndices<DESCRIPTOR>::xy]
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3XYY) * (aNeq3XXY + aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3XYY) * (aNeq3XXY - aNeq3XYY));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      else if constexpr (DESCRIPTOR::d == 3) {
        V aNeq3XXY = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3XXZ = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xz] + u[2] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XZZ = V(2) * u[2] * aHybrid[TensorIndices<DESCRIPTOR>::xz] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::zz];
        V aNeq3YYZ = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::yz] + u[2] * aHybrid[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3YZZ = V(2) * u[2] * aHybrid[TensorIndices<DESCRIPTOR>::yz] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::zz];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XZ = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2YZ = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2ZZ = descriptors::c<DESCRIPTOR>(iPop,2)*descriptors::c<DESCRIPTOR>(iPop,2) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3XXZ = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3XZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3YYZ = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3YZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,1);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * aHybrid[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * aHybrid[TensorIndices<DESCRIPTOR>::yy] + hermite2ZZ * aHybrid[TensorIndices<DESCRIPTOR>::zz])
                                                          + V(1.0) * reciSpeedOfSoundQc * (hermite2XY * aHybrid[TensorIndices<DESCRIPTOR>::xy] + hermite2XZ * aHybrid[TensorIndices<DESCRIPTOR>::xz] + hermite2YZ * aHybrid[TensorIndices<DESCRIPTOR>::yz])
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3YZZ) * (aNeq3XXY + aNeq3YZZ)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XZZ + hermite3XYY) * (aNeq3XZZ + aNeq3XYY)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3YYZ + hermite3XXZ) * (aNeq3YYZ + aNeq3XXZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3YZZ) * (aNeq3XXY - aNeq3YZZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XZZ - hermite3XYY) * (aNeq3XZZ - aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3YYZ - hermite3XXZ) * (aNeq3YYZ - aNeq3XXZ));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
};

struct ExternalRhoHRR {
  using parameters = typename meta::list<descriptors::OMEGA, collision::HYBRID, collision::HYBRID_RHO>;

  static std::string getName() {
    return "ExternalRhoHRR";
  }


  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      using namespace olb::util::tensorIndices3D;
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V hybridPar = parameters.template get<collision::HYBRID>();
      const V complementHybridPar = V(1) - hybridPar;
      const V hybridParRho = parameters.template get<collision::HYBRID_RHO>();
      const V complementHybridParRho = V(1) - hybridParRho;
      V kinVisc = (V(1)/omega - V(0.5))/descriptors::invCs2<V,DESCRIPTOR>();
      V invPreFactor = -V(1)/(omega * kinVisc * descriptors::invCs2<V,DESCRIPTOR>());
      auto strainRateTensorFD = cell.template getField<descriptors::TENSOR>();

      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeStress(cell, pi);

      rho *= hybridParRho;
      rho += complementHybridParRho * cell.template getField<descriptors::DENSITY>();

      V aHybrid[util::TensorVal<DESCRIPTOR>::n] {V(0)};
      for(int i = 0; i < util::TensorVal<DESCRIPTOR>::n; i++) {
        aHybrid[i] = hybridPar * pi[i] + complementHybridPar * invPreFactor * V(2) * kinVisc * rho * strainRateTensorFD[i];
      }

      if constexpr (DESCRIPTOR::d == 2) {
        V aNeq3XXY = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::yy];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * aHybrid[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * aHybrid[TensorIndices<DESCRIPTOR>::yy])
                                                          + V(1.0) * reciSpeedOfSoundQc * hermite2XY * aHybrid[TensorIndices<DESCRIPTOR>::xy]
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3XYY) * (aNeq3XXY + aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3XYY) * (aNeq3XXY - aNeq3XYY));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      else if constexpr (DESCRIPTOR::d == 3) {
        V aNeq3XXY = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XYY = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::xy] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3XXZ = V(2) * u[0] * aHybrid[TensorIndices<DESCRIPTOR>::xz] + u[2] * aHybrid[TensorIndices<DESCRIPTOR>::xx];
        V aNeq3XZZ = V(2) * u[2] * aHybrid[TensorIndices<DESCRIPTOR>::xz] + u[0] * aHybrid[TensorIndices<DESCRIPTOR>::zz];
        V aNeq3YYZ = V(2) * u[1] * aHybrid[TensorIndices<DESCRIPTOR>::yz] + u[2] * aHybrid[TensorIndices<DESCRIPTOR>::yy];
        V aNeq3YZZ = V(2) * u[2] * aHybrid[TensorIndices<DESCRIPTOR>::yz] + u[1] * aHybrid[TensorIndices<DESCRIPTOR>::zz];

        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          //Hermite polynomes  https://doi.org/10.1017/S0022112005008153
          V hermite2XX = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,0) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XY = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite2YY = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,1) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();
          V hermite2XZ = descriptors::c<DESCRIPTOR>(iPop,0)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2YZ = descriptors::c<DESCRIPTOR>(iPop,1)*descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite2ZZ = descriptors::c<DESCRIPTOR>(iPop,2)*descriptors::c<DESCRIPTOR>(iPop,2) - V{1}/descriptors::invCs2<V,DESCRIPTOR>();

          V hermite3XXY = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,1);
          V hermite3XYY = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3XXZ = hermite2XX * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3XZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,0);
          V hermite3YYZ = hermite2YY * descriptors::c<DESCRIPTOR>(iPop,2);
          V hermite3YZZ = hermite2ZZ * descriptors::c<DESCRIPTOR>(iPop,1);

          V reciSpeedOfSoundQc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
          V reciSpeedOfSoundHc = descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();

          V fEq = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);

          V fNeq = descriptors::t<V,DESCRIPTOR>(iPop) * ( V(0.5) * reciSpeedOfSoundQc * (hermite2XX * aHybrid[TensorIndices<DESCRIPTOR>::xx] + hermite2YY * aHybrid[TensorIndices<DESCRIPTOR>::yy] + hermite2ZZ * aHybrid[TensorIndices<DESCRIPTOR>::zz])
                                                          + V(1.0) * reciSpeedOfSoundQc * (hermite2XY * aHybrid[TensorIndices<DESCRIPTOR>::xy] + hermite2XZ * aHybrid[TensorIndices<DESCRIPTOR>::xz] + hermite2YZ * aHybrid[TensorIndices<DESCRIPTOR>::yz])
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XXY + hermite3YZZ) * (aNeq3XXY + aNeq3YZZ)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3XZZ + hermite3XYY) * (aNeq3XZZ + aNeq3XYY)
                                                          + V(1./2.) * reciSpeedOfSoundHc * (hermite3YYZ + hermite3XXZ) * (aNeq3YYZ + aNeq3XXZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XXY - hermite3YZZ) * (aNeq3XXY - aNeq3YZZ)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3XZZ - hermite3XYY) * (aNeq3XZZ - aNeq3XYY)
                                                          + V(1./6.) * reciSpeedOfSoundHc * (hermite3YYZ - hermite3XXZ) * (aNeq3YYZ - aNeq3XXZ));

          cell[iPop] = fEq + (V(1) - omega) * fNeq;
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
};

}

namespace forcing {

/// Dynamics combination rule implementing the forcing scheme by Guo et al. for 3rd order RLB
template <template <typename> typename Forced = momenta::Forced>
struct GuoThirdOrder {
  static std::string getName() {
    return "GuoThirdOrderForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V aNeqFst[DESCRIPTOR::d] {};
      for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
        V fneq = cell[iPop] - equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);
        for( int iD = 0; iD<DESCRIPTOR::d; iD++) {
          aNeqFst[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * fneq;
        }
      }
      CollisionO().apply(cell, parameters);
      const auto force = cell.template getField<descriptors::FORCE>();
      if ( util::normSqr<V,DESCRIPTOR::d>(force) != V(0) ) {
        const V omega = parameters.template get<descriptors::OMEGA>();
        lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
        for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
          V thirdOrderTerm{};
          for( int iD = 0; iD<DESCRIPTOR::d; iD++) {
            thirdOrderTerm += descriptors::c<DESCRIPTOR>(iPop,iD) * aNeqFst[iD];
          }
          cell[iPop] += (V(1) - omega) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>() * thirdOrderTerm;
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

/// Dynamics combination rule implementing the forcing scheme by Guo et al. for 3rd order RLB for simulation with wall model
template <template <typename> typename Forced = momenta::Forced>
struct WFGuoThirdOrder {
  static std::string getName() {
    return "WFGuoThirdOrderForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename Forced<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,Forced<MOMENTA>>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename Forced<MOMENTA>::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,Forced<MOMENTA>,EQUILIBRIUM>;

    static constexpr bool is_vectorizable = false;

    static_assert(COLLISION::parameters::template contains<descriptors::OMEGA>(),
                  "COLLISION must be parametrized using relaxation frequency OMEGA");

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V aNeqFst[DESCRIPTOR::d] {};
      for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
        V fneq = cell[iPop] - equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u);
        for( int iD = 0; iD<DESCRIPTOR::d; iD++) {
          aNeqFst[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * fneq;
        }
      }
      CollisionO().apply(cell, parameters);
      const V hybridPar = cell.template getField<collision::HYBRID>();
      if ( hybridPar != V(0) ) {
        const auto force = cell.template getField<descriptors::FORCE>();
        if ( util::normSqr<V,DESCRIPTOR::d>(force) != V(0) ) {
          const V omega = parameters.template get<descriptors::OMEGA>();
          lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
          for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++) {
            V thirdOrderTerm{};
            for( int iD = 0; iD<DESCRIPTOR::d; iD++) {
              thirdOrderTerm += descriptors::c<DESCRIPTOR>(iPop,iD) * aNeqFst[iD];
            }
            cell[iPop] += (V(1) - omega) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>() * thirdOrderTerm;
          }
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

/// PostProcessor for FDM strain rate tensor calculation
class FDMstrainRateTensorPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCell;

    int getPriority() const {
      return 0;
    }

    template <typename CELL, typename V = typename CELL::value_t>
    void apply(CELL& cell) any_platform{
      using DESCRIPTOR = typename CELL::descriptor_t;
      if constexpr (DESCRIPTOR::d == 2) {
        Vector<V,6> strainRate = strainRateTensorFDM2D(cell);
        cell.template setField<descriptors::TENSOR>(strainRate);
      }
      else if constexpr (DESCRIPTOR::d == 3) {
        Vector<V,6> strainRate = strainRateTensorFDM3D(cell);
        cell.template setField<descriptors::TENSOR>(strainRate);
      }
    }
  };
}

#endif
