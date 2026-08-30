/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Anas Selmi, Adrian Kummerlaender
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
#ifndef THERMAL_CREEP_COUPLING
#define THERMAL_CREEP_COUPLING

#include "../boundary/postprocessor/bouzidiSlipVelocityPostProcessor3D.h"
#include "core/operator.h"


namespace olb {

namespace descriptors {

struct BOUZIDI_SLIP_CREEP   : public FIELD_BASE<0,  0, 1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_SLIP_CREEP>()>(0.);
  }
};

}

// simple function to return closed lattice velocity to tangent to check for wall
template <typename V, typename DESCRIPTOR>
int findClosedLatticeVector(Vector<V, DESCRIPTOR::d> tangent) any_platform
{

  int iPopClosest = 0;
  V   cos         = 0.;
  tangent         = tangent / norm(tangent);

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto c = descriptors::c<DESCRIPTOR>(iPop);
    c      = c / norm(c);
    if (util::abs(tangent * c) > cos) {
      cos         = util::abs(tangent * c);
      iPopClosest = iPop;
    }
  }
  if (tangent * descriptors::c<DESCRIPTOR>(iPopClosest) < 0) {
    iPopClosest = descriptors::opposite<DESCRIPTOR>(iPopClosest);
  }
  return iPopClosest;
};

// thermal creep coupling: Implement velocity creep from temperature gradient on boundaries

struct ThermalCreepBouzidiCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct R : public descriptors::FIELD_BASE<1> {};
  struct MU : public descriptors::FIELD_BASE<1> {};
  struct CONVERSION_FACTOR_PRESSURE : public descriptors::FIELD_BASE<1> {};
  struct CHAR_PHYS_PRESSURE : public descriptors::FIELD_BASE<1> {};
  struct CONVERSION_FACTOR_TEMPERATURE : public descriptors::FIELD_BASE<1> {};
  struct N_XYZ : public descriptors::FIELD_BASE<0, 1, 0> {};

  using parameters =
      meta::list<R, MU, CONVERSION_FACTOR_PRESSURE, CHAR_PHYS_PRESSURE,
                 CONVERSION_FACTOR_TEMPERATURE,
                 descriptors::CONVERSION_FACTOR_VELOCITY,
                 descriptors::CONVERSION_FACTOR_LENGTH, N_XYZ>;

  int getPriority() const { return 0; }

  template <typename CELL, typename V = typename CELL::value_t>
  static V interpolateTemperatureOnPoint(CELL&        cell,
                                         Vector<V, 3> distance) any_platform
  {
    static_assert(CELL::descriptor_t::d == 3);
    Vector<int, 3> floorV(util::floor(distance[0]), util::floor(distance[1]),
                          util::floor(distance[2]));
    Vector<Vector<int, 3>, 8> surroundingPoints(floorV);
    // surroundingPoints[0] contains point at floor x, floor y, floor z
    surroundingPoints[1][0] += 1; // x+1, y, z
    surroundingPoints[2][1] += 1; // x, y+1, z
    surroundingPoints[3][2] += 1; // x, y, z+1
    surroundingPoints[4][0] += 1;
    surroundingPoints[4][1] += 1; // x+1, y+1,z
    surroundingPoints[5][0] += 1;
    surroundingPoints[5][2] += 1; // x+1, y, z+1
    surroundingPoints[6][1] += 1;
    surroundingPoints[6][2] += 1; // x, y+1, z+1
    surroundingPoints[7][0] += 1;
    surroundingPoints[7][1] += 1;
    surroundingPoints[7][2] += 1; // x+1, y+1, z+1
    V temperature = 0.;
    for (int nP = 0; nP < 8; nP++) {
      auto dist = distance - surroundingPoints[nP];
      //V volume = util::abs(dist[0] * dist[1] * dist[2]);
      V volume = (V(1) - util::abs(dist[0])) * (V(1) - util::abs(dist[1])) *
                 (V(1) - util::abs(dist[2]));
      V TnP = cell.neighbor(surroundingPoints[nP]).computeRho();
      temperature += TnP * volume;
    }
    return temperature;
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR =
        typename CELLS::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    // Set advection velocity for ADE lattice
    Vector<V, DESCRIPTOR::d> u {};
    cellNSE.computeU(u.data());
    cellADE.template setField<descriptors::VELOCITY>(u);

    auto normal = cellNSE.template getField<descriptors::NORMAL>();
    auto bouzidiVel =
        cellNSE.template getFieldPointer<descriptors::BOUZIDI_VELOCITY>();
    //auto creepPopulations = cellNSE.template getFieldPointer<descriptors::BOUZIDI_SLIP_CREEP>();
    auto thermal_creep =
        cellNSE.template getFieldPointer<descriptors::AVERAGE_VELOCITY>();
    //auto bouzidiSlip = cellNSE.template getFieldPointer<BouzidiSlipVelocityPostProcessor::BOUZIDI_SLIP>();
    auto q = cellNSE.template getField<descriptors::BOUZIDI_DISTANCE>();
    V    R_specific = parameters.template get<R>();
    V    mu         = parameters.template get<MU>();
    V    conversionPressure =
        parameters.template get<CONVERSION_FACTOR_PRESSURE>();
    V charPhysPressure = parameters.template get<CHAR_PHYS_PRESSURE>();
    V conversionTemperature =
        parameters.template get<CONVERSION_FACTOR_TEMPERATURE>();
    V conversionLength =
        parameters.template get<descriptors::CONVERSION_FACTOR_LENGTH>();
    V conversionVelocity =
        parameters.template get<descriptors::CONVERSION_FACTOR_VELOCITY>();
    //auto N_xyz = parameters.template get<N_XYZ>();

    // Only compute thermal creep for cells on boundary, i.e. normal vector is non zero
    if (norm(normal) != 0) {
      normal = normal / norm(normal);
      Vector<V, DESCRIPTOR::d> tangent1(0.);
      Vector<V, DESCRIPTOR::d> tangent2(0.);
      Vector<V, DESCRIPTOR::d> xyz(1., 0., 0.);

      // Calculate tangent to boundary
      if (std::abs(xyz * normal) != norm(xyz) * norm(normal)) {
        tangent1 = xyz - normal * (xyz * normal);
        tangent1 /= norm(tangent1);
      }
      else {
        tangent1[0] = 0;
        tangent1[1] = 1;
        tangent1[2] = 0;
      }

      tangent2 = crossProduct(normal, tangent1);
      tangent2 = tangent2 / norm(tangent2);
      Vector<V, DESCRIPTOR::d> u_creep_phys(V(0.));
      Vector<V, DESCRIPTOR::d> u_creep_latt(V(0.));
      V                        T_neighbor1 = 0;
      V                        T_neighbor2 = 0;
      V                        rho         = cellNSE.computeRho();
      V p = (rho - 1.0) / descriptors::invCs2<V, DESCRIPTOR>() *
                conversionPressure +
            charPhysPressure;

      assert(p!=0 && "Pressure must be non-zero");

      int iPopClosest1 = findClosedLatticeVector<V, DESCRIPTOR>(tangent1);
      //const auto iPopClosest1BouzidiWall = cellNSE.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPopClosest1);
      //auto cClosest1    = descriptors::c<DESCRIPTOR>(iPopClosest1);
      int iPopClosest2 = findClosedLatticeVector<V, DESCRIPTOR>(tangent2);
      //const auto iPopClosest2BouzidiWall = cellNSE.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPopClosest2);
      //auto cClosest2 = descriptors::c<DESCRIPTOR>(iPopClosest2);

      if (q[iPopClosest1] < 0 &&
          q[descriptors::opposite<DESCRIPTOR>(iPopClosest1)] <
              0) { //(iPopClosest1BouzidiWall < 0)
        T_neighbor1  = interpolateTemperatureOnPoint(cellADE, tangent1);
        T_neighbor2  = interpolateTemperatureOnPoint(cellADE, -1 * tangent1);
        T_neighbor1  = T_neighbor1 * conversionTemperature;
        T_neighbor2  = T_neighbor2 * conversionTemperature;
        u_creep_phys = tangent1 * V(0.75) * (mu * R_specific / p) *
                       (T_neighbor1 - T_neighbor2) / (2 * conversionLength);
      }
      if (q[iPopClosest2] < 0 &&
          q[descriptors::opposite<DESCRIPTOR>(iPopClosest2)] < 0) {
        T_neighbor1 = interpolateTemperatureOnPoint(cellADE, tangent2);
        T_neighbor2 = interpolateTemperatureOnPoint(cellADE, -1 * tangent2);
        T_neighbor1 = T_neighbor1 * conversionTemperature;
        T_neighbor2 = T_neighbor2 * conversionTemperature;
        u_creep_phys += tangent2 * V(0.75) * (mu * R_specific / p) *
                        (T_neighbor1 - T_neighbor2) / (2 * conversionLength);
      }

      thermal_creep = u_creep_phys;
      u_creep_latt  = u_creep_phys / conversionVelocity;

      // Loop over Lattice directions
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        const int  iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
        const auto opp_bouzidi_dist =
            cellNSE.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(
                iPop_opposite);
        if (opp_bouzidi_dist >= 0) {
          const auto c = descriptors::c<DESCRIPTOR>(iPop_opposite);
          //V          vel_creep_coeff = c * (u_creep_latt + u);
          V vel_creep_coeff = c * u_creep_latt;
          //creepPopulations[iPop_opposite] = vel_creep_coeff;
          bouzidiVel[iPop_opposite] += vel_creep_coeff;
        }
      }
    }
    // Write new Bouzidi velocity on cells
    applyBouzidiVelocity(cellNSE);
  }
};

} // namespace olb
#endif
