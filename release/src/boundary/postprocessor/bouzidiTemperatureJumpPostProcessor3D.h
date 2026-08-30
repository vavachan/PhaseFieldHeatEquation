
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2023, Anas Selmi, Shota Ito, Adrian Kummerländer, Mathias J. Krause
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

// This file contains the implementation for temperature jump boundary condition
// used in slip and moderate transition flow regimes (high Kn number)
// See example tempJumpTest for use
// It follows the equation below to compute the temperature of the fluid on the wall:
// T_f = c dT/dn
// To use it, first use setBouzidiBoundary<T,DESCRIPTOR,TemperatureJumpPostProcessor>();
// Then setBouzidiTempJump() with desired parameters
// The parameters need to be set by sLattice.setParameter as in example

#ifndef SET_TEMPERATURE_JUMP_POST_PROCESSOR_3D_H
#define SET_TEMPERATURE_JUMP_POST_PROCESSOR_3D_H
#include "bouzidiSlipVelocityPostProcessor3D.h"
namespace olb {

namespace descriptors {

struct BOUZIDI_WALL_TEMP : public descriptors::FIELD_BASE<0, 0, 1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue()
  {
    return Vector<value_type<T>,
                  DESCRIPTOR::template size<BOUZIDI_WALL_TEMP>()>(0);
  }
};

struct BOUZIDI_JUMP_TEMP : public descriptors::FIELD_BASE<0, 0, 1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue()
  {
    return Vector<value_type<T>,
                  DESCRIPTOR::template size<BOUZIDI_JUMP_TEMP>()>(0);
  }
};

} // namespace descriptors

template <typename CELL, typename V = typename CELL::value_t>
void applyBouzidiTemp(CELL& x_b) any_platform
{
  using DESCRIPTOR = typename CELL::descriptor_t;
  const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
  const auto phi_d =
      x_b.template getFieldPointer<descriptors::BOUZIDI_JUMP_TEMP>();
  V f = V {1 / 3}; // D2Q5
  if (descriptors::q<DESCRIPTOR>() == 7) {
    f = V {0.25}; // D3Q7
  }
  for (int iPop = 1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
    // update missing population if valid bouzidi distance
    if (q[iPop] > V {0}) {
      const auto c             = descriptors::c<DESCRIPTOR>(iPop);
      const int  iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      auto       x_s           = x_b.neighbor(c); // solid side neighbor
      auto       x_f           = x_b.neighbor(descriptors::c<DESCRIPTOR>(
          iPop_opposite)); // fluid side neighbor opposite to the missing population
      auto       source_d =
          phi_d[iPop] * f; // source term set by the dirichlet condition
      auto t_i    = descriptors::t<V, DESCRIPTOR>(iPop);
      auto t_iopp = descriptors::t<V, DESCRIPTOR>(iPop_opposite);

      x_b[iPop_opposite] =
          (q[iPop] <= V {0.5}) // cut is closer to the fluid cell
              * (V {-2} * q[iPop] * (x_s[iPop] + t_i) +
                 (V {2} * q[iPop] - V {1}) * (x_b[iPop] + t_i) + source_d) +
          (q[iPop] > V {0.5}) // cut is closer to the solid cell
              * (V {-1} / (V {2} * q[iPop]) * (x_s[iPop] + t_i) +
                 (V {1} - V {1} / (V {2} * q[iPop])) *
                     (x_f[iPop_opposite] + t_iopp) +
                 (V {1} / (V {2} * q[iPop])) * source_d) -
          t_iopp;
    }
    // if intersection point is on the cell then fall back to full-way bounce back
    else if (q[iPop] == V {0}) {
      const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      auto      source_d      = phi_d[iPop] * f;
      auto      t_i           = descriptors::t<V, DESCRIPTOR>(iPop);
      auto      t_iopp        = descriptors::t<V, DESCRIPTOR>(iPop_opposite);
      x_b[iPop_opposite]      = -(x_b[iPop] + t_i) + source_d - t_iopp;
    }
  }
};

//======================================================================
// ======== Bouzidi Temperature Jump 3D Post Processor ======//
//======================================================================
class TemperatureJumpPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  static V interpolateDensityOnPoint(CELL& cell, Vector<V, 3> distance) any_platform
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
    V rho = 0.;
    for (int nP = 0; nP < 8; nP++) {
      auto dist = distance - surroundingPoints[nP];
      //V volume = util::abs(dist[0] * dist[1] * dist[2]);
      V volume = (V(1) - util::abs(dist[0])) * (V(1) - util::abs(dist[1])) *
                 (V(1) - util::abs(dist[2]));
      rho += volume * cell.neighbor(surroundingPoints[nP]).computeRho();
    }
    return rho;
  }

  using parameters =
      meta::list<descriptors::LAMBDA, descriptors::CONVERSION_FACTOR_LENGTH,
                 descriptors::CHAR_LENGTH,
                 descriptors::SCALAR // SCALAR is C = (2*gamma/(gamma+1))*1/Pr
                 >;

  template <typename CELL, typename PARAMETERS,
            typename V = typename CELL::value_t>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {

    // get parameters
    using DESCRIPTOR = typename CELL::descriptor_t;

    const V lambda = parameters.template get<descriptors::LAMBDA>();
    const V conversionFactorLength =
        parameters.template get<descriptors::CONVERSION_FACTOR_LENGTH>();
    const V    charLength = parameters.template get<descriptors::CHAR_LENGTH>();
    const V    scalar     = parameters.template get<descriptors::SCALAR>();
    const auto phi_d =
        cell.template getFieldPointer<descriptors::BOUZIDI_WALL_TEMP>();
    auto bouzidiTemp =
        cell.template getFieldPointer<descriptors::BOUZIDI_JUMP_TEMP>();

    Vector<V, DESCRIPTOR::d> normal;
    normal = cell.template getField<descriptors::NORMAL>();

    // do nothing if the norm of the normal vector is 0
    V normal_norm = util::norm<DESCRIPTOR::d>(normal);
    if (normal_norm != V {0.}) {
      normal = normal / normal_norm;
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {

      const int  iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      const auto opp_bouzidi_dist =
          cell.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(
              iPop_opposite);

      if ((opp_bouzidi_dist >= 0) && (normal_norm != 0)) {
        V                        delta = 2;
        Vector<V, DESCRIPTOR::d> distance_bc(
            opp_bouzidi_dist * descriptors::c<DESCRIPTOR>(iPop_opposite));
        Vector<V, DESCRIPTOR::d> distance = delta * normal + distance_bc;
        if (norm(distance) != V {0.}) {
          V T_1 = interpolateDensityOnPoint(cell, distance);
          V T_2 = interpolateDensityOnPoint(cell, 2 * delta * normal + distance_bc);
          // compute Jump temperature using gradient along normal direction
          V lambda_eff = lambda;
          // First degree gradient
          V Dx = delta * conversionFactorLength;
          //V T_f = (lambda_eff*scalar*T_1 + Dx*phi_d[iPop_opposite]) / (Dx + scalar*lambda_eff) ;
          V T_f = (scalar * lambda * (4 * T_1 - T_2) +
                   2 * Dx * phi_d[iPop_opposite]) /
                  (2 * Dx + 3 * scalar * lambda);

          bouzidiTemp[iPop_opposite] = T_f;
        }
      }
    } // for iPop

    applyBouzidiTemp(cell);
  }
};

} // namespace olb
#endif
