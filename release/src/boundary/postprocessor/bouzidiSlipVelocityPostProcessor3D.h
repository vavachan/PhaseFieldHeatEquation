
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2023, Fedor Bukrev, Anas Selmi, Adrian Kummerländer, Mathias J. Krause
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

// This file contains the implementation for General slip boundary condition
// used in slip and moderate transition flow regimes (high Kn number)
// See example microPoiseuille flow for use
// It follows the equation below to compute the relative velocity on the wall:
// U_slip = lambda dU/dn
// To use it, first use setBouzidiBoundary<T,DESCRIPTOR,BouzidiSlipVelocityPostProcessor>();
// Then setBouzidiGeneralSlipVelocity() with desired parameters
// The parameters need to be set by sLattice.setParameter as in example

#ifndef SET_GENERAL_SLIP_VELOCITY_POST_PROCESSOR_3D_H
#define SET_GENERAL_SLIP_VELOCITY_POST_PROCESSOR_3D_H

namespace olb {

namespace descriptors {
struct FD_DEG               : public TYPED_FIELD_BASE<int, 1,0,0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<FD_DEG>()>(1);
  }
};
struct BOUZIDI_TUNER : public descriptors::FIELD_BASE<1,0,0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_TUNER>()>(0.);
  }
};
struct LAMBDA               : public FIELD_BASE<1,  0, 0> { };
struct CHAR_LENGTH          : public FIELD_BASE<1,  0, 0> { };
struct CONVERSION_FACTOR_VELOCITY: public FIELD_BASE<1,  0, 0> { };
struct CONVERSION_FACTOR_LENGTH  : public FIELD_BASE<1,  0, 0> { };

}

template <typename CELL, typename V = typename CELL::value_t,
          typename DESCRIPTOR = typename CELL::descriptor_t>
void applyBouzidiVelocity(CELL& x_b) any_platform
{
  const auto q = x_b.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
  const auto veloCoeff =
      x_b.template getFieldPointer<descriptors::BOUZIDI_VELOCITY>();
  //const auto material = x_b.template getField<descriptors::MATERIAL>();
  for (int iPop = 1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
    // update missing population if valid bouzidi distance
    if (q[iPop] > 0) {
      const auto c             = descriptors::c<DESCRIPTOR>(iPop);
      const int  iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      auto       x_s           = x_b.neighbor(c); // solid side neighbor
      auto       x_f           = x_b.neighbor(descriptors::c<DESCRIPTOR>(
          iPop_opposite)); // fluid side neighbor opposite to the missing population
      auto veloTerm = veloCoeff[iPop] * (descriptors::t<V, DESCRIPTOR>(iPop)) *
                      (descriptors::invCs2<V, DESCRIPTOR>());

      x_b[iPop_opposite] =
          (q[iPop] <= V {0.5}) // cut is closer to the fluid cell
              * (V {2} * q[iPop] * x_s[iPop] +
                 (V {1} - V {2} * q[iPop]) * x_b[iPop] - V {2} * veloTerm) +
          (q[iPop] > V {0.5}) // cut is closer to the solid cell
              * (V {0.5} / q[iPop] * x_s[iPop] +
                 V {0.5} * (V {2} * q[iPop] - V {1}) / q[iPop] *
                     x_f[iPop_opposite] -
                 V {1} / q[iPop] * veloTerm);
    }
    // if intersection point is on the cell then fall back to full-way bounce back
    else if (q[iPop] == V {0}) {
      const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      auto veloTerm = veloCoeff[iPop] * (descriptors::t<V, DESCRIPTOR>(iPop)) *
                      (descriptors::invCs2<V, DESCRIPTOR>());
      x_b[iPop_opposite] = x_b[iPop] - V {2} * veloTerm;
    }
  }
};

//======================================================================
// ======== Bouzidi General Slip Velocity 3D Post Processor ======//
//======================================================================

/* Post processor: BouzidiSlipVelocityPostProcessor
 * This class calculate the slip velocity on a cell from
 * neighbor velocities in normal direction.
 */
class BouzidiSlipVelocityPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct TANGENT : public descriptors::FIELD_BASE<0, 1, 0> {};
  struct BOUZIDI_SLIP : public descriptors::FIELD_BASE<0, 0, 1> {};

  int getPriority() const { return -1; }

  using parameters = meta::list<descriptors::BOUZIDI_TUNER, descriptors::LAMBDA,
                                descriptors::CONVERSION_FACTOR_VELOCITY,
                                descriptors::CONVERSION_FACTOR_LENGTH,
                                descriptors::CHAR_LENGTH, descriptors::FD_DEG>;

  template <typename CELL, typename PARAMETERS,
            typename V = typename CELL::value_t>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {

    // get parameters
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto bouzidiVel =
        cell.template getFieldPointer<descriptors::BOUZIDI_VELOCITY>();
    const V tuner  = parameters.template get<descriptors::BOUZIDI_TUNER>();
    const V lambda = parameters.template get<descriptors::LAMBDA>();
    const V conversionFactorVelocity =
        parameters.template get<descriptors::CONVERSION_FACTOR_VELOCITY>();
    const V conversionFactorLength =
        parameters.template get<descriptors::CONVERSION_FACTOR_LENGTH>();
    const V   charLength = parameters.template get<descriptors::CHAR_LENGTH>();
    const int fd_degree  = parameters.template get<
         descriptors::
             FD_DEG>(); // degree of the finite difference gradient calculation

    //lbm::computeEquilibriumSecondOrder() must be done for every population direction
    Vector<V, DESCRIPTOR::d> normal;
    Vector<V, DESCRIPTOR::d> tangent;

    Vector<V, DESCRIPTOR::d> u0_latt {};
    Vector<V, DESCRIPTOR::d> u0_phys {};
    cell.computeU(u0_latt.data());
    u0_phys = u0_latt * conversionFactorVelocity;

    // do nothing if the norm of the normal vector is 0

    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {

      const int  iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
      const auto opp_bouzidi_dist =
          cell.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(
              iPop_opposite);
      normal = cell.template getField<
          descriptors::NORMAL>(); // normal vector with length dx
      V normal_norm = util::norm<DESCRIPTOR::d>(normal);

      if (normal_norm != V {0.}) {
        normal = normal / normal_norm;
      }

      Vector<V, DESCRIPTOR::d> u1_latt {};
      Vector<V, DESCRIPTOR::d> u2_latt {};
      Vector<V, DESCRIPTOR::d> u3_latt {};
      Vector<V, DESCRIPTOR::d> u1_phys {};
      Vector<V, DESCRIPTOR::d> u2_phys {};
      Vector<V, DESCRIPTOR::d> u3_phys {};

      Vector<V, DESCRIPTOR::d> us_phys(V(0.));
      Vector<V, DESCRIPTOR::d> us_t_phys(V(0.));
      // delta is distance to go along the normal where to take value for gradient computation
      V delta = 2;

      if ((opp_bouzidi_dist >= 0) && (normal_norm != 0)) {

        // interpolate velocity at points delta, 2*delta and 3*delta along normal vector
        //Vector<V,3> coeff {-2.5, 4, -1.5}; // coefficients for 3 point gradient computation
        Vector<V, DESCRIPTOR::d> distance_bc(
            opp_bouzidi_dist * descriptors::c<DESCRIPTOR>(iPop_opposite));
        Vector<V, DESCRIPTOR::d> distance = delta * normal + distance_bc;
        if (norm(distance) != V {0.}) {
          interpolate3d<3>(cell, u1_latt, distance, [](auto cell, V res[3]){
            cell.computeU(res);
          });
          u1_phys  = u1_latt * conversionFactorVelocity;
          distance = V(2) * delta * normal + distance_bc;
          //distance = distance_bc + normal * charLength/(4*conversionFactorLength);
          interpolate3d<3>(cell, u2_latt, distance, [](auto cell, V res[3]){
            cell.computeU(res);
          });
          u2_phys = u2_latt * conversionFactorVelocity;
          V sign  = 1;
          if (u1_phys > u2_phys) {
            sign = -1;
          }

          //distance = 3*delta * normal + distance_bc;
          //interpolateVelocityOnPoint(cell, distance, u3_latt);
          //u3_phys = u3_latt * conversionFactorVelocity;
          // Only compute slip if neighbor velocity is non zero
          if (norm(u1_phys) > 0) {
            // compute tangent velocity to boundary from neighbor velocity
            tangent = u1_phys - normal * (u1_phys[0] * normal[0] +
                                          u1_phys[1] * normal[1] +
                                          u1_phys[2] * normal[2]);
            tangent /= norm(tangent);
            if (tangent[1] < 1e-3) {
              tangent[1] = 0;
            }
            if (tangent[2] < 1e-3) {
              tangent[2] = 0;
            }
            tangent /= norm(tangent);

            // compute slip velocity using gradient along normal direction
            auto kn         = lambda / charLength;
            V    lambda_eff = lambda;
            if (kn > 0) {
              lambda_eff = 1 / (1 / lambda + tuner / charLength);
              // tuner = 0 -> usual mean free path, tuner = 1 -> effective mean free path
            }

            if (fd_degree == 0) {
              us_phys = tuner * u1_phys;
            }
            else if (fd_degree == 1) {
              us_phys =
                  sign * lambda_eff /
                  (lambda_eff + delta * conversionFactorLength * norm(normal)) *
                  u1_phys;
            }
            else if (fd_degree == 2) {
              us_phys = sign * lambda_eff /
                        (1.5 * lambda_eff +
                         delta * conversionFactorLength * norm(normal)) *
                        (2 * u1_phys - 0.5 * u2_phys);
              //us_phys = sign * lambda_eff / (1.75*lambda_eff + delta*conversionFactorLength*norm(normal)) * (2.*u1_phys - 0.25*u2_phys);
            }
            //us_phys = lambda_eff * (2.*u1_phys - (7./4.)*u0_phys - 0.25*u2_phys) / (delta*conversionFactorLength); // u slip physical units
            //us_phys = -(1.)*lambda_eff / (1.75*lambda_eff + delta*conversionFactorLength*norm(normal)) * (2.*u1_phys - 0.25*u2_phys);
            //us_phys = (-1.)*lambda_eff / (lambda_eff + delta*conversionFactorLength*norm(normal)) * u1_phys;
            //us_phys = lambda_eff * (u2_phys - u1_phys)/(delta*conversionFactorLength*norm(normal));
            //us_phys = lambda_eff * (-2.5*u1_phys + 4*u2_phys - 1.5*u3_phys) / (delta*conversionFactorLength*norm(normal));

            if (kn >= 10) {
              us_phys = u1_phys;
            }
            // project u_s on the boundary-tangential vector
            us_t_phys = us_phys - (us_phys * normal) * normal;
            //us_t_phys = (us_phys*tangent)*tangent;
          }

          // Conversion from physical units to lattice units
          Vector<V, DESCRIPTOR::d> us_t_lattice =
              us_t_phys / conversionFactorVelocity;
          cell.template setField<descriptors::VELOCITY2>(us_t_phys);
          const auto c              = descriptors::c<DESCRIPTOR>(iPop_opposite);
          //V          rho            = interpolateDensityOnPoint(cell, 2 * normal);
          V          vel_slip_coeff = c * us_t_lattice;
          bouzidiVel[iPop_opposite] = vel_slip_coeff;
        }
      }
    } // for iPop

    applyBouzidiVelocity(cell);
    //BouzidiVelocityPostProcessor::apply(cell)
    // write Fields for output
    cell.template setField<descriptors::NORMAL>(normal);
    cell.template setField<TANGENT>(tangent);
  }
};

} // namespace olb
#endif
