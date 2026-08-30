/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerländer
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

#ifndef KNUDSEN_VELOCITY_POST_PROCESSOR_H
#define KNUDSEN_VELOCITY_POST_PROCESSOR_H

namespace olb {
//======================================================================
// ======== Knudsen General Slip PostProcessor ======//
//======================================================================
// interpolation of the momentum along the surface normal on the chosen distance from the boundary cell
  template <typename CELL, typename V = typename CELL::value_t, typename DESCRIPTOR = typename CELL::descriptor_t>
  Vector<V,DESCRIPTOR::d> interpolateMomentum(CELL& cell, Vector<V,DESCRIPTOR::d> distance) any_platform{
    Vector<V,DESCRIPTOR::d> floorV;
    Vector<V,DESCRIPTOR::d> momentum;
    for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
      floorV[iD] = util::floor(distance[iD]);
    }
    Vector<Vector<V,DESCRIPTOR::d>,(DESCRIPTOR::d==2)*4+(DESCRIPTOR::d==3)*8> surroundingPoints(floorV);
    surroundingPoints[1][0] += 1.;
    surroundingPoints[2][1] += 1.;
    surroundingPoints[3][0] += 1.; surroundingPoints[3][1] += 1.;
    if( DESCRIPTOR::d == 3) {
      surroundingPoints[4][2] += 1.;
      surroundingPoints[5][0] += 1.; surroundingPoints[5][2] += 1.;
      surroundingPoints[6][1] += 1.; surroundingPoints[6][2] += 1.;
      surroundingPoints[7][0] += 1.; surroundingPoints[7][1] += 1.; surroundingPoints[7][2] += 1.;
    }

    V sumWeight = V(0);
    for (auto point : surroundingPoints) {
      const Vector<V,DESCRIPTOR::d> dist = distance - point;
      V weight = V(1);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        weight *=(V(1.) - util::abs(dist[iD]));
      }
      auto jNP = cell.neighbor(point).template getField<descriptors::VELOCITY>();
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        if(!std::isnan(jNP[iD])) {
          momentum[iD] += jNP[iD]*weight;
          sumWeight += (iD==0)*weight;
        }
      }
    }
    momentum *= V(1)/sumWeight;
    return momentum;
  };

  template<bool thermalCreep = false>
  class KnudsenVelocityPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
    struct SLIPCOEFF             : public descriptors::FIELD_BASE<1> { };
    struct CREEPCOEFF             : public descriptors::FIELD_BASE<1> { };
    using parameters = meta::list<SLIPCOEFF,CREEPCOEFF>;

    int getPriority() const {
      return -10;
    }

    template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform{

      using DESCRIPTOR = typename CELL::descriptor_t;
      auto normal = cell.template getField<descriptors::NORMAL>();
      if(util::norm<DESCRIPTOR::d>(normal) != V(0)) {
        const V slipCoeff = parameters.template get<SLIPCOEFF>();
        const V creepCoeff = parameters.template get<CREEPCOEFF>();
        const V dX = V(2);
        V veloCoeff[DESCRIPTOR::q] {V(0)};
        Vector<V,DESCRIPTOR::d> uC {V(0)};
        if(thermalCreep) {
          V rho = cell.computeRho();
          uC = cell.template getField<descriptors::TEMPGRADIENT>();
          uC *= rho*creepCoeff;
        }
        for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
          const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
          const auto opp_bouzidi_dist = cell.template getFieldComponent<descriptors::BOUZIDI_DISTANCE>(iPop_opposite);
          if (opp_bouzidi_dist >= 0) {
            const auto c = descriptors::c<DESCRIPTOR>(iPop_opposite);
            auto dist0 = opp_bouzidi_dist * c + normal*dX;
            auto dist1 = opp_bouzidi_dist * c + normal*V(2)*dX;
            auto jSample0 = interpolateMomentum(cell, dist0);
            auto jSample1 = interpolateMomentum(cell, dist1);
            V jSampleNormal0 = jSample0*normal;
            V jSampleNormal1 = jSample1*normal;
            auto jSampleTang0 = jSample0 - jSampleNormal0*normal;
            auto jSampleTang1 = jSample1 - jSampleNormal1*normal;
            auto jSlip = V(1)/(1. + (1./dX)*3./2.*slipCoeff)*(slipCoeff/dX*0.5*(4.*jSampleTang0 - jSampleTang1) + uC);
            veloCoeff[iPop_opposite] =  c*jSlip;
          }
        }

        // Bouzidi step
        const auto q = cell.template getFieldPointer<descriptors::BOUZIDI_DISTANCE>();
        for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
          // update missing population if valid bouzidi distance
          if (q[iPop] > 0) {
            const auto c = descriptors::c<DESCRIPTOR>(iPop);
            const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
            auto x_s = cell.neighbor(c);                                             // solid side neighbor
            auto x_f = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop_opposite));     // fluid side neighbor opposite to the missing population
            auto veloTerm = veloCoeff[iPop] * (descriptors::t<V,DESCRIPTOR>(iPop)) * (descriptors::invCs2<V,DESCRIPTOR>());

            cell[iPop_opposite] = (q[iPop] <= V{0.5}) // cut is closer to the fluid cell
                               * (V{2} * q[iPop] * x_s[iPop] + (V{1} - V{2} * q[iPop]) * cell[iPop] - V{2} * veloTerm)
                               + (q[iPop] >  V{0.5}) // cut is closer to the solid cell
                               * (V{0.5} / q[iPop] * x_s[iPop] + V{0.5} * (V{2} * q[iPop] - V{1}) / q[iPop] * x_f[iPop_opposite] - V{1}/q[iPop] * veloTerm);
          }
          // if intersection point is on the cell then fall back to full-way bounce back
          else if (q[iPop] == V{0}) {
            const int iPop_opposite = descriptors::opposite<DESCRIPTOR>(iPop);
            auto veloTerm = veloCoeff[iPop] * (descriptors::t<V,DESCRIPTOR>(iPop)) * (descriptors::invCs2<V,DESCRIPTOR>());
            cell[iPop_opposite] = cell[iPop] - V{2} * veloTerm;
          }
        }
      }
    }
  };
}
#endif
