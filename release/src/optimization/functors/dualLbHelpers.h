/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
#ifndef DUAL_LB_HELPERS_H
#define DUAL_LB_HELPERS_H

#include "dynamics/lbm.h"

namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

/// This structure forwards the calls to the appropriate helper class
template<typename T, typename DESCRIPTOR>
struct dualLbHelpers {

  static T equilibrium(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    T eq = olb::equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u) + descriptors::t<T,DESCRIPTOR>(iPop);
    T sum = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      sum += (u[iD] - descriptors::c<DESCRIPTOR>(jPop,iD))*(u[iD] - descriptors::c<DESCRIPTOR>(iPop,iD));
    }
    return (descriptors::invCs2<T,DESCRIPTOR>()*sum + T(1))/rho*eq;
  }

  static T equilibrium2(int iPop, int jPop, T rho, const T u[DESCRIPTOR::d])
  {
    T u_c = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u_c += u[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
    T sum = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      sum += (u[iD]*u[iD] + 2.*descriptors::c<DESCRIPTOR>(jPop,iD)*(descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) )
            * descriptors::invCs2<T,DESCRIPTOR>()*0.5
         + (u_c*descriptors::c<DESCRIPTOR>(iPop,iD)*(2.*descriptors::c<DESCRIPTOR>(jPop,iD)-u[iD]))
            * descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*0.5;
    }
    return descriptors::t<T,DESCRIPTOR>(iPop)*(1.+sum);
  }

  /// Computation of hydrodynamic variables
  // identical computation to standard lbm, but different arguments
  static void computeRhoU(T const* cell, T& rho, T u[DESCRIPTOR::d])
  {

    rho = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
    }
  }

  static void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d])
  {
    computeRhoU(&cell[0], rho, u);
  }

};  // struct dualLbHelpers

template <typename DESCRIPTOR>
struct dualLbMomentaHelpers
{
  /// derive rho by population component
  template <concepts::MinimalCell CELL, typename V=typename CELL::value_t>
  static V dRhoDf(CELL& cell, int jPop) any_platform
  {
    return V{1};
  }

  /// derive momentum j by population component
  template <concepts::MinimalCell CELL, typename DJDF, typename V=typename CELL::value_t>
  static void dJDf(CELL& cell, DJDF& djdf, int jPop) any_platform
  {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      djdf[iD] = descriptors::c<DESCRIPTOR>(jPop,iD);
    }
  }

  /// derive velocity u by population component
  template <concepts::MinimalCell CELL, typename DUDF, typename V=typename CELL::value_t>
  static void dUDf(CELL& cell, DUDF& dudf, int jPop) any_platform
  {
    V rho, u[DESCRIPTOR::d];
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      dudf[iD] = (descriptors::c<DESCRIPTOR>(jPop,iD) - u[iD]) / rho;
    }
  }

  /// derive stress pi by population component
  template <concepts::MinimalCell CELL, typename DPIDF, typename V=typename CELL::value_t>
  static void dPiDf(CELL& cell, DPIDF& dpidf, int jPop) any_platform
  {
    V rho, u[DESCRIPTOR::d];
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);

    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        dpidf[iPi] = descriptors::c<DESCRIPTOR>(jPop,iAlpha) * descriptors::c<DESCRIPTOR>(jPop,iBeta)
          - descriptors::c<DESCRIPTOR>(jPop,iAlpha) * u[iBeta]
          - u[iAlpha] * (descriptors::c<DESCRIPTOR>(jPop,iBeta) - u[iBeta]);
        if (iAlpha==iBeta) {
          dpidf[iPi] -= V{1} / descriptors::invCs2<V,DESCRIPTOR>();
        }
        ++iPi;
      }
    }
  }
};

} // namespace opti

} // namespace olb


#endif
