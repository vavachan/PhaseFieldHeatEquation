/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Dennis Teutscher, Alexander Schulz
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

//This file contains the Partial Slip Boundary
//This is an onLattice boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef PARTIAL_SLIP_2D_HH
#define PARTIAL_SLIP_2D_HH

#include "partialSlip2D.h"

namespace olb {

namespace boundary {

template <int NX, int NY>
struct PartialSlipO2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return -1;
  }

  using parameters = typename meta::list<descriptors::TUNER>;


  template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
  void apply(CELL& x_b, PARAMETERS& parameters) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V tuner = parameters.template get<descriptors::TUNER>();
    int reflectionPop[DESCRIPTOR::q];
    int mirrorDirection0;
    int mirrorDirection1;
    int mult = 2 / (NX*NX + NY*NY);
    reflectionPop[0] =0;
    for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      reflectionPop[iPop] = 0;
      // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
      int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*NX + descriptors::c<DESCRIPTOR>(iPop,1)*NY;
      if ( scalarProduct < 0) {
        // bounce back for the case discreteNormalX = discreteNormalY = 1, that is mult=1
        if (mult == 1) {
          mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
          mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        }
        else {
          mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*NX;
          mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*NY;
        }

        // run through all lattice directions and look for match of direction
        for (int i = 1; i < DESCRIPTOR::q; i++) {
          if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
              && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1) {
            reflectionPop[iPop] = i;
            break;
          }
        }
      }
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      if (reflectionPop[iPop]!=0) {
        //do reflection
        x_b[iPop] = tuner*x_b[reflectionPop[iPop]];
      }
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q/2 ; ++iPop) {
      V provv = x_b[descriptors::opposite<DESCRIPTOR>(iPop)];
      x_b[descriptors::opposite<DESCRIPTOR>(iPop)] += (1.-tuner)*x_b[iPop];
      x_b[iPop] += (1.-tuner)*provv;
    }
  }

};

}

}

#endif
