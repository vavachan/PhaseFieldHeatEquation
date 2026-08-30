/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Dennis Teutscher, Alexander Schulz
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

//This file contains the Slip Boundary
//This is an onLattice boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SLIP_HH
#define SLIP_HH

#include "slip3D.h"

namespace olb {

namespace boundary {

template <int NX, int NY, int NZ>
struct FullSlipO3D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& x_b) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    int reflectionPop[DESCRIPTOR::q];
    int mirrorDirection0;
    int mirrorDirection1;
    int mirrorDirection2;
    int mult = 2 / (NX*NX + NY*NY + NZ*NZ);
    reflectionPop[0] =0;
    for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
      reflectionPop[iPop] = 0;
      // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
      int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*NX + descriptors::c<DESCRIPTOR>(iPop,1)*NY + descriptors::c<DESCRIPTOR>(iPop,2)*NZ;
      if ( scalarProduct < 0) {
        // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
        if (mult == 0) {
          mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
          mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
          mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
        }
        else {
          mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*NX;
          mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*NY;
          mirrorDirection2 = descriptors::c<DESCRIPTOR>(iPop,2) - mult*scalarProduct*NZ;
        }

        // run through all lattice directions and look for match of direction
        for (int i = 1; i < DESCRIPTOR::q; i++) {
          if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
              && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1
              && descriptors::c<DESCRIPTOR>(i,2)==mirrorDirection2) {
            reflectionPop[iPop] = i;
            break;
          }
        }
      }
    }
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      if (reflectionPop[iPop]!=0) {
        //do reflection
        x_b[iPop] = x_b[reflectionPop[iPop]];
      }
    }
  }

};

}

}

#endif
