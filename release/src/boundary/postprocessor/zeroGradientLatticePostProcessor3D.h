/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Fedor Bukreev, Adrian Kummerländer
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

#ifndef ZERO_GRADIENT_LATTICE_POST_PROCESSOR_3D_H
#define ZERO_GRADIENT_LATTICE_POST_PROCESSOR_3D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "descriptor/descriptor.h"
#include "utilities/omath.h"


namespace olb {

//======================================================================
// ======== Zero Gradient Boundary for AD 3D ======//
//======================================================================
template<typename T, typename DESCRIPTOR>
struct ZeroGradientLatticePostProcessor3D
{
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

template <typename CELL>
void apply(CELL& cell) any_platform{
  using V = typename CELL::value_t;
  for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++){
    const auto c = descriptors::c<DESCRIPTOR>(iPop);
    V v = V{0.};
    V fluidN = cell.neighbor({c[0], c[1], c[2]}).template getField<descriptors::NEIGHBOR>();
    V fluidNN = cell.neighbor({2*c[0], 2*c[1], 2*c[2]}).template getField<descriptors::NEIGHBOR>();
    v = fluidN * cell.neighbor({c[0], c[1], c[2]})[iPop] + fluidN * fluidNN * cell.neighbor({2*c[0], 2*c[1], 2*c[2]})[iPop];
    v *= (fluidN - V{0.5} * fluidN * fluidNN);
    cell[iPop] = v;
  }
}
};
}

#endif
