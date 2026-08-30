/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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
 *  Boston, MA  02110-1201, USA.
*/

#ifndef LATTICE_PHYS_EXTERNAL_SCALAR_2D_HH
#define LATTICE_PHYS_EXTERNAL_SCALAR_2D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysExternalScalar2D.h"
#include "superBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "indicator/superIndicatorF2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF2D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR, typename FIELD>
SuperLatticePhysExternalScalar2D<T,DESCRIPTOR,FIELD>::SuperLatticePhysExternalScalar2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, T convFactorToPhysUnits, std::string name)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = name;
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalScalar2D<T, DESCRIPTOR, FIELD>(
        this->_sLattice.getBlock(iC), convFactorToPhysUnits)
    );
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticePhysExternalScalar2D<T,DESCRIPTOR,FIELD>::BlockLatticePhysExternalScalar2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice, T convFactorToPhysUnits, std::string name)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice, 1),
    _convFactorToPhysUnits(convFactorToPhysUnits)
{
  this->getName() = name;
}

template <typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticePhysExternalScalar2D<T,DESCRIPTOR,FIELD>::operator()(
  T output[], const int input[])
{
  output[0] = this->_blockLattice.get( input[0], input[1] ).template getField<FIELD>() * _convFactorToPhysUnits;
  return true;
}

}
#endif
