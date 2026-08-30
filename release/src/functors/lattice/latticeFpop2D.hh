/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                     Albert Mink
 *                2025 Shota Ito
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

#ifndef LATTICE_FPOP_2D_HH
#define LATTICE_FPOP_2D_HH

#include "latticeFpop2D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeFpop2D<T, DESCRIPTOR>::SuperLatticeFpop2D(
  SuperLattice<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::q)
{
  this->getName() = "fPop";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeFpop2D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeFpop2D<T, DESCRIPTOR>::BlockLatticeFpop2D(BlockLattice<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::q)
{
  this->getName() = "fPop";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeFpop2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    output[iPop] =
      this->_blockLattice.get(input[0], input[1])[iPop]
      + descriptors::t<T,DESCRIPTOR>(iPop);
  }
  return true;
}

}
#endif
