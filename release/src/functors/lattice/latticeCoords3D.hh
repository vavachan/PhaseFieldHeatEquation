/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_COORDS_3D_HH
#define LATTICE_COORDS_3D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeCoords3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF3D.h"
#include "blockBaseF3D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeCoords3D<T,DESCRIPTOR>::SuperLatticeCoords3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& sGeometry)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 2)
{
  this->getName() = "coordinates";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCoords3D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC),
                                sGeometry.getBlock(iC)) );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeCoords3D<T,DESCRIPTOR>::BlockLatticeCoords3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const BlockGeometry<T,3>& iGeometry)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3), _iGeometry(iGeometry)
{
  this->getName() = "coordinates";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeCoords3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _iGeometry.getPhysR(output, input);
  return true;
}

}
#endif
