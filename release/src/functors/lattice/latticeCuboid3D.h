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
 *  Boston, MA  02110-1301, USA.
*/

#ifndef LATTICE_CUBOID_3D_H
#define LATTICE_CUBOID_3D_H

#include "superBaseF3D.h"
#include "superCalcF3D.h"

namespace olb {

/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeCuboid3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeCuboid3D(SuperLattice<T,DESCRIPTOR>& sLattice);
};

/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeCuboid3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  int _iC;
  const Cuboid3D<T>& _cuboid;

public:
  BlockLatticeCuboid3D(BlockLattice<T,DESCRIPTOR>& blockLattice, int iC, const Cuboid3D<T>& cuboid);

  bool operator() (T output[], const int input[]) override;

};

}
#endif
