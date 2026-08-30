/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerländer, Christoph Gaul, Mathias J. Krause
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

#ifndef PHYS_WALL_SHEAR_STRESS_ON_SURFACE_3D_H
#define PHYS_WALL_SHEAR_STRESS_ON_SURFACE_3D_H

#include<cmath>
#include<vector>

#include "analyticalF.h"
#include "functors/lattice/blockBaseF3D.h"
#include "functors/lattice/superBaseF3D.h"
#include "io/stlReader.h"
#include "core/unitConverter.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class SuperLatticeStress3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeStress3D(SuperLattice<T,DESCRIPTOR>& sLattice);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeStress3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeStress3D(BlockLattice<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

template<typename T, typename DESCRIPTOR>
SuperLatticeStress3D<T, DESCRIPTOR>::SuperLatticeStress3D(
  SuperLattice<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 6)
{
  this->getName() = "stress";
  int maxC = this->_sLattice.getLoadBalancer().size();
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeStress3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeStress3D<T, DESCRIPTOR>::BlockLatticeStress3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 6)
{
  this->getName() = "stress";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeStress3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeStress(output);
  return true;
}

template <typename T, typename DESCRIPTOR>
class PhysWallShearStressOnSurface3D final: public AnalyticalF<3,T,T> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  AnalyticalF3D<T,T>& _densityF;
  AnalyticalF3D<T,T>& _stressF;
  STLreader<T>& _stlReader;
  T _physFactor;

public:
  PhysWallShearStressOnSurface3D(const UnitConverter<T,DESCRIPTOR>& converter,
                   AnalyticalF3D<T,T>& densityF,
                   AnalyticalF3D<T,T>& stressF,
                   STLreader<T>& stlReader);

  bool operator() (T output[], const T physR[]) override;

};
} // namespace olb

#endif
