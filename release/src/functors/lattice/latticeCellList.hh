/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef LATTICE_CELL_LIST_HH
#define LATTICE_CELL_LIST_HH

namespace olb {

template<typename T, typename DESCRIPTOR, typename U>
SuperLatticeCellList<T,DESCRIPTOR,U>::SuperLatticeCellList(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  std::vector<LatticeR<4>>& cellList)
  : SuperLatticeF<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cellListF";
  auto& loadBalancer = this->_sLattice.getLoadBalancer();
  int maxC = loadBalancer.size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCellList<T,DESCRIPTOR>(
      this->_sLattice.getBlock(iC),
      loadBalancer.glob(iC), cellList ));
  }
}

template<typename T, typename DESCRIPTOR, typename U>
BlockLatticeCellList<T,DESCRIPTOR,U>::BlockLatticeCellList(
  BlockLattice<T,DESCRIPTOR>& blockLattice, int globiC,
  std::vector<LatticeR<4>>& cellList)
  : BlockLatticeF<T,DESCRIPTOR>(blockLattice, 1),
    _blockData(static_cast<BlockStructure<DESCRIPTOR>>(blockLattice))
{
  this->getName() = "cellListF";
  for (auto latticeRlist : cellList){
    int globiCList = latticeRlist[0];
    if (globiCList==globiC){
      LatticeR<DESCRIPTOR::d> latticeR(
        latticeRlist[1],
        latticeRlist[2],
        latticeRlist[3]
      );
      auto& dataPoint = _blockData.get(latticeR);
      dataPoint = true;
    }
  }
}

template<typename T, typename DESCRIPTOR, typename U>
bool BlockLatticeCellList<T, DESCRIPTOR, U>::operator()(T output[], const int input[])
{
  if constexpr (std::is_same_v<U,bool>){
    _blockData(output,input);
    return bool(output[0]);
  }
}

}
#endif
