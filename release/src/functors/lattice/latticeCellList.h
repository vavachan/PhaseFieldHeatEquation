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

#ifndef LATTICE_CELL_LIST_H
#define LATTICE_CELL_LIST_H

namespace olb {

template <typename T, typename DESCRIPTOR, typename U=bool>
class SuperLatticeCellList final : public SuperLatticeF<T,DESCRIPTOR> {
public:
  SuperLatticeCellList(SuperLattice<T,DESCRIPTOR>& sLattice,
                       std::vector<LatticeR<4>>& cellList);
};

template <typename T, typename DESCRIPTOR, typename U=bool>
class BlockLatticeCellList final : public BlockLatticeF<T,DESCRIPTOR> {
private:
  BlockData<DESCRIPTOR::d,T,U> _blockData;
public:
  BlockLatticeCellList(BlockLattice<T,DESCRIPTOR>& blockLattice, int globiC,
                      std::vector<LatticeR<4>>& cellList);
  bool operator() (T output[], const int input[]) override;

};

}
#endif
