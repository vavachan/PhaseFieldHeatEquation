/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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

/** \file
 * Classes to provide access to a reacting species - may be a LB density, or an external field.
 *  -- generic implementation
 */
#ifndef REACTING_SPECIES_2D_HH
#define REACTING_SPECIES_2D_HH

namespace olb {


///////////////////////////////////// class ReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::ReactingSpecies2D(T stoichioCoeff)
  : _stoichioCoeff(stoichioCoeff)
{}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
T ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::getStoichioCoeff()
{
  return _stoichioCoeff;
}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
T ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::getField(BlockStructureD<2>* blockStructure, int iX, int iY)
{
  return (static_cast<IMPL*>(this))->getField(blockStructure, iX, iY);
}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
T ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::getSource(BlockStructureD<2>* blockStructure, int iX, int iY)
{
  return static_cast<BlockLattice<T,DESCRIPTOR>*>(blockStructure)->get(iX, iY).template getFieldPointer<SOURCE>()[0];
}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
void ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::incrementSource(BlockStructureD<2>* blockStructure, T val, int iX, int iY)
{
  static_cast<BlockLattice<T,DESCRIPTOR>*>(blockStructure)->get(iX, iY).template setField<SOURCE>(val + getSource(blockStructure, iX, iY));
}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
void ReactingSpecies2D<T,DESCRIPTOR,SOURCE,IMPL>::resetSource(BlockStructureD<2>* blockStructure, int iX, int iY)
{
  static_cast<BlockLattice<T,DESCRIPTOR>*>(blockStructure)->get(iX, iY).template setField<SOURCE>(T());
}


///////////////////////////////////// class FiniteDifferenceReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE, typename FIELD>
FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,SOURCE,FIELD>::
FiniteDifferenceReactingSpecies2D(T stoichioCoeff, std::size_t& iT)
  : ReactingSpecies2D<T,DESCRIPTOR,SOURCE,FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,SOURCE,FIELD>>(stoichioCoeff),
    _iT(iT)
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
  static_assert(DESCRIPTOR::template size<SOURCE>() == 1, "SOURCE must have size 1.");
}

template <typename T, typename DESCRIPTOR, typename SOURCE, typename FIELD>
T FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,SOURCE,FIELD>::getField(BlockStructureD<2>* blockStructure, int iX, int iY)
{
  return *fd::accessNew<T,FIELD>(static_cast<BlockLattice<T,DESCRIPTOR>*>(blockStructure)->get(iX, iY), this->_iT);
}


///////////////////////////////////// class LatticeBoltzmannReactingSpecies2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, typename SOURCE>
LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>::
LatticeBoltzmannReactingSpecies2D(T stoichioCoeff)
  : ReactingSpecies2D<T,DESCRIPTOR,SOURCE,LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>>(stoichioCoeff)
{}

template <typename T, typename DESCRIPTOR,typename SOURCE>
T LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>::getField(BlockStructureD<2>* blockStructure, int iX, int iY)
{
  return static_cast<BlockLattice<T,DESCRIPTOR>*>(blockStructure)->get(iX, iY).computeRho();
}


}  // namespace olb

#endif
