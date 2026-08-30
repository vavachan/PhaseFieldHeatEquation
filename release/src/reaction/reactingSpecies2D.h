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
 *  -- header file
 */
#ifndef REACTING_SPECIES_2D_H
#define REACTING_SPECIES_2D_H

namespace olb {

/*
 * This class provides a unified set of reading-writing methods to reacting specie fields.
 * The methods work the same to the user regardless of whether the specie field is solved through explicit Finite-Difference
 * (thereby being defined as an external field), or through Lattice-Boltzmann (thereby being defined as the 0-th momenum of
 * Lattice-Boltzmann one-particle density functions). This is achieved through CRTP.
 * The Finite-Difference implementation is flexible enough to allow the external field being defined on both the main lattice,
 * and on a coupled lattice. Conversely, the Lattice-Boltzmann implementation allows only a coupled lattice.
 */
template <typename T, typename DESCRIPTOR, typename SOURCE, typename IMPL>
class ReactingSpecies2D {
public:
  ReactingSpecies2D(T stoichioCoeff);
  // Get the stoichiometric coefficient
  T getStoichioCoeff();
  // Read the field value at a input lattice site within the input BlockStructure.
  T getField(BlockStructureD<2>* blockStructure, int iX, int iY);
  // Read the source value at a input lattice site within the input BlockStructure.
  T getSource(BlockStructureD<2>* blockStructure, int iX, int iY);
  // Increment the source value at a input lattice site within the input BlockStructure.
  void incrementSource(BlockStructureD<2>* blockStructure, T val, int iX, int iY);
  // Resets the source value at a input lattice site within the input BlockStructure.
  void resetSource(BlockStructureD<2>* blockStructure, int iX, int iY);
protected:
  // Stoichiometric coefficient. If negative, the species is a reagent; if positive, a product.
  T _stoichioCoeff;
};

template <typename T, typename DESCRIPTOR, typename SOURCE, typename FIELD>
class FiniteDifferenceReactingSpecies2D final : public ReactingSpecies2D<T,DESCRIPTOR,SOURCE,FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,SOURCE,FIELD>> {
public:
  FiniteDifferenceReactingSpecies2D(T stoichioCoeff, std::size_t& iT);
  // Read the field value at a input lattice site within the input BlockLattice.
  T getField(BlockStructureD<2>* blockStructure, int iX, int iY);
protected:
  std::size_t& _iT;
};

template <typename T, typename DESCRIPTOR, typename SOURCE>
class LatticeBoltzmannReactingSpecies2D final : public ReactingSpecies2D<T,DESCRIPTOR,SOURCE,LatticeBoltzmannReactingSpecies2D<T,DESCRIPTOR,SOURCE>> {
public:
  LatticeBoltzmannReactingSpecies2D(T stoichioCoeff);
  // Read the field value at a input lattice site within the input BlockLattice.
  T getField(BlockStructureD<2>* blockStructure, int iX, int iY);
};

};  // namespace olb

#endif
