/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
 *                2024-2025 Danial Khazaeipoul
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

#ifndef FREE_SURFACE_POST_PROCESSOR_3D_H
#define FREE_SURFACE_POST_PROCESSOR_3D_H

#include "dynamics/freeSurfaceHelpers.h"
#include "core/postProcessing.h"
#include "core/blockLattice.h"
#include "core/superLattice.h"

namespace olb {

/**
 * Free Surface Processor 1
 * Mass Flow
 * Cleans up leftover flags from the previous simulation step.
 * This post processor is responsible for the calculation of exchange mass with the help of the distribution functions.
 * This whole step should be included in the collideAndStream step, though heavy modification of openlb would be necessary
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassFlowPostProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 2-3
 * Interface Reconstruction
 * Replaces incoming DFs by calculating equilibrium functions and using the laplace pressure to include surface tension.
 * Marks cells which may be changed at the last step.
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceInterfaceReconstructionPostProcessor3D {
public:
   using parameters = meta::list<
    FreeSurface::DROP_ISOLATED_CELLS,
    FreeSurface::TRANSITION,
    FreeSurface::LONELY_THRESHOLD,
    FreeSurface::HAS_SURFACE_TENSION,
    FreeSurface::SURFACE_TENSION_PARAMETER
  >;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 2;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& params) any_platform;

};

/*
 * Free Surface Processor 4
 * ToFluid
 * Converts cells to interface from gas if a neighbouring cell was converted to a fluid cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToFluidCellConversionPostProcessor3D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 3;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 5
 * ToGas
 * Converts cells to interface from fluid if a neighbouring cell was converted to a gas cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToGasCellConversionPostProcessor3D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 4;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 6
 * Calculates mass excess from the cell type conversions and distributes them to neighbouring interface cells
 * Keeps mass local if no neighbour exists until an interface reappears at this position
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassExcessPostProcessor3D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 5;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 7
 * Finishes up left over cell conversions and prepares the state for the next simulation step
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceFinalizeConversionPostProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 6;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/*
* Setup helper
*/
template<typename T, typename DESCRIPTOR>
class FreeSurface3DSetup {
public:
private:
  SuperLattice<T, DESCRIPTOR>& _sLattice;

  // SuperPostProcessors
  // Corresponding to the local block processors
public:
  FreeSurface3DSetup(SuperLattice<T, DESCRIPTOR>& sLattice);

  void addPostProcessor();
};

}
#endif
