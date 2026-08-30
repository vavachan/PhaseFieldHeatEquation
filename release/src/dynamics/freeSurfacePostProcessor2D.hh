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

#ifndef FREE_SURFACE_POST_PROCESSOR_2D_HH
#define FREE_SURFACE_POST_PROCESSOR_2D_HH

#include "freeSurfacePostProcessor2D.h"
#include "core/blockLattice.h"

namespace olb {

/**
 * Free Surface Processor 1: Mass Flow
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL>
void FreeSurfaceMassFlowPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {

  using namespace olb::FreeSurface;

  // Reset CELL_FLAGS and TEMP_MASS_EXCHANGE here, as they are needed in the final post-processor step
  setCellFlags(cell, FreeSurface::Flags::None);
  Vector<T, DESCRIPTOR::q> zero_mass_exchange{};
  cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>(zero_mass_exchange);

  // Skip if the current cell is not an interface cell.
  if (!isCellType(cell, FreeSurface::Type::Interface)) {
    return;
  }

  // Check if the current interface cell has gas or fluid neighbours
  FreeSurface::NeighbourInfo nbrInfo = getNeighbourInfo(cell);
  bool hasNoGasNeighbours = !nbrInfo.has_gas_neighbours;
  bool hasNoFluidNeighbours = !nbrInfo.has_fluid_neighbours;
  bool isHealthyInterfaceCell = isHealthyInterface(cell);

  // The cell has only interface cells as neighbours
  if (hasNoGasNeighbours && hasNoFluidNeighbours) {
    hasNoGasNeighbours = false;
    hasNoFluidNeighbours = false;
    isHealthyInterfaceCell = true;
  }

  // The notation used here differs from that in N. Thuerey's dissertation (2007), specifically in Section 4.1.
  // This is because the free-surface post-processors are scheduled for execution during the PostStream stage,
  // which occurs after the collideAndStream step. However, mass exchange computations require access to
  // post-collision distribution functions at this stage.
  // Since the post-collision distribution functions are not directly available here, we rely on the fact that,
  // for the current cell, the post-collision f_(i) is equal to that of the corresponding neighbouring cell
  // post-stream, i.e., cell[i](*) = nbrCell[i].
  // Keep in mind that each cell has DESCRIPTOR::q neighbours, and each neighbor contains DESCRIPTOR::q
  // distribution functions.
  T mass_exchange = T(0);
  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    // Skip if the neighbour is a gas cell, i.e., no mass exchange
    if (isCellType(nbrCell, FreeSurface::Type::Gas)) { continue; }

    // Mass exchange with a liquid neighbour, i.e., computed uisng the difference between incoming and outgoing DFs
    if (isCellType(nbrCell, FreeSurface::Type::Fluid)) {
      mass_exchange += cell[descriptors::opposite<DESCRIPTOR>(iPop)] - nbrCell[iPop];
      continue;
    }

    // Check if the neighbour cell has gas or fluid neighbours
    FreeSurface::NeighbourInfo nbrNeighbourInfo = getNeighbourInfo(nbrCell);
    bool nbrHasNoGasNeighbours = !nbrNeighbourInfo.has_gas_neighbours;
    bool nbrHasNoFluidNeighbours = !nbrNeighbourInfo.has_fluid_neighbours;
    bool nbrIsHealthyInterfaceCell = isHealthyInterface(nbrCell);

    // Neighbour cell has only interface cells as neighbours
    if (nbrHasNoGasNeighbours && nbrHasNoFluidNeighbours) {
      nbrHasNoGasNeighbours = false;
      nbrHasNoFluidNeighbours = false;
      nbrIsHealthyInterfaceCell = true;
    }

    T pdf_difference = T(0);
    if (isCellType(nbrCell, FreeSurface::Type::Interface)) {
      // Simple mass exchange if both current cell and neighbour cell are healthy interfaces
      if (isHealthyInterfaceCell && nbrIsHealthyInterfaceCell) {
        pdf_difference = cell[descriptors::opposite<DESCRIPTOR>(iPop)] - nbrCell[iPop];
      }
      else {
        // A healthy interface cell with a neighbour that has no gas neighbours itself,
        // or an interface cell with no liquid neighbours, but its neighbor has a liquid neighbour.
        // make this cell empty, ref. N. Thuerey dissertation, 2007.
        if ((isHealthyInterfaceCell && nbrHasNoGasNeighbours) || (hasNoFluidNeighbours && !nbrHasNoFluidNeighbours)) {
          pdf_difference = -nbrCell[iPop];
        }
        else {
          // A healthy interface cell with a neighbour that has no fluid neighbours itself,
          // or an interface cell with no gas neighbours, but its neighbour has a gas neighbour.
          // make the neighbour cell empty, ref. N. Thuerey dissertation, 2007.
          if ((isHealthyInterfaceCell && nbrHasNoFluidNeighbours) || (hasNoGasNeighbours && !nbrHasNoGasNeighbours)) {
            pdf_difference = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
          }
          else {
            // An interface cell with no liquid neighbors, whose neighbor also has no liquid neighbors,
            // or an interface cell with no gas neighbors, whose neighbor also has no liquid neighbors.
            // ref. N. Thuerey dissertation, 2007.
            if ((hasNoFluidNeighbours && nbrHasNoFluidNeighbours) || (hasNoGasNeighbours && nbrHasNoGasNeighbours)) {
              pdf_difference = cell[descriptors::opposite<DESCRIPTOR>(iPop)] - nbrCell[iPop];
            }
          }
        }
      }
    }

    const T epsilon_average = T(0.5) * (getClampedEpsilon(cell) + getClampedEpsilon(nbrCell));
    mass_exchange += pdf_difference * epsilon_average;
  }

  // Update the mass of interface cell
  const auto mass = cell.template getField<FreeSurface::MASS>();
  cell.template setField<FreeSurface::MASS>(mass + mass_exchange);
}

/**
 * Free Surface Processor 2-3: Interface Reconstruction
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL, typename PARAMETERS>
void FreeSurfaceInterfaceReconstructionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell, PARAMETERS& params) {

  using namespace olb::FreeSurface;

  const bool drop_isolated_cells = params.template get<FreeSurface::DROP_ISOLATED_CELLS>();
  const auto transition          = params.template get<FreeSurface::TRANSITION>();
  const auto lonely_threshold    = params.template get<FreeSurface::LONELY_THRESHOLD>();
  const bool has_surface_tension = params.template get<FreeSurface::HAS_SURFACE_TENSION>();
  const auto surface_tension     = params.template get<FreeSurface::SURFACE_TENSION_PARAMETER>();

  if (!isCellType(cell, FreeSurface::Type::Interface)) {
    return;
  }

  FreeSurface::NeighbourInfo nbrInfo = getNeighbourInfo(cell);

  // Adjust gas density, using Laplace pressure but without support for the bubble model
  T deltaRho = T(0);
  if (has_surface_tension && nbrInfo.has_gas_neighbours) {
    const T curvature = computeCurvature(cell);
    deltaRho = T(6) * surface_tension * curvature;
  }
  const T bubble_density = T(1);
  const T gas_density = bubble_density - deltaRho;

  // Replace the distribution functions of the gas cell using the updated density.
  // A race condition is not expected here since only the distribution functions
  // of an interface cell are modified, while only the distribution functions of
  // a neighboring gas cell are read.
  Vector<T, DESCRIPTOR::q> dfs{};
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    // Distribution functions of an interface cell are reconstructed only for the missing directions,
    // i.e., incoming DFs from gas cells, using the pressure anti bounce back boundary condition.
    // Note: the notation used here differs from that in N. Thuerey's dissertation (2007), refer to the
    // explanation given above.
    if (isCellType(nbrCell, FreeSurface::Type::Gas)) {
      constexpr bool usePrecise = true;
      Vector<T, DESCRIPTOR::d> prev_velocity = cell.template getField<FreeSurface::PREVIOUS_VELOCITY>();

      // (1) Precise formulation of the pressure anti bounce back boundary condition: feq_(i) + feq_(-i) - nbr(i)
      // This is mathematically identical to the regular formulation, but the linear term of the 2nd order
      // equilibrium equation is manually eliminated due to its cancellation during summation, i.e., c_(-i).u = -(c_(i).u).
      // This formulation offers consistent floating-point rounding behaviour during arithmetic operations.
      if constexpr (usePrecise) {
        const T uSqr = util::normSqr<T, DESCRIPTOR::d>(prev_velocity);
        Vector<T, DESCRIPTOR::d> direction{};
        for (int i = 0; i < DESCRIPTOR::d; ++i) { direction[i] = descriptors::c<DESCRIPTOR>(iPop, i); }
        const T c_u = util::dotProduct(direction, prev_velocity);

        dfs[descriptors::opposite<DESCRIPTOR>(iPop)] =
            T(2) * gas_density * descriptors::t<T, DESCRIPTOR>(iPop) * (T(1) + T(4.5) * c_u * c_u - T(1.5) * uSqr)
          - T(2) * descriptors::t<T, DESCRIPTOR>(iPop)
          - nbrCell[iPop];
      }
      // (2) Regular formulation of the pressure anti bounce back boundary condition: feq_(i) + feq_(-i) - nbr(i)
      // This formulation might suffer from inconsistent floating-point rounding behaviour during arithmetic
      // operations, where the expression c_(-i).u = -(c_(i).u) could become untrue.
      else {
        T u[DESCRIPTOR::d] = {prev_velocity[0], prev_velocity[1]};
        dfs[descriptors::opposite<DESCRIPTOR>(iPop)] =
            equilibrium<DESCRIPTOR>::secondOrder(iPop, gas_density, u)
          + equilibrium<DESCRIPTOR>::secondOrder(descriptors::opposite<DESCRIPTOR>(iPop), gas_density, u)
          - nbrCell[iPop];
      }
    }
    else {
      dfs[descriptors::opposite<DESCRIPTOR>(iPop)] =
        cell[descriptors::opposite<DESCRIPTOR>(iPop)];
    }
  }

  // Update interafce cell populations using the computed dfs in previous step.
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) { cell[iPop] = dfs[iPop]; }

  // Using the updated mass and density, find and set appropriate flag on the current interface cell
  const auto rho = cell.computeRho();
  const auto mass = cell.template getField<FreeSurface::MASS>();

  // An interface cell is flagged as toGas if its volume fraction is below the threshold
  if (mass < -transition * rho) {
    setCellFlags(cell, FreeSurface::Flags::ToGas);
    return;
  }

  // An interface cell is flagged as toFluid if its volume fraction is above the (1.0 + threshold)
  if (mass > (T(1) + transition) * rho) {
    setCellFlags(cell, FreeSurface::Flags::ToFluid);
    return;
  }

  // An interface cell with no fluid neighbours is flagged as toGas:
  // (1) If its volume fraction is below the (lonely threshold) and has interface neighbours (not isolated)
  // (2) If it has no interface neighbours (isolated) and drop_isolated_cells is set
  if (!nbrInfo.has_fluid_neighbours) {
    if (mass < lonely_threshold * rho && nbrInfo.interface_neighbours != 0) {
      setCellFlags(cell, FreeSurface::Flags::ToGas);
      return;
    }

    if (nbrInfo.interface_neighbours == 0 && drop_isolated_cells) {
      setCellFlags(cell, FreeSurface::Flags::ToGas);
      return;
    }
  }

  // An interface cell with no gas neighbours is flagged as toFluid:
  // (1) If its volume fraction is above the (1.0 - lonely threshold) and has interface neighbours (not isolated)
  // (2) If it has no interface neighbours (isolated) and drop_isolated_cells is set
  if (!nbrInfo.has_gas_neighbours) {
    if (mass > (T(1) - lonely_threshold) * rho && nbrInfo.interface_neighbours != 0) {
      setCellFlags(cell, FreeSurface::Flags::ToFluid);
      return;
    }

    if (nbrInfo.interface_neighbours == 0 && drop_isolated_cells) {
      setCellFlags(cell, FreeSurface::Flags::ToFluid);
      return;
    }
  }
}

/*
 * Free Surface Processor 4: ToFluid
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL>
void FreeSurfaceToFluidCellConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;

  // Convert a gas cell to interface if it is in the neighborhood of a toFluid neighbour cell
  if (isCellType(cell, FreeSurface::Type::Gas) && hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)) {
    // Initialize distribution functions for a cell converted to interface from gas, using the
    // equilibrium refilling method, i.e., the average velocity and density of neighboring cells.
    setCellFlags(cell, FreeSurface::Flags::NewInterface);
    T rho_average = T(0);
    T u_average[DESCRIPTOR::d] = {T(0), T(0)};
    std::uint8_t count = std::uint8_t(0);

    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

      // Only consider neighbouring cells which are fluid or interface cells.
      // Note: newly converted interface cells are not considered here, i.e., CELL_TYPE is not updated yet.
      if (isCellType(nbrCell, FreeSurface::Type::Fluid) || isCellType(nbrCell, FreeSurface::Type::Interface)) {
        T rho_tmp = T(0);
        T u_tmp[DESCRIPTOR::d] = {T(0), T(0)};
        ++count;
        nbrCell.computeRhoU(rho_tmp, u_tmp);
        rho_average += rho_tmp;
        for (std::size_t i = 0; i < DESCRIPTOR::d; ++i) { u_average[i] += u_tmp[i]; }
      }
    }

    if (count > std::uint8_t(0)) {
      rho_average /= static_cast<T>(count);
      for (std::size_t i = 0; i < DESCRIPTOR::d; ++i) { u_average[i] /= static_cast<T>(count); }
    }
    else {
      // If no valid neighbouring cells are found
      rho_average = T(1);
      for (std::size_t i = 0; i < DESCRIPTOR::d; ++i) { u_average[i] = T(0); }
    }

    cell.iniEquilibrium(rho_average, u_average);
  }

  // If an interface cell with a ToGas flag has a neighbouring ToFluid cell, unset the ToGas flag
  if (hasCellFlags(cell, FreeSurface::Flags::ToGas)) {
    if (hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)) {
      setCellFlags(cell, FreeSurface::Flags::None);
    }
  }
}

/**
 * Free Surface Processor 5: ToGas
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL>
void FreeSurfaceToGasCellConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;

  if (isCellType(cell, FreeSurface::Type::Fluid)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

      // Convert a liquid cell to interface if it is in the neighborhood of a toGas cell
      if (hasCellFlags(nbrCell, FreeSurface::Flags::ToGas)) {
        setCellFlags(cell, FreeSurface::Flags::NewInterface);
        T rho = cell.computeRho();
        cell.template setField<FreeSurface::MASS>(rho);

        // Cell already flagged as NewInterface, skip the remaining neighbours.
        break;
      }
    }
  }
}

/**
 * Free Surface Processor 6: Mass Excess
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL>
void FreeSurfaceMassExcessPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  using namespace olb::FreeSurface;

  // Note that EPSILON cannot be set in this operator because it is needed for the normal computation,
  // thus MASS is set here and EPSILON is set by the next operator.
  // Skip if the current cell is not an interface cell.
  if (!isCellType(cell, FreeSurface::Type::Interface)) {
    return;
  }

  T mass_excess = T(0);
  // If an interface cell is flagged as toGas,
  if (hasCellFlags(cell, FreeSurface::Flags::ToGas)) {
    // Get mass of the current interface cell
    const T mass = cell.template getField<FreeSurface::MASS>();
    mass_excess = mass;
    cell.template setField<FreeSurface::MASS>(T(0));
  }
  // If an interface cell is flagged as toFluid,
  else if (hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
    // Get density and mass of the current interface cell
    const T rho = cell.computeRho();
    const T mass = cell.template getField<FreeSurface::MASS>();
    mass_excess = mass - rho;
    cell.template setField<FreeSurface::MASS>(rho);
  }
  else { return; }

  // Redistribute the excess mass among the neighbouring interface cells, note that
  // if no neighbouring interface cells are found, the excess mass is lost or gained
  // by the current cell.
  // TODO: Thuerey Paper says we can't use new interface cells or flagged cells
  std::uint8_t oldIntefaceNbrs = std::uint8_t(0);
  std::uint8_t newIntefaceNbrs = std::uint8_t(0);
  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    // Note that CELL_TYPE is not updated yet, and if no flag is set for this neighbour cell, it is an old interface cell
    // This means that this neighbour cell will remain an interface cell for the next time step
    if (isCellType(nbrCell, FreeSurface::Type::Interface) && hasCellFlags(nbrCell, FreeSurface::Flags::None)) {
      ++oldIntefaceNbrs;
    }
    // If the neighbour cell is flagged as Newinterface, it is a new interface cell
    else if (hasCellFlags(nbrCell, FreeSurface::Flags::NewInterface)) {
      ++newIntefaceNbrs;
    }
  }

  // Return if no interface cell is available to distribute the excess mass, i.e., the mass is lost or gained.
  const std::uint8_t intefaceNbrs = newIntefaceNbrs + oldIntefaceNbrs;
  if (intefaceNbrs == std::uint8_t(0)) { return; }

  // Distribution of excess mass among the neighbouring interface cells, weighted or evenly.
  bool enableAllInterfaces = true;
  constexpr bool weightedDistribution = true;
  Vector<T, DESCRIPTOR::q> weights{};
  T weightsSum = T(0);

  // (1) Distribute excess mass among the neighbouring interface cells, weighted by normal vector.
  if constexpr (weightedDistribution) {
    // Compute the weights for the mass excess distribution using normal vector of the current interface cell.
    auto normal = computeInterfaceNormal(cell);
    computeMassExcessWeights(cell, normal, enableAllInterfaces, weights);

    // Calculate sum of the computed weights
    for (const auto& weight : weights) { weightsSum += weight; }

    if (util::fabs(weightsSum) < zeroThreshold<T>()) {
      // (a) If no old interface cell is available in normal direction,
      //     fall back to weighted all interface cells distribution model.
      if (!enableAllInterfaces) {
        enableAllInterfaces = true;
        weightsSum = T(0);
        computeMassExcessWeights(cell, normal, enableAllInterfaces, weights);
        for (const auto& weight : weights) { weightsSum += weight; }

        // (b) If no interface cell is available in normal direction,
        //     fall back to evenly all interface cells distribution model.
        if (util::fabs(weightsSum) < zeroThreshold<T>()) {
          for (auto& weight : weights) { weight = T(1); }
          weightsSum = T(intefaceNbrs);
        }
      }
      // (c) If all interface cells are selected but none is available in normal direction,
      //     fall back to evenly all interface cells distribution model.
      else {
        for (auto& weight : weights) { weight = T(1); }
        weightsSum = T(intefaceNbrs);
      }
    }
  }
  // (2) Distribute excess mass among the neighbouring interface cells, evenly.
  else {
    // (a) Use old neighbouring interface cells to distribute excess mass.
    if (!enableAllInterfaces && oldIntefaceNbrs > std::uint8_t(0)) {
      for (auto& weight : weights) { weight = T(1); }
      weightsSum = T(oldIntefaceNbrs);
    }
    // (b) Use all neighbouring interface cells to distribute excess mass.
    else {
      enableAllInterfaces = true;
      for (auto& weight : weights) { weight = T(1); }
      weightsSum = T(intefaceNbrs);
    }
  }

  Vector<T, DESCRIPTOR::q> mass_exchange{};
  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if (isCellType(nbrCell, FreeSurface::Type::Interface) && hasCellFlags(nbrCell, FreeSurface::Flags::None)) {
      mass_exchange[iPop] = mass_excess * weights[iPop] / weightsSum;
    }
    else if (hasCellFlags(nbrCell, FreeSurface::Flags::NewInterface) && enableAllInterfaces) {
      mass_exchange[iPop] = mass_excess * weights[iPop] / weightsSum;
    }
    else {
      mass_exchange[iPop] = T(0);
    }
  }

  cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>(mass_exchange);
}

/**
 * Free Surface Processor 7: Finalize Conversion
 */
template<typename T, typename DESCRIPTOR>
template<typename CELL>
void FreeSurfaceFinalizeConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {

  using namespace olb::FreeSurface;

  // Convert flagged cells to appropriate cell types
  if (hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
    setCellType(cell, FreeSurface::Type::Fluid);
    cell.template setField<FreeSurface::EPSILON>(T(1));
  }
  else if (hasCellFlags(cell, FreeSurface::Flags::ToGas)) {
    setCellType(cell, FreeSurface::Type::Gas);
    cell.template setField<FreeSurface::EPSILON>(T(0));
  }
  else if (hasCellFlags(cell, FreeSurface::Flags::NewInterface)) {
    setCellType(cell, FreeSurface::Type::Interface);
  }
  else {
    // if (hasCellFlags(cell, FreeSurface::Flags::None) == true),
    // no further action is needed, as no conversion is triggered.
  }

  // Mass excess is distributed to old and new interface cells only,
  // thus skip the non-interface cells.
  if (!isCellType(cell, FreeSurface::Type::Interface)) {
    return;
  }

  // Collection of mass excess in a pulling step for the interface cells
  T collected_excess = T(0);
  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if (hasCellFlags(nbrCell, FreeSurface::Flags::ToFluid | FreeSurface::Flags::ToGas)) {
      auto tempMassExchange = nbrCell.template getField<FreeSurface::TEMP_MASS_EXCHANGE>();
      collected_excess += tempMassExchange[descriptors::opposite<DESCRIPTOR>(iPop)];
    }
  }

  T mass_tmp = cell.template getField<FreeSurface::MASS>();
  mass_tmp += collected_excess;
  T rho = T(0);
  T u_tmp[DESCRIPTOR::d] = {T(0), T(0)};
  cell.computeRhoU(rho, u_tmp);
  Vector<T,DESCRIPTOR::d> u_vel{u_tmp[0], u_tmp[1]};
  cell.template setField<FreeSurface::EPSILON>(mass_tmp / rho);
  cell.template setField<FreeSurface::MASS>(mass_tmp);
  cell.template setField<FreeSurface::PREVIOUS_VELOCITY>(u_vel);
}

// Setup
template<typename T, typename DESCRIPTOR>
FreeSurface2DSetup<T,DESCRIPTOR>::FreeSurface2DSetup(SuperLattice<T, DESCRIPTOR>& sLattice)
:
  _sLattice(sLattice)
{}

template<typename T, typename DESCRIPTOR>
void FreeSurface2DSetup<T,DESCRIPTOR>::addPostProcessor() {
  _sLattice.template addPostProcessor<FreeSurface::Stage0>(
    meta::id<FreeSurfaceMassFlowPostProcessor2D<T,DESCRIPTOR>>{});
  _sLattice.template addPostProcessor<FreeSurface::Stage1>(
    meta::id<FreeSurfaceInterfaceReconstructionPostProcessor2D<T,DESCRIPTOR>>{});
  _sLattice.template addPostProcessor<FreeSurface::Stage2>(
    meta::id<FreeSurfaceToFluidCellConversionPostProcessor2D<T,DESCRIPTOR>>{});
  _sLattice.template addPostProcessor<FreeSurface::Stage3>(
    meta::id<FreeSurfaceToGasCellConversionPostProcessor2D<T,DESCRIPTOR>>{});
  _sLattice.template addPostProcessor<FreeSurface::Stage4>(
    meta::id<FreeSurfaceMassExcessPostProcessor2D<T,DESCRIPTOR>>{});
  _sLattice.template addPostProcessor<FreeSurface::Stage5>(
    meta::id<FreeSurfaceFinalizeConversionPostProcessor2D<T,DESCRIPTOR>>{});

  {
    // Communicate fields: EPSILON, CELL_TYPE, POPULATION, and PREVIOUS_VELOCITY
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage0());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::EPSILON>();
    communicator.template requestField<FreeSurface::CELL_TYPE>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.template requestField<FreeSurface::PREVIOUS_VELOCITY>();
    communicator.exchangeRequests();
  }

  {
    // Communicate fields: MASS and PREVIOUS_VELOCITY
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage1());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::MASS>();
    communicator.template requestField<FreeSurface::PREVIOUS_VELOCITY>();
    communicator.exchangeRequests();
  }

  {
    // Communicate fields: CELL_FLAGS and POPULATION
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage2());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.exchangeRequests();
  }

  {
    // Communicate fields: CELL_FLAGS
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage3());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate fields: CELL_FLAGS
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage4());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate fields: TEMP_MASS_EXCHANGE
    auto& communicator = _sLattice.getCommunicator(FreeSurface::Stage5());
    communicator.requestOverlap(_sLattice.getOverlap());
    communicator.template requestField<FreeSurface::TEMP_MASS_EXCHANGE>();
    communicator.exchangeRequests();
  }

  // Add custom tasks to be executed post stream
  _sLattice.template addCustomTask<stage::PostStream>([&]() {
    _sLattice.executePostProcessors(FreeSurface::Stage0());
    _sLattice.executePostProcessors(FreeSurface::Stage1());
    _sLattice.executePostProcessors(FreeSurface::Stage2());
    _sLattice.executePostProcessors(FreeSurface::Stage3());
    _sLattice.executePostProcessors(FreeSurface::Stage4());
    _sLattice.executePostProcessors(FreeSurface::Stage5());
  });
}

}
#endif
