/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef FD_BOUNDARIES_2D_HH
#define FD_BOUNDARIES_2D_HH

#include "boundaryPostProcessors2D.h"

#include "core/util.h"
#include "utilities/finiteDifference2D.h"

#include "dynamics/dynamics.h"
#include "dynamics/lbm.h"

namespace olb {

///////////  StraightFdBoundaryProcessor2D ///////////////////////////////////

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightFdBoundaryProcessor2D<T, DESCRIPTOR, direction,
                                   orientation>::apply(CELL& cell)
{
  using namespace olb::util::tensorIndices2D;
  using V = typename CELL::value_t;

  V dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
  V rho, u[DESCRIPTOR::d];

  auto& dynamics = cell.getDynamics();

  cell.computeRhoU(rho, u);

  interpolateGradients<0>(cell, dx_u);
  interpolateGradients<1>(cell, dy_u);

  V dx_ux = dx_u[0];
  V dy_ux = dy_u[0];
  V dx_uy = dx_u[1];
  V dy_uy = dy_u[1];
  V omega =
      dynamics.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
  V sToPi = -rho / descriptors::invCs2<V, DESCRIPTOR>() / omega;
  V pi[util::TensorVal<DESCRIPTOR>::n];
  pi[xx] = (V)2 * dx_ux * sToPi;
  pi[yy] = (V)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  V fEq[DESCRIPTOR::q] { };
  dynamics.computeEquilibrium(cell, rho, u, fEq);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = fEq[iPop] + equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
  }
}

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <int deriveDirection, typename CELL, typename V>
void StraightFdBoundaryProcessor2D<T, DESCRIPTOR, direction, orientation>::
    interpolateGradients(CELL& cell, V velDeriv[DESCRIPTOR::d]) const
{
  fd::DirectedGradients2D<V, DESCRIPTOR, direction, orientation,
                          direction ==
                              deriveDirection>::interpolateVector(velDeriv,
                                                                  cell);
}

////////  StraightConvectionBoundaryProcessor2D ////////////////////////////////

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightConvectionBoundaryProcessor2D<T, DESCRIPTOR, direction,
                                           orientation>::initialize(CELL& cell)
{
  constexpr auto missing =
      util::populationsContributingToVelocity<DESCRIPTOR, direction,
                                              -orientation>();
  auto prevCell = cell.template getFieldPointer<PREV_CELL>();
  for (unsigned i = 0; i < missing.size(); ++i) {
    prevCell[i] = cell[missing[i]];
  }
}

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightConvectionBoundaryProcessor2D<T, DESCRIPTOR, direction,
                                           orientation>::apply(CELL& cell)
{
  using V = typename CELL::value_t;
  constexpr auto missing =
      util::populationsContributingToVelocity<DESCRIPTOR, direction,
                                              -orientation>();

  auto prevCell = cell.template getField<PREV_CELL>();

  for (unsigned i = 0; i < missing.size(); ++i) {
    cell[missing[i]] = prevCell[i];
  }

  V rho0, u0[2];
  V rho1, u1[2];
  V rho2, u2[2];

  cell.computeRhoU(rho0, u0);

  static_assert(direction == 0 || direction == 1
                ,"Direction must be one of 0 or 1 in 2D");
  if constexpr (direction == 0) {
    cell.neighbor({-orientation, 0}).computeRhoU(rho1, u1);
    cell.neighbor({-orientation * 2, 0}).computeRhoU(rho2, u2);
  }
  else if constexpr (direction == 1) {
    cell.neighbor({0, -orientation}).computeRhoU(rho1, u1);
    cell.neighbor({0, -orientation * 2}).computeRhoU(rho2, u2);
  }

  V uDelta[2];
  V uAverage = rho0 * u0[direction];

  uDelta[0] =
      -uAverage * 0.5 * (3 * rho0 * u0[0] - 4 * rho1 * u1[0] + rho2 * u2[0]);
  uDelta[1] =
      -uAverage * 0.5 * (3 * rho0 * u0[1] - 4 * rho1 * u1[1] + rho2 * u2[1]);

  for (unsigned i = 0; i < missing.size(); ++i) {
    auto iPop = missing[i];
    prevCell[i] =
        cell[iPop] + descriptors::invCs2<V, DESCRIPTOR>() *
                         descriptors::t<V, DESCRIPTOR>(iPop) *
                         (uDelta[0] * descriptors::c<DESCRIPTOR>(iPop, 0) +
                          uDelta[1] * descriptors::c<DESCRIPTOR>(iPop, 1));
  }

  cell.template setField<PREV_CELL>(prevCell);
}

////////  SlipBoundaryProcessor2D ////////////////////////////////

template <typename T, typename DESCRIPTOR>
SlipBoundaryProcessor2D<T, DESCRIPTOR>::SlipBoundaryProcessor2D(
    int x0_, int x1_, int y0_, int y1_, int discreteNormalX,
    int discreteNormalY)
    : x0(x0_)
    , x1(x1_)
    , y0(y0_)
    , y1(y1_)
{
  this->getName() = "SlipBoundaryProcessor2D";
  OLB_PRECONDITION(x0 == x1 || y0 == y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop, 0) * discreteNormalX +
            descriptors::c<DESCRIPTOR>(iPop, 1) * discreteNormalY <
        0) {
      // std::cout << "-----" <<s td::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX * discreteNormalY == 0) {
        mult = 2;
      }
      mirrorDirection0 =
          (descriptors::c<DESCRIPTOR>(iPop, 0) -
           mult *
               (descriptors::c<DESCRIPTOR>(iPop, 0) * discreteNormalX +
                descriptors::c<DESCRIPTOR>(iPop, 1) * discreteNormalY) *
               discreteNormalX);
      mirrorDirection1 =
          (descriptors::c<DESCRIPTOR>(iPop, 1) -
           mult *
               (descriptors::c<DESCRIPTOR>(iPop, 0) * discreteNormalX +
                descriptors::c<DESCRIPTOR>(iPop, 1) * discreteNormalY) *
               discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q;
           reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop], 0) ==
                mirrorDirection0 &&
            descriptors::c<DESCRIPTOR>(reflectionPop[iPop], 1) ==
                mirrorDirection1) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template <typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T, DESCRIPTOR>::processSubDomain(
    BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_,
    int y1_)
{
  int newX0, newX1, newY0, newY1;
  if (util::intersect(x0, x1, y0, y1, x0_, x1_, y0_, y1_, newX0, newX1, newY0,
                      newY1)) {

    int iX;
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (iX = newX0; iX <= newX1; ++iX) {
      for (int iY = newY0; iY <= newY1; ++iY) {
        for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
          if (reflectionPop[iPop] != 0) {
            //do reflection
            blockLattice.get(iX, iY)[iPop] =
                blockLattice.get(iX, iY)[reflectionPop[iPop]];
          }
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T, DESCRIPTOR>::process(
    BlockLattice<T, DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  SlipBoundaryProcessorGenerator2D ////////////////////////////////

template <typename T, typename DESCRIPTOR>
SlipBoundaryProcessorGenerator2D<
    T, DESCRIPTOR>::SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_,
                                                     int y1_,
                                                     int discreteNormalX_,
                                                     int discreteNormalY_)
    : PostProcessorGenerator2D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_)
    , discreteNormalX(discreteNormalX_)
    , discreteNormalY(discreteNormalY_)
{}

template <typename T, typename DESCRIPTOR>
PostProcessor2D<T, DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>::generate() const
{
  return new SlipBoundaryProcessor2D<T, DESCRIPTOR>(
      this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template <typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T, DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>::clone() const
{
  return new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(
      this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

/////////// OuterVelocityCornerProcessor2D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, int xNormal, int yNormal>
template <concepts::DynamicCell CELL>
void OuterVelocityCornerProcessor2D<T, DESCRIPTOR, xNormal, yNormal>::apply(
    CELL& cell)
{
  using namespace olb::util::tensorIndices2D;
  using V = typename CELL::value_t;

  V rho10 = cell.neighbor({-1 * xNormal, -0 * yNormal}).computeRho();
  V rho01 = cell.neighbor({-0 * xNormal, -1 * yNormal}).computeRho();

  V rho20 = cell.neighbor({-2 * xNormal, -0 * yNormal}).computeRho();
  V rho02 = cell.neighbor({-0 * xNormal, -2 * yNormal}).computeRho();

  V rho = (V)2 / (V)3 * (rho01 + rho10) - (V)1 / (V)6 * (rho02 + rho20);

  V dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
  fd::DirectedGradients2D<V, DESCRIPTOR, 0, xNormal, true>::interpolateVector(
      dx_u, cell);
  fd::DirectedGradients2D<V, DESCRIPTOR, 1, yNormal, true>::interpolateVector(
      dy_u, cell);
  V dx_ux = dx_u[0];
  V dy_ux = dy_u[0];
  V dx_uy = dx_u[1];
  V dy_uy = dy_u[1];

  auto& dynamics = cell.getDynamics();
  V     omega =
      dynamics.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());

  V sToPi = -rho / descriptors::invCs2<V, DESCRIPTOR>() / omega;
  V pi[util::TensorVal<DESCRIPTOR>::n];
  pi[xx] = (V)2 * dx_ux * sToPi;
  pi[yy] = (V)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  V u[DESCRIPTOR::d];
  cell.computeU(u);
  V fEq[DESCRIPTOR::q] { };
  dynamics.computeEquilibrium(cell, rho, u, fEq);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = fEq[iPop] + equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
  }
}

}  // namespace olb


#endif
