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

#ifndef BOUNDARY_POST_PROCESSORS_3D_HH
#define BOUNDARY_POST_PROCESSORS_3D_HH

#include "boundaryPostProcessors3D.h"

#include "utilities/finiteDifference3D.h"
#include "core/util.h"

#include "dynamics/dynamics.h"
#include "dynamics/lbm.h"

namespace olb {

////////  PlaneFdBoundaryProcessor3D ///////////////////////////////////

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::apply(CELL& cell)
{
  using V = typename CELL::value_t;
  using namespace olb::util::tensorIndices3D;

  V dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d], dz_u[DESCRIPTOR::d];
  V rho, u[DESCRIPTOR::d];

  auto& dynamics = cell.getDynamics();

  cell.computeRhoU(rho,u);

  interpolateGradients<0>(cell, dx_u);
  interpolateGradients<1>(cell, dy_u);
  interpolateGradients<2>(cell, dz_u);

  V dx_ux = dx_u[0];
  V dy_ux = dy_u[0];
  V dz_ux = dz_u[0];
  V dx_uy = dx_u[1];
  V dy_uy = dy_u[1];
  V dz_uy = dz_u[1];
  V dx_uz = dx_u[2];
  V dy_uz = dy_u[2];
  V dz_uz = dz_u[2];
  V omega = dynamics.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
  V sToPi = - rho / descriptors::invCs2<V,DESCRIPTOR>() / omega;
  V pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (V)2 * dx_ux * sToPi;
  pi[yy] = (V)2 * dy_uy * sToPi;
  pi[zz] = (V)2 * dz_uz * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;
  pi[xz] = (dx_uz + dz_ux) * sToPi;
  pi[yz] = (dy_uz + dz_uy) * sToPi;

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
void PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::interpolateGradients(
  CELL& cell, V velDeriv[DESCRIPTOR::d]) const
{
  fd::DirectedGradients3D<V,DESCRIPTOR,direction,orientation,deriveDirection,direction==deriveDirection>
    ::interpolateVector(velDeriv, cell);
}

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::initialize(CELL& cell)
{
  constexpr auto missing = util::populationsContributingToVelocity<DESCRIPTOR,direction,-orientation>();
  auto prevCell = cell.template getFieldPointer<PREV_CELL>();
  for (unsigned i=0; i < missing.size(); ++i) {
    prevCell[i] = cell[missing[i]];
  }
}

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::apply(CELL& cell)
{
  using V = typename CELL::value_t;
  constexpr auto missing = util::populationsContributingToVelocity<DESCRIPTOR,direction,-orientation>();

  auto prevCell = cell.template getField<PREV_CELL>();

  for (unsigned i=0; i < missing.size(); ++i) {
    cell[missing[i]] = prevCell[i];
  }

  V rho0, u0[3];
  V rho1, u1[3];
  V rho2, u2[3];

  cell.computeRhoU(rho0, u0);

  static_assert(direction == 0 || direction == 1 || direction ==2,
                "Direction must be one of 0, 1 or 2");
  if constexpr (direction == 0) {
    cell.neighbor({-orientation  ,0,0}).computeRhoU(rho1, u1);
    cell.neighbor({-orientation*2,0,0}).computeRhoU(rho2, u2);
  }
  else if constexpr (direction == 1) {
    cell.neighbor({0,-orientation  ,0}).computeRhoU(rho1, u1);
    cell.neighbor({0,-orientation*2,0}).computeRhoU(rho2, u2);
  }
  else if constexpr (direction == 2) {
    cell.neighbor({0,0,-orientation  }).computeRhoU(rho1, u1);
    cell.neighbor({0,0,-orientation*2}).computeRhoU(rho2, u2);
  }

  V uDelta[3];
  V uAverage = rho0*u0[direction];
  //if (uAv!=nullptr) {
  //  rho0 = V{1};
  //  rho1 = V{1};
  //  rho2 = V{1};
  //  uAverage = *uAv * rho0;
  //}

  uDelta[0] = -uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
  uDelta[1] = -uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);
  uDelta[2] = -uAverage*0.5*(3*rho0*u0[2]-4*rho1*u1[2]+rho2*u2[2]);

  for (unsigned i=0; i < missing.size(); ++i) {
    auto iPop = missing[i];
    prevCell[i] = cell[iPop]
                + descriptors::invCs2<V,DESCRIPTOR>()*descriptors::t<V,DESCRIPTOR>(iPop)
                  * ( uDelta[0]*descriptors::c<DESCRIPTOR>(iPop,0)
                    + uDelta[1]*descriptors::c<DESCRIPTOR>(iPop,1)
                    + uDelta[2]*descriptors::c<DESCRIPTOR>(iPop,2));
  }

  cell.template setField<PREV_CELL>(prevCell);
}

////////  OuterVelocityEdgeProcessor3D ///////////////////////////////////

template <typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
template <concepts::DynamicCell CELL>
void OuterVelocityEdgeProcessor3D<T,DESCRIPTOR,plane,normal1,normal2>::apply(CELL& cell)
{
  using V = typename CELL::value_t;
  using namespace olb::util::tensorIndices3D;

  constexpr auto direction1 = (plane+1)%3;
  constexpr auto direction2 = (plane+2)%3;

  auto& dynamics = cell.getDynamics();

  V rho10 = getNeighborRho(cell, 1,0);
  V rho01 = getNeighborRho(cell, 0,1);
  V rho20 = getNeighborRho(cell, 2,0);
  V rho02 = getNeighborRho(cell, 0,2);
  V rho = (V)2/(V)3*(rho01+rho10)-(V)1/(V)6*(rho02+rho20);

  V dA_uB_[3][3];

  interpolateGradients<plane,0>           (cell, dA_uB_[0]);
  interpolateGradients<direction1,normal1>(cell, dA_uB_[1]);
  interpolateGradients<direction2,normal2>(cell, dA_uB_[2]);

  V dA_uB[3][3];
  for (int iBeta=0; iBeta<3; ++iBeta) {
    dA_uB[plane][iBeta]      = dA_uB_[0][iBeta];
    dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
    dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
  }
  V omega = dynamics.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
  V sToPi = - rho / descriptors::invCs2<V,DESCRIPTOR>() / omega;
  V pi[util::TensorVal<DESCRIPTOR>::n];
  pi[xx] = (V)2 * dA_uB[0][0] * sToPi;
  pi[yy] = (V)2 * dA_uB[1][1] * sToPi;
  pi[zz] = (V)2 * dA_uB[2][2] * sToPi;
  pi[xy] = (dA_uB[0][1]+dA_uB[1][0]) * sToPi;
  pi[xz] = (dA_uB[0][2]+dA_uB[2][0]) * sToPi;
  pi[yz] = (dA_uB[1][2]+dA_uB[2][1]) * sToPi;

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

template <typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
template <typename CELL, typename V>
V OuterVelocityEdgeProcessor3D<T,DESCRIPTOR,plane,normal1,normal2>
  ::getNeighborRho(CELL& cell, int step1, int step2)
{
  constexpr auto direction1 = (plane+1)%3;
  constexpr auto direction2 = (plane+2)%3;
  int coords[3] { };
  coords[direction1] = -normal1*step1;
  coords[direction2] = -normal2*step2;
  return cell.neighbor(coords).computeRho();
}

template <typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
template <int deriveDirection, int orientation, typename CELL, typename V>
void OuterVelocityEdgeProcessor3D<T,DESCRIPTOR,plane,normal1,normal2>
  ::interpolateGradients(CELL& cell, V velDeriv[DESCRIPTOR::d]) const
{
  fd::DirectedGradients3D<V,DESCRIPTOR,deriveDirection,orientation,deriveDirection,deriveDirection!=plane>
    ::interpolateVector(velDeriv, cell);
}

/////////// OuterVelocityCornerProcessor3D /////////////////////////////////////

template <typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
template <concepts::DynamicCell CELL>
void OuterVelocityCornerProcessor3D<T,DESCRIPTOR,xNormal,yNormal,zNormal>::apply(CELL& cell)
{
  using V = typename CELL::value_t;
  using namespace olb::util::tensorIndices3D;

  auto& dynamics = cell.getDynamics();

  V rho100 = cell.neighbor({-1*xNormal, -0*yNormal, -0*zNormal}).computeRho();
  V rho010 = cell.neighbor({-0*xNormal, -1*yNormal, -0*zNormal}).computeRho();
  V rho001 = cell.neighbor({-0*xNormal, -0*yNormal, -1*zNormal}).computeRho();
  V rho200 = cell.neighbor({-2*xNormal, -0*yNormal, -0*zNormal}).computeRho();
  V rho020 = cell.neighbor({-0*xNormal, -2*yNormal, -0*zNormal}).computeRho();
  V rho002 = cell.neighbor({-0*xNormal, -0*yNormal, -2*zNormal}).computeRho();

  V rho = (V)4/(V)9 * (rho001 + rho010 + rho100) - (V)1/(V)9 * (rho002 + rho020 + rho200);

  V dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d], dz_u[DESCRIPTOR::d];
  fd::DirectedGradients3D<V, DESCRIPTOR, 0, xNormal, 0, true>::interpolateVector(dx_u, cell);
  fd::DirectedGradients3D<V, DESCRIPTOR, 1, yNormal, 0, true>::interpolateVector(dy_u, cell);
  fd::DirectedGradients3D<V, DESCRIPTOR, 2, zNormal, 0, true>::interpolateVector(dz_u, cell);

  V dx_ux = dx_u[0];
  V dy_ux = dy_u[0];
  V dz_ux = dz_u[0];
  V dx_uy = dx_u[1];
  V dy_uy = dy_u[1];
  V dz_uy = dz_u[1];
  V dx_uz = dx_u[2];
  V dy_uz = dy_u[2];
  V dz_uz = dz_u[2];
  V omega = dynamics.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
  V sToPi = - rho / descriptors::invCs2<V,DESCRIPTOR>() / omega;
  V pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (V)2 * dx_ux * sToPi;
  pi[yy] = (V)2 * dy_uy * sToPi;
  pi[zz] = (V)2 * dz_uz * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;
  pi[xz] = (dx_uz + dz_ux) * sToPi;
  pi[yz] = (dy_uz + dz_uy) * sToPi;

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

////////  SlipBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessor3D<T,DESCRIPTOR>::
SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  this->getName() = "SlipBoundaryProcessor3D";
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  int mirrorDirection0;
  int mirrorDirection1;
  int mirrorDirection2;
  int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
    int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ;
    if ( scalarProduct < 0) {
      // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
      if (mult == 0) {
        mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
        mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
      }
      else {
        mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*discreteNormalX;
        mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*discreteNormalY;
        mirrorDirection2 = descriptors::c<DESCRIPTOR>(iPop,2) - mult*scalarProduct*discreteNormalZ;
      }

      // run through all lattice directions and look for match of direction
      for (int i = 1; i < DESCRIPTOR::q; i++) {
        if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
            && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1
            && descriptors::c<DESCRIPTOR>(i,2)==mirrorDirection2) {
          reflectionPop[iPop] = i;
          break;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=x0; iX<=x1; ++iX) {
      for (int iY=y0; iY<=y1; ++iY) {
        for (int iZ=z0; iZ<=z1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (reflectionPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] = blockLattice.get(iX,iY,iZ)[reflectionPop[iPop]];
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  SlipBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new SlipBoundaryProcessor3D<T,DESCRIPTOR>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  PartialSlipBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessor3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), tuner(tuner_)
{
  this->getName() = "PartialSlipBoundaryProcessor3D";
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ < 0) {
      //std::cout << "-----" <<std::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mirrorDirection2;
      int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);

      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalX);

      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalY);
      mirrorDirection2 = (descriptors::c<DESCRIPTOR>(iPop,2) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalZ);

      // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
      if (mult == 0) {
        mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
        mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
      }

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],2)==mirrorDirection2) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (reflectionPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] = tuner*blockLattice.get(iX,iY,iZ)[reflectionPop[iPop]];
            }
          }
          for (int iPop = 1; iPop < DESCRIPTOR::q/2 ; ++iPop) {
            T provv = blockLattice.get(iX,iY,iZ)[descriptors::opposite<DESCRIPTOR>(iPop)];
            blockLattice.get(iX,iY,iZ)[descriptors::opposite<DESCRIPTOR>(iPop)] +=
              (1.-tuner)*blockLattice.get(iX,iY,iZ)[iPop];
            blockLattice.get(iX,iY,iZ)[iPop] += (1.-tuner)*provv;
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  PartialSlipBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessorGenerator3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), tuner(tuner_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  FreeEnergyWallProcessor3D ////////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
template <concepts::DynamicCell CELL, typename PARAMETERS>
void FreeEnergyWallProcessor3D<T,DESCRIPTOR,NORMAL_X,NORMAL_Y,NORMAL_Z>::apply(CELL& cell, PARAMETERS& vars) {

  T addend = cell.template getField<descriptors::ADDEND>();

  T rhoBulk = cell.neighbor({-NORMAL_X, -NORMAL_Y, -NORMAL_Z}).computeRho();
  T rhoTmp = 0.;

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    rhoTmp += cell[iPop];
  }

  T rho = rhoBulk + addend;
  rho -= rhoTmp;

  cell[0] = rho - 1.;

}


////////  FreeEnergyChemPotBoundaryProcessor3D ////////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
template <concepts::DynamicCell CELL>
void FreeEnergyChemPotBoundaryProcessor3DA<T,DESCRIPTOR,NORMAL_X,NORMAL_Y,NORMAL_Z>::apply(CELL& cell) {

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

  T rho0 = cell.computeRho();
  T rho1 = cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).computeRho();

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.template getField<descriptors::CHEM_POTENTIAL>() + (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>());
}

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
template <concepts::DynamicCell CELL>
void FreeEnergyChemPotBoundaryProcessor3DB<T,DESCRIPTOR,NORMAL_X,NORMAL_Y,NORMAL_Z>::apply(CELL& cell) {

  cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());
}


////////  FreeEnergyConvectiveProcessor3D ////////////////////////////////
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
template <concepts::DynamicCell CELL>
void FreeEnergyConvectiveProcessor3D<T,DESCRIPTOR,NORMAL_X,NORMAL_Y,NORMAL_Z>::apply(CELL& cell) {

  T rho, rho0, rho1, u[3];

  rho0 = cell.computeRho();

  cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).computeRhoU(rho1, u);

  T uPerp = 0;

  Vector<T, 3> normalVec({NORMAL_X,NORMAL_Y,NORMAL_Z});

  if (normalVec[2] == 0) {
    if (normalVec[1] == 0) {
      if (normalVec[0] < 0) {
        uPerp = -u[0];
      } else {
        uPerp = u[0];
      }
    } else if (normalVec[0] == 0) {
      if (normalVec[1] < 0) {
        uPerp = -u[1];
      } else {
        uPerp = u[1];
      }
    } else {
      uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1]);
    }
  } else if (normalVec[1] == 0) {
    if (normalVec[0] == 0) {
      if (normalVec[2] < 0) {
        uPerp = -u[2];
      } else {
        uPerp = u[2];
      }
    } else {
      uPerp = util::sqrt(u[0] * u[0] + u[2] * u[2]);
    }
  } else if (normalVec[0] == 0) {
    uPerp = util::sqrt(u[1] * u[1] + u[2] * u[2]);
  } else {
    uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  }

  rho = (rho0 + uPerp * rho1) / (1. + uPerp);
  cell.defineRho(rho);
}


}  // namespace olb

#endif
