/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_LBM_CSE_H
#define DYNAMICS_LBM_CSE_H


#ifndef DISABLE_CSE

#include "lbm.h"
#include "descriptor/descriptor.h"

namespace olb {


template <typename... FIELDS>
struct lbm<descriptors::D2Q5<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
j[0] = -V{1}*cell[1] + V{1}*cell[3];
j[1] = -V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = V{1}/x0;
rho = x0;
u[0] = -x1*(cell[1] - cell[3]);
u[1] = -x1*(cell[2] - cell[4]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
j[0] = -V{1}*cell[1] + V{1}*cell[3];
j[1] = -V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{0.333333333333333} - V{0.333333333333333}*rho;
pi[0] = -rho*u[0]*u[0] + x0 + V{1}*cell[1] + V{1}*cell[3];
pi[1] = -rho*u[0]*u[1];
pi[2] = -rho*u[1]*u[1] + x0 + V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = cell[1] - cell[3];
auto x2 = V{1}/x0;
auto x3 = x1*x2;
auto x4 = cell[2] - cell[4];
auto x5 = V{0.333333333333333}*cell[0];
rho = x0;
u[0] = -x3;
u[1] = -x2*x4;
pi[0] = -x2*x1*x1 - x5 + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4];
pi[1] = -x3*x4;
pi[2] = -x2*x4*x4 - x5 - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = cell[2] - cell[4];
auto x6 = -x5;
auto x7 = x6*x6;
auto x8 = cell[1] - cell[3];
auto x9 = -x8;
auto x15 = x4*(x9*x9) + V{-1};
auto x16 = V{1} / (x2);
auto x17 = x16*(V{3}*cell[1] - V{3}*cell[3]);
auto x18 = V{3}*x3;
auto x19 = x8*x8;
auto x20 = x5*x5;
auto x21 = x20*x4;
auto x22 = V{3}*cell[2] - V{3}*cell[4];
fEq[0] = -V{0.333333333333333}*x1*(x15 + x4*x7) + V{-0.333333333333333};
fEq[1] = V{0.166666666666667}*x1*(x17 + x18*x19 - x21 + V{1}) + V{-0.166666666666667};
fEq[2] = -V{0.166666666666667}*x1*(x15 - x16*x22 - x18*x7) + V{-0.166666666666667};
fEq[3] = V{0.166666666666667}*x1*(-x17 + V{3}*x19*x3 - x21 + V{1}) + V{-0.166666666666667};
fEq[4] = -V{0.166666666666667}*x1*(x16*x22 - x18*x20 + x19*x4 + V{-1}) + V{-0.166666666666667};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = x3 + V{-1};
auto x5 = V{0.166666666666667}*rho;
auto x6 = V{3}*u[0];
auto x7 = V{3}*u[1];
auto x8 = V{3}*x2;
fNeq[0] = V{0.333333333333333}*rho*(x1 + x4) + cell[0] + V{0.333333333333333};
fNeq[1] = -x5*(V{3}*x0 - x4 - x6) + cell[1] + V{0.166666666666667};
fNeq[2] = x5*(x1 + x7 - x8 + V{-1}) + cell[2] + V{0.166666666666667};
fNeq[3] = -x5*(V{3}*x0 - x3 + x6 + V{1}) + cell[3] + V{0.166666666666667};
fNeq[4] = -x5*(-x1 + x7 + x8 + V{1}) + cell[4] + V{0.166666666666667};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{1.5}*x1;
auto x3 = cell[1] - cell[3];
auto x4 = x3*x3;
auto x5 = x2*x4;
auto x6 = cell[2] - cell[4];
auto x7 = x6*x6;
auto x8 = x2*x7;
auto x9 = x8 + V{-1};
auto x15 = V{0.166666666666667}*cell[0] + V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667};
auto x16 = V{1} / (x0);
auto x17 = x16*(V{3}*cell[1] - V{3}*cell[3]);
auto x18 = V{3}*x1;
auto x19 = x16*(V{3}*cell[2] - V{3}*cell[4]);
auto x20 = x18*x7;
fNeq[0] = (x5 + x9)*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}) + cell[0] + V{0.333333333333333};
fNeq[1] = -x15*(x17 + x18*x4 - x8 + V{1}) + cell[1] + V{0.166666666666667};
fNeq[2] = x15*(V{1.5}*x1*x4 - x19 - x20 + V{-1}) + cell[2] + V{0.166666666666667};
fNeq[3] = -x15*(V{3}*x1*x4 - x17 - x9) + cell[3] + V{0.166666666666667};
fNeq[4] = x15*(x19 - x20 + x5 + V{-1}) + cell[4] + V{0.166666666666667};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = u[0]*u[0];
auto x7 = V{1.5}*x6;
auto x8 = u[1]*u[1];
auto x9 = V{1.5}*x8;
auto x10 = x9 + V{-1};
auto x11 = V{0.166666666666667}*omega;
auto x12 = V{3}*u[0];
auto x13 = V{3}*u[1];
auto x14 = V{3}*x8;
cell[0] = -V{0.333333333333333}*omega*(rho*(x10 + x7) + V{1}) - x5*cell[0];
cell[1] = x11*(rho*(-x10 - x12 + V{3}*x6) + V{-1}) - x5*cell[1];
cell[2] = -x11*(rho*(x13 - x14 + x7 + V{-1}) + V{1}) - x5*cell[2];
cell[3] = x11*(rho*(x12 + V{3}*x6 - x9 + V{1}) + V{-1}) - x5*cell[3];
cell[4] = x11*(rho*(x13 + x14 - x7 + V{1}) + V{-1}) - x5*cell[4];
return x6 + x8;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = V{3}*u[0];
auto x7 = V{0.166666666666667}*omega;
auto x8 = V{3}*u[1];
cell[0] = V{0.333333333333333}*omega*(rho + V{-1}) - x5*cell[0];
cell[1] = -x5*cell[1] - x7*(rho*(x6 + V{-1}) + V{1});
cell[2] = -x5*cell[2] - x7*(rho*(x8 + V{-1}) + V{1});
cell[3] = -x5*cell[3] + x7*(rho*(x6 + V{1}) + V{-1});
cell[4] = -x5*cell[4] + x7*(rho*(x8 + V{1}) + V{-1});
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x5 = j[0]*j[0];
auto x6 = V{0.5}*x5;
auto x7 = j[1]*j[1];
auto x8 = V{0.5}*x7;
auto x9 = omega + V{-1};
auto x10 = V{0.25}*x7;
auto x11 = V{0.5}*j[0];
auto x12 = V{0.5}*pressure;
auto x13 = V{0.166666666666667} - x12;
auto x14 = V{0.25}*x5;
auto x15 = V{0.5}*j[1];
auto x16 = x12 + V{-0.166666666666667};
cell[0] = -omega*(-V{1}*pressure + x6 + x8 + V{0.333333333333333}) - x9*cell[0];
cell[1] = -omega*(x10 + x11 + x13 - x6) - x9*cell[1];
cell[2] = -omega*(x13 + x14 + x15 - x8) - x9*cell[2];
cell[3] = omega*(-x10 + x11 + x16 + x6) - x9*cell[3];
cell[4] = omega*(-x14 + x15 + x16 + x8) - x9*cell[4];
return x5 + x7;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = V{0.333333333333333}*rho;
auto x7 = u[0]*u[0];
auto x8 = V{1.5}*x7;
auto x9 = u[1]*u[1];
auto x10 = V{1.5}*x9;
auto x11 = x10 + V{-1};
auto x12 = x11 + x8;
auto x13 = V{0.166666666666667}*rho;
auto x14 = V{3}*u[0];
auto x15 = -x11 - x14 + V{3}*x7;
auto x16 = V{3}*u[1];
auto x17 = V{3}*x9;
auto x18 = x16 - x17 + x8 + V{-1};
auto x19 = -x10 + x14 + V{3}*x7 + V{1};
auto x20 = x16 + x17 - x8 + V{1};
cell[0] = -ratioRho*x12*x6 - x5*(x12*x6 + cell[0] + V{0.333333333333333}) + V{-0.333333333333333};
cell[1] = V{0.166666666666667}*ratioRho*rho*x15 - x5*(-x13*x15 + cell[1] + V{0.166666666666667}) + V{-0.166666666666667};
cell[2] = -ratioRho*x13*x18 - x5*(x13*x18 + cell[2] + V{0.166666666666667}) + V{-0.166666666666667};
cell[3] = V{0.166666666666667}*ratioRho*rho*x19 - x5*(-x13*x19 + cell[3] + V{0.166666666666667}) + V{-0.166666666666667};
cell[4] = V{0.166666666666667}*ratioRho*rho*x20 - x5*(-x13*x20 + cell[4] + V{0.166666666666667}) + V{-0.166666666666667};
return x7 + x9;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = V{0.5}*cell[1];
auto x7 = V{3}*u[0];
auto x8 = x7 + V{-1};
auto x9 = V{0.0833333333333333}*rho;
auto x10 = x7 + V{1};
auto x11 = x10*x9;
auto x12 = V{0.166666666666667}*rho;
auto x13 = V{0.5}*cell[2];
auto x14 = V{3}*u[1];
auto x15 = x14 + V{-1};
auto x16 = x14 + V{1};
auto x17 = x16*x9;
cell[0] = V{0.333333333333333}*rho + V{-0.333333333333333};
cell[1] = -x12*x8 + x5*(-x11 - x6 - x8*x9 + V{0.5}*cell[3]) + V{-0.166666666666667};
cell[2] = -x12*x15 + x5*(-x13 - x15*x9 - x17 + V{0.5}*cell[4]) + V{-0.166666666666667};
cell[3] = V{0.166666666666667}*rho*x10 - x5*(-x11 - x6 - x8*x9 + V{0.5}*cell[3]) + V{-0.166666666666667};
cell[4] = V{0.166666666666667}*rho*x16 - x5*(-x13 - x15*x9 - x17 + V{0.5}*cell[4]) + V{-0.166666666666667};
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x5 = u[0]*u[0];
auto x6 = V{1.5}*x5;
auto x7 = u[1]*u[1];
auto x8 = V{1.5}*x7;
auto x9 = x8 + V{-1};
auto x10 = omega + V{-1};
auto x11 = V{3}*u[0];
auto x12 = x10*(V{0.5}*pi[0] - V{0.25}*pi[2]) + V{0.166666666666667};
auto x13 = V{0.25}*pi[0] - V{0.5}*pi[2];
auto x14 = V{0.166666666666667}*rho;
auto x15 = V{3}*u[1];
auto x16 = V{3}*x7;
cell[0] = -V{0.333333333333333}*rho*(x6 + x9) + V{0.5}*x10*(pi[0] + pi[2]) + V{-0.333333333333333};
cell[1] = V{0.166666666666667}*rho*(-x11 + V{3}*x5 - x9) - x12;
cell[2] = x10*x13 - x14*(x15 - x16 + x6 + V{-1}) + V{-0.166666666666667};
cell[3] = V{0.166666666666667}*rho*(x11 + V{3}*x5 - x8 + V{1}) - x12;
cell[4] = x10*x13 + x14*(x15 + x16 - x6 + V{1}) + V{-0.166666666666667};
return x5 + x7;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x5 = V{3}*newU[0];
auto x6 = V{3}*newU[1];
cell[0] = V{0.333333333333333}*newRho + V{-0.333333333333333};
cell[1] = -V{0.166666666666667}*newRho*(x5 + V{-1}) + V{-0.166666666666667};
cell[2] = -V{0.166666666666667}*newRho*(x6 + V{-1}) + V{-0.166666666666667};
cell[3] = V{0.166666666666667}*newRho*(x5 + V{1}) + V{-0.166666666666667};
cell[4] = V{0.166666666666667}*newRho*(x6 + V{1}) + V{-0.166666666666667};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x5 = oldU[0]*oldU[0];
auto x6 = V{1.5}*x5;
auto x7 = oldU[1]*oldU[1];
auto x8 = V{1.5}*x7;
auto x9 = x8 + V{-1};
auto x10 = newU[0]*newU[0];
auto x11 = V{1.5}*x10;
auto x12 = newU[1]*newU[1];
auto x13 = V{1.5}*x12;
auto x14 = x13 + V{-1};
auto x15 = V{0.166666666666667}*newRho;
auto x16 = V{3}*newU[0];
auto x17 = V{0.166666666666667}*oldRho;
auto x18 = V{3}*oldU[0];
auto x19 = V{3}*oldU[1];
auto x20 = V{3}*x7;
auto x21 = V{3}*newU[1];
auto x22 = V{3}*x12;
cell[0] = -V{0.333333333333333}*newRho*(x11 + x14) + V{0.333333333333333}*oldRho*(x6 + x9) + cell[0];
cell[1] = x15*(V{3}*x10 - x14 - x16) - x17*(-x18 + V{3}*x5 - x9) + cell[1];
cell[2] = -x15*(x11 + x21 - x22 + V{-1}) + x17*(x19 - x20 + x6 + V{-1}) + cell[2];
cell[3] = x15*(V{3}*x10 - x13 + x16 + V{1}) - x17*(x18 + V{3}*x5 - x8 + V{1}) + cell[3];
cell[4] = x15*(-x11 + x21 + x22 + V{1}) - x17*(x19 + x20 - x6 + V{1}) + cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x5 = V{0.5}*pi[0];
auto x6 = V{0.5}*pi[2];
auto x7 = u[0]*u[0];
auto x8 = V{1.5}*x7;
auto x9 = u[1]*u[1];
auto x10 = V{1.5}*x9;
auto x11 = x10 + V{-1};
auto x12 = V{0.166666666666667}*rho;
auto x13 = V{3}*u[0];
auto x14 = x5 - V{0.25}*pi[2] + V{-0.166666666666667};
auto x15 = V{0.25}*pi[0];
auto x16 = V{3}*u[1];
auto x17 = V{3}*x9;
cell[0] = -V{0.333333333333333}*rho*(x11 + x8) - x5 - x6 + V{-0.333333333333333};
cell[1] = x12*(-x11 - x13 + V{3}*x7) + x14;
cell[2] = -x12*(x16 - x17 + x8 + V{-1}) - x15 + V{0.5}*pi[2] + V{-0.166666666666667};
cell[3] = x12*(-x10 + x13 + V{3}*x7 + V{1}) + x14;
cell[4] = x12*(x16 + x17 - x8 + V{1}) - x15 + x6 + V{-0.166666666666667};

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
auto x1 = x0 + V{1};
auto x2 = cell[1] - cell[3];
auto x3 = cell[2] - cell[4];
auto x4 = x0 + V{1};
auto x5 = x4*(x2*force[1] + x3*force[0]);
auto x6 = x2*x3;
auto x7 = V{0.333333333333333}*cell[0];
auto x8 = V{1} / (x1);
auto x9 = x4*x8;
auto x10 = x3*x9*force[1] + x7 + x8*(x3*x3) + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4];
auto x11 = x2*x9*force[0] + x7 + x8*(x2*x2) - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4];
return x10*x10 + x11*x11 + (V{0.5}*x5 + V{1}*x6)*(V{1}*x5 + V{2}*x6)/((x1)*(x1));
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = cell[1] - cell[3];
auto x2 = x1*x1;
auto x3 = cell[2] - cell[4];
auto x4 = x3*x3;
auto x5 = V{0.333333333333333}*cell[0];
auto x6 = V{1}/x0;
auto x7 = x4*x6 + x5 + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4];
auto x8 = x2*x6 + x5 - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4];
return x7*x7 + x8*x8 + V{2}*x2*x4/((x0)*(x0));
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x5 = force[0]*u[0];
auto x6 = force[1]*u[1];
auto x7 = rho*(V{0.5}*omega + V{-1});
auto x8 = V{6}*u[0];
auto x9 = V{0.166666666666667}*force[0];
auto x10 = -V{0.5}*x6;
auto x11 = V{6}*u[1];
auto x12 = V{0.166666666666667}*force[1];
auto x13 = -V{0.5}*x5;
cell[0] = V{1}*x7*(x5 + x6) + cell[0];
cell[1] = -x7*(x10 + x9*(x8 + V{-3})) + cell[1];
cell[2] = -x7*(x12*(x11 + V{-3}) + x13) + cell[2];
cell[3] = -x7*(x10 + x9*(x8 + V{3})) + cell[3];
cell[4] = -x7*(x12*(x11 + V{3}) + x13) + cell[4];

}

template <typename CELL, typename RHO, typename NABLARHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addLiangForce(CELL& cell, RHO& rho, NABLARHO& nablarho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x5 = rho*force[0];
auto x6 = nablarho[0]*u[0];
auto x7 = V{0.25}*omega + V{-0.5};
auto x8 = rho*force[1];
auto x9 = nablarho[1]*u[1];
cell[0] = cell[0];
cell[1] = x7*(x5 - x6) + cell[1];
cell[2] = x7*(x8 - x9) + cell[2];
cell[3] = -x7*(x5 + x6) + cell[3];
cell[4] = -x7*(x8 + x9) + cell[4];

}

template <typename CELL, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addAllenCahnForce(CELL& cell, OMEGA& omega, FORCE& force) any_platform
{
auto x5 = V{0.25}*omega + V{-0.5};
auto x6 = x5*force[0];
auto x7 = x5*force[1];
cell[0] = cell[0];
cell[1] = x6 + cell[1];
cell[2] = x7 + cell[2];
cell[3] = -x6 + cell[3];
cell[4] = -x7 + cell[4];

}

template <typename CELL, typename OMEGA, typename SOURCE, typename V=typename CELL::value_t>
static void addAllenCahnSource(CELL& cell, OMEGA& omega, SOURCE& source) any_platform
{
auto x5 = V{0.166666666666667}*source;
cell[0] = V{0.333333333333333}*source + cell[0];
cell[1] = x5 + cell[1];
cell[2] = x5 + cell[2];
cell[3] = x5 + cell[3];
cell[4] = x5 + cell[4];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D2Q9<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = cell[1] - cell[5];
j[0] = -V{1}*x0 - V{1}*cell[2] - V{1}*cell[3] + V{1}*cell[6] + V{1}*cell[7];
j[1] = V{1}*x0 - V{1}*cell[3] - V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[2] + cell[3];
auto x1 = cell[7] + cell[8];
auto x2 = x0 + x1 + cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = cell[1] - cell[5];
auto x4 = V{1}/x2;
rho = x2;
u[0] = -x4*(x0 + x3 - cell[6] - cell[7]);
u[1] = x4*(x1 + x3 - cell[3] - cell[4]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[2] + cell[3];
auto x1 = cell[7] + cell[8];
auto x2 = cell[1] - cell[5];
rho = x0 + x1 + cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + V{1};
j[0] = -V{1}*x0 - V{1}*x2 + V{1}*cell[6] + V{1}*cell[7];
j[1] = V{1}*x1 + V{1}*x2 - V{1}*cell[3] - V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{1}*cell[1] + V{1}*cell[5];
auto x1 = -V{0.333333333333333}*rho + x0 + V{1}*cell[3] + V{1}*cell[7] + V{0.333333333333333};
pi[0] = -rho*u[0]*u[0] + x1 + V{1}*cell[2] + V{1}*cell[6];
pi[1] = -rho*u[0]*u[1] - x0 + V{1}*cell[3] + V{1}*cell[7];
pi[2] = -rho*u[1]*u[1] + x1 + V{1}*cell[4] + V{1}*cell[8];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[1] + cell[2];
auto x1 = cell[7] + cell[8];
auto x2 = x0 + x1 + cell[0] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = -cell[5];
auto x4 = x3 + cell[3];
auto x5 = x0 + x4 - cell[6] - cell[7];
auto x6 = V{1} / (x2);
auto x7 = V{1}*x6;
auto x8 = x1 + x3 + cell[1] - cell[3] - cell[4];
auto x9 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
rho = x2;
u[0] = -x5*x7;
u[1] = x7*x8;
pi[0] = -x7*x5*x5 + x9 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
pi[1] = V{1}*x4 + V{1}*x5*x6*x8 - V{1}*cell[1] + V{1}*cell[7];
pi[2] = -x7*x8*x8 + x9 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = cell[1] - cell[5];
auto x6 = -cell[4] + cell[8];
auto x7 = x5 + x6 - cell[3] + cell[7];
auto x8 = x7*x7;
auto x9 = x4*x8;
auto x10 = cell[2] - cell[6];
auto x11 = x10 + x5 + cell[3] - cell[7];
auto x12 = -x11;
auto x13 = x4*(x12*x12) + V{-1};
auto x23 = x13 + x9;
auto x24 = V{4.5}*x3;
auto x25 = x10 + x6 + V{2}*cell[1] - V{2}*cell[5];
auto x26 = x24*(x25*x25);
auto x27 = V{1} / (x2);
auto x28 = V{3}*cell[3];
auto x29 = V{3}*cell[7];
auto x30 = V{3}*cell[1] - V{3}*cell[5];
auto x31 = x28 - x29 + x30 + V{3}*cell[2] - V{3}*cell[6];
auto x32 = x27*x31;
auto x33 = x32 - x9 + V{1};
auto x34 = -x28 + x29 + x30 - V{3}*cell[4] + V{3}*cell[8];
auto x35 = x27*x34;
auto x36 = x11*x11;
auto x37 = x36*x4;
auto x38 = x35 - x37;
auto x39 = V{3}*x3;
auto x40 = x10 + V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8];
auto x41 = -x40;
auto x42 = x39*x8;
auto x43 = x32 + x9 + V{-1};
auto x44 = x37 + x43;
fEq[0] = -V{0.444444444444444}*x1*x23 + V{-0.444444444444444};
fEq[1] = V{0.0277777777777778}*x1*(x26 + x33 + x38) + V{-0.0277777777777778};
fEq[2] = V{0.111111111111111}*x1*(x33 + x36*x39) + V{-0.111111111111111};
fEq[3] = -V{0.0277777777777778}*(x1*(x23 - x24*x41*x41 - x27*x31 + x35) + V{1});
fEq[4] = -V{0.111111111111111}*x1*(x13 + x35 - x42) + V{-0.111111111111111};
fEq[5] = -V{0.0277777777777778}*x1*(-x26 + x35 + x44) + V{-0.0277777777777778};
fEq[6] = V{0.111111111111111}*x1*(V{3}*x3*x36 - x43) + V{-0.111111111111111};
fEq[7] = V{0.0277777777777778}*(-x1*(-x27*x34 - V{4.5}*x3*x40*x40 + x44) + V{-1});
fEq[8] = V{0.111111111111111}*x1*(x38 + x42 + V{1}) + V{-0.111111111111111};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = x3 + V{-1};
auto x5 = x1 + x4;
auto x6 = V{0.0277777777777778}*rho;
auto x7 = V{3}*u[1];
auto x8 = -x7;
auto x9 = u[0] - u[1];
auto x10 = x9*x9;
auto x11 = V{3}*u[0];
auto x12 = x11 + x5;
auto x13 = V{0.111111111111111}*rho;
auto x23 = u[0] + u[1];
auto x24 = V{4.5}*(x23*x23);
auto x25 = V{3}*x2;
auto x26 = -x1;
auto x27 = x11 - x3 + V{1};
auto x28 = x26 + x7;
fNeq[0] = V{0.444444444444444}*rho*x5 + cell[0] + V{0.444444444444444};
fNeq[1] = -x6*(V{4.5}*x10 - x12 - x8) + cell[1] + V{0.0277777777777778};
fNeq[2] = -x13*(V{3}*x0 - x11 - x4) + cell[2] + V{0.111111111111111};
fNeq[3] = x6*(x12 - x24 + x7) + cell[3] + V{0.0277777777777778};
fNeq[4] = x13*(x1 - x25 + x7 + V{-1}) + cell[4] + V{0.111111111111111};
fNeq[5] = -x6*(V{4.5}*x10 + x26 + x27 + x8) + cell[5] + V{0.0277777777777778};
fNeq[6] = -x13*(V{3}*x0 + x27) + cell[6] + V{0.111111111111111};
fNeq[7] = -x6*(x24 + x27 + x28) + cell[7] + V{0.0277777777777778};
fNeq[8] = -x13*(x25 + x28 + V{1}) + cell[8] + V{0.111111111111111};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{1.5}*x1;
auto x3 = cell[1] - cell[5];
auto x4 = -cell[4] + cell[8];
auto x5 = x3 + x4 - cell[3] + cell[7];
auto x6 = x5*x5;
auto x7 = x2*x6;
auto x8 = cell[2] - cell[6];
auto x9 = x3 + x8 + cell[3] - cell[7];
auto x10 = x9*x9;
auto x11 = x10*x2;
auto x12 = x11 + V{-1};
auto x13 = x12 + x7;
auto x23 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x24 = V{4.5}*x1;
auto x25 = x4 + x8 + V{2}*cell[1] - V{2}*cell[5];
auto x26 = x24*(x25*x25);
auto x27 = V{1} / (x0);
auto x28 = V{3}*cell[3];
auto x29 = V{3}*cell[7];
auto x30 = V{3}*cell[1] - V{3}*cell[5];
auto x31 = x27*(x28 - x29 + x30 + V{3}*cell[2] - V{3}*cell[6]);
auto x32 = x31 - x7 + V{1};
auto x33 = -x28 + x29 + x30 - V{3}*cell[4] + V{3}*cell[8];
auto x34 = x27*x33;
auto x35 = -x11;
auto x36 = x34 + x35;
auto x37 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x38 = V{3}*x1;
auto x39 = -x27*x33;
auto x40 = x8 + V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8];
auto x41 = -x40;
auto x42 = x38*x6;
auto x43 = x13 + x31;
fNeq[0] = x13*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + cell[0] + V{0.444444444444444};
fNeq[1] = -x23*(x26 + x32 + x36) + cell[1] + V{0.0277777777777778};
fNeq[2] = -x37*(x10*x38 + x32) + cell[2] + V{0.111111111111111};
fNeq[3] = -x23*(x24*(x41*x41) + x32 + x35 + x39) + cell[3] + V{0.0277777777777778};
fNeq[4] = x37*(x12 + x34 - x42) + cell[4] + V{0.111111111111111};
fNeq[5] = x23*(-x26 + x34 + x43) + cell[5] + V{0.0277777777777778};
fNeq[6] = -x37*(V{3}*x1*x10 - x31 - x7 + V{1}) + cell[6] + V{0.111111111111111};
fNeq[7] = x23*(-V{4.5}*x1*x40*x40 + x39 + x43) + cell[7] + V{0.0277777777777778};
fNeq[8] = -x37*(x36 + x42 + V{1}) + cell[8] + V{0.111111111111111};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = u[0]*u[0];
auto x11 = V{1.5}*x10;
auto x12 = u[1]*u[1];
auto x13 = V{1.5}*x12;
auto x14 = x13 + V{-1};
auto x15 = x11 + x14;
auto x16 = V{0.0277777777777778}*omega;
auto x17 = V{3}*u[1];
auto x18 = -x17;
auto x19 = u[0] - u[1];
auto x20 = x19*x19;
auto x21 = V{3}*u[0];
auto x22 = x15 + x21;
auto x23 = V{0.111111111111111}*omega;
auto x24 = u[0] + u[1];
auto x25 = V{4.5}*(x24*x24);
auto x26 = V{3}*x12;
auto x27 = -x11;
auto x28 = -x13 + x21 + V{1};
auto x29 = x17 + x27;
cell[0] = -V{0.444444444444444}*omega*(rho*x15 + V{1}) - x9*cell[0];
cell[1] = x16*(rho*(-x18 + V{4.5}*x20 - x22) + V{-1}) - x9*cell[1];
cell[2] = x23*(rho*(V{3}*x10 - x14 - x21) + V{-1}) - x9*cell[2];
cell[3] = -x16*(rho*(x17 + x22 - x25) + V{1}) - x9*cell[3];
cell[4] = -x23*(rho*(x11 + x17 - x26 + V{-1}) + V{1}) - x9*cell[4];
cell[5] = x16*(rho*(x18 + V{4.5}*x20 + x27 + x28) + V{-1}) - x9*cell[5];
cell[6] = x23*(rho*(V{3}*x10 + x28) + V{-1}) - x9*cell[6];
cell[7] = x16*(rho*(x25 + x28 + x29) + V{-1}) - x9*cell[7];
cell[8] = x23*(rho*(x26 + x29 + V{1}) + V{-1}) - x9*cell[8];
return x10 + x12;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = V{3}*u[0];
auto x11 = V{3}*u[1];
auto x12 = x11 + V{1};
auto x13 = V{0.0277777777777778}*omega;
auto x14 = x10 + V{-1};
auto x15 = V{0.111111111111111}*omega;
auto x16 = x10 + V{1};
cell[0] = V{0.444444444444444}*omega*(rho + V{-1}) - x9*cell[0];
cell[1] = x13*(rho*(-x10 + x12) + V{-1}) - x9*cell[1];
cell[2] = -x15*(rho*x14 + V{1}) - x9*cell[2];
cell[3] = -x13*(rho*(x11 + x14) + V{1}) - x9*cell[3];
cell[4] = -x15*(rho*(x11 + V{-1}) + V{1}) - x9*cell[4];
cell[5] = x13*(rho*(-x11 + x16) + V{-1}) - x9*cell[5];
cell[6] = x15*(rho*x16 + V{-1}) - x9*cell[6];
cell[7] = x13*(rho*(x11 + x16) + V{-1}) - x9*cell[7];
cell[8] = x15*(rho*x12 + V{-1}) - x9*cell[8];
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x9 = j[0]*j[0];
auto x10 = j[1]*j[1];
auto x11 = omega + V{-1};
auto x12 = j[0] - j[1];
auto x13 = -x12;
auto x14 = V{0.0833333333333333}*j[1];
auto x15 = V{0.0833333333333333}*j[0];
auto x16 = V{0.0416666666666667}*x9;
auto x17 = V{0.0416666666666667}*x10;
auto x18 = V{0.0833333333333333}*pressure;
auto x19 = x16 + x17 - x18 + V{0.0277777777777778};
auto x20 = x15 + x19;
auto x21 = V{0.166666666666667}*x10;
auto x22 = V{0.333333333333333}*j[0];
auto x23 = V{0.333333333333333}*x9;
auto x24 = V{0.333333333333333}*pressure;
auto x25 = V{0.111111111111111} - x24;
auto x26 = j[0] + j[1];
auto x27 = V{0.125}*(x26*x26);
auto x28 = V{0.166666666666667}*x9;
auto x29 = V{0.333333333333333}*j[1];
auto x30 = V{0.333333333333333}*x10;
auto x31 = x24 + V{-0.111111111111111};
cell[0] = -omega*(-V{1.33333333333333}*pressure + V{0.666666666666667}*x10 + V{0.666666666666667}*x9 + V{0.444444444444444}) - x11*cell[0];
cell[1] = -omega*(-x14 + x20 - V{0.125}*x13*x13) - x11*cell[1];
cell[2] = -omega*(x21 + x22 - x23 + x25) - x11*cell[2];
cell[3] = -omega*(x14 + x20 - x27) - x11*cell[3];
cell[4] = -omega*(x25 + x28 + x29 - x30) - x11*cell[4];
cell[5] = -omega*(x14 - x15 + x19 - V{0.125}*x12*x12) - x11*cell[5];
cell[6] = omega*(-x21 + x22 + x23 + x31) - x11*cell[6];
cell[7] = omega*(x14 + x15 - x16 - x17 + x18 + x27 + V{-0.0277777777777778}) - x11*cell[7];
cell[8] = omega*(-x28 + x29 + x30 + x31) - x11*cell[8];
return x10 + x9;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = V{0.444444444444444}*rho;
auto x11 = u[0]*u[0];
auto x12 = V{1.5}*x11;
auto x13 = u[1]*u[1];
auto x14 = V{1.5}*x13;
auto x15 = x14 + V{-1};
auto x16 = x12 + x15;
auto x17 = V{0.0277777777777778}*rho;
auto x18 = V{3}*u[1];
auto x19 = -x18;
auto x20 = u[0] - u[1];
auto x21 = x20*x20;
auto x22 = V{3}*u[0];
auto x23 = x16 + x22;
auto x24 = -x19 + V{4.5}*x21 - x23;
auto x25 = V{0.111111111111111}*rho;
auto x26 = V{3}*x11 - x15 - x22;
auto x27 = u[0] + u[1];
auto x28 = V{4.5}*(x27*x27);
auto x29 = x18 + x23 - x28;
auto x30 = V{3}*x13;
auto x31 = x12 + x18 - x30 + V{-1};
auto x32 = -x12;
auto x33 = -x14 + x22 + V{1};
auto x34 = x19 + V{4.5}*x21 + x32 + x33;
auto x35 = V{3}*x11 + x33;
auto x36 = x18 + x32;
auto x37 = x28 + x33 + x36;
auto x38 = x30 + x36 + V{1};
cell[0] = -ratioRho*x10*x16 - x9*(x10*x16 + cell[0] + V{0.444444444444444}) + V{-0.444444444444444};
cell[1] = V{0.0277777777777778}*ratioRho*rho*x24 - x9*(-x17*x24 + cell[1] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[2] = V{0.111111111111111}*ratioRho*rho*x26 - x9*(-x25*x26 + cell[2] + V{0.111111111111111}) + V{-0.111111111111111};
cell[3] = -ratioRho*x17*x29 - x9*(x17*x29 + cell[3] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[4] = -ratioRho*x25*x31 - x9*(x25*x31 + cell[4] + V{0.111111111111111}) + V{-0.111111111111111};
cell[5] = V{0.0277777777777778}*ratioRho*rho*x34 - x9*(-x17*x34 + cell[5] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[6] = V{0.111111111111111}*ratioRho*rho*x35 - x9*(-x25*x35 + cell[6] + V{0.111111111111111}) + V{-0.111111111111111};
cell[7] = V{0.0277777777777778}*ratioRho*rho*x37 - x9*(-x17*x37 + cell[7] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[8] = V{0.111111111111111}*ratioRho*rho*x38 - x9*(-x25*x38 + cell[8] + V{0.111111111111111}) + V{-0.111111111111111};
return x11 + x13;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = V{0.0833333333333333}*cell[8];
auto x11 = V{3}*u[1];
auto x12 = V{3}*u[0];
auto x13 = x12 + V{1};
auto x14 = -x11 + x13;
auto x15 = V{0.00462962962962963}*rho;
auto x16 = x11 + V{1};
auto x17 = -x12 + x16;
auto x18 = x11 + V{-1};
auto x19 = V{0.00925925925925926}*rho;
auto x20 = x12 + V{-1};
auto x21 = x13*x19 + x19*x20 + V{0.0833333333333333}*cell[2] - V{0.0833333333333333}*cell[6];
auto x22 = x9*(V{0.00925925925925926}*rho*x16 + V{0.00462962962962963}*rho*x17 + V{0.00925925925925926}*rho*x18 - x10 - x14*x15 - x21 - V{0.166666666666667}*cell[1] + V{0.0833333333333333}*cell[4] + V{0.166666666666667}*cell[5]);
auto x23 = V{0.0277777777777778}*rho;
auto x24 = V{0.333333333333333}*cell[1];
auto x25 = V{0.333333333333333}*cell[5];
auto x26 = x14*x19;
auto x27 = V{0.037037037037037}*rho;
auto x28 = x17*x19;
auto x29 = x12 + x16;
auto x30 = x12 + x18;
auto x31 = x19*x29 + x19*x30 + V{0.333333333333333}*cell[3] - V{0.333333333333333}*cell[7] + V{4.62592926927149e-18};
auto x32 = x9*(x13*x27 + x20*x27 + x24 - x25 + x26 - x28 + x31 + V{0.333333333333333}*cell[2] - V{0.333333333333333}*cell[6]);
auto x33 = V{0.111111111111111}*rho;
auto x34 = x9*(-x10 + x15*x29 + x15*x30 + x16*x19 + x18*x19 + x21 + V{0.166666666666667}*cell[3] + V{0.0833333333333333}*cell[4] - V{0.166666666666667}*cell[7] + V{2.31296463463574e-18});
auto x35 = x9*(x16*x27 + x18*x27 - x24 + x25 - x26 + x28 + x31 + V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[8]);
cell[0] = V{0.444444444444444}*rho + V{-0.444444444444444};
cell[1] = x17*x23 + x22 + V{-0.0277777777777778};
cell[2] = -x20*x33 - x32 + V{-0.111111111111111};
cell[3] = -x23*x30 - x34 + V{-0.0277777777777778};
cell[4] = -x18*x33 - x35 + V{-0.111111111111111};
cell[5] = V{0.0277777777777778}*rho*x14 - x22 + V{-0.0277777777777778};
cell[6] = x13*x33 + x32 + V{-0.111111111111111};
cell[7] = x23*x29 + x34 + V{-0.0277777777777778};
cell[8] = x16*x33 + x35 + V{-0.111111111111111};
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x9 = u[0]*u[0];
auto x10 = V{1.5}*x9;
auto x11 = u[1]*u[1];
auto x12 = V{1.5}*x11;
auto x13 = x12 + V{-1};
auto x14 = x10 + x13;
auto x15 = omega + V{-1};
auto x16 = V{3}*u[1];
auto x17 = -x16;
auto x18 = u[0] - u[1];
auto x19 = x18*x18;
auto x20 = V{3}*u[0];
auto x21 = x14 + x20;
auto x22 = V{0.25}*pi[1];
auto x23 = V{0.0833333333333333}*pi[0] + V{0.0833333333333333}*pi[2];
auto x24 = x15*(-x22 + x23) + V{0.0277777777777778};
auto x25 = x15*(V{0.333333333333333}*pi[0] - V{0.166666666666667}*pi[2]) + V{0.111111111111111};
auto x26 = u[0] + u[1];
auto x27 = V{4.5}*(x26*x26);
auto x28 = x15*(x22 + x23) + V{0.0277777777777778};
auto x29 = V{0.166666666666667}*pi[0] - V{0.333333333333333}*pi[2];
auto x30 = V{0.111111111111111}*rho;
auto x31 = V{3}*x11;
auto x32 = -x10;
auto x33 = -x12 + x20 + V{1};
auto x34 = x16 + x32;
cell[0] = -V{0.444444444444444}*rho*x14 + V{0.666666666666667}*x15*(pi[0] + pi[2]) + V{-0.444444444444444};
cell[1] = V{0.0277777777777778}*rho*(-x17 + V{4.5}*x19 - x21) - x24;
cell[2] = V{0.111111111111111}*rho*(-x13 - x20 + V{3}*x9) - x25;
cell[3] = -V{0.0277777777777778}*rho*(x16 + x21 - x27) - x28;
cell[4] = x15*x29 - x30*(x10 + x16 - x31 + V{-1}) + V{-0.111111111111111};
cell[5] = V{0.0277777777777778}*rho*(x17 + V{4.5}*x19 + x32 + x33) - x24;
cell[6] = V{0.111111111111111}*rho*(x33 + V{3}*x9) - x25;
cell[7] = V{0.0277777777777778}*rho*(x27 + x33 + x34) - x28;
cell[8] = x15*x29 + x30*(x31 + x34 + V{1}) + V{-0.111111111111111};
return x11 + x9;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x9 = V{3}*newU[0];
auto x10 = V{3}*newU[1];
auto x11 = x10 + V{1};
auto x12 = x9 + V{-1};
auto x13 = x9 + V{1};
cell[0] = V{0.444444444444444}*newRho + V{-0.444444444444444};
cell[1] = V{0.0277777777777778}*newRho*(x11 - x9) + V{-0.0277777777777778};
cell[2] = -V{0.111111111111111}*newRho*x12 + V{-0.111111111111111};
cell[3] = -V{0.0277777777777778}*newRho*(x10 + x12) + V{-0.0277777777777778};
cell[4] = -V{0.111111111111111}*newRho*(x10 + V{-1}) + V{-0.111111111111111};
cell[5] = V{0.0277777777777778}*newRho*(-x10 + x13) + V{-0.0277777777777778};
cell[6] = V{0.111111111111111}*newRho*x13 + V{-0.111111111111111};
cell[7] = V{0.0277777777777778}*newRho*(x10 + x13) + V{-0.0277777777777778};
cell[8] = V{0.111111111111111}*newRho*x11 + V{-0.111111111111111};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x9 = oldU[0]*oldU[0];
auto x10 = V{1.5}*x9;
auto x11 = oldU[1]*oldU[1];
auto x12 = V{1.5}*x11;
auto x13 = x12 + V{-1};
auto x14 = x10 + x13;
auto x15 = newU[0]*newU[0];
auto x16 = V{1.5}*x15;
auto x17 = newU[1]*newU[1];
auto x18 = V{1.5}*x17;
auto x19 = x18 + V{-1};
auto x20 = x16 + x19;
auto x21 = V{0.0277777777777778}*newRho;
auto x22 = V{3}*newU[1];
auto x23 = -x22;
auto x24 = newU[0] - newU[1];
auto x25 = x24*x24;
auto x26 = V{3}*newU[0];
auto x27 = x20 + x26;
auto x28 = V{0.0277777777777778}*oldRho;
auto x29 = V{3}*oldU[1];
auto x30 = -x29;
auto x31 = oldU[0] - oldU[1];
auto x32 = x31*x31;
auto x33 = V{3}*oldU[0];
auto x34 = x14 + x33;
auto x35 = V{0.111111111111111}*newRho;
auto x36 = V{0.111111111111111}*oldRho;
auto x37 = oldU[0] + oldU[1];
auto x38 = V{4.5}*(x37*x37);
auto x39 = newU[0] + newU[1];
auto x40 = V{4.5}*(x39*x39);
auto x41 = V{3}*x11;
auto x42 = V{3}*x17;
auto x43 = -x16;
auto x44 = -x18 + x26 + V{1};
auto x45 = -x10;
auto x46 = -x12 + x33 + V{1};
auto x47 = x22 + x43;
auto x48 = x29 + x45;
cell[0] = -V{0.444444444444444}*newRho*x20 + V{0.444444444444444}*oldRho*x14 + cell[0];
cell[1] = x21*(-x23 + V{4.5}*x25 - x27) - x28*(-x30 + V{4.5}*x32 - x34) + cell[1];
cell[2] = x35*(V{3}*x15 - x19 - x26) - x36*(-x13 - x33 + V{3}*x9) + cell[2];
cell[3] = -x21*(x22 + x27 - x40) + x28*(x29 + x34 - x38) + cell[3];
cell[4] = -x35*(x16 + x22 - x42 + V{-1}) + x36*(x10 + x29 - x41 + V{-1}) + cell[4];
cell[5] = x21*(x23 + V{4.5}*x25 + x43 + x44) - x28*(x30 + V{4.5}*x32 + x45 + x46) + cell[5];
cell[6] = x35*(V{3}*x15 + x44) - x36*(x46 + V{3}*x9) + cell[6];
cell[7] = x21*(x40 + x44 + x47) - x28*(x38 + x46 + x48) + cell[7];
cell[8] = x35*(x42 + x47 + V{1}) - x36*(x41 + x48 + V{1}) + cell[8];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x9 = u[0]*u[0];
auto x10 = V{1.5}*x9;
auto x11 = u[1]*u[1];
auto x12 = V{1.5}*x11;
auto x13 = x12 + V{-1};
auto x14 = x10 + x13;
auto x15 = V{0.0277777777777778}*rho;
auto x16 = V{3}*u[1];
auto x17 = -x16;
auto x18 = u[0] - u[1];
auto x19 = x18*x18;
auto x20 = V{3}*u[0];
auto x21 = x14 + x20;
auto x22 = V{0.25}*pi[1];
auto x23 = V{0.0833333333333333}*pi[0] + V{0.0833333333333333}*pi[2] + V{-0.0277777777777778};
auto x24 = -x22 + x23;
auto x25 = V{0.111111111111111}*rho;
auto x26 = V{0.333333333333333}*pi[0] - V{0.166666666666667}*pi[2] + V{-0.111111111111111};
auto x27 = u[0] + u[1];
auto x28 = V{4.5}*(x27*x27);
auto x29 = x22 + x23;
auto x30 = V{0.166666666666667}*pi[0];
auto x31 = V{3}*x11;
auto x32 = -x10;
auto x33 = -x12 + x20 + V{1};
auto x34 = x16 + x32;
cell[0] = -V{0.444444444444444}*rho*x14 - V{0.666666666666667}*pi[0] - V{0.666666666666667}*pi[2] + V{-0.444444444444444};
cell[1] = x15*(-x17 + V{4.5}*x19 - x21) + x24;
cell[2] = x25*(-x13 - x20 + V{3}*x9) + x26;
cell[3] = -x15*(x16 + x21 - x28) + x29;
cell[4] = -x25*(x10 + x16 - x31 + V{-1}) - x30 + V{0.333333333333333}*pi[2] + V{-0.111111111111111};
cell[5] = x15*(x17 + V{4.5}*x19 + x32 + x33) + x24;
cell[6] = x25*(x33 + V{3}*x9) + x26;
cell[7] = x15*(x28 + x33 + x34) + x29;
cell[8] = x25*(x31 + x34 + V{1}) - x30 + V{0.333333333333333}*pi[2] + V{-0.111111111111111};

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[7] + cell[8];
auto x1 = cell[2] + cell[3];
auto x2 = x0 + x1 + cell[0] + cell[1] + cell[4] + cell[5] + cell[6];
auto x3 = V{1} / (x2 + V{1});
auto x4 = cell[1] - cell[5];
auto x5 = x1 + x4 - cell[6] - cell[7];
auto x6 = x0 + x4 - cell[3] - cell[4];
auto x7 = x2 + V{1};
auto x8 = x7*(x5*force[1] - x6*force[0]);
auto x9 = x3*x7;
auto x10 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x11 = x10 - x3*x5*x5 - x5*x9*force[0] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
auto x12 = x10 - x3*x6*x6 + x6*x9*force[1] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];
return (V{1}*x3*x5*x6 - V{0.5}*x3*x8 - V{1}*cell[1] + V{1}*cell[3] - V{1}*cell[5] + V{1}*cell[7])*(V{2}*x3*x5*x6 - V{1}*x3*x8 - V{2}*cell[1] + V{2}*cell[3] - V{2}*cell[5] + V{2}*cell[7]) + x11*x11 + x12*x12;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = cell[1] + cell[8];
auto x1 = cell[2] + cell[3];
auto x2 = V{1} / (x0 + x1 + cell[0] + cell[4] + cell[5] + cell[6] + cell[7] + V{1});
auto x3 = -cell[5];
auto x4 = x3 + cell[7];
auto x5 = x0 + x4 - cell[3] - cell[4];
auto x6 = x1 + x3 + cell[1] - cell[6] - cell[7];
auto x7 = -x2*x5*x6 - x4 + cell[1] - cell[3];
auto x8 = V{1}*x2;
auto x9 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x10 = -x8*x6*x6 + x9 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
auto x11 = -x8*x5*x5 + x9 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];
return x10*x10 + x11*x11 + V{2}*(x7*x7);
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x9 = force[0]*u[0];
auto x10 = force[1]*u[1];
auto x11 = rho*(V{0.5}*omega + V{-1});
auto x12 = V{9}*u[0];
auto x13 = V{6}*u[1];
auto x14 = x13 + V{3};
auto x15 = V{9}*u[1];
auto x16 = V{6}*u[0];
auto x17 = V{0.0277777777777778}*x11;
auto x18 = x16 + V{-3};
auto x19 = V{0.111111111111111}*force[0];
auto x20 = -V{0.333333333333333}*x10;
auto x21 = x13 + V{-3};
auto x22 = V{0.111111111111111}*force[1];
auto x23 = -V{0.333333333333333}*x9;
auto x24 = x16 + V{3};
cell[0] = V{1.33333333333333}*x11*(x10 + x9) + cell[0];
cell[1] = -x17*((-x12 + x14)*force[1] - (x15 - x16 + V{3})*force[0]) + cell[1];
cell[2] = -x11*(x18*x19 + x20) + cell[2];
cell[3] = -x17*((x12 + x21)*force[1] + (x15 + x18)*force[0]) + cell[3];
cell[4] = -x11*(x21*x22 + x23) + cell[4];
cell[5] = -x17*((-x15 + x24)*force[0] - (x12 - x13 + V{3})*force[1]) + cell[5];
cell[6] = -x11*(x19*x24 + x20) + cell[6];
cell[7] = -x17*((x12 + x14)*force[1] + (x15 + x24)*force[0]) + cell[7];
cell[8] = -x11*(x14*x22 + x23) + cell[8];

}

template <typename CELL, typename RHO, typename NABLARHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addLiangForce(CELL& cell, RHO& rho, NABLARHO& nablarho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x9 = u[0] - u[1];
auto x10 = x9*nablarho[1];
auto x11 = rho*force[0];
auto x12 = rho*force[1];
auto x13 = x11 - x12;
auto x14 = V{0.5}*omega + V{-1};
auto x15 = V{0.0833333333333333}*x14;
auto x16 = nablarho[0]*u[0];
auto x17 = V{0.333333333333333}*x14;
auto x18 = u[0] + u[1];
auto x19 = x18*nablarho[0];
auto x20 = x18*nablarho[1];
auto x21 = x11 + x12;
auto x22 = nablarho[1]*u[1];
cell[0] = cell[0];
cell[1] = -x15*(-x10 - x13 + x9*nablarho[0]) + cell[1];
cell[2] = x17*(x11 - x16) + cell[2];
cell[3] = x15*(-x19 - x20 + x21) + cell[3];
cell[4] = x17*(x12 - x22) + cell[4];
cell[5] = -x15*(-x10 + x13 + x9*nablarho[0]) + cell[5];
cell[6] = -x17*(x11 + x16) + cell[6];
cell[7] = -x15*(x19 + x20 + x21) + cell[7];
cell[8] = -x17*(x12 + x22) + cell[8];

}

template <typename CELL, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addAllenCahnForce(CELL& cell, OMEGA& omega, FORCE& force) any_platform
{
auto x9 = V{0.5}*omega + V{-1};
auto x10 = V{0.0833333333333333}*x9;
auto x11 = x10*(force[0] - force[1]);
auto x12 = V{0.333333333333333}*x9;
auto x13 = x12*force[0];
auto x14 = x10*(force[0] + force[1]);
auto x15 = x12*force[1];
cell[0] = cell[0];
cell[1] = x11 + cell[1];
cell[2] = x13 + cell[2];
cell[3] = x14 + cell[3];
cell[4] = x15 + cell[4];
cell[5] = -x11 + cell[5];
cell[6] = -x13 + cell[6];
cell[7] = -x14 + cell[7];
cell[8] = -x15 + cell[8];

}

template <typename CELL, typename OMEGA, typename SOURCE, typename V=typename CELL::value_t>
static void addAllenCahnSource(CELL& cell, OMEGA& omega, SOURCE& source) any_platform
{
auto x9 = V{0.0277777777777778}*source;
auto x10 = V{0.111111111111111}*source;
cell[0] = V{0.444444444444444}*source + cell[0];
cell[1] = x9 + cell[1];
cell[2] = x10 + cell[2];
cell[3] = x9 + cell[3];
cell[4] = x10 + cell[4];
cell[5] = x9 + cell[5];
cell[6] = x10 + cell[6];
cell[7] = x9 + cell[7];
cell[8] = x10 + cell[8];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q7<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
j[0] = -V{1}*cell[1] + V{1}*cell[4];
j[1] = -V{1}*cell[2] + V{1}*cell[5];
j[2] = -V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{1}/x0;
rho = x0;
u[0] = -x1*(cell[1] - cell[4]);
u[1] = -x1*(cell[2] - cell[5]);
u[2] = -x1*(cell[3] - cell[6]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
j[0] = -V{1}*cell[1] + V{1}*cell[4];
j[1] = -V{1}*cell[2] + V{1}*cell[5];
j[2] = -V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{0.25} - V{0.25}*rho;
auto x1 = rho*u[0];
pi[0] = -rho*u[0]*u[0] + x0 + V{1}*cell[1] + V{1}*cell[4];
pi[1] = -x1*u[1];
pi[2] = -x1*u[2];
pi[3] = -rho*u[1]*u[1] + x0 + V{1}*cell[2] + V{1}*cell[5];
pi[4] = -rho*u[1]*u[2];
pi[5] = -rho*u[2]*u[2] + x0 + V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = cell[1] - cell[4];
auto x2 = V{1}/x0;
auto x3 = x1*x2;
auto x4 = cell[2] - cell[5];
auto x5 = x2*x4;
auto x6 = cell[3] - cell[6];
auto x7 = V{0.25}*cell[0];
auto x8 = x7 + V{0.25}*cell[3] + V{0.25}*cell[6];
auto x9 = V{0.25}*cell[2] + V{0.25}*cell[5];
auto x10 = V{0.25}*cell[1] + V{0.25}*cell[4];
rho = x0;
u[0] = -x3;
u[1] = -x5;
u[2] = -x2*x6;
pi[0] = -x2*x1*x1 - x8 - x9 + V{0.75}*cell[1] + V{0.75}*cell[4];
pi[1] = -x3*x4;
pi[2] = -x3*x6;
pi[3] = -x10 - x2*x4*x4 - x8 + V{0.75}*cell[2] + V{0.75}*cell[5];
pi[4] = -x5*x6;
pi[5] = -x10 - x2*x6*x6 - x7 - x9 + V{0.75}*cell[3] + V{0.75}*cell[6];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{2}*x3;
auto x5 = cell[1] - cell[4];
auto x6 = -x5;
auto x7 = x6*x6;
auto x8 = x4*x7;
auto x9 = cell[2] - cell[5];
auto x10 = -x9;
auto x11 = x10*x10;
auto x12 = x11*x4;
auto x13 = cell[3] - cell[6];
auto x21 = -x13;
auto x22 = x21*x21;
auto x23 = x22*x4;
auto x24 = x12 + x23 + V{-1};
auto x25 = V{1} / (x2);
auto x26 = V{4}*cell[1] - V{4}*cell[4];
auto x27 = V{6}*x3;
auto x28 = V{4}*cell[2] - V{4}*cell[5];
auto x29 = x8 + V{-1};
auto x30 = V{4}*cell[3] - V{4}*cell[6];
auto x31 = x9*x9;
auto x32 = x31*x4;
auto x33 = x5*x5;
auto x34 = x13*x13;
auto x35 = x34*x4 + V{-1};
auto x36 = x33*x4;
fEq[0] = -V{0.25}*x1*(x24 + x8) + V{-0.25};
fEq[1] = -V{0.125}*x1*(x24 - x25*x26 - x27*x7) + V{-0.125};
fEq[2] = -V{0.125}*x1*(-x11*x27 + x23 - x25*x28 + x29) + V{-0.125};
fEq[3] = -V{0.125}*x1*(x12 - x22*x27 - x25*x30 + x29) + V{-0.125};
fEq[4] = -V{0.125}*x1*(x25*x26 - x27*x33 + x32 + x35) + V{-0.125};
fEq[5] = -V{0.125}*x1*(x25*x28 - x27*x31 + x35 + x36) + V{-0.125};
fEq[6] = -V{0.125}*x1*(x25*x30 - x27*x34 + x32 + x36 + V{-1}) + V{-0.125};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{2}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{2}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{2}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = V{0.125}*rho;
auto x8 = V{4}*u[0];
auto x9 = V{6}*x0;
auto x10 = V{4}*u[1];
auto x11 = V{6}*x2;
auto x12 = x1 + V{-1};
auto x13 = V{4}*u[2];
auto x21 = V{6}*x4;
auto x22 = -x3;
auto x23 = V{1} - x5;
auto x24 = -x1;
fNeq[0] = V{0.25}*rho*(x1 + x6) + cell[0] + V{0.25};
fNeq[1] = x7*(x6 + x8 - x9) + cell[1] + V{0.125};
fNeq[2] = x7*(x10 - x11 + x12 + x5) + cell[2] + V{0.125};
fNeq[3] = x7*(x12 + x13 - x21 + x3) + cell[3] + V{0.125};
fNeq[4] = -x7*(x22 + x23 + x8 + x9) + cell[4] + V{0.125};
fNeq[5] = -x7*(x10 + x11 + x23 + x24) + cell[5] + V{0.125};
fNeq[6] = -x7*(x13 + x21 + x22 + x24 + V{1}) + cell[6] + V{0.125};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{2}*x1;
auto x3 = cell[1] - cell[4];
auto x4 = x3*x3;
auto x5 = x2*x4;
auto x6 = cell[2] - cell[5];
auto x7 = x6*x6;
auto x8 = x2*x7;
auto x9 = cell[3] - cell[6];
auto x10 = x9*x9;
auto x11 = x10*x2;
auto x12 = x11 + x8 + V{-1};
auto x13 = V{0.125}*cell[0] + V{0.125}*cell[1] + V{0.125}*cell[2] + V{0.125}*cell[3] + V{0.125}*cell[4] + V{0.125}*cell[5] + V{0.125}*cell[6] + V{0.125};
auto x21 = V{1} / (x0);
auto x22 = x21*(V{4}*cell[1] - V{4}*cell[4]);
auto x23 = V{6}*x1;
auto x24 = x23*x4;
auto x25 = -V{2}*x1*x7;
auto x26 = -V{2}*x1*x10 + V{1};
auto x27 = x21*(V{4}*cell[2] - V{4}*cell[5]);
auto x28 = x23*x7;
auto x29 = -V{2}*x1*x4;
auto x30 = x21*(V{4}*cell[3] - V{4}*cell[6]);
auto x31 = x10*x23;
auto x32 = x5 + V{-1};
fNeq[0] = (x12 + x5)*(V{0.25}*cell[0] + V{0.25}*cell[1] + V{0.25}*cell[2] + V{0.25}*cell[3] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.25}*cell[6] + V{0.25}) + cell[0] + V{0.25};
fNeq[1] = x13*(-x22 - x24 - x25 - x26) + cell[1] + V{0.125};
fNeq[2] = x13*(-x26 - x27 - x28 - x29) + cell[2] + V{0.125};
fNeq[3] = x13*(-x25 - x29 - x30 - x31 + V{-1}) + cell[3] + V{0.125};
fNeq[4] = x13*(x12 + x22 - x24) + cell[4] + V{0.125};
fNeq[5] = x13*(x11 + x27 - x28 + x32) + cell[5] + V{0.125};
fNeq[6] = x13*(x30 - x31 + x32 + x8) + cell[6] + V{0.125};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = u[0]*u[0];
auto x9 = V{2}*x8;
auto x10 = u[1]*u[1];
auto x11 = V{2}*x10;
auto x12 = u[2]*u[2];
auto x13 = V{2}*x12;
auto x14 = x11 + x13 + V{-1};
auto x15 = V{0.125}*omega;
auto x16 = V{4}*u[0];
auto x17 = V{6}*x8;
auto x18 = V{4}*u[1];
auto x19 = V{6}*x10;
auto x20 = x9 + V{-1};
auto x21 = V{4}*u[2];
auto x22 = V{6}*x12;
auto x23 = -x11;
auto x24 = V{1} - x13;
auto x25 = -x9;
cell[0] = -V{0.25}*omega*(rho*(x14 + x9) + V{1}) - x7*cell[0];
cell[1] = -x15*(rho*(x14 + x16 - x17) + V{1}) - x7*cell[1];
cell[2] = -x15*(rho*(x13 + x18 - x19 + x20) + V{1}) - x7*cell[2];
cell[3] = -x15*(rho*(x11 + x20 + x21 - x22) + V{1}) - x7*cell[3];
cell[4] = x15*(rho*(x16 + x17 + x23 + x24) + V{-1}) - x7*cell[4];
cell[5] = x15*(rho*(x18 + x19 + x24 + x25) + V{-1}) - x7*cell[5];
cell[6] = x15*(rho*(x21 + x22 + x23 + x25 + V{1}) + V{-1}) - x7*cell[6];
return x10 + x12 + x8;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = V{4}*u[0];
auto x9 = V{0.125}*omega;
auto x10 = V{4}*u[1];
auto x11 = V{4}*u[2];
cell[0] = V{0.25}*omega*(rho + V{-1}) - x7*cell[0];
cell[1] = -x7*cell[1] - x9*(rho*(x8 + V{-1}) + V{1});
cell[2] = -x7*cell[2] - x9*(rho*(x10 + V{-1}) + V{1});
cell[3] = -x7*cell[3] - x9*(rho*(x11 + V{-1}) + V{1});
cell[4] = -x7*cell[4] + x9*(rho*(x8 + V{1}) + V{-1});
cell[5] = -x7*cell[5] + x9*(rho*(x10 + V{1}) + V{-1});
cell[6] = -x7*cell[6] + x9*(rho*(x11 + V{1}) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x7 = j[0]*j[0];
auto x8 = j[1]*j[1];
auto x9 = j[2]*j[2];
auto x10 = omega + V{-1};
auto x11 = V{0.25}*x8;
auto x12 = V{0.5}*j[0];
auto x13 = V{0.75}*x7;
auto x14 = V{0.25}*x9;
auto x15 = V{0.5}*pressure;
auto x16 = -x15;
auto x17 = x14 + x16 + V{0.125};
auto x18 = V{0.25}*x7;
auto x19 = V{0.5}*j[1];
auto x20 = V{0.75}*x8;
auto x21 = V{0.5}*j[2];
auto x22 = V{0.75}*x9;
auto x23 = -x11;
auto x24 = -x14 + x15 + V{-0.125};
auto x25 = -x18;
cell[0] = -omega*(-V{1}*pressure + V{0.5}*x7 + V{0.5}*x8 + V{0.5}*x9 + V{0.25}) - x10*cell[0];
cell[1] = -omega*(x11 + x12 - x13 + x17) - x10*cell[1];
cell[2] = -omega*(x17 + x18 + x19 - x20) - x10*cell[2];
cell[3] = -omega*(x11 + x16 + x18 + x21 - x22 + V{0.125}) - x10*cell[3];
cell[4] = omega*(x12 + x13 + x23 + x24) - x10*cell[4];
cell[5] = omega*(x19 + x20 + x24 + x25) - x10*cell[5];
cell[6] = omega*(x15 + x21 + x22 + x23 + x25 + V{-0.125}) - x10*cell[6];
return x7 + x8 + x9;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = V{0.25}*rho;
auto x9 = u[0]*u[0];
auto x10 = V{2}*x9;
auto x11 = u[1]*u[1];
auto x12 = V{2}*x11;
auto x13 = u[2]*u[2];
auto x14 = V{2}*x13;
auto x15 = x12 + x14 + V{-1};
auto x16 = x10 + x15;
auto x17 = V{0.125}*rho;
auto x18 = V{4}*u[0];
auto x19 = V{6}*x9;
auto x20 = x15 + x18 - x19;
auto x21 = ratioRho*x17;
auto x22 = V{4}*u[1];
auto x23 = V{6}*x11;
auto x24 = x10 + V{-1};
auto x25 = x14 + x22 - x23 + x24;
auto x26 = V{4}*u[2];
auto x27 = V{6}*x13;
auto x28 = x12 + x24 + x26 - x27;
auto x29 = -x12;
auto x30 = V{1} - x14;
auto x31 = x18 + x19 + x29 + x30;
auto x32 = -x10;
auto x33 = x22 + x23 + x30 + x32;
auto x34 = x26 + x27 + x29 + x32 + V{1};
cell[0] = -ratioRho*x16*x8 - x7*(x16*x8 + cell[0] + V{0.25}) + V{-0.25};
cell[1] = -x20*x21 - x7*(x17*x20 + cell[1] + V{0.125}) + V{-0.125};
cell[2] = -x21*x25 - x7*(x17*x25 + cell[2] + V{0.125}) + V{-0.125};
cell[3] = -x21*x28 - x7*(x17*x28 + cell[3] + V{0.125}) + V{-0.125};
cell[4] = V{0.125}*ratioRho*rho*x31 - x7*(-x17*x31 + cell[4] + V{0.125}) + V{-0.125};
cell[5] = V{0.125}*ratioRho*rho*x33 - x7*(-x17*x33 + cell[5] + V{0.125}) + V{-0.125};
cell[6] = V{0.125}*ratioRho*rho*x34 - x7*(-x17*x34 + cell[6] + V{0.125}) + V{-0.125};
return x11 + x13 + x9;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = V{0.5}*cell[1];
auto x9 = V{4}*u[0];
auto x10 = x9 + V{-1};
auto x11 = V{0.0625}*rho;
auto x12 = x9 + V{1};
auto x13 = x11*x12;
auto x14 = V{0.125}*rho;
auto x15 = V{0.5}*cell[2];
auto x16 = V{4}*u[1];
auto x17 = x16 + V{-1};
auto x18 = x16 + V{1};
auto x19 = x11*x18;
auto x20 = V{0.5}*cell[3];
auto x21 = V{4}*u[2];
auto x22 = x21 + V{-1};
auto x23 = x21 + V{1};
auto x24 = x11*x23;
cell[0] = V{0.25}*rho + V{-0.25};
cell[1] = -x10*x14 + x7*(-x10*x11 - x13 - x8 + V{0.5}*cell[4]) + V{-0.125};
cell[2] = -x14*x17 + x7*(-x11*x17 - x15 - x19 + V{0.5}*cell[5]) + V{-0.125};
cell[3] = -x14*x22 + x7*(-x11*x22 - x20 - x24 + V{0.5}*cell[6]) + V{-0.125};
cell[4] = V{0.125}*rho*x12 - x7*(-x10*x11 - x13 - x8 + V{0.5}*cell[4]) + V{-0.125};
cell[5] = V{0.125}*rho*x18 - x7*(-x11*x17 - x15 - x19 + V{0.5}*cell[5]) + V{-0.125};
cell[6] = V{0.125}*rho*x23 - x7*(-x11*x22 - x20 - x24 + V{0.5}*cell[6]) + V{-0.125};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x7 = u[0]*u[0];
auto x8 = V{2}*x7;
auto x9 = u[1]*u[1];
auto x10 = V{2}*x9;
auto x11 = u[2]*u[2];
auto x12 = V{2}*x11;
auto x13 = x10 + x12 + V{-1};
auto x14 = omega + V{-1};
auto x15 = V{0.25}*pi[3];
auto x16 = V{0.25}*pi[5];
auto x17 = x15 + x16 - V{0.75}*pi[0];
auto x18 = V{0.125}*rho;
auto x19 = V{4}*u[0];
auto x20 = V{6}*x7;
auto x21 = V{0.25}*pi[0];
auto x22 = x16 + x21 - V{0.75}*pi[3];
auto x23 = V{4}*u[1];
auto x24 = V{6}*x9;
auto x25 = x8 + V{-1};
auto x26 = x15 + x21 - V{0.75}*pi[5];
auto x27 = V{4}*u[2];
auto x28 = V{6}*x11;
auto x29 = -x10;
auto x30 = V{1} - x12;
auto x31 = -x8;
cell[0] = -V{0.25}*rho*(x13 + x8) + V{0.5}*x14*(pi[0] + pi[3] + pi[5]) + V{-0.25};
cell[1] = x14*x17 - x18*(x13 + x19 - x20) + V{-0.125};
cell[2] = x14*x22 - x18*(x12 + x23 - x24 + x25) + V{-0.125};
cell[3] = x14*x26 - x18*(x10 + x25 + x27 - x28) + V{-0.125};
cell[4] = x14*x17 + x18*(x19 + x20 + x29 + x30) + V{-0.125};
cell[5] = x14*x22 + x18*(x23 + x24 + x30 + x31) + V{-0.125};
cell[6] = x14*x26 + x18*(x27 + x28 + x29 + x31 + V{1}) + V{-0.125};
return x11 + x7 + x9;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x7 = V{4}*newU[0];
auto x8 = V{4}*newU[1];
auto x9 = V{4}*newU[2];
cell[0] = V{0.25}*newRho + V{-0.25};
cell[1] = -V{0.125}*newRho*(x7 + V{-1}) + V{-0.125};
cell[2] = -V{0.125}*newRho*(x8 + V{-1}) + V{-0.125};
cell[3] = -V{0.125}*newRho*(x9 + V{-1}) + V{-0.125};
cell[4] = V{0.125}*newRho*(x7 + V{1}) + V{-0.125};
cell[5] = V{0.125}*newRho*(x8 + V{1}) + V{-0.125};
cell[6] = V{0.125}*newRho*(x9 + V{1}) + V{-0.125};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x7 = oldU[0]*oldU[0];
auto x8 = V{2}*x7;
auto x9 = oldU[1]*oldU[1];
auto x10 = V{2}*x9;
auto x11 = oldU[2]*oldU[2];
auto x12 = V{2}*x11;
auto x13 = x10 + x12 + V{-1};
auto x14 = newU[0]*newU[0];
auto x15 = V{2}*x14;
auto x16 = newU[1]*newU[1];
auto x17 = V{2}*x16;
auto x18 = newU[2]*newU[2];
auto x19 = V{2}*x18;
auto x20 = x17 + x19 + V{-1};
auto x21 = V{0.125}*oldRho;
auto x22 = V{4}*oldU[0];
auto x23 = V{6}*x7;
auto x24 = V{0.125}*newRho;
auto x25 = V{4}*newU[0];
auto x26 = V{6}*x14;
auto x27 = V{4}*oldU[1];
auto x28 = V{6}*x9;
auto x29 = x8 + V{-1};
auto x30 = V{4}*newU[1];
auto x31 = V{6}*x16;
auto x32 = x15 + V{-1};
auto x33 = V{4}*oldU[2];
auto x34 = V{6}*x11;
auto x35 = V{4}*newU[2];
auto x36 = V{6}*x18;
auto x37 = -x17;
auto x38 = V{1} - x19;
auto x39 = -x10;
auto x40 = V{1} - x12;
auto x41 = -x15;
auto x42 = -x8;
cell[0] = -V{0.25}*newRho*(x15 + x20) + V{0.25}*oldRho*(x13 + x8) + cell[0];
cell[1] = x21*(x13 + x22 - x23) - x24*(x20 + x25 - x26) + cell[1];
cell[2] = x21*(x12 + x27 - x28 + x29) - x24*(x19 + x30 - x31 + x32) + cell[2];
cell[3] = x21*(x10 + x29 + x33 - x34) - x24*(x17 + x32 + x35 - x36) + cell[3];
cell[4] = -x21*(x22 + x23 + x39 + x40) + x24*(x25 + x26 + x37 + x38) + cell[4];
cell[5] = -x21*(x27 + x28 + x40 + x42) + x24*(x30 + x31 + x38 + x41) + cell[5];
cell[6] = -x21*(x33 + x34 + x39 + x42 + V{1}) + x24*(x35 + x36 + x37 + x41 + V{1}) + cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x7 = u[0]*u[0];
auto x8 = V{2}*x7;
auto x9 = u[1]*u[1];
auto x10 = V{2}*x9;
auto x11 = u[2]*u[2];
auto x12 = V{2}*x11;
auto x13 = x10 + x12 + V{-1};
auto x14 = V{0.125}*rho;
auto x15 = V{4}*u[0];
auto x16 = V{6}*x7;
auto x17 = V{0.25}*pi[3];
auto x18 = V{0.25}*pi[5] + V{0.125};
auto x19 = x17 + x18 - V{0.75}*pi[0];
auto x20 = V{4}*u[1];
auto x21 = V{6}*x9;
auto x22 = x8 + V{-1};
auto x23 = V{0.25}*pi[0];
auto x24 = x18 + x23 - V{0.75}*pi[3];
auto x25 = V{4}*u[2];
auto x26 = V{6}*x11;
auto x27 = x17 + x23 - V{0.75}*pi[5] + V{0.125};
auto x28 = -x10;
auto x29 = V{1} - x12;
auto x30 = -x8;
cell[0] = -V{0.25}*rho*(x13 + x8) - V{0.5}*pi[0] - V{0.5}*pi[3] - V{0.5}*pi[5] + V{-0.25};
cell[1] = -x14*(x13 + x15 - x16) - x19;
cell[2] = -x14*(x12 + x20 - x21 + x22) - x24;
cell[3] = -x14*(x10 + x22 + x25 - x26) - x27;
cell[4] = V{0.125}*rho*(x15 + x16 + x28 + x29) - x19;
cell[5] = V{0.125}*rho*(x20 + x21 + x29 + x30) - x24;
cell[6] = V{0.125}*rho*(x25 + x26 + x28 + x30 + V{1}) - x27;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
auto x1 = x0 + V{1};
auto x2 = V{1} / ((x1)*(x1));
auto x3 = cell[1] - cell[4];
auto x4 = cell[3] - cell[6];
auto x5 = x0 + V{1};
auto x6 = V{0.5}*x5;
auto x7 = V{1}*x3;
auto x8 = x4*x7 + x6*(x3*force[2] + x4*force[0]);
auto x9 = cell[2] - cell[5];
auto x10 = x3*force[1] + x9*force[0];
auto x11 = V{1}*x5;
auto x12 = V{2}*x9;
auto x13 = x4*force[1] + x9*force[2];
auto x14 = V{1} / (x1);
auto x15 = x14*x5;
auto x16 = V{0.25}*cell[0];
auto x17 = x16 + V{0.25}*cell[1] + V{0.25}*cell[4];
auto x18 = V{0.25}*cell[2] + V{0.25}*cell[5];
auto x19 = x14*(x4*x4) + x15*x4*force[2] + x17 + x18 - V{0.75}*cell[3] - V{0.75}*cell[6];
auto x20 = V{0.25}*cell[3] + V{0.25}*cell[6];
auto x21 = x14*(x9*x9) + x15*x9*force[1] + x17 + x20 - V{0.75}*cell[2] - V{0.75}*cell[5];
auto x22 = x14*(x3*x3) + x15*x3*force[0] + x16 + x18 + x20 - V{0.75}*cell[1] - V{0.75}*cell[4];
return x2*(x10*x11 + x12*x3)*(x10*x6 + x7*x9) + x2*(x11*x13 + x12*x4)*(x13*x6 + V{1}*x4*x9) + 2*x2*(x8*x8) + x19*x19 + x21*x21 + x22*x22;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{2}/((x0)*(x0));
auto x2 = cell[1] - cell[4];
auto x3 = x2*x2;
auto x4 = cell[2] - cell[5];
auto x5 = x4*x4;
auto x6 = cell[3] - cell[6];
auto x7 = x6*x6;
auto x8 = V{1}/x0;
auto x9 = V{0.25}*cell[0];
auto x10 = x9 + V{0.25}*cell[1] + V{0.25}*cell[4];
auto x11 = V{0.25}*cell[2] + V{0.25}*cell[5];
auto x12 = x10 + x11 + x7*x8 - V{0.75}*cell[3] - V{0.75}*cell[6];
auto x13 = V{0.25}*cell[3] + V{0.25}*cell[6];
auto x14 = x10 + x13 + x5*x8 - V{0.75}*cell[2] - V{0.75}*cell[5];
auto x15 = x11 + x13 + x3*x8 + x9 - V{0.75}*cell[1] - V{0.75}*cell[4];
return x1*x3*x5 + x1*x3*x7 + x1*x5*x7 + x12*x12 + x14*x14 + x15*x15;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x7 = force[0]*u[0];
auto x8 = force[1]*u[1];
auto x9 = force[2]*u[2];
auto x10 = rho*(V{0.5}*omega + V{-1});
auto x11 = V{12}*u[0];
auto x12 = V{0.125}*force[0];
auto x13 = V{0.5}*x8;
auto x14 = V{0.5}*x9;
auto x15 = x13 + x14;
auto x16 = V{12}*u[1];
auto x17 = V{0.125}*force[1];
auto x18 = V{0.5}*x7;
auto x19 = x14 + x18;
auto x20 = V{12}*u[2];
auto x21 = V{0.125}*force[2];
auto x22 = x13 + x18;
cell[0] = V{1}*x10*(x7 + x8 + x9) + cell[0];
cell[1] = x10*(-x12*(x11 + V{-4}) + x15) + cell[1];
cell[2] = x10*(-x17*(x16 + V{-4}) + x19) + cell[2];
cell[3] = x10*(-x21*(x20 + V{-4}) + x22) + cell[3];
cell[4] = x10*(-x12*(x11 + V{4}) + x15) + cell[4];
cell[5] = x10*(-x17*(x16 + V{4}) + x19) + cell[5];
cell[6] = x10*(-x21*(x20 + V{4}) + x22) + cell[6];

}

template <typename CELL, typename RHO, typename NABLARHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addLiangForce(CELL& cell, RHO& rho, NABLARHO& nablarho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x7 = rho*force[0];
auto x8 = nablarho[0]*u[0];
auto x9 = V{0.25}*omega + V{-0.5};
auto x10 = rho*force[1];
auto x11 = nablarho[1]*u[1];
auto x12 = rho*force[2];
auto x13 = nablarho[2]*u[2];
cell[0] = cell[0];
cell[1] = x9*(x7 - x8) + cell[1];
cell[2] = x9*(x10 - x11) + cell[2];
cell[3] = x9*(x12 - x13) + cell[3];
cell[4] = -x9*(x7 + x8) + cell[4];
cell[5] = -x9*(x10 + x11) + cell[5];
cell[6] = -x9*(x12 + x13) + cell[6];

}

template <typename CELL, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addAllenCahnForce(CELL& cell, OMEGA& omega, FORCE& force) any_platform
{
auto x7 = V{0.25}*omega + V{-0.5};
auto x8 = x7*force[0];
auto x9 = x7*force[1];
auto x10 = x7*force[2];
cell[0] = cell[0];
cell[1] = x8 + cell[1];
cell[2] = x9 + cell[2];
cell[3] = x10 + cell[3];
cell[4] = -x8 + cell[4];
cell[5] = -x9 + cell[5];
cell[6] = -x10 + cell[6];

}

template <typename CELL, typename OMEGA, typename SOURCE, typename V=typename CELL::value_t>
static void addAllenCahnSource(CELL& cell, OMEGA& omega, SOURCE& source) any_platform
{
auto x7 = V{0.125}*source;
cell[0] = V{0.25}*source + cell[0];
cell[1] = x7 + cell[1];
cell[2] = x7 + cell[2];
cell[3] = x7 + cell[3];
cell[4] = x7 + cell[4];
cell[5] = x7 + cell[5];
cell[6] = x7 + cell[6];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q19<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = cell[13] - cell[4];
auto x1 = cell[15] - cell[6];
auto x2 = cell[17] - cell[8];
j[0] = V{1}*x0 + V{1}*x1 + V{1}*cell[10] + V{1}*cell[14] + V{1}*cell[16] - V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x0 + V{1}*x2 + V{1}*cell[11] - V{1}*cell[14] + V{1}*cell[18] - V{1}*cell[2] + V{1}*cell[5] - V{1}*cell[9];
j[2] = V{1}*x1 + V{1}*x2 + V{1}*cell[12] - V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + V{1}*cell[7] + V{1}*cell[9];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[10] + cell[14] + cell[16];
auto x1 = cell[11] + cell[18] + cell[5];
auto x2 = cell[12] + cell[7] + cell[9];
auto x3 = x0 + x1 + x2 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = cell[13] - cell[4];
auto x5 = cell[15] - cell[6];
auto x6 = V{1}/x3;
auto x7 = cell[17] - cell[8];
rho = x3;
u[0] = x6*(x0 + x4 + x5 - cell[1] - cell[5] - cell[7]);
u[1] = x6*(x1 + x4 + x7 - cell[14] - cell[2] - cell[9]);
u[2] = x6*(x2 + x5 + x7 - cell[16] - cell[18] - cell[3]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[10] + cell[14] + cell[16];
auto x1 = cell[11] + cell[18] + cell[5];
auto x2 = cell[12] + cell[7] + cell[9];
auto x3 = cell[13] - cell[4];
auto x4 = cell[15] - cell[6];
auto x5 = cell[17] - cell[8];
rho = x0 + x1 + x2 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
j[0] = V{1}*x0 + V{1}*x3 + V{1}*x4 - V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x1 + V{1}*x3 + V{1}*x5 - V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9];
j[2] = V{1}*x2 + V{1}*x4 + V{1}*x5 - V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = -V{0.333333333333333}*rho;
auto x1 = V{1}*cell[14] + V{1}*cell[5];
auto x2 = x0 + x1 + V{1}*cell[13] + V{1}*cell[4] + V{0.333333333333333};
auto x3 = V{1}*cell[16] + V{1}*cell[7];
auto x4 = x3 + V{1}*cell[15] + V{1}*cell[6];
auto x5 = rho*u[0];
auto x6 = V{1}*cell[18] + V{1}*cell[9];
auto x7 = x6 + V{1}*cell[17] + V{1}*cell[8];
pi[0] = -rho*u[0]*u[0] + x2 + x4 + V{1}*cell[10] + V{1}*cell[1];
pi[1] = -x1 - x5*u[1] + V{1}*cell[13] + V{1}*cell[4];
pi[2] = -x3 - x5*u[2] + V{1}*cell[15] + V{1}*cell[6];
pi[3] = -rho*u[1]*u[1] + x2 + x7 + V{1}*cell[11] + V{1}*cell[2];
pi[4] = -rho*u[1]*u[2] - x6 + V{1}*cell[17] + V{1}*cell[8];
pi[5] = -rho*u[2]*u[2] + x0 + x4 + x7 + V{1}*cell[12] + V{1}*cell[3] + V{0.333333333333333};

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[10] + cell[13] + cell[15];
auto x1 = cell[11] + cell[17] + cell[5];
auto x2 = cell[12] + cell[7] + cell[9];
auto x3 = x0 + x1 + x2 + cell[0] + cell[14] + cell[16] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = -cell[4];
auto x5 = x4 + cell[14];
auto x6 = -cell[6];
auto x7 = x6 + cell[16];
auto x8 = x0 + x5 + x7 - cell[1] - cell[5] - cell[7];
auto x9 = V{1} / (x3);
auto x10 = V{1}*x9;
auto x11 = -cell[8];
auto x12 = x11 + cell[18];
auto x13 = x1 + x12 + x4 + cell[13] - cell[14] - cell[2] - cell[9];
auto x14 = x11 + x2 + x6 + cell[15] - cell[16] + cell[17] - cell[18] - cell[3];
auto x15 = V{0.333333333333333}*cell[0];
auto x16 = x15 + V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x17 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x18 = x8*x9;
auto x19 = V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9];
rho = x3;
u[0] = x10*x8;
u[1] = x10*x13;
u[2] = x10*x14;
pi[0] = -x10*x8*x8 - x16 - x17 + V{0.666666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
pi[1] = -V{1}*x13*x18 - V{1}*x5 + V{1}*cell[13] - V{1}*cell[5];
pi[2] = -V{1}*x14*x18 - V{1}*x7 + V{1}*cell[15] - V{1}*cell[7];
pi[3] = -x10*x13*x13 - x16 - x19 + V{0.666666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
pi[4] = -V{1}*x12 - V{1}*x13*x14*x9 + V{1}*cell[17] - V{1}*cell[9];
pi[5] = -x10*x14*x14 - x15 - x17 - x19 + V{0.666666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[10] + cell[14];
auto x1 = cell[12] + cell[7];
auto x2 = x0 + x1 + cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9];
auto x3 = x2 + V{1};
auto x4 = x2 + V{1};
auto x5 = V{1} / ((x4)*(x4));
auto x6 = V{1.5}*x5;
auto x7 = cell[13] - cell[4];
auto x8 = cell[15] - cell[6];
auto x9 = x7 + x8;
auto x10 = -cell[1];
auto x11 = cell[16] - cell[7];
auto x12 = x10 + x11;
auto x13 = x0 - cell[5];
auto x14 = x12 + x13 + x9;
auto x15 = x14*x14;
auto x16 = x15*x6;
auto x17 = cell[17] - cell[8];
auto x18 = x17 + x7;
auto x19 = cell[18] - cell[9];
auto x20 = -cell[2];
auto x21 = x20 + cell[11] - cell[14] + cell[5];
auto x22 = x18 + x19 + x21;
auto x23 = x22*x22;
auto x24 = x23*x6;
auto x25 = x17 + x8;
auto x45 = -cell[3];
auto x46 = -cell[18] + cell[9];
auto x47 = x45 + x46;
auto x48 = x1 - cell[16];
auto x49 = x25 + x47 + x48;
auto x50 = x49*x49;
auto x51 = x50*x6;
auto x52 = x24 + x51 + V{-1};
auto x53 = x16 + x52;
auto x54 = V{1} / (x4);
auto x55 = V{3}*cell[14];
auto x56 = V{3}*cell[16];
auto x57 = V{3}*cell[5];
auto x58 = V{3}*cell[7];
auto x59 = V{3}*cell[13] - V{3}*cell[4];
auto x60 = V{3}*cell[15] - V{3}*cell[6];
auto x61 = x54*(x55 + x56 - x57 - x58 + x59 + x60 + V{3}*cell[10] - V{3}*cell[1]);
auto x62 = V{3}*x5;
auto x63 = x15*x62;
auto x64 = V{3}*cell[18];
auto x65 = V{3}*cell[9];
auto x66 = V{3}*cell[17] - V{3}*cell[8];
auto x67 = x54*(-x55 + x57 + x59 + x64 - x65 + x66 + V{3}*cell[11] - V{3}*cell[2]);
auto x68 = x23*x62;
auto x69 = x16 + V{-1};
auto x70 = x54*(-x56 + x58 + x60 - x64 + x65 + x66 + V{3}*cell[12] - V{3}*cell[3]);
auto x71 = x50*x62;
auto x72 = V{4.5}*x5;
auto x73 = x12 + cell[10];
auto x74 = x19 + x20 + x25 + x73 + cell[11] + V{2}*cell[13] - V{2}*cell[4];
auto x75 = x72*(x74*x74);
auto x76 = x53 + x61;
auto x77 = -x67;
auto x78 = -cell[17] + cell[8];
auto x79 = x46 + x73 + x78 + x8 - cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5];
auto x80 = -x79;
auto x81 = x10 + x13;
auto x82 = x18 + x47 + x81 + cell[12] + V{2}*cell[15] - V{2}*cell[6];
auto x83 = x72*(x82*x82);
auto x84 = -x70;
auto x85 = x7 - cell[12] + cell[3];
auto x86 = x19 + x78 + x81 + x85 + V{2}*cell[16] - V{2}*cell[7];
auto x87 = -x86;
auto x88 = x21 + x45 + x48 + x9 + V{2}*cell[17] - V{2}*cell[8];
auto x89 = x72*(x88*x88);
auto x90 = x53 + x67;
auto x91 = x11 + x21 + x85 - cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9];
auto x92 = -x91;
auto x93 = -x24;
auto x94 = V{1} - x51;
auto x95 = x93 + x94;
auto x96 = x61 + x95;
auto x97 = -x16;
auto x98 = x67 + x97;
auto x99 = x70 + x97;
auto x100 = -x61;
auto x101 = x53 + x70;
fEq[0] = -V{0.333333333333333}*x3*x53 + V{-0.333333333333333};
fEq[1] = -V{0.0555555555555556}*x3*(x52 + x61 - x63) + V{-0.0555555555555556};
fEq[2] = -V{0.0555555555555556}*x3*(x51 + x67 - x68 + x69) + V{-0.0555555555555556};
fEq[3] = -V{0.0555555555555556}*x3*(x24 + x69 + x70 - x71) + V{-0.0555555555555556};
fEq[4] = -V{0.0277777777777778}*x3*(x67 - x75 + x76) + V{-0.0277777777777778};
fEq[5] = -V{0.0277777777777778}*(x3*(-x72*x80*x80 + x76 + x77) + V{1});
fEq[6] = -V{0.0277777777777778}*x3*(x70 + x76 - x83) + V{-0.0277777777777778};
fEq[7] = -V{0.0277777777777778}*(x3*(-x72*x87*x87 + x76 + x84) + V{1});
fEq[8] = -V{0.0277777777777778}*x3*(x70 - x89 + x90) + V{-0.0277777777777778};
fEq[9] = -V{0.0277777777777778}*(x3*(-x72*x92*x92 + x84 + x90) + V{1});
fEq[10] = V{0.0555555555555556}*x3*(x63 + x96) + V{-0.0555555555555556};
fEq[11] = V{0.0555555555555556}*x3*(x68 + x94 + x98) + V{-0.0555555555555556};
fEq[12] = V{0.0555555555555556}*x3*(x71 + x93 + x99 + V{1}) + V{-0.0555555555555556};
fEq[13] = V{0.0277777777777778}*x3*(x75 + x96 + x98) + V{-0.0277777777777778};
fEq[14] = -V{0.0277777777777778}*(x3*(x100 - x72*x79*x79 + x90) + V{1});
fEq[15] = V{0.0277777777777778}*x3*(x83 + x96 + x99) + V{-0.0277777777777778};
fEq[16] = -V{0.0277777777777778}*(x3*(x100 + x101 - x72*x86*x86) + V{1});
fEq[17] = V{0.0277777777777778}*x3*(x70 + x89 + x95 + x98) + V{-0.0277777777777778};
fEq[18] = -V{0.0277777777777778}*(x3*(x101 - x72*x91*x91 + x77) + V{1});

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{1.5}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = x1 + x6;
auto x8 = V{0.0555555555555556}*rho;
auto x9 = V{3}*u[0];
auto x10 = V{3}*x0;
auto x11 = V{3}*u[1];
auto x12 = V{3}*x2;
auto x13 = x1 + V{-1};
auto x14 = V{3}*u[2];
auto x15 = V{3}*x4;
auto x16 = V{0.0277777777777778}*rho;
auto x17 = u[0] + u[1];
auto x18 = V{4.5}*(x17*x17);
auto x19 = x7 + x9;
auto x20 = -x11;
auto x21 = u[0] - u[1];
auto x22 = -V{4.5}*x21*x21;
auto x23 = u[0] + u[2];
auto x24 = V{4.5}*(x23*x23);
auto x25 = -x14;
auto x45 = -u[2];
auto x46 = x45 + u[0];
auto x47 = -V{4.5}*x46*x46;
auto x48 = u[1] + u[2];
auto x49 = V{4.5}*(x48*x48);
auto x50 = x11 + x7;
auto x51 = x45 + u[1];
auto x52 = -V{4.5}*x51*x51;
auto x53 = -x3;
auto x54 = V{1} - x5;
auto x55 = x53 + x54;
auto x56 = x55 + x9;
auto x57 = -x1;
auto x58 = x11 + x57;
auto x59 = x14 + x57;
auto x60 = -x9;
auto x61 = x14 + x7;
fNeq[0] = V{0.333333333333333}*rho*x7 + cell[0] + V{0.333333333333333};
fNeq[1] = x8*(-x10 + x6 + x9) + cell[1] + V{0.0555555555555556};
fNeq[2] = x8*(x11 - x12 + x13 + x5) + cell[2] + V{0.0555555555555556};
fNeq[3] = x8*(x13 + x14 - x15 + x3) + cell[3] + V{0.0555555555555556};
fNeq[4] = x16*(x11 - x18 + x19) + cell[4] + V{0.0277777777777778};
fNeq[5] = x16*(x19 + x20 + x22) + cell[5] + V{0.0277777777777778};
fNeq[6] = x16*(x14 + x19 - x24) + cell[6] + V{0.0277777777777778};
fNeq[7] = x16*(x19 + x25 + x47) + cell[7] + V{0.0277777777777778};
fNeq[8] = x16*(x14 - x49 + x50) + cell[8] + V{0.0277777777777778};
fNeq[9] = x16*(x25 + x50 + x52) + cell[9] + V{0.0277777777777778};
fNeq[10] = -x8*(x10 + x56) + cell[10] + V{0.0555555555555556};
fNeq[11] = -x8*(x12 + x54 + x58) + cell[11] + V{0.0555555555555556};
fNeq[12] = -x8*(x15 + x53 + x59 + V{1}) + cell[12] + V{0.0555555555555556};
fNeq[13] = -x16*(x18 + x56 + x58) + cell[13] + V{0.0277777777777778};
fNeq[14] = x16*(x22 + x50 + x60) + cell[14] + V{0.0277777777777778};
fNeq[15] = -x16*(x24 + x56 + x59) + cell[15] + V{0.0277777777777778};
fNeq[16] = x16*(x47 + x60 + x61) + cell[16] + V{0.0277777777777778};
fNeq[17] = -x16*(x14 + x49 + x55 + x58) + cell[17] + V{0.0277777777777778};
fNeq[18] = x16*(x20 + x52 + x61) + cell[18] + V{0.0277777777777778};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[10] + cell[14];
auto x1 = cell[12] + cell[7];
auto x2 = x0 + x1 + cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = cell[13] - cell[4];
auto x6 = cell[15] - cell[6];
auto x7 = x5 + x6;
auto x8 = -cell[1];
auto x9 = cell[16] - cell[7];
auto x10 = x8 + x9;
auto x11 = x0 - cell[5];
auto x12 = x10 + x11 + x7;
auto x13 = x12*x12;
auto x14 = x13*x4;
auto x15 = cell[17] - cell[8];
auto x16 = x15 + x5;
auto x17 = cell[18] - cell[9];
auto x18 = -cell[2];
auto x19 = x18 + cell[11] - cell[14] + cell[5];
auto x20 = x16 + x17 + x19;
auto x21 = x20*x20;
auto x22 = x21*x4;
auto x23 = x15 + x6;
auto x24 = -cell[3];
auto x25 = -cell[18] + cell[9];
auto x45 = x24 + x25;
auto x46 = x1 - cell[16];
auto x47 = x23 + x45 + x46;
auto x48 = x47*x47;
auto x49 = x4*x48;
auto x50 = x22 + x49 + V{-1};
auto x51 = x14 + x50;
auto x52 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x53 = V{1} / (x2);
auto x54 = V{3}*cell[14];
auto x55 = V{3}*cell[16];
auto x56 = V{3}*cell[5];
auto x57 = V{3}*cell[7];
auto x58 = V{3}*cell[13] - V{3}*cell[4];
auto x59 = V{3}*cell[15] - V{3}*cell[6];
auto x60 = x53*(x54 + x55 - x56 - x57 + x58 + x59 + V{3}*cell[10] - V{3}*cell[1]);
auto x61 = V{3}*x3;
auto x62 = x13*x61;
auto x63 = V{3}*cell[18];
auto x64 = V{3}*cell[9];
auto x65 = V{3}*cell[17] - V{3}*cell[8];
auto x66 = x53*(-x54 + x56 + x58 + x63 - x64 + x65 + V{3}*cell[11] - V{3}*cell[2]);
auto x67 = x21*x61;
auto x68 = x14 + V{-1};
auto x69 = x53*(-x55 + x57 + x59 - x63 + x64 + x65 + V{3}*cell[12] - V{3}*cell[3]);
auto x70 = x48*x61;
auto x71 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x72 = V{4.5}*x3;
auto x73 = x10 + cell[10];
auto x74 = x17 + x18 + x23 + x73 + cell[11] + V{2}*cell[13] - V{2}*cell[4];
auto x75 = x72*(x74*x74);
auto x76 = x51 + x60;
auto x77 = -x66;
auto x78 = -cell[17] + cell[8];
auto x79 = x25 + x6 + x73 + x78 - cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5];
auto x80 = -x72*x79*x79;
auto x81 = x11 + x8;
auto x82 = x16 + x45 + x81 + cell[12] + V{2}*cell[15] - V{2}*cell[6];
auto x83 = x72*(x82*x82);
auto x84 = -x69;
auto x85 = x5 - cell[12] + cell[3];
auto x86 = x17 + x78 + x81 + x85 + V{2}*cell[16] - V{2}*cell[7];
auto x87 = -x72*x86*x86;
auto x88 = x19 + x24 + x46 + x7 + V{2}*cell[17] - V{2}*cell[8];
auto x89 = x72*(x88*x88);
auto x90 = x51 + x66;
auto x91 = x19 + x85 + x9 - cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9];
auto x92 = -x72*x91*x91;
auto x93 = -x22;
auto x94 = V{1} - x49;
auto x95 = x93 + x94;
auto x96 = x60 + x95;
auto x97 = -x14;
auto x98 = x66 + x97;
auto x99 = x69 + x97;
auto x100 = -x60;
auto x101 = x51 + x69;
fNeq[0] = x51*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + cell[0] + V{0.333333333333333};
fNeq[1] = x52*(x50 + x60 - x62) + cell[1] + V{0.0555555555555556};
fNeq[2] = x52*(x49 + x66 - x67 + x68) + cell[2] + V{0.0555555555555556};
fNeq[3] = x52*(x22 + x68 + x69 - x70) + cell[3] + V{0.0555555555555556};
fNeq[4] = x71*(x66 - x75 + x76) + cell[4] + V{0.0277777777777778};
fNeq[5] = x71*(x76 + x77 + x80) + cell[5] + V{0.0277777777777778};
fNeq[6] = x71*(x69 + x76 - x83) + cell[6] + V{0.0277777777777778};
fNeq[7] = x71*(x76 + x84 + x87) + cell[7] + V{0.0277777777777778};
fNeq[8] = x71*(x69 - x89 + x90) + cell[8] + V{0.0277777777777778};
fNeq[9] = x71*(x84 + x90 + x92) + cell[9] + V{0.0277777777777778};
fNeq[10] = -x52*(x62 + x96) + cell[10] + V{0.0555555555555556};
fNeq[11] = -x52*(x67 + x94 + x98) + cell[11] + V{0.0555555555555556};
fNeq[12] = -x52*(x70 + x93 + x99 + V{1}) + cell[12] + V{0.0555555555555556};
fNeq[13] = -x71*(x75 + x96 + x98) + cell[13] + V{0.0277777777777778};
fNeq[14] = x71*(x100 + x80 + x90) + cell[14] + V{0.0277777777777778};
fNeq[15] = -x71*(x83 + x96 + x99) + cell[15] + V{0.0277777777777778};
fNeq[16] = x71*(x100 + x101 + x87) + cell[16] + V{0.0277777777777778};
fNeq[17] = -x71*(x69 + x89 + x95 + x98) + cell[17] + V{0.0277777777777778};
fNeq[18] = x71*(x101 + x77 + x92) + cell[18] + V{0.0277777777777778};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = u[0]*u[0];
auto x21 = V{1.5}*x20;
auto x22 = u[1]*u[1];
auto x23 = V{1.5}*x22;
auto x24 = u[2]*u[2];
auto x25 = V{1.5}*x24;
auto x26 = x23 + x25 + V{-1};
auto x27 = x21 + x26;
auto x28 = V{0.0555555555555556}*omega;
auto x29 = V{3}*u[0];
auto x30 = V{3}*x20;
auto x31 = V{3}*u[1];
auto x32 = V{3}*x22;
auto x33 = x21 + V{-1};
auto x34 = V{3}*u[2];
auto x35 = V{3}*x24;
auto x36 = V{0.0277777777777778}*omega;
auto x37 = u[0] + u[1];
auto x38 = V{4.5}*(x37*x37);
auto x39 = x27 + x29;
auto x40 = -x31;
auto x41 = u[0] - u[1];
auto x42 = -x41;
auto x43 = u[0] + u[2];
auto x44 = V{4.5}*(x43*x43);
auto x45 = -x34;
auto x46 = -u[2];
auto x47 = x46 + u[0];
auto x48 = -x47;
auto x49 = u[1] + u[2];
auto x50 = V{4.5}*(x49*x49);
auto x51 = x27 + x31;
auto x52 = x46 + u[1];
auto x53 = -x52;
auto x54 = -x23;
auto x55 = V{1} - x25;
auto x56 = x54 + x55;
auto x57 = x29 + x56;
auto x58 = -x21;
auto x59 = x31 + x58;
auto x60 = x34 + x58;
auto x61 = -x29;
auto x62 = x27 + x34;
cell[0] = -V{0.333333333333333}*omega*(rho*x27 + V{1}) - x19*cell[0];
cell[1] = -x19*cell[1] - x28*(rho*(x26 + x29 - x30) + V{1});
cell[2] = -x19*cell[2] - x28*(rho*(x25 + x31 - x32 + x33) + V{1});
cell[3] = -x19*cell[3] - x28*(rho*(x23 + x33 + x34 - x35) + V{1});
cell[4] = -x19*cell[4] - x36*(rho*(x31 - x38 + x39) + V{1});
cell[5] = -x19*cell[5] - x36*(rho*(x39 + x40 - V{4.5}*x42*x42) + V{1});
cell[6] = -x19*cell[6] - x36*(rho*(x34 + x39 - x44) + V{1});
cell[7] = -x19*cell[7] - x36*(rho*(x39 + x45 - V{4.5}*x48*x48) + V{1});
cell[8] = -x19*cell[8] - x36*(rho*(x34 - x50 + x51) + V{1});
cell[9] = -x19*cell[9] - x36*(rho*(x45 + x51 - V{4.5}*x53*x53) + V{1});
cell[10] = -x19*cell[10] + x28*(rho*(x30 + x57) + V{-1});
cell[11] = -x19*cell[11] + x28*(rho*(x32 + x55 + x59) + V{-1});
cell[12] = -x19*cell[12] + x28*(rho*(x35 + x54 + x60 + V{1}) + V{-1});
cell[13] = -x19*cell[13] + x36*(rho*(x38 + x57 + x59) + V{-1});
cell[14] = -x19*cell[14] - x36*(rho*(x51 + x61 - V{4.5}*x41*x41) + V{1});
cell[15] = -x19*cell[15] + x36*(rho*(x44 + x57 + x60) + V{-1});
cell[16] = -x19*cell[16] - x36*(rho*(x61 + x62 - V{4.5}*x47*x47) + V{1});
cell[17] = -x19*cell[17] + x36*(rho*(x34 + x50 + x56 + x59) + V{-1});
cell[18] = -x19*cell[18] - x36*(rho*(x40 + x62 - V{4.5}*x52*x52) + V{1});
return x20 + x22 + x24;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = V{3}*u[0];
auto x21 = x20 + V{-1};
auto x22 = V{0.0555555555555556}*omega;
auto x23 = V{3}*u[1];
auto x24 = x23 + V{-1};
auto x25 = V{3}*u[2];
auto x26 = V{0.0277777777777778}*omega;
auto x27 = -x20;
auto x28 = x23 + V{1};
auto x29 = x25 + V{1};
auto x30 = -x23;
auto x31 = x20 + V{1};
auto x32 = -x25;
cell[0] = V{0.333333333333333}*omega*(rho + V{-1}) - x19*cell[0];
cell[1] = -x19*cell[1] - x22*(rho*x21 + V{1});
cell[2] = -x19*cell[2] - x22*(rho*x24 + V{1});
cell[3] = -x19*cell[3] - x22*(rho*(x25 + V{-1}) + V{1});
cell[4] = -x19*cell[4] - x26*(rho*(x21 + x23) + V{1});
cell[5] = -x19*cell[5] + x26*(rho*(x27 + x28) + V{-1});
cell[6] = -x19*cell[6] - x26*(rho*(x21 + x25) + V{1});
cell[7] = -x19*cell[7] + x26*(rho*(x27 + x29) + V{-1});
cell[8] = -x19*cell[8] - x26*(rho*(x24 + x25) + V{1});
cell[9] = -x19*cell[9] + x26*(rho*(x29 + x30) + V{-1});
cell[10] = -x19*cell[10] + x22*(rho*x31 + V{-1});
cell[11] = -x19*cell[11] + x22*(rho*x28 + V{-1});
cell[12] = -x19*cell[12] + x22*(rho*x29 + V{-1});
cell[13] = -x19*cell[13] + x26*(rho*(x23 + x31) + V{-1});
cell[14] = -x19*cell[14] + x26*(rho*(x30 + x31) + V{-1});
cell[15] = -x19*cell[15] + x26*(rho*(x25 + x31) + V{-1});
cell[16] = -x19*cell[16] + x26*(rho*(x31 + x32) + V{-1});
cell[17] = -x19*cell[17] + x26*(rho*(x25 + x28) + V{-1});
cell[18] = -x19*cell[18] + x26*(rho*(x28 + x32) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x19 = j[0]*j[0];
auto x20 = j[1]*j[1];
auto x21 = j[2]*j[2];
auto x22 = omega + V{-1};
auto x23 = V{0.0833333333333333}*x20;
auto x24 = V{0.166666666666667}*j[0];
auto x25 = V{0.166666666666667}*x19;
auto x26 = V{0.0833333333333333}*x21;
auto x27 = V{0.166666666666667}*pressure;
auto x28 = -x27;
auto x29 = x26 + x28 + V{0.0555555555555556};
auto x30 = V{0.0833333333333333}*x19;
auto x31 = V{0.166666666666667}*j[1];
auto x32 = V{0.166666666666667}*x20;
auto x33 = V{0.166666666666667}*j[2];
auto x34 = V{0.166666666666667}*x21;
auto x35 = j[0] + j[1];
auto x36 = V{0.125}*(x35*x35);
auto x37 = V{0.0833333333333333}*j[0];
auto x38 = V{0.0833333333333333}*j[1];
auto x39 = x37 + x38;
auto x40 = V{0.0416666666666667}*x19;
auto x41 = V{0.0416666666666667}*x20;
auto x42 = V{0.0416666666666667}*x21;
auto x43 = V{0.0833333333333333}*pressure;
auto x44 = x40 + x41 + x42 - x43 + V{0.0277777777777778};
auto x45 = j[0] - j[1];
auto x46 = -x45;
auto x47 = -x38;
auto x48 = x37 + x44;
auto x49 = V{0.0833333333333333}*j[2];
auto x50 = j[0] + j[2];
auto x51 = V{0.125}*(x50*x50);
auto x52 = -j[2];
auto x53 = x52 + j[0];
auto x54 = -x53;
auto x55 = -x49;
auto x56 = j[1] + j[2];
auto x57 = V{0.125}*(x56*x56);
auto x58 = x38 + x44;
auto x59 = x52 + j[1];
auto x60 = -x59;
auto x61 = -x23;
auto x62 = -x26 + x27 + V{-0.0555555555555556};
auto x63 = -x30;
auto x64 = -x40 - x41 - x42 + x43 + V{-0.0277777777777778};
auto x65 = -x37;
auto x66 = x49 + x64;
auto x67 = x44 + x49;
cell[0] = -omega*(-V{1}*pressure + V{0.5}*x19 + V{0.5}*x20 + V{0.5}*x21 + V{0.333333333333333}) - x22*cell[0];
cell[1] = -omega*(x23 + x24 - x25 + x29) - x22*cell[1];
cell[2] = -omega*(x29 + x30 + x31 - x32) - x22*cell[2];
cell[3] = -omega*(x23 + x28 + x30 + x33 - x34 + V{0.0555555555555556}) - x22*cell[3];
cell[4] = -omega*(-x36 + x39 + x44) - x22*cell[4];
cell[5] = -omega*(x47 + x48 - V{0.125}*x46*x46) - x22*cell[5];
cell[6] = -omega*(x48 + x49 - x51) - x22*cell[6];
cell[7] = -omega*(x48 + x55 - V{0.125}*x54*x54) - x22*cell[7];
cell[8] = -omega*(x49 - x57 + x58) - x22*cell[8];
cell[9] = -omega*(x55 + x58 - V{0.125}*x60*x60) - x22*cell[9];
cell[10] = omega*(x24 + x25 + x61 + x62) - x22*cell[10];
cell[11] = omega*(x31 + x32 + x62 + x63) - x22*cell[11];
cell[12] = omega*(x27 + x33 + x34 + x61 + x63 + V{-0.0555555555555556}) - x22*cell[12];
cell[13] = omega*(x36 + x39 + x64) - x22*cell[13];
cell[14] = -omega*(x58 + x65 - V{0.125}*x45*x45) - x22*cell[14];
cell[15] = omega*(x37 + x51 + x66) - x22*cell[15];
cell[16] = -omega*(x65 + x67 - V{0.125}*x53*x53) - x22*cell[16];
cell[17] = omega*(x38 + x57 + x66) - x22*cell[17];
cell[18] = -omega*(x47 + x67 - V{0.125}*x59*x59) - x22*cell[18];
return x19 + x20 + x21;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = V{0.333333333333333}*rho;
auto x21 = u[0]*u[0];
auto x22 = V{1.5}*x21;
auto x23 = u[1]*u[1];
auto x24 = V{1.5}*x23;
auto x25 = u[2]*u[2];
auto x26 = V{1.5}*x25;
auto x27 = x24 + x26 + V{-1};
auto x28 = x22 + x27;
auto x29 = V{0.0555555555555556}*rho;
auto x30 = V{3}*u[0];
auto x31 = V{3}*x21;
auto x32 = x27 + x30 - x31;
auto x33 = ratioRho*x29;
auto x34 = V{3}*u[1];
auto x35 = V{3}*x23;
auto x36 = x22 + V{-1};
auto x37 = x26 + x34 - x35 + x36;
auto x38 = V{3}*u[2];
auto x39 = V{3}*x25;
auto x40 = x24 + x36 + x38 - x39;
auto x41 = V{0.0277777777777778}*rho;
auto x42 = u[0] + u[1];
auto x43 = V{4.5}*(x42*x42);
auto x44 = x28 + x30;
auto x45 = x34 - x43 + x44;
auto x46 = ratioRho*x41;
auto x47 = -x34;
auto x48 = u[0] - u[1];
auto x49 = -x48;
auto x50 = x44 + x47 - V{4.5}*x49*x49;
auto x51 = u[0] + u[2];
auto x52 = V{4.5}*(x51*x51);
auto x53 = x38 + x44 - x52;
auto x54 = -x38;
auto x55 = -u[2];
auto x56 = x55 + u[0];
auto x57 = -x56;
auto x58 = x44 + x54 - V{4.5}*x57*x57;
auto x59 = u[1] + u[2];
auto x60 = V{4.5}*(x59*x59);
auto x61 = x28 + x34;
auto x62 = x38 - x60 + x61;
auto x63 = x55 + u[1];
auto x64 = -x63;
auto x65 = x54 + x61 - V{4.5}*x64*x64;
auto x66 = -x24;
auto x67 = V{1} - x26;
auto x68 = x66 + x67;
auto x69 = x30 + x68;
auto x70 = x31 + x69;
auto x71 = -x22;
auto x72 = x34 + x71;
auto x73 = x35 + x67 + x72;
auto x74 = x38 + x71;
auto x75 = x39 + x66 + x74 + V{1};
auto x76 = x43 + x69 + x72;
auto x77 = -x30;
auto x78 = x61 + x77 - V{4.5}*x48*x48;
auto x79 = x52 + x69 + x74;
auto x80 = x28 + x38;
auto x81 = x77 + x80 - V{4.5}*x56*x56;
auto x82 = x38 + x60 + x68 + x72;
auto x83 = x47 + x80 - V{4.5}*x63*x63;
cell[0] = -ratioRho*x20*x28 - x19*(x20*x28 + cell[0] + V{0.333333333333333}) + V{-0.333333333333333};
cell[1] = -x19*(x29*x32 + cell[1] + V{0.0555555555555556}) - x32*x33 + V{-0.0555555555555556};
cell[2] = -x19*(x29*x37 + cell[2] + V{0.0555555555555556}) - x33*x37 + V{-0.0555555555555556};
cell[3] = -x19*(x29*x40 + cell[3] + V{0.0555555555555556}) - x33*x40 + V{-0.0555555555555556};
cell[4] = -x19*(x41*x45 + cell[4] + V{0.0277777777777778}) - x45*x46 + V{-0.0277777777777778};
cell[5] = -x19*(x41*x50 + cell[5] + V{0.0277777777777778}) - x46*x50 + V{-0.0277777777777778};
cell[6] = -x19*(x41*x53 + cell[6] + V{0.0277777777777778}) - x46*x53 + V{-0.0277777777777778};
cell[7] = -x19*(x41*x58 + cell[7] + V{0.0277777777777778}) - x46*x58 + V{-0.0277777777777778};
cell[8] = -x19*(x41*x62 + cell[8] + V{0.0277777777777778}) - x46*x62 + V{-0.0277777777777778};
cell[9] = -x19*(x41*x65 + cell[9] + V{0.0277777777777778}) - x46*x65 + V{-0.0277777777777778};
cell[10] = V{0.0555555555555556}*ratioRho*rho*x70 - x19*(-x29*x70 + cell[10] + V{0.0555555555555556}) + V{-0.0555555555555556};
cell[11] = V{0.0555555555555556}*ratioRho*rho*x73 - x19*(-x29*x73 + cell[11] + V{0.0555555555555556}) + V{-0.0555555555555556};
cell[12] = V{0.0555555555555556}*ratioRho*rho*x75 - x19*(-x29*x75 + cell[12] + V{0.0555555555555556}) + V{-0.0555555555555556};
cell[13] = V{0.0277777777777778}*ratioRho*rho*x76 - x19*(-x41*x76 + cell[13] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[14] = -x19*(x41*x78 + cell[14] + V{0.0277777777777778}) - x46*x78 + V{-0.0277777777777778};
cell[15] = V{0.0277777777777778}*ratioRho*rho*x79 - x19*(-x41*x79 + cell[15] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[16] = -x19*(x41*x81 + cell[16] + V{0.0277777777777778}) - x46*x81 + V{-0.0277777777777778};
cell[17] = V{0.0277777777777778}*ratioRho*rho*x82 - x19*(-x41*x82 + cell[17] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[18] = -x19*(x41*x83 + cell[18] + V{0.0277777777777778}) - x46*x83 + V{-0.0277777777777778};
return x21 + x23 + x25;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = V{3}*u[0];
auto x21 = x20 + V{1};
auto x22 = V{0.00925925925925926}*rho;
auto x23 = x20 + V{-1};
auto x24 = V{0.166666666666667}*cell[14];
auto x25 = V{3}*u[1];
auto x26 = -x25;
auto x27 = x21 + x26;
auto x28 = V{0.00462962962962963}*rho;
auto x29 = -x20;
auto x30 = x25 + V{1};
auto x31 = x29 + x30;
auto x32 = x28*x31;
auto x33 = -x24 + x27*x28 - x32 + V{0.166666666666667}*cell[5];
auto x34 = V{0.166666666666667}*cell[16];
auto x35 = V{3}*u[2];
auto x36 = -x35;
auto x37 = x21 + x36;
auto x38 = x35 + V{1};
auto x39 = x29 + x38;
auto x40 = x28*x39;
auto x41 = x28*x37 - x34 - x40 + V{0.166666666666667}*cell[7];
auto x42 = V{0.166666666666667}*cell[4];
auto x43 = V{0.166666666666667}*cell[13];
auto x44 = x21 + x25;
auto x45 = x28*x44;
auto x46 = x23 + x25;
auto x47 = x28*x46 + x42 - x43 + x45;
auto x48 = V{0.166666666666667}*cell[6];
auto x49 = V{0.166666666666667}*cell[15];
auto x50 = x21 + x35;
auto x51 = x28*x50;
auto x52 = x23 + x35;
auto x53 = x28*x52 + x48 - x49 + x51;
auto x54 = x19*(x21*x22 + x22*x23 + x33 + x41 + x47 + x53 - V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{4.62592926927149e-18});
auto x55 = V{0.0555555555555556}*rho;
auto x56 = V{0.166666666666667}*cell[2];
auto x57 = -x46;
auto x58 = x25 + V{-1};
auto x59 = -x58;
auto x60 = x22*x30;
auto x61 = V{0.166666666666667}*cell[8];
auto x62 = x35 + x58;
auto x63 = x30 + x35;
auto x64 = x28*x63;
auto x65 = -x28*x62 - x61 - x64 + V{0.166666666666667}*cell[17];
auto x66 = V{0.166666666666667}*cell[18];
auto x67 = V{0.166666666666667}*cell[9];
auto x68 = x26 + x38;
auto x69 = x28*x68;
auto x70 = x30 + x36;
auto x71 = x28*x70;
auto x72 = x66 - x67 + x69 - x71;
auto x73 = V{0.166666666666667}*cell[3];
auto x74 = -x52;
auto x75 = x35 + V{-1};
auto x76 = -x75;
auto x77 = x22*x38;
auto x78 = -x66 + x67 - x69 + x71;
auto x79 = V{0.0833333333333333}*cell[2];
auto x80 = V{0.0833333333333333}*cell[11];
auto x81 = -x80;
auto x82 = x28*x30;
auto x83 = x28*x58;
auto x84 = V{0.0833333333333333}*cell[7];
auto x85 = V{0.0833333333333333}*cell[16];
auto x86 = V{0.00231481481481481}*rho;
auto x87 = x37*x86;
auto x88 = x39*x86;
auto x89 = x84 - x85 + x87 - x88;
auto x90 = x79 + x81 + x82 + x83 + x89;
auto x91 = V{0.0833333333333333}*cell[1];
auto x92 = -V{0.0833333333333333}*cell[10];
auto x93 = x21*x28;
auto x94 = x23*x28;
auto x95 = V{0.0833333333333333}*cell[6];
auto x96 = V{0.0833333333333333}*cell[15];
auto x97 = x50*x86;
auto x98 = x52*x86;
auto x99 = x95 - x96 + x97 + x98;
auto x100 = x91 + x92 + x93 + x94 + x99 + V{2.31296463463574e-18};
auto x101 = V{0.0833333333333333}*cell[8];
auto x102 = V{0.0833333333333333}*cell[17];
auto x103 = x63*x86;
auto x104 = x62*x86;
auto x105 = x101 - x102 + x103 + x104;
auto x106 = V{0.0833333333333333}*cell[9];
auto x107 = V{0.0833333333333333}*cell[18];
auto x108 = x70*x86;
auto x109 = x68*x86;
auto x110 = x106 - x107 + x108 - x109;
auto x111 = x19*(x100 + x105 + x110 + x47 + x90);
auto x112 = V{0.0277777777777778}*rho;
auto x113 = -x79 + x80 - x82 + x89;
auto x114 = -x106 + x107 - x108 + x109;
auto x115 = -x101 + x102 - x103 - x104;
auto x116 = x19*(x100 + x113 + x114 + x115 + x33 - x83);
auto x117 = V{0.0833333333333333}*cell[14];
auto x118 = x31*x86;
auto x119 = -x117 - x118 + x27*x86 + V{0.0833333333333333}*cell[5];
auto x120 = V{0.0833333333333333}*cell[3];
auto x121 = V{0.0833333333333333}*cell[12];
auto x122 = x28*x38;
auto x123 = x28*x75;
auto x124 = V{0.0833333333333333}*cell[4];
auto x125 = V{0.0833333333333333}*cell[13];
auto x126 = x44*x86;
auto x127 = x124 - x125 + x126 + x46*x86;
auto x128 = x120 - x121 + x122 + x123 + x127;
auto x129 = x91 + x92 + x93 + x94 + V{2.31296463463574e-18};
auto x130 = x19*(x105 + x114 + x119 + x128 + x129 + x53);
auto x131 = -x120;
auto x132 = -x122;
auto x133 = x119 + x121 + x131 + x132;
auto x134 = -x123 + x127;
auto x135 = x19*(x110 + x115 + x129 + x133 + x134 + x41);
auto x136 = -x95 + x96 - x97;
auto x137 = -V{0.00231481481481481}*rho*x27 + x117 + x118 - V{0.0833333333333333}*cell[5];
auto x138 = x19*(-x121 - x131 - x132 - x134 - x136 - x137 - x78 - x90 + x98);
auto x139 = x28*x62 + x61 + x64 - V{0.166666666666667}*cell[17];
cell[0] = V{0.333333333333333}*rho + V{-0.333333333333333};
cell[1] = -x23*x55 - x54 + V{-0.0555555555555556};
cell[2] = x19*(x22*x59 + x28*x57 + x33 - x42 + x43 - x45 - x56 - x60 + x65 + x72 + V{0.166666666666667}*cell[11]) - x55*x58 + V{-0.0555555555555556};
cell[3] = x19*(x22*x76 + x28*x74 + x41 - x48 + x49 - x51 + x65 - x73 - x77 + x78 + V{0.166666666666667}*cell[12]) - x55*x75 + V{-0.0555555555555556};
cell[4] = -x111 - x112*x46 + V{-0.0277777777777778};
cell[5] = V{0.0277777777777778}*rho*x31 - x116 + V{-0.0277777777777778};
cell[6] = -x112*x52 - x130 + V{-0.0277777777777778};
cell[7] = V{0.0277777777777778}*rho*x39 - x135 + V{-0.0277777777777778};
cell[8] = -x112*x62 + x19*(x113 - x124 + x125 - x126 + x133 + x136 + x28*x59 + x28*x76 + x57*x86 + x65 + x74*x86) + V{-0.0277777777777778};
cell[9] = x112*x68 + x138 + V{-0.0277777777777778};
cell[10] = x21*x55 + x54 + V{-0.0555555555555556};
cell[11] = V{0.0555555555555556}*rho*x30 - x19*(V{0.00462962962962963}*rho*x27 - x139 - x22*x58 - x24 - x32 - x47 - x56 - x60 - x78 + V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[5]) + V{-0.0555555555555556};
cell[12] = V{0.0555555555555556}*rho*x38 - x19*(V{0.00462962962962963}*rho*x37 - x139 - x22*x75 - x34 - x40 - x53 - x72 - x73 - x77 + V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[7]) + V{-0.0555555555555556};
cell[13] = x111 + x112*x44 + V{-0.0277777777777778};
cell[14] = x112*x27 + x116 + V{-0.0277777777777778};
cell[15] = x112*x50 + x130 + V{-0.0277777777777778};
cell[16] = x112*x37 + x135 + V{-0.0277777777777778};
cell[17] = V{0.0277777777777778}*rho*x63 - x19*(-x128 - x137 - x139 - x79 - x81 - x82 - x83 + x84 - x85 + x87 - x88 - x99) + V{-0.0277777777777778};
cell[18] = V{0.0277777777777778}*rho*x70 - x138 + V{-0.0277777777777778};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x19 = u[0]*u[0];
auto x20 = V{1.5}*x19;
auto x21 = u[1]*u[1];
auto x22 = V{1.5}*x21;
auto x23 = u[2]*u[2];
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = omega + V{-1};
auto x28 = V{0.0833333333333333}*pi[3];
auto x29 = V{0.0833333333333333}*pi[5];
auto x30 = x28 + x29 - V{0.166666666666667}*pi[0];
auto x31 = V{0.0555555555555556}*rho;
auto x32 = V{3}*u[0];
auto x33 = V{3}*x19;
auto x34 = V{0.0833333333333333}*pi[0];
auto x35 = x29 + x34 - V{0.166666666666667}*pi[3];
auto x36 = V{3}*u[1];
auto x37 = V{3}*x21;
auto x38 = x20 + V{-1};
auto x39 = x28 + x34 - V{0.166666666666667}*pi[5];
auto x40 = V{3}*u[2];
auto x41 = V{3}*x23;
auto x42 = V{0.0277777777777778}*rho;
auto x43 = u[0] + u[1];
auto x44 = V{4.5}*(x43*x43);
auto x45 = x26 + x32;
auto x46 = V{0.25}*pi[1];
auto x47 = V{0.0833333333333333}*pi[0];
auto x48 = V{0.0833333333333333}*pi[3];
auto x49 = x47 + x48 - V{0.0416666666666667}*pi[5];
auto x50 = x27*(x46 + x49) + V{0.0277777777777778};
auto x51 = -x36;
auto x52 = u[0] - u[1];
auto x53 = -x52;
auto x54 = x27*(-x46 + x49) + V{0.0277777777777778};
auto x55 = u[0] + u[2];
auto x56 = V{4.5}*(x55*x55);
auto x57 = V{0.25}*pi[2];
auto x58 = V{0.0833333333333333}*pi[5];
auto x59 = x47 + x58 - V{0.0416666666666667}*pi[3];
auto x60 = x27*(x57 + x59) + V{0.0277777777777778};
auto x61 = -x40;
auto x62 = -u[2];
auto x63 = x62 + u[0];
auto x64 = -x63;
auto x65 = x27*(-x57 + x59) + V{0.0277777777777778};
auto x66 = u[1] + u[2];
auto x67 = V{4.5}*(x66*x66);
auto x68 = x26 + x36;
auto x69 = V{0.25}*pi[4];
auto x70 = V{0.0416666666666667}*pi[0];
auto x71 = x27*(x48 + x58 + x69 - x70) + V{0.0277777777777778};
auto x72 = x62 + u[1];
auto x73 = -x72;
auto x74 = -x27*(-x48 - x58 + x69 + x70) + V{0.0277777777777778};
auto x75 = -x22;
auto x76 = V{1} - x24;
auto x77 = x75 + x76;
auto x78 = x32 + x77;
auto x79 = -x20;
auto x80 = x36 + x79;
auto x81 = x40 + x79;
auto x82 = -x32;
auto x83 = x26 + x40;
cell[0] = -V{0.333333333333333}*rho*x26 + V{0.5}*x27*(pi[0] + pi[3] + pi[5]) + V{-0.333333333333333};
cell[1] = x27*x30 - x31*(x25 + x32 - x33) + V{-0.0555555555555556};
cell[2] = x27*x35 - x31*(x24 + x36 - x37 + x38) + V{-0.0555555555555556};
cell[3] = x27*x39 - x31*(x22 + x38 + x40 - x41) + V{-0.0555555555555556};
cell[4] = -x42*(x36 - x44 + x45) - x50;
cell[5] = -x42*(x45 + x51 - V{4.5}*x53*x53) - x54;
cell[6] = -x42*(x40 + x45 - x56) - x60;
cell[7] = -x42*(x45 + x61 - V{4.5}*x64*x64) - x65;
cell[8] = -x42*(x40 - x67 + x68) - x71;
cell[9] = -x42*(x61 + x68 - V{4.5}*x73*x73) - x74;
cell[10] = x27*x30 + x31*(x33 + x78) + V{-0.0555555555555556};
cell[11] = x27*x35 + x31*(x37 + x76 + x80) + V{-0.0555555555555556};
cell[12] = x27*x39 + x31*(x41 + x75 + x81 + V{1}) + V{-0.0555555555555556};
cell[13] = V{0.0277777777777778}*rho*(x44 + x78 + x80) - x50;
cell[14] = -x42*(x68 + x82 - V{4.5}*x52*x52) - x54;
cell[15] = V{0.0277777777777778}*rho*(x56 + x78 + x81) - x60;
cell[16] = -x42*(x82 + x83 - V{4.5}*x63*x63) - x65;
cell[17] = V{0.0277777777777778}*rho*(x40 + x67 + x77 + x80) - x71;
cell[18] = -x42*(x51 + x83 - V{4.5}*x72*x72) - x74;
return x19 + x21 + x23;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x19 = V{3}*newU[0];
auto x20 = x19 + V{-1};
auto x21 = V{3}*newU[1];
auto x22 = x21 + V{-1};
auto x23 = V{3}*newU[2];
auto x24 = -x19;
auto x25 = x21 + V{1};
auto x26 = x23 + V{1};
auto x27 = -x21;
auto x28 = x19 + V{1};
auto x29 = -x23;
cell[0] = V{0.333333333333333}*newRho + V{-0.333333333333333};
cell[1] = -V{0.0555555555555556}*newRho*x20 + V{-0.0555555555555556};
cell[2] = -V{0.0555555555555556}*newRho*x22 + V{-0.0555555555555556};
cell[3] = -V{0.0555555555555556}*newRho*(x23 + V{-1}) + V{-0.0555555555555556};
cell[4] = -V{0.0277777777777778}*newRho*(x20 + x21) + V{-0.0277777777777778};
cell[5] = V{0.0277777777777778}*newRho*(x24 + x25) + V{-0.0277777777777778};
cell[6] = -V{0.0277777777777778}*newRho*(x20 + x23) + V{-0.0277777777777778};
cell[7] = V{0.0277777777777778}*newRho*(x24 + x26) + V{-0.0277777777777778};
cell[8] = -V{0.0277777777777778}*newRho*(x22 + x23) + V{-0.0277777777777778};
cell[9] = V{0.0277777777777778}*newRho*(x26 + x27) + V{-0.0277777777777778};
cell[10] = V{0.0555555555555556}*newRho*x28 + V{-0.0555555555555556};
cell[11] = V{0.0555555555555556}*newRho*x25 + V{-0.0555555555555556};
cell[12] = V{0.0555555555555556}*newRho*x26 + V{-0.0555555555555556};
cell[13] = V{0.0277777777777778}*newRho*(x21 + x28) + V{-0.0277777777777778};
cell[14] = V{0.0277777777777778}*newRho*(x27 + x28) + V{-0.0277777777777778};
cell[15] = V{0.0277777777777778}*newRho*(x23 + x28) + V{-0.0277777777777778};
cell[16] = V{0.0277777777777778}*newRho*(x28 + x29) + V{-0.0277777777777778};
cell[17] = V{0.0277777777777778}*newRho*(x23 + x25) + V{-0.0277777777777778};
cell[18] = V{0.0277777777777778}*newRho*(x25 + x29) + V{-0.0277777777777778};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x19 = oldU[0]*oldU[0];
auto x20 = V{1.5}*x19;
auto x21 = oldU[1]*oldU[1];
auto x22 = V{1.5}*x21;
auto x23 = oldU[2]*oldU[2];
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = newU[0]*newU[0];
auto x28 = V{1.5}*x27;
auto x29 = newU[1]*newU[1];
auto x30 = V{1.5}*x29;
auto x31 = newU[2]*newU[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = V{0.0555555555555556}*oldRho;
auto x36 = V{3}*oldU[0];
auto x37 = V{3}*x19;
auto x38 = V{0.0555555555555556}*newRho;
auto x39 = V{3}*newU[0];
auto x40 = V{3}*x27;
auto x41 = V{3}*oldU[1];
auto x42 = V{3}*x21;
auto x43 = x20 + V{-1};
auto x44 = V{3}*newU[1];
auto x45 = V{3}*x29;
auto x46 = x28 + V{-1};
auto x47 = V{3}*oldU[2];
auto x48 = V{3}*x23;
auto x49 = V{3}*newU[2];
auto x50 = V{3}*x31;
auto x51 = V{0.0277777777777778}*oldRho;
auto x52 = oldU[0] + oldU[1];
auto x53 = V{4.5}*(x52*x52);
auto x54 = x26 + x36;
auto x55 = V{0.0277777777777778}*newRho;
auto x56 = newU[0] + newU[1];
auto x57 = V{4.5}*(x56*x56);
auto x58 = x34 + x39;
auto x59 = -x41;
auto x60 = oldU[0] - oldU[1];
auto x61 = -V{4.5}*x60*x60;
auto x62 = -x44;
auto x63 = newU[0] - newU[1];
auto x64 = -V{4.5}*x63*x63;
auto x65 = oldU[0] + oldU[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = newU[0] + newU[2];
auto x68 = V{4.5}*(x67*x67);
auto x69 = -x47;
auto x70 = -oldU[2];
auto x71 = x70 + oldU[0];
auto x72 = -V{4.5}*x71*x71;
auto x73 = -x49;
auto x74 = -newU[2];
auto x75 = x74 + newU[0];
auto x76 = -V{4.5}*x75*x75;
auto x77 = oldU[1] + oldU[2];
auto x78 = V{4.5}*(x77*x77);
auto x79 = x26 + x41;
auto x80 = newU[1] + newU[2];
auto x81 = V{4.5}*(x80*x80);
auto x82 = x34 + x44;
auto x83 = x70 + oldU[1];
auto x84 = -V{4.5}*x83*x83;
auto x85 = x74 + newU[1];
auto x86 = -V{4.5}*x85*x85;
auto x87 = -x30;
auto x88 = V{1} - x32;
auto x89 = x87 + x88;
auto x90 = x39 + x89;
auto x91 = -x22;
auto x92 = V{1} - x24;
auto x93 = x91 + x92;
auto x94 = x36 + x93;
auto x95 = -x28;
auto x96 = x44 + x95;
auto x97 = -x20;
auto x98 = x41 + x97;
auto x99 = x49 + x95;
auto x100 = x47 + x97;
auto x101 = -x36;
auto x102 = -x39;
auto x103 = x26 + x47;
auto x104 = x34 + x49;
cell[0] = -V{0.333333333333333}*newRho*x34 + V{0.333333333333333}*oldRho*x26 + cell[0];
cell[1] = x35*(x25 + x36 - x37) - x38*(x33 + x39 - x40) + cell[1];
cell[2] = x35*(x24 + x41 - x42 + x43) - x38*(x32 + x44 - x45 + x46) + cell[2];
cell[3] = x35*(x22 + x43 + x47 - x48) - x38*(x30 + x46 + x49 - x50) + cell[3];
cell[4] = x51*(x41 - x53 + x54) - x55*(x44 - x57 + x58) + cell[4];
cell[5] = x51*(x54 + x59 + x61) - x55*(x58 + x62 + x64) + cell[5];
cell[6] = x51*(x47 + x54 - x66) - x55*(x49 + x58 - x68) + cell[6];
cell[7] = x51*(x54 + x69 + x72) - x55*(x58 + x73 + x76) + cell[7];
cell[8] = x51*(x47 - x78 + x79) - x55*(x49 - x81 + x82) + cell[8];
cell[9] = x51*(x69 + x79 + x84) - x55*(x73 + x82 + x86) + cell[9];
cell[10] = -x35*(x37 + x94) + x38*(x40 + x90) + cell[10];
cell[11] = -x35*(x42 + x92 + x98) + x38*(x45 + x88 + x96) + cell[11];
cell[12] = -x35*(x100 + x48 + x91 + V{1}) + x38*(x50 + x87 + x99 + V{1}) + cell[12];
cell[13] = -x51*(x53 + x94 + x98) + x55*(x57 + x90 + x96) + cell[13];
cell[14] = x51*(x101 + x61 + x79) - x55*(x102 + x64 + x82) + cell[14];
cell[15] = -x51*(x100 + x66 + x94) + x55*(x68 + x90 + x99) + cell[15];
cell[16] = x51*(x101 + x103 + x72) - x55*(x102 + x104 + x76) + cell[16];
cell[17] = -x51*(x47 + x78 + x93 + x98) + x55*(x49 + x81 + x89 + x96) + cell[17];
cell[18] = x51*(x103 + x59 + x84) - x55*(x104 + x62 + x86) + cell[18];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x19 = u[0]*u[0];
auto x20 = V{1.5}*x19;
auto x21 = u[1]*u[1];
auto x22 = V{1.5}*x21;
auto x23 = u[2]*u[2];
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = V{0.0555555555555556}*rho;
auto x28 = V{3}*u[0];
auto x29 = V{3}*x19;
auto x30 = V{0.0833333333333333}*pi[3];
auto x31 = V{0.0833333333333333}*pi[5] + V{0.0555555555555556};
auto x32 = x30 + x31 - V{0.166666666666667}*pi[0];
auto x33 = V{3}*u[1];
auto x34 = V{3}*x21;
auto x35 = x20 + V{-1};
auto x36 = V{0.0833333333333333}*pi[0];
auto x37 = x31 + x36 - V{0.166666666666667}*pi[3];
auto x38 = V{3}*u[2];
auto x39 = V{3}*x23;
auto x40 = x30 + x36 - V{0.166666666666667}*pi[5] + V{0.0555555555555556};
auto x41 = V{0.25}*pi[1];
auto x42 = V{0.0277777777777778}*rho;
auto x43 = u[0] + u[1];
auto x44 = V{4.5}*(x43*x43);
auto x45 = x26 + x28;
auto x46 = V{0.0416666666666667}*pi[5];
auto x47 = -V{0.0833333333333333}*pi[3];
auto x48 = V{0.0277777777777778} - V{0.0833333333333333}*pi[0];
auto x49 = x46 + x47 + x48;
auto x50 = -x33;
auto x51 = u[0] - u[1];
auto x52 = -x51;
auto x53 = x41 + x49;
auto x54 = V{0.25}*pi[2];
auto x55 = u[0] + u[2];
auto x56 = V{4.5}*(x55*x55);
auto x57 = V{0.0416666666666667}*pi[3];
auto x58 = -V{0.0833333333333333}*pi[5];
auto x59 = x48 + x57 + x58;
auto x60 = -x38;
auto x61 = -u[2];
auto x62 = x61 + u[0];
auto x63 = -x62;
auto x64 = x54 + x59;
auto x65 = V{0.25}*pi[4];
auto x66 = u[1] + u[2];
auto x67 = V{4.5}*(x66*x66);
auto x68 = x26 + x33;
auto x69 = V{0.0416666666666667}*pi[0];
auto x70 = x47 + x58 + x69 + V{0.0277777777777778};
auto x71 = x61 + u[1];
auto x72 = -x71;
auto x73 = x65 + x70;
auto x74 = -x22;
auto x75 = V{1} - x24;
auto x76 = x74 + x75;
auto x77 = x28 + x76;
auto x78 = -x20;
auto x79 = x33 + x78;
auto x80 = x38 + x78;
auto x81 = V{0.0833333333333333}*pi[3];
auto x82 = V{0.0833333333333333}*pi[0] + V{-0.0277777777777778};
auto x83 = -x28;
auto x84 = V{0.0833333333333333}*pi[5];
auto x85 = x26 + x38;
cell[0] = -V{0.333333333333333}*rho*x26 - V{0.5}*pi[0] - V{0.5}*pi[3] - V{0.5}*pi[5] + V{-0.333333333333333};
cell[1] = -x27*(x25 + x28 - x29) - x32;
cell[2] = -x27*(x24 + x33 - x34 + x35) - x37;
cell[3] = -x27*(x22 + x35 + x38 - x39) - x40;
cell[4] = x41 - x42*(x33 - x44 + x45) - x49;
cell[5] = -x42*(x45 + x50 - V{4.5}*x52*x52) - x53;
cell[6] = -x42*(x38 + x45 - x56) + x54 - x59;
cell[7] = -x42*(x45 + x60 - V{4.5}*x63*x63) - x64;
cell[8] = -x42*(x38 - x67 + x68) + x65 - x70;
cell[9] = -x42*(x60 + x68 - V{4.5}*x72*x72) - x73;
cell[10] = V{0.0555555555555556}*rho*(x29 + x77) - x32;
cell[11] = V{0.0555555555555556}*rho*(x34 + x75 + x79) - x37;
cell[12] = V{0.0555555555555556}*rho*(x39 + x74 + x80 + V{1}) - x40;
cell[13] = x41 + x42*(x44 + x77 + x79) - x46 + x81 + x82;
cell[14] = -x42*(x68 + x83 - V{4.5}*x51*x51) - x53;
cell[15] = x42*(x56 + x77 + x80) + x54 - x57 + x82 + x84;
cell[16] = -x42*(x83 + x85 - V{4.5}*x62*x62) - x64;
cell[17] = x42*(x38 + x67 + x76 + x79) + x65 - x69 + x81 + x84 + V{-0.0277777777777778};
cell[18] = -x42*(x50 + x85 - V{4.5}*x71*x71) - x73;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[15] - cell[6];
auto x1 = cell[13] - cell[4];
auto x2 = cell[10] + cell[14] + cell[16];
auto x3 = x0 + x1 + x2 - cell[1] - cell[5] - cell[7];
auto x4 = cell[17] - cell[8];
auto x5 = cell[12] + cell[7] + cell[9];
auto x6 = x0 + x4 + x5 - cell[16] - cell[18] - cell[3];
auto x7 = cell[11] + cell[18] + cell[5];
auto x8 = x2 + x5 + x7 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8];
auto x9 = V{1} / (x8 + V{1});
auto x10 = x9*(x8 + V{1});
auto x11 = V{0.5}*x10;
auto x12 = x3*x9;
auto x13 = V{1}*x12;
auto x14 = x11*(x3*force[2] + x6*force[0]) - x13*x6 + V{1}*cell[15] - V{1}*cell[16] + V{1}*cell[6] - V{1}*cell[7];
auto x15 = x1 + x4 + x7 - cell[14] - cell[2] - cell[9];
auto x16 = x15*force[0] + x3*force[1];
auto x17 = V{1}*x10;
auto x18 = V{2}*x15;
auto x19 = x15*force[2] + x6*force[1];
auto x20 = x6*x9;
auto x21 = -V{0.333333333333333}*cell[0];
auto x22 = x21 - V{0.333333333333333}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x23 = -V{0.333333333333333}*cell[11] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x24 = x10*x3*force[0] + x22 + x23 - x9*x3*x3 + V{0.666666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x25 = -V{0.333333333333333}*cell[10] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x26 = x10*x15*force[1] + x22 + x25 - x9*x15*x15 + V{0.666666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x27 = x10*x6*force[2] + x21 + x23 + x25 - x9*x6*x6 + V{0.666666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return (x11*x16 - x13*x15 + V{1}*cell[13] - V{1}*cell[14] + V{1}*cell[4] - V{1}*cell[5])*(-x12*x18 + x16*x17 + V{2}*cell[13] - V{2}*cell[14] + V{2}*cell[4] - V{2}*cell[5]) + (x11*x19 - V{1}*x15*x20 + V{1}*cell[17] - V{1}*cell[18] + V{1}*cell[8] - V{1}*cell[9])*(x17*x19 - x18*x20 + V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[8] - V{2}*cell[9]) + 2*(x14*x14) + x24*x24 + x26*x26 + x27*x27;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = -cell[4];
auto x1 = -cell[8];
auto x2 = x1 + cell[18];
auto x3 = cell[11] + cell[17] + cell[5];
auto x4 = x0 + x2 + x3 + cell[13] - cell[14] - cell[2] - cell[9];
auto x5 = cell[10] + cell[13] + cell[15];
auto x6 = cell[12] + cell[7] + cell[9];
auto x7 = V{1} / (x3 + x5 + x6 + cell[0] + cell[14] + cell[16] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1});
auto x8 = x0 + cell[14];
auto x9 = -cell[6];
auto x10 = x9 + cell[16];
auto x11 = x10 + x5 + x8 - cell[1] - cell[5] - cell[7];
auto x12 = x11*x7;
auto x13 = x12*x4 + x8 - cell[13] + cell[5];
auto x14 = x1 + x6 + x9 + cell[15] - cell[16] + cell[17] - cell[18] - cell[3];
auto x15 = x10 + x12*x14 - cell[15] + cell[7];
auto x16 = x14*x4*x7 + x2 - cell[17] + cell[9];
auto x17 = V{1}*x7;
auto x18 = V{0.333333333333333}*cell[0];
auto x19 = x18 + V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9];
auto x20 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x21 = x17*(x14*x14) + x19 + x20 - V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5];
auto x22 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x23 = x17*(x4*x4) + x19 + x22 - V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7];
auto x24 = x17*(x11*x11) + x18 + x20 + x22 - V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9];
return V{2}*(x13*x13) + V{2}*(x15*x15) + V{2}*(x16*x16) + x21*x21 + x23*x23 + x24*x24;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x19 = force[0]*u[0];
auto x20 = force[1]*u[1];
auto x21 = force[2]*u[2];
auto x22 = rho*(V{0.5}*omega + V{-1});
auto x23 = V{6}*u[0];
auto x24 = x23 + V{-3};
auto x25 = V{0.0555555555555556}*force[0];
auto x26 = V{0.166666666666667}*x20;
auto x27 = V{0.166666666666667}*x21;
auto x28 = x26 + x27;
auto x29 = V{6}*u[1];
auto x30 = x29 + V{-3};
auto x31 = V{0.0555555555555556}*force[1];
auto x32 = V{0.166666666666667}*x19;
auto x33 = x27 + x32;
auto x34 = V{6}*u[2];
auto x35 = x34 + V{-3};
auto x36 = V{0.0555555555555556}*force[2];
auto x37 = x26 + x32;
auto x38 = V{9}*u[1];
auto x39 = V{0.0277777777777778}*force[0];
auto x40 = V{9}*u[0];
auto x41 = V{0.0277777777777778}*force[1];
auto x42 = V{0.0833333333333333}*x21;
auto x43 = -x42;
auto x44 = V{3} - x23;
auto x45 = -x40;
auto x46 = x29 + V{3};
auto x47 = V{9}*u[2];
auto x48 = V{0.0277777777777778}*force[2];
auto x49 = V{0.0833333333333333}*x20;
auto x50 = -x49;
auto x51 = x34 + V{3};
auto x52 = V{0.0833333333333333}*x19;
auto x53 = -x52;
auto x54 = V{3} - x29;
auto x55 = -x38;
auto x56 = x23 + V{3};
auto x57 = V{3} - x34;
auto x58 = -x47;
cell[0] = V{1}*x22*(x19 + x20 + x21) + cell[0];
cell[1] = x22*(-x24*x25 + x28) + cell[1];
cell[2] = x22*(-x30*x31 + x33) + cell[2];
cell[3] = x22*(-x35*x36 + x37) + cell[3];
cell[4] = -x22*(x39*(x24 + x38) + x41*(x30 + x40) + x43) + cell[4];
cell[5] = -x22*(-x39*(x38 + x44) - x42 + V{0.0277777777777778}*(x45 + x46)*force[1]) + cell[5];
cell[6] = -x22*(x39*(x24 + x47) + x48*(x35 + x40) + x50) + cell[6];
cell[7] = -x22*(-x39*(x44 + x47) - x49 + V{0.0277777777777778}*(x45 + x51)*force[2]) + cell[7];
cell[8] = -x22*(x41*(x30 + x47) + x48*(x35 + x38) + x53) + cell[8];
cell[9] = -x22*(-x41*(x47 + x54) - x52 + V{0.0277777777777778}*(x51 + x55)*force[2]) + cell[9];
cell[10] = x22*(-x25*x56 + x28) + cell[10];
cell[11] = x22*(-x31*x46 + x33) + cell[11];
cell[12] = x22*(-x36*x51 + x37) + cell[12];
cell[13] = -x22*(x39*(x38 + x56) + x41*(x40 + x46) + x43) + cell[13];
cell[14] = -x22*(-x41*(x40 + x54) - x42 + V{0.0277777777777778}*(x55 + x56)*force[0]) + cell[14];
cell[15] = -x22*(x39*(x47 + x56) + x48*(x40 + x51) + x50) + cell[15];
cell[16] = -x22*(-x48*(x40 + x57) - x49 + V{0.0277777777777778}*(x56 + x58)*force[0]) + cell[16];
cell[17] = -x22*(x41*(x46 + x47) + x48*(x38 + x51) + x53) + cell[17];
cell[18] = -x22*(-x48*(x38 + x57) - x52 + V{0.0277777777777778}*(x46 + x58)*force[1]) + cell[18];

}

template <typename CELL, typename RHO, typename NABLARHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addLiangForce(CELL& cell, RHO& rho, NABLARHO& nablarho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x19 = rho*force[0];
auto x20 = nablarho[0]*u[0];
auto x21 = V{0.5}*omega + V{-1};
auto x22 = V{0.166666666666667}*x21;
auto x23 = rho*force[1];
auto x24 = nablarho[1]*u[1];
auto x25 = rho*force[2];
auto x26 = nablarho[2]*u[2];
auto x27 = u[0] + u[1];
auto x28 = x27*nablarho[0];
auto x29 = x27*nablarho[1];
auto x30 = x19 + x23;
auto x31 = V{0.0833333333333333}*x21;
auto x32 = u[0] - u[1];
auto x33 = x32*nablarho[1];
auto x34 = x19 - x23;
auto x35 = u[0] + u[2];
auto x36 = x35*nablarho[0];
auto x37 = x35*nablarho[2];
auto x38 = x19 + x25;
auto x39 = -u[2];
auto x40 = x39 + u[0];
auto x41 = x40*nablarho[2];
auto x42 = -x25;
auto x43 = x19 + x42;
auto x44 = u[1] + u[2];
auto x45 = x44*nablarho[1];
auto x46 = x44*nablarho[2];
auto x47 = x23 + x25;
auto x48 = x39 + u[1];
auto x49 = x48*nablarho[2];
auto x50 = x23 + x42;
cell[0] = cell[0];
cell[1] = x22*(x19 - x20) + cell[1];
cell[2] = x22*(x23 - x24) + cell[2];
cell[3] = x22*(x25 - x26) + cell[3];
cell[4] = x31*(-x28 - x29 + x30) + cell[4];
cell[5] = -x31*(x32*nablarho[0] - x33 - x34) + cell[5];
cell[6] = x31*(-x36 - x37 + x38) + cell[6];
cell[7] = -x31*(x40*nablarho[0] - x41 - x43) + cell[7];
cell[8] = x31*(-x45 - x46 + x47) + cell[8];
cell[9] = -x31*(x48*nablarho[1] - x49 - x50) + cell[9];
cell[10] = -x22*(x19 + x20) + cell[10];
cell[11] = -x22*(x23 + x24) + cell[11];
cell[12] = -x22*(x25 + x26) + cell[12];
cell[13] = -x31*(x28 + x29 + x30) + cell[13];
cell[14] = -x31*(x32*nablarho[0] - x33 + x34) + cell[14];
cell[15] = -x31*(x36 + x37 + x38) + cell[15];
cell[16] = -x31*(x40*nablarho[0] - x41 + x43) + cell[16];
cell[17] = -x31*(x45 + x46 + x47) + cell[17];
cell[18] = -x31*(x48*nablarho[1] - x49 + x50) + cell[18];

}

template <typename CELL, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addAllenCahnForce(CELL& cell, OMEGA& omega, FORCE& force) any_platform
{
auto x19 = V{0.5}*omega + V{-1};
auto x20 = V{0.166666666666667}*x19;
auto x21 = x20*force[0];
auto x22 = x20*force[1];
auto x23 = x20*force[2];
auto x24 = V{0.0833333333333333}*x19;
auto x25 = x24*(force[0] + force[1]);
auto x26 = x24*(force[0] - force[1]);
auto x27 = x24*(force[0] + force[2]);
auto x28 = -force[2];
auto x29 = x24*(x28 + force[0]);
auto x30 = x24*(force[1] + force[2]);
auto x31 = x24*(x28 + force[1]);
cell[0] = cell[0];
cell[1] = x21 + cell[1];
cell[2] = x22 + cell[2];
cell[3] = x23 + cell[3];
cell[4] = x25 + cell[4];
cell[5] = x26 + cell[5];
cell[6] = x27 + cell[6];
cell[7] = x29 + cell[7];
cell[8] = x30 + cell[8];
cell[9] = x31 + cell[9];
cell[10] = -x21 + cell[10];
cell[11] = -x22 + cell[11];
cell[12] = -x23 + cell[12];
cell[13] = -x25 + cell[13];
cell[14] = -x26 + cell[14];
cell[15] = -x27 + cell[15];
cell[16] = -x29 + cell[16];
cell[17] = -x30 + cell[17];
cell[18] = -x31 + cell[18];

}

template <typename CELL, typename OMEGA, typename SOURCE, typename V=typename CELL::value_t>
static void addAllenCahnSource(CELL& cell, OMEGA& omega, SOURCE& source) any_platform
{
auto x19 = V{0.0555555555555556}*source;
auto x20 = V{0.0277777777777778}*source;
cell[0] = V{0.333333333333333}*source + cell[0];
cell[1] = x19 + cell[1];
cell[2] = x19 + cell[2];
cell[3] = x19 + cell[3];
cell[4] = x20 + cell[4];
cell[5] = x20 + cell[5];
cell[6] = x20 + cell[6];
cell[7] = x20 + cell[7];
cell[8] = x20 + cell[8];
cell[9] = x20 + cell[9];
cell[10] = x19 + cell[10];
cell[11] = x19 + cell[11];
cell[12] = x19 + cell[12];
cell[13] = x20 + cell[13];
cell[14] = x20 + cell[14];
cell[15] = x20 + cell[15];
cell[16] = x20 + cell[16];
cell[17] = x20 + cell[17];
cell[18] = x20 + cell[18];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q27<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[21] + cell[22] + cell[23] + cell[24] + cell[25] + cell[26] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = -cell[23];
auto x1 = x0 + cell[10] + cell[11] - cell[17] - cell[24] + cell[4];
auto x2 = cell[12] - cell[19] - cell[25] + cell[6];
auto x3 = -cell[13] - cell[21] + cell[26] + cell[8];
j[0] = -V{1}*x1 - V{1}*x2 - V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[18] - V{1}*cell[1] + V{1}*cell[20] + V{1}*cell[26] - V{1}*cell[5] - V{1}*cell[7];
j[1] = -V{1}*x1 - V{1}*x3 + V{1}*cell[12] + V{1}*cell[15] - V{1}*cell[18] + V{1}*cell[22] - V{1}*cell[25] - V{1}*cell[2] + V{1}*cell[5] - V{1}*cell[9];
j[2] = -V{1}*x0 - V{1}*x2 - V{1}*x3 - V{1}*cell[10] + V{1}*cell[11] + V{1}*cell[16] - V{1}*cell[20] - V{1}*cell[22] - V{1}*cell[24] - V{1}*cell[3] + V{1}*cell[7] + V{1}*cell[9];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[13] + cell[1] + cell[5] + cell[7];
auto x1 = cell[18] + cell[25] + cell[2] + cell[9];
auto x2 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x3 = x0 + x1 + x2 + cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = -cell[23];
auto x5 = x4 + cell[10] + cell[11] - cell[17] - cell[24] + cell[4];
auto x6 = cell[12] - cell[19] - cell[25] + cell[6];
auto x7 = V{1}/x3;
auto x8 = -cell[13] - cell[21] + cell[26] + cell[8];
rho = x3;
u[0] = -x7*(x0 + x5 + x6 - cell[14] - cell[18] - cell[20] - cell[26]);
u[1] = -x7*(x1 + x5 + x8 - cell[12] - cell[15] - cell[22] - cell[5]);
u[2] = -x7*(x2 + x4 + x6 + x8 - cell[11] - cell[16] - cell[7] - cell[9]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[13] + cell[1] + cell[5] + cell[7];
auto x1 = cell[18] + cell[25] + cell[2] + cell[9];
auto x2 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x3 = -cell[23];
auto x4 = x3 + cell[10] + cell[11] - cell[17] - cell[24] + cell[4];
auto x5 = cell[12] - cell[19] - cell[25] + cell[6];
auto x6 = -cell[13] - cell[21] + cell[26] + cell[8];
rho = x0 + x1 + x2 + cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + V{1};
j[0] = -V{1}*x0 - V{1}*x4 - V{1}*x5 + V{1}*cell[14] + V{1}*cell[18] + V{1}*cell[20] + V{1}*cell[26];
j[1] = -V{1}*x1 - V{1}*x4 - V{1}*x6 + V{1}*cell[12] + V{1}*cell[15] + V{1}*cell[22] + V{1}*cell[5];
j[2] = -V{1}*x2 - V{1}*x3 - V{1}*x5 - V{1}*x6 + V{1}*cell[11] + V{1}*cell[16] + V{1}*cell[7] + V{1}*cell[9];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{1}*cell[13];
auto x1 = V{1}*cell[26];
auto x2 = V{1}*cell[10];
auto x3 = V{1}*cell[23];
auto x4 = -V{0.333333333333333}*rho;
auto x5 = V{1}*cell[12];
auto x6 = V{1}*cell[25];
auto x7 = x5 + x6;
auto x8 = x7 + V{1}*cell[18] + V{1}*cell[5];
auto x9 = x0 + x1 + x2 + x3 + x4 + x8 + V{1}*cell[17] + V{1}*cell[4] + V{0.333333333333333};
auto x10 = V{1}*cell[11];
auto x11 = V{1}*cell[24];
auto x12 = x10 + x11;
auto x13 = x12 + V{1}*cell[20] + V{1}*cell[7];
auto x14 = x13 + V{1}*cell[19] + V{1}*cell[6];
auto x15 = rho*u[0];
auto x16 = -V{1}*cell[10];
auto x17 = -V{1}*cell[23];
auto x18 = x0 + x1 + x16 + x17;
auto x19 = V{1}*cell[9];
auto x20 = V{1}*cell[22];
auto x21 = x12 + x19 + x20;
auto x22 = V{1}*cell[21] + V{1}*cell[8];
pi[0] = -rho*u[0]*u[0] + x14 + x9 + V{1}*cell[14] + V{1}*cell[1];
pi[1] = x10 + x11 - x15*u[1] - x18 - x8 + V{1}*cell[17] + V{1}*cell[4];
pi[2] = -x13 - x15*u[2] - x18 + x5 + x6 + V{1}*cell[19] + V{1}*cell[6];
pi[3] = -rho*u[1]*u[1] + x21 + x22 + x9 + V{1}*cell[15] + V{1}*cell[2];
pi[4] = -rho*u[1]*u[2] + x0 + x1 - x16 - x17 - x21 - x7 + V{1}*cell[21] + V{1}*cell[8];
pi[5] = -rho*u[2]*u[2] + x0 + x1 + x14 + x19 + x2 + x20 + x22 + x3 + x4 + x7 + V{1}*cell[16] + V{1}*cell[3] + V{0.333333333333333};

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[12] + cell[25];
auto x1 = x0 + cell[18];
auto x2 = cell[11] + cell[24];
auto x3 = x2 + cell[20];
auto x4 = cell[22] + cell[9];
auto x5 = cell[13] + cell[1];
auto x6 = cell[26] + cell[2] + cell[8];
auto x7 = cell[10] + cell[3];
auto x8 = x1 + x3 + x4 + x5 + x6 + x7 + cell[0] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[4] + cell[5] + cell[6] + cell[7] + V{1};
auto x9 = -cell[17];
auto x10 = -cell[24];
auto x11 = x10 + x9 + cell[5];
auto x12 = -cell[19];
auto x13 = -cell[25];
auto x14 = x12 + x13 + cell[7];
auto x15 = -cell[23];
auto x16 = x15 - cell[26];
auto x17 = cell[10] + cell[11] + cell[4];
auto x18 = cell[12] + cell[6];
auto x19 = x11 + x14 + x16 + x17 + x18 + x5 - cell[14] - cell[18] - cell[20];
auto x20 = V{1} / (x8);
auto x21 = V{1}*x20;
auto x22 = -cell[12];
auto x23 = -cell[13] - cell[21];
auto x24 = x15 + x23;
auto x25 = x10 + x17 + x22 + x24 + x6 + x9 - cell[15] + cell[18] - cell[22] + cell[25] - cell[5] + cell[9];
auto x26 = -cell[11];
auto x27 = x12 + x13 + x18 + x24 + x26 + x7 - cell[16] + cell[20] + cell[22] + cell[24] + cell[26] - cell[7] + cell[8] - cell[9];
auto x28 = V{0.666666666666667}*cell[10];
auto x29 = V{0.666666666666667}*cell[11];
auto x30 = V{0.666666666666667}*cell[12];
auto x31 = V{0.666666666666667}*cell[13];
auto x32 = V{0.666666666666667}*cell[23];
auto x33 = V{0.666666666666667}*cell[24];
auto x44 = V{0.666666666666667}*cell[25];
auto x45 = V{0.666666666666667}*cell[26];
auto x46 = -V{0.333333333333333}*cell[0];
auto x47 = x28 + x29 + x30 + x31 + x32 + x33 + x44 + x45 + x46 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x48 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x49 = x19*x20;
auto x50 = -cell[10];
auto x51 = x15 + x50 + cell[13] + cell[26];
auto x52 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
rho = x8;
u[0] = -x19*x21;
u[1] = -x21*x25;
u[2] = -x21*x27;
pi[0] = -x21*x19*x19 + x47 + x48 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
pi[1] = -V{1}*x1 - V{1}*x11 - V{1}*x25*x49 - V{1}*x26 - V{1}*x51 + V{1}*cell[4];
pi[2] = -V{1}*x14 - V{1}*x22 - V{1}*x27*x49 - V{1}*x3 - V{1}*x51 + V{1}*cell[6];
pi[3] = -x21*x25*x25 + x47 + x52 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
pi[4] = -V{1}*x0 - V{1}*x16 - V{1}*x2 - V{1}*x20*x25*x27 - V{1}*x23 - V{1}*x4 - V{1}*x50 + V{1}*cell[8];
pi[5] = -x21*x27*x27 + x28 + x29 + x30 + x31 + x32 + x33 + x44 + x45 + x46 + x48 + x52 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[8] + cell[9];
auto x1 = cell[15] + cell[21];
auto x2 = cell[11] + cell[26];
auto x3 = x0 + x1 + x2 + cell[0] + cell[10] + cell[12] + cell[13] + cell[14] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[22] + cell[23] + cell[24] + cell[25] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7];
auto x4 = x3 + V{1};
auto x5 = x3 + V{1};
auto x6 = V{1} / ((x5)*(x5));
auto x7 = V{1.5}*x6;
auto x8 = -cell[17];
auto x9 = -cell[19];
auto x10 = x8 + x9 + cell[4] + cell[6];
auto x11 = -cell[24];
auto x12 = cell[10] - cell[23];
auto x13 = x11 + x12 + cell[11];
auto x14 = cell[12] - cell[25];
auto x15 = -cell[14];
auto x16 = -cell[20] + cell[7];
auto x17 = x15 + x16 + cell[1];
auto x18 = -cell[18] + cell[5];
auto x19 = cell[13] - cell[26];
auto x20 = x10 + x13 + x14 + x17 + x18 + x19;
auto x21 = -x20;
auto x22 = x21*x21;
auto x23 = x22*x7;
auto x24 = -cell[13];
auto x25 = -cell[21];
auto x26 = x24 + x25 + cell[26] + cell[8];
auto x27 = -cell[22];
auto x28 = x8 + cell[4];
auto x29 = x27 + x28 + cell[9];
auto x30 = -cell[15];
auto x31 = x30 + cell[18] + cell[2] - cell[5];
auto x32 = -cell[12] + cell[25];
auto x33 = x13 + x26 + x29 + x31 + x32;
auto x61 = -x33;
auto x62 = x61*x61;
auto x63 = x62*x7;
auto x64 = -cell[16];
auto x65 = x14 + x64 + cell[3];
auto x66 = -cell[9];
auto x67 = x9 + cell[6];
auto x68 = x66 + x67 + cell[22];
auto x69 = cell[20] - cell[7];
auto x70 = x12 - cell[11] + cell[24];
auto x71 = x26 + x65 + x68 + x69 + x70;
auto x72 = -x71;
auto x73 = x72*x72;
auto x74 = x7*x73;
auto x75 = x63 + x74 + V{-1};
auto x76 = x23 + x75;
auto x77 = V{1} / (x5);
auto x78 = V{3}*cell[5];
auto x79 = V{3}*cell[7];
auto x80 = V{3}*cell[13];
auto x81 = V{3}*cell[18];
auto x82 = V{3}*cell[20];
auto x83 = V{3}*cell[26];
auto x84 = V{3}*cell[10];
auto x85 = V{3}*cell[11];
auto x86 = -V{3}*cell[23];
auto x87 = V{3}*cell[24];
auto x88 = x84 + x85 + x86 - x87 - V{3}*cell[17] + V{3}*cell[4];
auto x89 = V{3}*cell[12];
auto x90 = V{3}*cell[25];
auto x91 = x89 - x90 - V{3}*cell[19] + V{3}*cell[6];
auto x92 = x78 + x79 + x80 - x81 - x82 - x83 + x88 + x91 - V{3}*cell[14] + V{3}*cell[1];
auto x93 = -x77*x92;
auto x94 = V{3}*x6;
auto x95 = V{3}*cell[9];
auto x96 = V{3}*cell[22];
auto x97 = -x80 + x83 - V{3}*cell[21] + V{3}*cell[8];
auto x98 = -x78 + x81 + x88 - x89 + x90 + x95 - x96 + x97 - V{3}*cell[15] + V{3}*cell[2];
auto x99 = -x77*x98;
auto x100 = x23 + V{-1};
auto x101 = -x79 + x82 + x84 - x85 + x86 + x87 + x91 - x95 + x96 + x97 - V{3}*cell[16] + V{3}*cell[3];
auto x102 = -x101*x77;
auto x103 = V{4.5}*x6;
auto x104 = V{2}*cell[11] - V{2}*cell[24];
auto x105 = V{2}*cell[10] - V{2}*cell[23];
auto x106 = x105 + x25;
auto x107 = x30 + cell[2];
auto x108 = x107 - V{2}*cell[17] + V{2}*cell[4];
auto x109 = x0 + x104 + x106 + x108 + x17 + x27 + x67;
auto x110 = -x109;
auto x111 = x76 + x93;
auto x112 = x111 + x99;
auto x113 = V{2}*cell[25];
auto x114 = V{2}*cell[12];
auto x115 = -x113 + x114;
auto x116 = V{2}*cell[26];
auto x117 = V{2}*cell[13];
auto x118 = -x116 + x117 - cell[8];
auto x119 = V{2}*cell[18];
auto x120 = V{2}*cell[5];
auto x121 = -x119 + x120 - cell[2];
auto x122 = x1 + x115 + x118 + x121 + x17 + x68;
auto x123 = -x99;
auto x124 = x111 + x123;
auto x125 = x15 + cell[1];
auto x126 = x125 + x18;
auto x127 = x64 + cell[3];
auto x128 = -V{2}*cell[19] + V{2}*cell[6];
auto x129 = x127 + x128;
auto x130 = x106 + x115 + x126 + x129 + x28 + x66 + cell[22] + cell[8];
auto x131 = -x130;
auto x132 = -x102;
auto x133 = -cell[3];
auto x134 = x104 + x133 + cell[16];
auto x135 = V{2}*cell[20];
auto x136 = V{2}*cell[7];
auto x137 = -x135 + x136;
auto x138 = x118 + x126 + x134 + x137 + x29 + cell[21];
auto x139 = -V{2}*cell[21] + V{2}*cell[8];
auto x140 = x127 + x139;
auto x141 = x10 + x105 + x116 - x117 + x140 + x31 + x69;
auto x142 = -x141;
auto x143 = x76 + x99;
auto x144 = x102 + x143;
auto x145 = V{2}*cell[22];
auto x146 = V{2}*cell[9];
auto x147 = -x145 + x146;
auto x148 = x113 - x114 + x134 + x147 + x16 + x28 + x31 + cell[19] - cell[6];
auto x149 = x108 + x125;
auto x150 = x11 + x128 + x139 + x149 + x2 + x24 + x65 + V{3}*cell[10] - V{3}*cell[23];
auto x151 = -x150;
auto x152 = x12 + x133 + x137 + x147 + x149 + x19 + x32 + V{3}*cell[11] + cell[16] - V{3}*cell[24];
auto x153 = x121 + x125 + x129 + x145 - x146 + x19 + x70 + V{3}*cell[12] + cell[15] - V{3}*cell[25];
auto x154 = x77*x92;
auto x155 = x107 + x119 - x120 + x135 - x136 + x140 + x32 + x70 - V{3}*cell[13] + cell[14] - cell[1] + V{3}*cell[26];
auto x156 = -x155;
auto x157 = x101*x77;
auto x158 = x33*x33;
auto x159 = x158*x7;
auto x160 = x71*x71;
auto x161 = x160*x7 + V{-1};
auto x162 = x159 + x161;
auto x163 = x77*x98;
auto x164 = x20*x20;
auto x165 = x164*x7;
auto x166 = x163 + x165;
auto x167 = x157 + x162 + x166;
auto x168 = x154 + x162;
auto x169 = x157 + x165;
auto x170 = x166 + x168;
auto x171 = -x93;
auto x172 = -x122;
auto x173 = x168 + x169;
auto x174 = -x138;
auto x175 = x102 + x76;
auto x176 = -x148;
auto x177 = -x152;
auto x178 = -x153;
fEq[0] = -V{0.296296296296296}*x4*x76 + V{-0.296296296296296};
fEq[1] = -V{0.0740740740740741}*x4*(-x22*x94 + x75 + x93) + V{-0.0740740740740741};
fEq[2] = -V{0.0740740740740741}*x4*(x100 - x62*x94 + x74 + x99) + V{-0.0740740740740741};
fEq[3] = -V{0.0740740740740741}*x4*(x100 + x102 + x63 - x73*x94) + V{-0.0740740740740741};
fEq[4] = -V{0.0185185185185185}*(x4*(-x103*x110*x110 + x112) + V{1});
fEq[5] = -V{0.0185185185185185}*(x4*(-x103*x122*x122 + x124) + V{1});
fEq[6] = -V{0.0185185185185185}*(x4*(x102 - x103*x131*x131 + x111) + V{1});
fEq[7] = -V{0.0185185185185185}*(x4*(-x103*x138*x138 + x111 + x132) + V{1});
fEq[8] = -V{0.0185185185185185}*(x4*(-x103*x142*x142 + x144) + V{1});
fEq[9] = -V{0.0185185185185185}*(x4*(-x103*x148*x148 + x132 + x143) + V{1});
fEq[10] = -V{0.00462962962962963}*(x4*(x102 - x103*x151*x151 + x112) + V{1});
fEq[11] = -V{0.00462962962962963}*(x4*(-x103*x152*x152 + x112 + x132) + V{1});
fEq[12] = -V{0.00462962962962963}*(x4*(x102 - x103*x153*x153 + x124) + V{1});
fEq[13] = -V{0.00462962962962963}*(x4*(-x103*x156*x156 - x154 + x167) + V{1});
fEq[14] = -V{0.0740740740740741}*x4*(-x164*x94 + x168) + V{-0.0740740740740741};
fEq[15] = -V{0.0740740740740741}*x4*(-x158*x94 + x161 + x166) + V{-0.0740740740740741};
fEq[16] = -V{0.0740740740740741}*x4*(x159 - x160*x94 + x169 + V{-1}) + V{-0.0740740740740741};
fEq[17] = -V{0.0185185185185185}*(x4*(-x103*x109*x109 + x170) + V{1});
fEq[18] = -V{0.0185185185185185}*(x4*(-x103*x172*x172 + x143 + x171) + V{1});
fEq[19] = -V{0.0185185185185185}*(x4*(-x103*x130*x130 + x173) + V{1});
fEq[20] = -V{0.0185185185185185}*(x4*(-x103*x174*x174 + x171 + x175) + V{1});
fEq[21] = -V{0.0185185185185185}*(x4*(-x103*x141*x141 + x167) + V{1});
fEq[22] = -V{0.0185185185185185}*(x4*(-x103*x176*x176 + x123 + x175) + V{1});
fEq[23] = -V{0.00462962962962963}*(x4*(-x103*x150*x150 + x157 + x170) + V{1});
fEq[24] = -V{0.00462962962962963}*(x4*(-x103*x177*x177 - x157 + x170) + V{1});
fEq[25] = -V{0.00462962962962963}*(x4*(-x103*x178*x178 - x163 + x173) + V{1});
fEq[26] = -V{0.00462962962962963}*(x4*(-x103*x155*x155 + x144 + x171) + V{1});

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{1.5}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = x1 + x6;
auto x8 = V{0.0740740740740741}*rho;
auto x9 = V{3}*u[0];
auto x10 = V{3}*x0;
auto x11 = V{3}*u[1];
auto x12 = V{3}*x2;
auto x13 = x1 + V{-1};
auto x14 = V{3}*u[2];
auto x15 = V{3}*x4;
auto x16 = V{0.0185185185185185}*rho;
auto x17 = u[0] + u[1];
auto x18 = V{4.5}*(x17*x17);
auto x19 = x7 + x9;
auto x20 = x11 + x19;
auto x21 = u[0] - u[1];
auto x22 = -V{4.5}*x21*x21;
auto x23 = -x11;
auto x24 = x19 + x23;
auto x25 = u[0] + u[2];
auto x26 = V{4.5}*(x25*x25);
auto x27 = -x14;
auto x28 = -u[2];
auto x29 = x28 + u[0];
auto x30 = -V{4.5}*x29*x29;
auto x31 = u[1] + u[2];
auto x32 = V{4.5}*(x31*x31);
auto x33 = x11 + x7;
auto x61 = x14 + x33;
auto x62 = x28 + u[1];
auto x63 = -V{4.5}*x62*x62;
auto x64 = V{0.00462962962962963}*rho;
auto x65 = x17 + u[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = x17 + x28;
auto x68 = V{4.5}*(x67*x67);
auto x69 = x21 + u[2];
auto x70 = V{4.5}*(x69*x69);
auto x71 = x31 - u[0];
auto x72 = -V{4.5}*x71*x71;
auto x73 = -x3;
auto x74 = V{1} - x5;
auto x75 = x73 + x74;
auto x76 = x75 + x9;
auto x77 = -x1;
auto x78 = x11 + x77;
auto x79 = x14 + x77;
auto x80 = x76 + x78;
auto x81 = -x9;
auto x82 = x76 + x79;
auto x83 = x14 + x7;
fNeq[0] = V{0.296296296296296}*rho*x7 + cell[0] + V{0.296296296296296};
fNeq[1] = x8*(-x10 + x6 + x9) + cell[1] + V{0.0740740740740741};
fNeq[2] = x8*(x11 - x12 + x13 + x5) + cell[2] + V{0.0740740740740741};
fNeq[3] = x8*(x13 + x14 - x15 + x3) + cell[3] + V{0.0740740740740741};
fNeq[4] = x16*(-x18 + x20) + cell[4] + V{0.0185185185185185};
fNeq[5] = x16*(x22 + x24) + cell[5] + V{0.0185185185185185};
fNeq[6] = x16*(x14 + x19 - x26) + cell[6] + V{0.0185185185185185};
fNeq[7] = x16*(x19 + x27 + x30) + cell[7] + V{0.0185185185185185};
fNeq[8] = x16*(-x32 + x61) + cell[8] + V{0.0185185185185185};
fNeq[9] = x16*(x27 + x33 + x63) + cell[9] + V{0.0185185185185185};
fNeq[10] = x64*(x14 + x20 - x66) + cell[10] + V{0.00462962962962963};
fNeq[11] = x64*(x20 + x27 - x68) + cell[11] + V{0.00462962962962963};
fNeq[12] = x64*(x14 + x24 - x70) + cell[12] + V{0.00462962962962963};
fNeq[13] = -x64*(-x24 - x27 - x72) + cell[13] + V{0.00462962962962963};
fNeq[14] = -x8*(x10 + x76) + cell[14] + V{0.0740740740740741};
fNeq[15] = -x8*(x12 + x74 + x78) + cell[15] + V{0.0740740740740741};
fNeq[16] = -x8*(x15 + x73 + x79 + V{1}) + cell[16] + V{0.0740740740740741};
fNeq[17] = -x16*(x18 + x80) + cell[17] + V{0.0185185185185185};
fNeq[18] = x16*(x22 + x33 + x81) + cell[18] + V{0.0185185185185185};
fNeq[19] = -x16*(x26 + x82) + cell[19] + V{0.0185185185185185};
fNeq[20] = x16*(x30 + x81 + x83) + cell[20] + V{0.0185185185185185};
fNeq[21] = -x16*(x14 + x32 + x75 + x78) + cell[21] + V{0.0185185185185185};
fNeq[22] = x16*(x23 + x63 + x83) + cell[22] + V{0.0185185185185185};
fNeq[23] = -x64*(x14 + x66 + x80) + cell[23] + V{0.00462962962962963};
fNeq[24] = -x64*(x27 + x68 + x80) + cell[24] + V{0.00462962962962963};
fNeq[25] = -x64*(x23 + x70 + x82) + cell[25] + V{0.00462962962962963};
fNeq[26] = x64*(x61 + x72 + x81) + cell[26] + V{0.00462962962962963};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[15] + cell[21];
auto x1 = cell[12] + cell[26];
auto x2 = x0 + x1 + cell[0] + cell[10] + cell[11] + cell[13] + cell[14] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[22] + cell[23] + cell[24] + cell[25] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = -cell[19];
auto x6 = -cell[25];
auto x7 = x5 + x6 + cell[12] + cell[6];
auto x8 = -cell[14];
auto x9 = cell[11] - cell[24];
auto x10 = x8 + x9 + cell[1];
auto x11 = -cell[18] + cell[5];
auto x12 = -cell[17];
auto x13 = x12 + cell[4];
auto x14 = -cell[20] + cell[7];
auto x15 = x13 + x14;
auto x16 = cell[10] - cell[23];
auto x17 = x16 + cell[13] - cell[26];
auto x18 = x10 + x11 + x15 + x17 + x7;
auto x19 = x18*x18;
auto x20 = x19*x4;
auto x21 = -cell[21];
auto x22 = x12 + x21 + cell[4] + cell[8];
auto x23 = -cell[13];
auto x24 = x16 + x23 + cell[26];
auto x25 = -cell[15];
auto x26 = -cell[22] + cell[9];
auto x27 = x25 + x26 + cell[2];
auto x28 = cell[18] - cell[5];
auto x29 = -cell[12] + cell[25];
auto x30 = x22 + x24 + x27 + x28 + x29 + x9;
auto x31 = x30*x30;
auto x32 = x31*x4;
auto x33 = x21 + cell[8];
auto x61 = -cell[16];
auto x62 = cell[22] - cell[9];
auto x63 = x61 + x62 + cell[3];
auto x64 = -cell[11] + cell[24];
auto x65 = cell[20] - cell[7];
auto x66 = x24 + x33 + x63 + x64 + x65 + x7;
auto x67 = x66*x66;
auto x68 = x4*x67;
auto x69 = x32 + x68 + V{-1};
auto x70 = x20 + x69;
auto x71 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x72 = V{3}*x3;
auto x73 = x19*x72;
auto x74 = V{1} / (x2);
auto x75 = V{3}*cell[5];
auto x76 = V{3}*cell[7];
auto x77 = V{3}*cell[13];
auto x78 = V{3}*cell[18];
auto x79 = V{3}*cell[20];
auto x80 = V{3}*cell[26];
auto x81 = V{3}*cell[10];
auto x82 = V{3}*cell[11];
auto x83 = -V{3}*cell[23];
auto x84 = V{3}*cell[24];
auto x85 = x81 + x82 + x83 - x84 - V{3}*cell[17] + V{3}*cell[4];
auto x86 = V{3}*cell[12];
auto x87 = V{3}*cell[25];
auto x88 = x86 - x87 - V{3}*cell[19] + V{3}*cell[6];
auto x89 = x74*(x75 + x76 + x77 - x78 - x79 - x80 + x85 + x88 - V{3}*cell[14] + V{3}*cell[1]);
auto x90 = -V{1.5}*x3*x31;
auto x91 = -V{1.5}*x3*x67 + V{1};
auto x92 = x90 + x91;
auto x93 = x89 + x92;
auto x94 = x31*x72;
auto x95 = V{3}*cell[9];
auto x96 = V{3}*cell[22];
auto x97 = -x77 + x80 - V{3}*cell[21] + V{3}*cell[8];
auto x98 = x74*(-x75 + x78 + x85 - x86 + x87 + x95 - x96 + x97 - V{3}*cell[15] + V{3}*cell[2]);
auto x99 = -V{1.5}*x19*x3;
auto x100 = x98 + x99;
auto x101 = x67*x72;
auto x102 = x74*(-x76 + x79 + x81 - x82 + x83 + x84 + x88 - x95 + x96 + x97 - V{3}*cell[16] + V{3}*cell[3]);
auto x103 = x102 + x99;
auto x104 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x105 = V{4.5}*x3;
auto x106 = V{2}*cell[11] - V{2}*cell[24];
auto x107 = -V{2}*cell[17] + V{2}*cell[4];
auto x108 = -V{2}*cell[23];
auto x109 = V{2}*cell[10];
auto x110 = x8 + cell[1];
auto x111 = x108 + x109 + x110;
auto x112 = x5 + cell[6];
auto x113 = x112 + x14;
auto x114 = x106 + x107 + x111 + x113 + x27 + x33;
auto x115 = -x114;
auto x116 = x100 + x93;
auto x117 = V{2}*cell[25];
auto x118 = V{2}*cell[12];
auto x119 = -x117 + x118;
auto x120 = V{2}*cell[26];
auto x121 = V{2}*cell[13];
auto x122 = x110 - x120 + x121 - cell[8];
auto x123 = V{2}*cell[18];
auto x124 = V{2}*cell[5];
auto x125 = -x123 + x124 - cell[2];
auto x126 = -x0 - x113 - x119 - x122 - x125 - x62;
auto x127 = -x105*x126*x126;
auto x128 = -x89;
auto x129 = x128 + x70;
auto x130 = x129 + x98;
auto x131 = -V{2}*cell[19] + V{2}*cell[6];
auto x132 = x11 + x111 + x119 + x131 + x22 + x63;
auto x133 = -x132;
auto x134 = x103 + x93;
auto x135 = -cell[3];
auto x136 = x106 + x135 + cell[16];
auto x137 = V{2}*cell[20];
auto x138 = V{2}*cell[7];
auto x139 = -x137 + x138;
auto x140 = -x11 - x122 - x13 - x136 - x139 - x26 - cell[21];
auto x141 = -x105*x140*x140;
auto x142 = x102 + x129;
auto x143 = x25 + cell[2];
auto x144 = x143 + x28;
auto x145 = x61 + cell[3];
auto x146 = x145 - V{2}*cell[21] + V{2}*cell[8];
auto x147 = x108 + x109 + x112 + x120 - x121 + x13 + x144 + x146 + x65;
auto x148 = -x147;
auto x149 = x100 + x102 + x92;
auto x150 = V{2}*cell[22];
auto x151 = V{2}*cell[9];
auto x152 = -x150 + x151;
auto x153 = -x117 + x118 - x136 - x144 - x15 - x152 - cell[19] + cell[6];
auto x154 = -x105*x153*x153;
auto x155 = -x98;
auto x156 = x155 + x70;
auto x157 = x102 + x156;
auto x158 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x159 = x107 + x143;
auto x160 = x1 + x10 + x131 + x146 + x159 + x23 + x6 + V{3}*cell[10] - V{3}*cell[23];
auto x161 = -x160;
auto x162 = x110 + x17;
auto x163 = -x135 - x139 - x152 - x159 - x162 - x29 - V{3}*cell[11] - cell[16] + V{3}*cell[24];
auto x164 = -x105*x163*x163;
auto x165 = -x125 - x131 - x145 - x150 + x151 - x162 - x64 - V{3}*cell[12] - cell[15] + V{3}*cell[25];
auto x166 = x105*(x165*x165);
auto x167 = -x123 + x124 - x137 + x138 - x143 - x146 - x16 - x29 - x64 + V{3}*cell[13] - cell[14] + cell[1] - V{3}*cell[26];
auto x168 = x105*(x167*x167);
auto x169 = x20 + V{-1};
auto x170 = x70 + x89;
auto x171 = x170 + x98;
auto x172 = -x102;
auto x173 = x70 + x98;
fNeq[0] = x70*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + cell[0] + V{0.296296296296296};
fNeq[1] = x71*(-x73 - x93) + cell[1] + V{0.0740740740740741};
fNeq[2] = x71*(-x100 - x91 - x94) + cell[2] + V{0.0740740740740741};
fNeq[3] = x71*(-x101 - x103 - x90 + V{-1}) + cell[3] + V{0.0740740740740741};
fNeq[4] = -x104*(x105*(x115*x115) + x116) + cell[4] + V{0.0185185185185185};
fNeq[5] = x104*(x127 + x130) + cell[5] + V{0.0185185185185185};
fNeq[6] = -x104*(x105*(x133*x133) + x134) + cell[6] + V{0.0185185185185185};
fNeq[7] = x104*(x141 + x142) + cell[7] + V{0.0185185185185185};
fNeq[8] = -x104*(x105*(x148*x148) + x149) + cell[8] + V{0.0185185185185185};
fNeq[9] = x104*(x154 + x157) + cell[9] + V{0.0185185185185185};
fNeq[10] = -x158*(x102 + x105*(x161*x161) + x116) + cell[10] + V{0.00462962962962963};
fNeq[11] = x158*(x142 + x155 + x164) + cell[11] + V{0.00462962962962963};
fNeq[12] = x158*(-x134 - x155 - x166) + cell[12] + V{0.00462962962962963};
fNeq[13] = x158*(x102 + x130 - x168) + cell[13] + V{0.00462962962962963};
fNeq[14] = x71*(x69 - x73 + x89) + cell[14] + V{0.0740740740740741};
fNeq[15] = x71*(x169 + x68 - x94 + x98) + cell[15] + V{0.0740740740740741};
fNeq[16] = x71*(-x101 + x102 + x169 + x32) + cell[16] + V{0.0740740740740741};
fNeq[17] = x104*(-x105*x114*x114 + x171) + cell[17] + V{0.0185185185185185};
fNeq[18] = x104*(x127 + x156 + x89) + cell[18] + V{0.0185185185185185};
fNeq[19] = x104*(x102 - x105*x132*x132 + x170) + cell[19] + V{0.0185185185185185};
fNeq[20] = x104*(x141 + x170 + x172) + cell[20] + V{0.0185185185185185};
fNeq[21] = x104*(x102 - x105*x147*x147 + x173) + cell[21] + V{0.0185185185185185};
fNeq[22] = x104*(x154 + x172 + x173) + cell[22] + V{0.0185185185185185};
fNeq[23] = x158*(x102 - x105*x160*x160 + x171) + cell[23] + V{0.00462962962962963};
fNeq[24] = x158*(x164 + x171 + x172) + cell[24] + V{0.00462962962962963};
fNeq[25] = x158*(x157 - x166 + x89) + cell[25] + V{0.00462962962962963};
fNeq[26] = x158*(-x128 - x149 - x168) + cell[26] + V{0.00462962962962963};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = u[0]*u[0];
auto x29 = V{1.5}*x28;
auto x30 = u[1]*u[1];
auto x31 = V{1.5}*x30;
auto x32 = u[2]*u[2];
auto x33 = V{1.5}*x32;
auto x34 = x31 + x33 + V{-1};
auto x35 = x29 + x34;
auto x36 = V{0.0740740740740741}*omega;
auto x37 = V{3}*u[0];
auto x38 = V{3}*x28;
auto x39 = V{3}*u[1];
auto x40 = V{3}*x30;
auto x41 = x29 + V{-1};
auto x42 = V{3}*u[2];
auto x43 = V{3}*x32;
auto x44 = V{0.0185185185185185}*omega;
auto x45 = u[0] + u[1];
auto x46 = V{4.5}*(x45*x45);
auto x47 = x35 + x37;
auto x48 = x39 + x47;
auto x49 = u[0] - u[1];
auto x50 = -x49;
auto x51 = -x39;
auto x52 = x47 + x51;
auto x53 = u[0] + u[2];
auto x54 = V{4.5}*(x53*x53);
auto x55 = -x42;
auto x56 = -u[2];
auto x57 = x56 + u[0];
auto x58 = -x57;
auto x59 = u[1] + u[2];
auto x60 = V{4.5}*(x59*x59);
auto x61 = x35 + x39;
auto x62 = x42 + x61;
auto x63 = x56 + u[1];
auto x64 = -x63;
auto x65 = V{0.00462962962962963}*omega;
auto x66 = x45 + u[2];
auto x67 = V{4.5}*(x66*x66);
auto x68 = x45 + x56;
auto x69 = -x68;
auto x70 = x49 + u[2];
auto x71 = -x70;
auto x72 = x59 - u[0];
auto x73 = -x31;
auto x74 = V{1} - x33;
auto x75 = x73 + x74;
auto x76 = x37 + x75;
auto x77 = -x29;
auto x78 = x39 + x77;
auto x79 = x42 + x77;
auto x80 = x76 + x78;
auto x81 = -x37;
auto x82 = x76 + x79;
auto x83 = x35 + x42;
auto x84 = -x72;
cell[0] = -V{0.296296296296296}*omega*(rho*x35 + V{1}) - x27*cell[0];
cell[1] = -x27*cell[1] - x36*(rho*(x34 + x37 - x38) + V{1});
cell[2] = -x27*cell[2] - x36*(rho*(x33 + x39 - x40 + x41) + V{1});
cell[3] = -x27*cell[3] - x36*(rho*(x31 + x41 + x42 - x43) + V{1});
cell[4] = -x27*cell[4] - x44*(rho*(-x46 + x48) + V{1});
cell[5] = -x27*cell[5] - x44*(rho*(x52 - V{4.5}*x50*x50) + V{1});
cell[6] = -x27*cell[6] - x44*(rho*(x42 + x47 - x54) + V{1});
cell[7] = -x27*cell[7] - x44*(rho*(x47 + x55 - V{4.5}*x58*x58) + V{1});
cell[8] = -x27*cell[8] - x44*(rho*(-x60 + x62) + V{1});
cell[9] = -x27*cell[9] - x44*(rho*(x55 + x61 - V{4.5}*x64*x64) + V{1});
cell[10] = -x27*cell[10] - x65*(rho*(x42 + x48 - x67) + V{1});
cell[11] = -x27*cell[11] - x65*(rho*(x48 + x55 - V{4.5}*x69*x69) + V{1});
cell[12] = -x27*cell[12] - x65*(rho*(x42 + x52 - V{4.5}*x71*x71) + V{1});
cell[13] = -x27*cell[13] + x65*(-rho*(x52 + x55 - V{4.5}*x72*x72) + V{-1});
cell[14] = -x27*cell[14] + x36*(rho*(x38 + x76) + V{-1});
cell[15] = -x27*cell[15] + x36*(rho*(x40 + x74 + x78) + V{-1});
cell[16] = -x27*cell[16] + x36*(rho*(x43 + x73 + x79 + V{1}) + V{-1});
cell[17] = -x27*cell[17] + x44*(rho*(x46 + x80) + V{-1});
cell[18] = -x27*cell[18] - x44*(rho*(x61 + x81 - V{4.5}*x49*x49) + V{1});
cell[19] = -x27*cell[19] + x44*(rho*(x54 + x82) + V{-1});
cell[20] = -x27*cell[20] - x44*(rho*(x81 + x83 - V{4.5}*x57*x57) + V{1});
cell[21] = -x27*cell[21] + x44*(rho*(x42 + x60 + x75 + x78) + V{-1});
cell[22] = -x27*cell[22] - x44*(rho*(x51 + x83 - V{4.5}*x63*x63) + V{1});
cell[23] = -x27*cell[23] + x65*(rho*(x42 + x67 + x80) + V{-1});
cell[24] = -x27*cell[24] + x65*(rho*(x55 + x80 + V{4.5}*(x68*x68)) + V{-1});
cell[25] = -x27*cell[25] + x65*(rho*(x51 + x82 + V{4.5}*(x70*x70)) + V{-1});
cell[26] = -x27*cell[26] - x65*(rho*(x62 + x81 - V{4.5}*x84*x84) + V{1});
return x28 + x30 + x32;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = V{3}*u[0];
auto x29 = x28 + V{-1};
auto x30 = V{0.0740740740740741}*omega;
auto x31 = V{3}*u[1];
auto x32 = x31 + V{-1};
auto x33 = V{3}*u[2];
auto x34 = x29 + x31;
auto x35 = V{0.0185185185185185}*omega;
auto x36 = -x28;
auto x37 = x31 + V{1};
auto x38 = x36 + x37;
auto x39 = x29 + x33;
auto x40 = x33 + V{1};
auto x41 = -x31;
auto x42 = V{0.00462962962962963}*omega;
auto x43 = -x33;
auto x44 = x28 + V{1};
auto x45 = x31 + x44;
auto x46 = x41 + x44;
cell[0] = V{0.296296296296296}*omega*(rho + V{-1}) - x27*cell[0];
cell[1] = -x27*cell[1] - x30*(rho*x29 + V{1});
cell[2] = -x27*cell[2] - x30*(rho*x32 + V{1});
cell[3] = -x27*cell[3] - x30*(rho*(x33 + V{-1}) + V{1});
cell[4] = -x27*cell[4] - x35*(rho*x34 + V{1});
cell[5] = -x27*cell[5] + x35*(rho*x38 + V{-1});
cell[6] = -x27*cell[6] - x35*(rho*x39 + V{1});
cell[7] = -x27*cell[7] + x35*(rho*(x36 + x40) + V{-1});
cell[8] = -x27*cell[8] - x35*(rho*(x32 + x33) + V{1});
cell[9] = -x27*cell[9] + x35*(rho*(x40 + x41) + V{-1});
cell[10] = -x27*cell[10] - x42*(rho*(x33 + x34) + V{1});
cell[11] = -x27*cell[11] - x42*(rho*(x34 + x43) + V{1});
cell[12] = -x27*cell[12] - x42*(rho*(x39 + x41) + V{1});
cell[13] = -x27*cell[13] + x42*(rho*(x33 + x38) + V{-1});
cell[14] = -x27*cell[14] + x30*(rho*x44 + V{-1});
cell[15] = -x27*cell[15] + x30*(rho*x37 + V{-1});
cell[16] = -x27*cell[16] + x30*(rho*x40 + V{-1});
cell[17] = -x27*cell[17] + x35*(rho*x45 + V{-1});
cell[18] = -x27*cell[18] + x35*(rho*x46 + V{-1});
cell[19] = -x27*cell[19] + x35*(rho*(x33 + x44) + V{-1});
cell[20] = -x27*cell[20] + x35*(rho*(x43 + x44) + V{-1});
cell[21] = -x27*cell[21] + x35*(rho*(x33 + x37) + V{-1});
cell[22] = -x27*cell[22] + x35*(rho*(x37 + x43) + V{-1});
cell[23] = -x27*cell[23] + x42*(rho*(x33 + x45) + V{-1});
cell[24] = -x27*cell[24] + x42*(rho*(x43 + x45) + V{-1});
cell[25] = -x27*cell[25] + x42*(rho*(x33 + x46) + V{-1});
cell[26] = -x27*cell[26] + x42*(rho*(x43 + x46) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x27 = j[0]*j[0];
auto x28 = j[1]*j[1];
auto x29 = j[2]*j[2];
auto x30 = omega + V{-1};
auto x31 = V{0.111111111111111}*x28;
auto x32 = V{0.222222222222222}*j[0];
auto x33 = V{0.222222222222222}*x27;
auto x34 = V{0.111111111111111}*x29;
auto x35 = V{0.222222222222222}*pressure;
auto x36 = -x35;
auto x37 = x34 + x36 + V{0.0740740740740741};
auto x38 = V{0.111111111111111}*x27;
auto x39 = V{0.222222222222222}*j[1];
auto x40 = V{0.222222222222222}*x28;
auto x41 = V{0.222222222222222}*j[2];
auto x42 = V{0.222222222222222}*x29;
auto x43 = j[0] + j[1];
auto x44 = V{0.0833333333333333}*(x43*x43);
auto x45 = V{0.0555555555555556}*j[0];
auto x46 = V{0.0555555555555556}*j[1];
auto x47 = x45 + x46;
auto x48 = V{0.0277777777777778}*x27;
auto x49 = V{0.0277777777777778}*x28;
auto x50 = V{0.0277777777777778}*x29;
auto x51 = V{0.0555555555555556}*pressure;
auto x52 = x48 + x49 + x50 - x51 + V{0.0185185185185185};
auto x53 = -x46;
auto x54 = j[0] - j[1];
auto x55 = -x54;
auto x56 = x45 + x52;
auto x57 = V{0.0555555555555556}*j[2];
auto x58 = j[0] + j[2];
auto x59 = V{0.0833333333333333}*(x58*x58);
auto x60 = -x57;
auto x61 = -j[2];
auto x62 = x61 + j[0];
auto x63 = -x62;
auto x64 = j[1] + j[2];
auto x65 = V{0.0833333333333333}*(x64*x64);
auto x66 = x46 + x52;
auto x67 = x61 + j[1];
auto x68 = -x67;
auto x69 = x43 + j[2];
auto x70 = V{0.0208333333333333}*(x69*x69);
auto x71 = V{0.0138888888888889}*j[0];
auto x72 = V{0.0138888888888889}*j[1];
auto x73 = V{0.0138888888888889}*j[2];
auto x74 = x71 + x72 + x73;
auto x75 = V{0.00694444444444444}*x27;
auto x76 = V{0.00694444444444444}*x28;
auto x77 = V{0.00694444444444444}*x29;
auto x78 = V{0.0138888888888889}*pressure;
auto x79 = x75 + x76 + x77 - x78 + V{0.00462962962962963};
auto x80 = x43 + x61;
auto x81 = -x80;
auto x82 = -x73;
auto x83 = x71 + x79;
auto x84 = x82 + x83;
auto x85 = x54 + j[2];
auto x86 = -x85;
auto x87 = -x72;
auto x88 = x73 + x87;
auto x89 = x64 - j[0];
auto x90 = -x31;
auto x91 = -x34 + x35 + V{-0.0740740740740741};
auto x92 = -x38;
auto x93 = -x48 - x49 - x50 + x51 + V{-0.0185185185185185};
auto x94 = -x45;
auto x95 = x57 + x93;
auto x96 = x52 + x57;
auto x97 = -x71 + x79;
auto x98 = x72 + x97;
auto x99 = -x89;
cell[0] = -omega*(-V{0.888888888888889}*pressure + V{0.444444444444444}*x27 + V{0.444444444444444}*x28 + V{0.444444444444444}*x29 + V{0.296296296296296}) - x30*cell[0];
cell[1] = -omega*(x31 + x32 - x33 + x37) - x30*cell[1];
cell[2] = -omega*(x37 + x38 + x39 - x40) - x30*cell[2];
cell[3] = -omega*(x31 + x36 + x38 + x41 - x42 + V{0.0740740740740741}) - x30*cell[3];
cell[4] = -omega*(-x44 + x47 + x52) - x30*cell[4];
cell[5] = -omega*(x53 + x56 - V{0.0833333333333333}*x55*x55) - x30*cell[5];
cell[6] = -omega*(x56 + x57 - x59) - x30*cell[6];
cell[7] = -omega*(x56 + x60 - V{0.0833333333333333}*x63*x63) - x30*cell[7];
cell[8] = -omega*(x57 - x65 + x66) - x30*cell[8];
cell[9] = -omega*(x60 + x66 - V{0.0833333333333333}*x68*x68) - x30*cell[9];
cell[10] = -omega*(-x70 + x74 + x79) - x30*cell[10];
cell[11] = -omega*(x72 + x84 - V{0.0208333333333333}*x81*x81) - x30*cell[11];
cell[12] = -omega*(x83 + x88 - V{0.0208333333333333}*x86*x86) - x30*cell[12];
cell[13] = -omega*(x84 + x87 - V{0.0208333333333333}*x89*x89) - x30*cell[13];
cell[14] = omega*(x32 + x33 + x90 + x91) - x30*cell[14];
cell[15] = omega*(x39 + x40 + x91 + x92) - x30*cell[15];
cell[16] = omega*(x35 + x41 + x42 + x90 + x92 + V{-0.0740740740740741}) - x30*cell[16];
cell[17] = omega*(x44 + x47 + x93) - x30*cell[17];
cell[18] = -omega*(x66 + x94 - V{0.0833333333333333}*x54*x54) - x30*cell[18];
cell[19] = omega*(x45 + x59 + x95) - x30*cell[19];
cell[20] = -omega*(x94 + x96 - V{0.0833333333333333}*x62*x62) - x30*cell[20];
cell[21] = omega*(x46 + x65 + x95) - x30*cell[21];
cell[22] = -omega*(x53 + x96 - V{0.0833333333333333}*x67*x67) - x30*cell[22];
cell[23] = omega*(x70 + x74 - x75 - x76 - x77 + x78 + V{-0.00462962962962963}) - x30*cell[23];
cell[24] = -omega*(x88 + x97 - V{0.0208333333333333}*x80*x80) - x30*cell[24];
cell[25] = -omega*(x82 + x98 - V{0.0208333333333333}*x85*x85) - x30*cell[25];
cell[26] = -omega*(x73 + x98 - V{0.0208333333333333}*x99*x99) - x30*cell[26];
return x27 + x28 + x29;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = V{0.296296296296296}*rho;
auto x29 = u[0]*u[0];
auto x30 = V{1.5}*x29;
auto x31 = u[1]*u[1];
auto x32 = V{1.5}*x31;
auto x33 = u[2]*u[2];
auto x34 = V{1.5}*x33;
auto x35 = x32 + x34 + V{-1};
auto x36 = x30 + x35;
auto x37 = V{0.0740740740740741}*rho;
auto x38 = V{3}*u[0];
auto x39 = V{3}*x29;
auto x40 = x35 + x38 - x39;
auto x41 = ratioRho*x37;
auto x42 = V{3}*u[1];
auto x43 = V{3}*x31;
auto x44 = x30 + V{-1};
auto x45 = x34 + x42 - x43 + x44;
auto x46 = V{3}*u[2];
auto x47 = V{3}*x33;
auto x48 = x32 + x44 + x46 - x47;
auto x49 = V{0.0185185185185185}*rho;
auto x50 = u[0] + u[1];
auto x51 = V{4.5}*(x50*x50);
auto x52 = x36 + x38;
auto x53 = x42 + x52;
auto x54 = -x51 + x53;
auto x55 = ratioRho*x49;
auto x56 = u[0] - u[1];
auto x57 = -x56;
auto x58 = -x42;
auto x59 = x52 + x58;
auto x60 = x59 - V{4.5}*x57*x57;
auto x61 = u[0] + u[2];
auto x62 = V{4.5}*(x61*x61);
auto x63 = x46 + x52 - x62;
auto x64 = -x46;
auto x65 = -u[2];
auto x66 = x65 + u[0];
auto x67 = -x66;
auto x68 = x52 + x64 - V{4.5}*x67*x67;
auto x69 = u[1] + u[2];
auto x70 = V{4.5}*(x69*x69);
auto x71 = x36 + x42;
auto x72 = x46 + x71;
auto x73 = -x70 + x72;
auto x74 = x65 + u[1];
auto x75 = -x74;
auto x76 = x64 + x71 - V{4.5}*x75*x75;
auto x77 = V{0.00462962962962963}*rho;
auto x78 = x50 + u[2];
auto x79 = V{4.5}*(x78*x78);
auto x80 = x46 + x53 - x79;
auto x81 = ratioRho*x77;
auto x82 = x50 + x65;
auto x83 = -x82;
auto x84 = x53 + x64 - V{4.5}*x83*x83;
auto x85 = x56 + u[2];
auto x86 = -x85;
auto x87 = x46 + x59 - V{4.5}*x86*x86;
auto x88 = x69 - u[0];
auto x89 = -x59 - x64 + V{4.5}*(x88*x88);
auto x90 = -x32;
auto x91 = V{1} - x34;
auto x92 = x90 + x91;
auto x93 = x38 + x92;
auto x94 = x39 + x93;
auto x95 = -x30;
auto x96 = x42 + x95;
auto x97 = x43 + x91 + x96;
auto x98 = x46 + x95;
auto x99 = x47 + x90 + x98 + V{1};
auto x100 = x93 + x96;
auto x101 = x100 + x51;
auto x102 = -x38;
auto x103 = x102 + x71 - V{4.5}*x56*x56;
auto x104 = x93 + x98;
auto x105 = x104 + x62;
auto x106 = x36 + x46;
auto x107 = x102 + x106 - V{4.5}*x66*x66;
auto x108 = x46 + x70 + x92 + x96;
auto x109 = x106 + x58 - V{4.5}*x74*x74;
auto x110 = x100 + x46 + x79;
auto x111 = x100 + x64 + V{4.5}*(x82*x82);
auto x112 = x104 + x58 + V{4.5}*(x85*x85);
auto x113 = -x88;
auto x114 = x102 + x72 - V{4.5}*x113*x113;
cell[0] = -ratioRho*x28*x36 - x27*(x28*x36 + cell[0] + V{0.296296296296296}) + V{-0.296296296296296};
cell[1] = -x27*(x37*x40 + cell[1] + V{0.0740740740740741}) - x40*x41 + V{-0.0740740740740741};
cell[2] = -x27*(x37*x45 + cell[2] + V{0.0740740740740741}) - x41*x45 + V{-0.0740740740740741};
cell[3] = -x27*(x37*x48 + cell[3] + V{0.0740740740740741}) - x41*x48 + V{-0.0740740740740741};
cell[4] = -x27*(x49*x54 + cell[4] + V{0.0185185185185185}) - x54*x55 + V{-0.0185185185185185};
cell[5] = -x27*(x49*x60 + cell[5] + V{0.0185185185185185}) - x55*x60 + V{-0.0185185185185185};
cell[6] = -x27*(x49*x63 + cell[6] + V{0.0185185185185185}) - x55*x63 + V{-0.0185185185185185};
cell[7] = -x27*(x49*x68 + cell[7] + V{0.0185185185185185}) - x55*x68 + V{-0.0185185185185185};
cell[8] = -x27*(x49*x73 + cell[8] + V{0.0185185185185185}) - x55*x73 + V{-0.0185185185185185};
cell[9] = -x27*(x49*x76 + cell[9] + V{0.0185185185185185}) - x55*x76 + V{-0.0185185185185185};
cell[10] = -x27*(x77*x80 + cell[10] + V{0.00462962962962963}) - x80*x81 + V{-0.00462962962962963};
cell[11] = -x27*(x77*x84 + cell[11] + V{0.00462962962962963}) - x81*x84 + V{-0.00462962962962963};
cell[12] = -x27*(x77*x87 + cell[12] + V{0.00462962962962963}) - x81*x87 + V{-0.00462962962962963};
cell[13] = V{0.00462962962962963}*ratioRho*rho*x89 - x27*(-x77*x89 + cell[13] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[14] = V{0.0740740740740741}*ratioRho*rho*x94 - x27*(-x37*x94 + cell[14] + V{0.0740740740740741}) + V{-0.0740740740740741};
cell[15] = V{0.0740740740740741}*ratioRho*rho*x97 - x27*(-x37*x97 + cell[15] + V{0.0740740740740741}) + V{-0.0740740740740741};
cell[16] = V{0.0740740740740741}*ratioRho*rho*x99 - x27*(-x37*x99 + cell[16] + V{0.0740740740740741}) + V{-0.0740740740740741};
cell[17] = V{0.0185185185185185}*ratioRho*rho*x101 - x27*(-x101*x49 + cell[17] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[18] = -x103*x55 - x27*(x103*x49 + cell[18] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[19] = V{0.0185185185185185}*ratioRho*rho*x105 - x27*(-x105*x49 + cell[19] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[20] = -x107*x55 - x27*(x107*x49 + cell[20] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[21] = V{0.0185185185185185}*ratioRho*rho*x108 - x27*(-x108*x49 + cell[21] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[22] = -x109*x55 - x27*(x109*x49 + cell[22] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[23] = V{0.00462962962962963}*ratioRho*rho*x110 - x27*(-x110*x77 + cell[23] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[24] = V{0.00462962962962963}*ratioRho*rho*x111 - x27*(-x111*x77 + cell[24] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[25] = V{0.00462962962962963}*ratioRho*rho*x112 - x27*(-x112*x77 + cell[25] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[26] = -x114*x81 - x27*(x114*x77 + cell[26] + V{0.00462962962962963}) + V{-0.00462962962962963};
return x29 + x31 + x33;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = V{0.222222222222222}*cell[23];
auto x29 = V{0.222222222222222}*cell[24];
auto x30 = V{0.222222222222222}*cell[10];
auto x31 = -x30;
auto x32 = V{0.222222222222222}*cell[11];
auto x33 = -x32;
auto x34 = V{3}*u[2];
auto x35 = V{3}*u[1];
auto x36 = V{3}*u[0];
auto x37 = x36 + V{-1};
auto x38 = x35 + x37;
auto x39 = x34 + x38;
auto x40 = -x39;
auto x41 = V{0.00102880658436214}*rho;
auto x42 = x40*x41;
auto x43 = x36 + V{1};
auto x44 = x35 + x43;
auto x45 = x34 + x44;
auto x46 = x41*x45;
auto x47 = -x46;
auto x48 = -x34;
auto x49 = x44 + x48;
auto x50 = x41*x49;
auto x51 = -x50;
auto x52 = V{0.222222222222222}*cell[1];
auto x53 = -x37;
auto x54 = V{0.0164609053497942}*rho;
auto x55 = x43*x54;
auto x56 = V{0.222222222222222}*cell[4];
auto x57 = x38 + x48;
auto x58 = -x57;
auto x59 = x41*x58;
auto x60 = -x38;
auto x61 = V{0.00411522633744856}*rho;
auto x62 = x44*x61;
auto x63 = -x56 + x59 + x60*x61 - x62 + V{0.222222222222222}*cell[17];
auto x64 = V{0.222222222222222}*cell[25];
auto x65 = V{0.222222222222222}*cell[6];
auto x66 = V{0.222222222222222}*cell[12];
auto x67 = -x66;
auto x68 = -x35;
auto x69 = x34 + x37;
auto x70 = x68 + x69;
auto x71 = -x70;
auto x72 = x41*x71;
auto x73 = -x69;
auto x74 = x34 + x43;
auto x75 = x68 + x74;
auto x76 = x41*x75;
auto x77 = -x76;
auto x78 = x61*x74;
auto x79 = x61*x73 + x64 - x65 + x67 + x72 + x77 - x78 + V{0.222222222222222}*cell[19];
auto x80 = V{0.222222222222222}*cell[18];
auto x81 = V{0.222222222222222}*cell[26];
auto x82 = V{0.222222222222222}*cell[5];
auto x83 = V{0.222222222222222}*cell[13];
auto x84 = -x83;
auto x85 = -x36;
auto x86 = x35 + V{1};
auto x87 = x85 + x86;
auto x88 = x34 + x87;
auto x89 = x41*x88;
auto x90 = x61*x87;
auto x91 = x43 + x68;
auto x92 = x48 + x91;
auto x93 = x41*x92;
auto x94 = -x93;
auto x95 = x61*x91;
auto x96 = x80 + x81 - x82 + x84 + x89 + x90 + x94 - x95;
auto x97 = V{0.222222222222222}*cell[20];
auto x98 = V{0.222222222222222}*cell[7];
auto x99 = x34 + V{1};
auto x100 = x85 + x99;
auto x101 = x100*x61;
auto x102 = x43 + x48;
auto x103 = x102*x61;
auto x104 = x101 - x103 + x97 - x98;
auto x105 = V{0.0740740740740741}*rho;
auto x106 = V{0.222222222222222}*cell[2];
auto x107 = x35 + V{-1};
auto x108 = -x107;
auto x109 = x54*x86;
auto x110 = V{0.222222222222222}*cell[22];
auto x111 = V{0.222222222222222}*cell[9];
auto x112 = -x64;
auto x113 = x68 + x99;
auto x114 = x113*x61;
auto x115 = x48 + x86;
auto x116 = x115*x61;
auto x117 = x110 - x111 + x112 + x114 - x116 + x29 + x33 + x51 + x66 + x76;
auto x118 = -x81;
auto x119 = -x89;
auto x120 = x118 + x119 - x80 + x82 + x83 - x90 + x93 + x95;
auto x121 = V{0.222222222222222}*cell[8];
auto x122 = x107 + x34;
auto x123 = -x122;
auto x124 = x34 + x86;
auto x125 = x124*x61;
auto x126 = -x121 + x123*x61 - x125 + x28 + x31 + x42 + x47 + V{0.222222222222222}*cell[21];
auto x127 = V{0.222222222222222}*cell[3];
auto x128 = x34 + V{-1};
auto x129 = -x128;
auto x130 = x54*x99;
auto x131 = -x29;
auto x132 = -x101 + x103 + x131 + x32 + x50 - x97 + x98;
auto x133 = -x110 + x111 - x114 + x116;
auto x134 = V{0.0555555555555556}*cell[15];
auto x135 = V{0.0555555555555556}*cell[19];
auto x136 = V{0.0555555555555556}*cell[2];
auto x137 = -x136;
auto x138 = V{0.0555555555555556}*cell[6];
auto x139 = -x138;
auto x140 = x41*x74;
auto x141 = -x140;
auto x142 = x61*x86;
auto x143 = -x142;
auto x144 = V{0.111111111111111}*cell[4];
auto x145 = V{0.111111111111111}*cell[11];
auto x146 = V{0.00051440329218107}*rho;
auto x147 = V{0.00205761316872428}*rho;
auto x148 = x146*x49;
auto x149 = x147*x44;
auto x150 = V{0.111111111111111}*cell[23];
auto x151 = V{0.111111111111111}*cell[10];
auto x152 = -x151;
auto x153 = x146*x40;
auto x154 = x146*x45;
auto x155 = -x154;
auto x156 = x108*x61 + x150 + x152 + x153 + x155 + x41*x73;
auto x157 = V{0.0555555555555556}*cell[20];
auto x158 = V{0.0555555555555556}*cell[7];
auto x159 = x100*x41;
auto x160 = x102*x41;
auto x161 = x157 - x158 + x159 - x160;
auto x162 = V{0.0555555555555556}*cell[21];
auto x163 = V{0.0555555555555556}*cell[8];
auto x164 = x124*x41;
auto x165 = x162 - x163 - x164;
auto x166 = V{0.0555555555555556}*cell[22];
auto x167 = V{0.0555555555555556}*cell[9];
auto x168 = x113*x41;
auto x169 = x115*x41;
auto x170 = x166 - x167 + x168 - x169;
auto x171 = x165 + x170;
auto x172 = V{0.0555555555555556}*cell[1];
auto x173 = x43*x61;
auto x174 = x123*x41 - x172 - x173 + x53*x61 + V{0.0555555555555556}*cell[14] + V{1.5419764230905e-18};
auto x175 = V{0.0185185185185185}*rho;
auto x176 = x107*x61;
auto x177 = -x157 + x158 - x159 + x160;
auto x178 = V{0.111111111111111}*cell[13];
auto x179 = V{0.111111111111111}*cell[26];
auto x180 = -x179;
auto x181 = x146*x92;
auto x182 = x146*x88;
auto x183 = -x182;
auto x184 = x134 + x137 + x143 + x178 + x180 + x181 + x183;
auto x185 = x41*x69;
auto x186 = -x135 + x138 + x140 + x185;
auto x187 = -V{0.0555555555555556}*cell[14];
auto x188 = x37*x61;
auto x189 = x122*x41;
auto x190 = x172 + x173 + x187 + x188 - x189 + V{-1.5419764230905e-18};
auto x191 = V{0.111111111111111}*cell[12];
auto x192 = V{0.111111111111111}*cell[25];
auto x193 = x146*x75;
auto x194 = x146*x70;
auto x195 = x191 - x192 + x193 + x194;
auto x196 = x27*(V{0.00205761316872428}*rho*x87 - x147*x91 - x171 + x176 - x177 - x184 - x186 - x190 - x195 + V{0.111111111111111}*cell[18] - V{0.111111111111111}*cell[5]);
auto x197 = V{0.111111111111111}*cell[6];
auto x198 = x147*x74;
auto x199 = V{0.0555555555555556}*cell[4];
auto x200 = x41*x44;
auto x201 = x129*x61 - x199 - x200 + x41*x60 + V{0.0555555555555556}*cell[17];
auto x202 = V{0.0555555555555556}*cell[18];
auto x203 = V{0.0555555555555556}*cell[5];
auto x204 = -x203;
auto x205 = x41*x87;
auto x206 = x41*x91;
auto x207 = -x206;
auto x208 = -x191 + x192 - x193 + x202 + x204 + x205 + x207;
auto x209 = V{0.0555555555555556}*cell[16];
auto x210 = V{0.0555555555555556}*cell[3];
auto x211 = -x210;
auto x212 = x61*x99;
auto x213 = -x212;
auto x214 = -x166 + x167 - x168 + x169;
auto x215 = x165 + x209 + x211 + x213 + x214;
auto x216 = -V{0.111111111111111}*cell[24];
auto x217 = x146*x57;
auto x218 = -x202 + x203 - x205 + x206;
auto x219 = x128*x61;
auto x220 = x199 + x200 + x38*x41 - V{0.0555555555555556}*cell[17];
auto x221 = -x219 + x220;
auto x222 = x27*(V{0.00205761316872428}*rho*x100 - x102*x147 - x145 - x148 - x178 - x180 - x181 - x183 - x190 - x215 - x216 - x217 - x218 - x221 + V{0.111111111111111}*cell[20] - V{0.111111111111111}*cell[7]);
auto x223 = V{0.111111111111111}*cell[8];
auto x224 = x124*x147;
auto x225 = x135 + x139 + x141 + x177 + x209 + x211 + x213;
auto x226 = -x134;
auto x227 = x136 + x142 + x145 + x148 + x176 + x216 + x217 + x226;
auto x228 = x27*(V{0.00205761316872428}*rho*x113 - x115*x147 + x185 + x194 - x208 - x221 - x225 - x227 + V{0.111111111111111}*cell[22] - V{0.111111111111111}*cell[9]);
auto x229 = V{0.0277777777777778}*cell[21];
auto x230 = V{0.0277777777777778}*cell[8];
auto x231 = V{0.0416666666666667}*cell[10];
auto x232 = V{0.000192901234567901}*rho;
auto x233 = V{6.43004115226337e-05}*rho;
auto x234 = x233*x71;
auto x235 = x233*x58;
auto x236 = x232*x45;
auto x237 = x124*x146;
auto x238 = V{0.0138888888888889}*cell[14];
auto x239 = V{0.0138888888888889}*cell[24];
auto x240 = V{0.0138888888888889}*cell[25];
auto x241 = V{0.0138888888888889}*cell[1];
auto x242 = -x241;
auto x243 = V{0.0138888888888889}*cell[11];
auto x244 = -x243;
auto x245 = V{0.0138888888888889}*cell[12];
auto x246 = -x245;
auto x247 = x233*x49;
auto x248 = -x247;
auto x249 = x233*x75;
auto x250 = -x249;
auto x251 = x41*x43;
auto x252 = -x251;
auto x253 = x238 + x239 + x240 + x242 + x244 + x246 + x248 + x250 + x252 + V{3.85494105772624e-19};
auto x254 = V{0.0138888888888889}*cell[13];
auto x255 = V{0.0138888888888889}*cell[16];
auto x256 = V{0.0138888888888889}*cell[3];
auto x257 = -x256;
auto x258 = V{0.0138888888888889}*cell[26];
auto x259 = -x258;
auto x260 = x233*x92;
auto x261 = x233*x88;
auto x262 = -x261;
auto x263 = x41*x99;
auto x264 = -x263;
auto x265 = x254 + x255 + x257 + x259 + x260 + x262 + x264;
auto x266 = V{0.0277777777777778}*cell[4];
auto x267 = x41*x53;
auto x268 = x108*x41;
auto x269 = x146*x44;
auto x270 = V{0.0138888888888889}*cell[15];
auto x271 = V{0.0138888888888889}*cell[2];
auto x272 = x41*x86;
auto x273 = x270 - x271 - x272;
auto x274 = x146*x60 - x266 + x267 + x268 - x269 + x273 + V{0.0277777777777778}*cell[17];
auto x275 = V{0.0277777777777778}*cell[6];
auto x276 = x129*x41;
auto x277 = x146*x74;
auto x278 = x146*x73 - x275 + x276 - x277 + V{0.0277777777777778}*cell[19];
auto x279 = V{0.00462962962962963}*rho;
auto x280 = V{0.0416666666666667}*cell[11];
auto x281 = x232*x49;
auto x282 = -x255;
auto x283 = V{0.0277777777777778}*cell[7];
auto x284 = x102*x146;
auto x285 = x100*x146 + x256 + x263 + x282 - x283 - x284 + V{0.0277777777777778}*cell[20];
auto x286 = V{0.0277777777777778}*cell[22];
auto x287 = V{0.0277777777777778}*cell[9];
auto x288 = x113*x146;
auto x289 = x115*x146;
auto x290 = x286 - x287 + x288 - x289;
auto x291 = -x254;
auto x292 = -x260;
auto x293 = -x240 + x245 + x249 + x258 + x261 + x291 + x292;
auto x294 = V{0.0138888888888889}*cell[23];
auto x295 = V{0.0138888888888889}*cell[10];
auto x296 = x233*x45;
auto x297 = x233*x40 + x238 + x242 + x252 + x294 - x295 - x296 + V{3.85494105772624e-19};
auto x298 = V{0.0416666666666667}*cell[12];
auto x299 = x232*x75;
auto x300 = -x270;
auto x301 = -x239 + x243 + x247 + x271 + x272 + x300;
auto x302 = V{0.0277777777777778}*cell[5];
auto x303 = x146*x91;
auto x304 = x146*x87 - x302 - x303 + V{0.0277777777777778}*cell[18];
auto x305 = -x286 + x287 - x288 + x289;
auto x306 = x37*x41;
auto x307 = x107*x41;
auto x308 = x128*x41;
auto x309 = x122*x146 - x229 + x230 + x237 + x307 + x308;
auto x310 = -x294;
auto x311 = x233*x39;
auto x312 = x233*x57;
auto x313 = x295 + x296 + x310 + x311 - x312;
auto x314 = x233*x70;
auto x315 = x271 + x272 + x300 - x314;
auto x316 = x27*(x232*x88 - x232*x92 + x253 + x285 + x304 - x306 + x309 + x313 + x315 - V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[26]);
auto x317 = -V{0.222222222222222}*cell[23];
auto x318 = x39*x41;
auto x319 = x41*x57;
auto x320 = x30 + x317 + x318 + x319 + x38*x61 + x46 + x56 + x62 - V{0.222222222222222}*cell[17];
auto x321 = x41*x70;
auto x322 = x321 + x61*x69 + x65 + x78 - V{0.222222222222222}*cell[19];
auto x323 = x121 + x122*x61 + x125 - V{0.222222222222222}*cell[21];
auto x324 = -V{0.111111111111111}*cell[23];
auto x325 = x146*x39;
auto x326 = x151 + x154 + x186 + x324 + x325;
auto x327 = -x162 + x163 + x164 + x172 + x173 + x187 + x188 + x189 + V{-1.5419764230905e-18};
auto x328 = -x209 + x210 + x212 + x219 + x220;
auto x329 = -V{0.0138888888888889}*cell[14];
auto x330 = x146*x38 + x241 + x251 + x266 + x269 + x306 + x329 - V{0.0277777777777778}*cell[17] + V{-3.85494105772624e-19};
auto x331 = x146*x69 + x256 + x263 + x275 + x277 + x282 - V{0.0277777777777778}*cell[19];
cell[0] = V{0.296296296296296}*rho + V{-0.296296296296296};
cell[1] = -x105*x37 + x27*(x104 + x28 + x29 + x31 + x33 + x42 + x47 + x51 - x52 + x53*x54 - x55 + x63 + x79 + x96 + V{0.222222222222222}*cell[14] + V{6.16790569236198e-18}) + V{-0.0740740740740741};
cell[2] = -x105*x107 + x27*(-x106 + x108*x54 - x109 + x117 + x120 + x126 + x63 - x72 + V{0.222222222222222}*cell[15]) + V{-0.0740740740740741};
cell[3] = -x105*x128 + x27*(x118 + x119 + x126 - x127 + x129*x54 - x130 + x132 + x133 - x59 + x79 + x83 + x93 + V{0.222222222222222}*cell[16]) + V{-0.0740740740740741};
cell[4] = -x175*x38 + x27*(x134 + x135 + x137 + x139 + x141 + x143 - x144 - x145 + x146*x58 + x147*x60 - x148 - x149 + x156 + x161 + x171 + x174 + V{0.111111111111111}*cell[17] + V{0.111111111111111}*cell[24]) + V{-0.0185185185185185};
cell[5] = x175*x87 + x196 + V{-0.0185185185185185};
cell[6] = -x175*x69 + x27*(x146*x71 + x147*x73 + x150 + x152 + x153 + x155 + x174 - x197 - x198 + x201 + x208 + x215 + V{0.111111111111111}*cell[19]) + V{-0.0185185185185185};
cell[7] = x100*x175 + x222 + V{-0.0185185185185185};
cell[8] = -x122*x175 + x27*(x123*x147 + x156 + x184 + x201 + x218 - x223 - x224 + x225 + V{0.111111111111111}*cell[21]) + V{-0.0185185185185185};
cell[9] = x113*x175 + x228 + V{-0.0185185185185185};
cell[10] = x27*(x123*x146 + x229 - x230 - x231 + x232*x40 + x234 + x235 - x236 - x237 + x253 + x265 + x274 + x278 + V{0.0416666666666667}*cell[23]) - x279*x39 + V{-0.00462962962962963};
cell[11] = x27*(x232*x58 - x234 + x274 - x276 - x280 - x281 + x285 + x290 + x293 + x297 + V{0.0416666666666667}*cell[24]) - x279*x57 + V{-0.00462962962962963};
cell[12] = x27*(x232*x71 - x235 + x255 + x257 + x258 + x261 + x264 + x267 - x268 + x278 + x291 + x292 + x297 - x298 - x299 + x301 + x304 + x305 + V{0.0416666666666667}*cell[25]) - x279*x70 + V{-0.00462962962962963};
cell[13] = x279*x88 + x316 + V{-0.00462962962962963};
cell[14] = V{0.0740740740740741}*rho*x43 - x27*(-x112 - x120 - x132 - x320 - x322 - x37*x54 - x52 - x55 - x66 - x76 + V{0.222222222222222}*cell[14] + V{6.16790569236198e-18}) + V{-0.0740740740740741};
cell[15] = V{0.0740740740740741}*rho*x86 - x27*(-x106 - x107*x54 - x109 - x131 - x133 - x32 - x320 + x321 - x323 - x50 - x64 - x67 - x77 - x96 + V{0.222222222222222}*cell[15]) + V{-0.0740740740740741};
cell[16] = V{0.0740740740740741}*rho*x99 - x27*(-x104 - x117 - x127 - x128*x54 - x130 - x30 - x317 - x318 + x319 - x322 - x323 - x46 - x81 - x84 - x89 - x94 + V{0.222222222222222}*cell[16]) + V{-0.0740740740740741};
cell[17] = V{0.0185185185185185}*rho*x44 - x27*(-x144 - x147*x38 - x149 - x177 - x214 - x227 - x326 - x327 + V{0.111111111111111}*cell[17]) + V{-0.0185185185185185};
cell[18] = V{0.0185185185185185}*rho*x91 - x196 + V{-0.0185185185185185};
cell[19] = V{0.0185185185185185}*rho*x74 - x27*(-x147*x69 - x151 - x154 - x170 - x195 - x197 - x198 - x218 - x324 - x325 - x327 - x328 + V{0.111111111111111}*cell[19]) + V{-0.0185185185185185};
cell[20] = V{0.0185185185185185}*rho*x102 - x222 + V{-0.0185185185185185};
cell[21] = V{0.0185185185185185}*rho*x124 - x27*(-x122*x147 - x136 - x142 - x161 - x176 + x178 - x179 + x181 - x182 - x202 - x204 - x205 - x207 - x223 - x224 - x226 - x326 - x328 + V{0.111111111111111}*cell[21]) + V{-0.0185185185185185};
cell[22] = V{0.0185185185185185}*rho*x115 - x228 + V{-0.0185185185185185};
cell[23] = V{0.00462962962962963}*rho*x45 - x27*(-x231 - x232*x39 - x236 - x293 - x301 - x309 - x312 - x314 - x330 - x331 + V{0.0416666666666667}*cell[23]) + V{-0.00462962962962963};
cell[24] = V{0.00462962962962963}*rho*x49 - x27*(V{0.00051440329218107}*rho*x100 - x232*x57 - x240 - x246 - x250 - x265 - x280 - x281 - x283 - x284 - x295 - x296 - x305 - x307 + x308 - x310 - x311 - x315 - x330 + V{0.0277777777777778}*cell[20] + V{0.0416666666666667}*cell[24]) + V{-0.00462962962962963};
cell[25] = V{0.00462962962962963}*rho*x75 - x27*(V{0.00051440329218107}*rho*x87 - x232*x70 - x239 - x241 - x244 - x248 - x251 - x254 - x259 - x260 - x262 - x273 - x290 - x298 - x299 - x302 - x303 - x306 + x307 - x308 - x313 - x329 - x331 + V{0.0277777777777778}*cell[18] + V{0.0416666666666667}*cell[25] + V{3.85494105772624e-19}) + V{-0.00462962962962963};
cell[26] = V{0.00462962962962963}*rho*x92 - x316 + V{-0.00462962962962963};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x27 = u[0]*u[0];
auto x28 = V{1.5}*x27;
auto x29 = u[1]*u[1];
auto x30 = V{1.5}*x29;
auto x31 = u[2]*u[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = omega + V{-1};
auto x36 = V{0.111111111111111}*pi[3];
auto x37 = V{0.111111111111111}*pi[5];
auto x38 = x36 + x37 - V{0.222222222222222}*pi[0];
auto x39 = V{0.0740740740740741}*rho;
auto x40 = V{3}*u[0];
auto x41 = V{3}*x27;
auto x42 = V{0.111111111111111}*pi[0];
auto x43 = x37 + x42 - V{0.222222222222222}*pi[3];
auto x44 = V{3}*u[1];
auto x45 = V{3}*x29;
auto x46 = x28 + V{-1};
auto x47 = x36 + x42 - V{0.222222222222222}*pi[5];
auto x48 = V{3}*u[2];
auto x49 = V{3}*x31;
auto x50 = V{0.0185185185185185}*rho;
auto x51 = u[0] + u[1];
auto x52 = V{4.5}*(x51*x51);
auto x53 = x34 + x40;
auto x54 = x44 + x53;
auto x55 = V{0.166666666666667}*pi[1];
auto x56 = V{0.0555555555555556}*pi[0];
auto x57 = V{0.0555555555555556}*pi[3];
auto x58 = x56 + x57 - V{0.0277777777777778}*pi[5];
auto x59 = x35*(x55 + x58) + V{0.0185185185185185};
auto x60 = u[0] - u[1];
auto x61 = -x60;
auto x62 = -x44;
auto x63 = x53 + x62;
auto x64 = x35*(-x55 + x58) + V{0.0185185185185185};
auto x65 = u[0] + u[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = V{0.166666666666667}*pi[2];
auto x68 = V{0.0555555555555556}*pi[5];
auto x69 = x56 + x68 - V{0.0277777777777778}*pi[3];
auto x70 = x35*(x67 + x69) + V{0.0185185185185185};
auto x71 = -x48;
auto x72 = -u[2];
auto x73 = x72 + u[0];
auto x74 = -x73;
auto x75 = x35*(-x67 + x69) + V{0.0185185185185185};
auto x76 = u[1] + u[2];
auto x77 = V{4.5}*(x76*x76);
auto x78 = x34 + x44;
auto x79 = x48 + x78;
auto x80 = V{0.166666666666667}*pi[4];
auto x81 = V{0.0277777777777778}*pi[0];
auto x82 = x35*(x57 + x68 + x80 - x81) + V{0.0185185185185185};
auto x83 = x72 + u[1];
auto x84 = -x83;
auto x85 = -x35*(-x57 - x68 + x80 + x81) + V{0.0185185185185185};
auto x86 = V{0.00462962962962963}*rho;
auto x87 = x51 + u[2];
auto x88 = V{4.5}*(x87*x87);
auto x89 = V{0.0416666666666667}*pi[2];
auto x90 = V{0.0416666666666667}*pi[4];
auto x91 = V{0.0416666666666667}*pi[1];
auto x92 = V{0.0138888888888889}*pi[0];
auto x93 = V{0.0138888888888889}*pi[3];
auto x94 = V{0.0138888888888889}*pi[5];
auto x95 = x91 + x92 + x93 + x94;
auto x96 = x35*(x89 + x90 + x95) + V{0.00462962962962963};
auto x97 = x51 + x72;
auto x98 = -x97;
auto x99 = -x89;
auto x100 = -x90;
auto x101 = x35*(x100 + x95 + x99) + V{0.00462962962962963};
auto x102 = x60 + u[2];
auto x103 = -x102;
auto x104 = -x91 + x92 + x93 + x94;
auto x105 = x35*(x100 + x104 + x89) + V{0.00462962962962963};
auto x106 = x76 - u[0];
auto x107 = x35*(x104 + x90 + x99) + V{0.00462962962962963};
auto x108 = -x30;
auto x109 = V{1} - x32;
auto x110 = x108 + x109;
auto x111 = x110 + x40;
auto x112 = -x28;
auto x113 = x112 + x44;
auto x114 = x112 + x48;
auto x115 = x111 + x113;
auto x116 = -x40;
auto x117 = x111 + x114;
auto x118 = x34 + x48;
auto x119 = -x106;
cell[0] = -V{0.296296296296296}*rho*x34 + V{0.444444444444444}*x35*(pi[0] + pi[3] + pi[5]) + V{-0.296296296296296};
cell[1] = x35*x38 - x39*(x33 + x40 - x41) + V{-0.0740740740740741};
cell[2] = x35*x43 - x39*(x32 + x44 - x45 + x46) + V{-0.0740740740740741};
cell[3] = x35*x47 - x39*(x30 + x46 + x48 - x49) + V{-0.0740740740740741};
cell[4] = -x50*(-x52 + x54) - x59;
cell[5] = -x50*(x63 - V{4.5}*x61*x61) - x64;
cell[6] = -x50*(x48 + x53 - x66) - x70;
cell[7] = -x50*(x53 + x71 - V{4.5}*x74*x74) - x75;
cell[8] = -x50*(-x77 + x79) - x82;
cell[9] = -x50*(x71 + x78 - V{4.5}*x84*x84) - x85;
cell[10] = -x86*(x48 + x54 - x88) - x96;
cell[11] = -x101 - x86*(x54 + x71 - V{4.5}*x98*x98);
cell[12] = -x105 - x86*(x48 + x63 - V{4.5}*x103*x103);
cell[13] = -V{0.00462962962962963}*rho*(x63 + x71 - V{4.5}*x106*x106) - x107;
cell[14] = x35*x38 + x39*(x111 + x41) + V{-0.0740740740740741};
cell[15] = x35*x43 + x39*(x109 + x113 + x45) + V{-0.0740740740740741};
cell[16] = x35*x47 + x39*(x108 + x114 + x49 + V{1}) + V{-0.0740740740740741};
cell[17] = V{0.0185185185185185}*rho*(x115 + x52) - x59;
cell[18] = -x50*(x116 + x78 - V{4.5}*x60*x60) - x64;
cell[19] = V{0.0185185185185185}*rho*(x117 + x66) - x70;
cell[20] = -x50*(x116 + x118 - V{4.5}*x73*x73) - x75;
cell[21] = V{0.0185185185185185}*rho*(x110 + x113 + x48 + x77) - x82;
cell[22] = -x50*(x118 + x62 - V{4.5}*x83*x83) - x85;
cell[23] = V{0.00462962962962963}*rho*(x115 + x48 + x88) - x96;
cell[24] = V{0.00462962962962963}*rho*(x115 + x71 + V{4.5}*(x97*x97)) - x101;
cell[25] = V{0.00462962962962963}*rho*(x117 + x62 + V{4.5}*(x102*x102)) - x105;
cell[26] = -x107 - x86*(x116 + x79 - V{4.5}*x119*x119);
return x27 + x29 + x31;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x27 = V{3}*newU[0];
auto x28 = x27 + V{-1};
auto x29 = V{3}*newU[1];
auto x30 = x29 + V{-1};
auto x31 = V{3}*newU[2];
auto x32 = x28 + x29;
auto x33 = -x27;
auto x34 = x29 + V{1};
auto x35 = x33 + x34;
auto x36 = x28 + x31;
auto x37 = x31 + V{1};
auto x38 = -x29;
auto x39 = -x31;
auto x40 = x27 + V{1};
auto x41 = x29 + x40;
auto x42 = x38 + x40;
cell[0] = V{0.296296296296296}*newRho + V{-0.296296296296296};
cell[1] = -V{0.0740740740740741}*newRho*x28 + V{-0.0740740740740741};
cell[2] = -V{0.0740740740740741}*newRho*x30 + V{-0.0740740740740741};
cell[3] = -V{0.0740740740740741}*newRho*(x31 + V{-1}) + V{-0.0740740740740741};
cell[4] = -V{0.0185185185185185}*newRho*x32 + V{-0.0185185185185185};
cell[5] = V{0.0185185185185185}*newRho*x35 + V{-0.0185185185185185};
cell[6] = -V{0.0185185185185185}*newRho*x36 + V{-0.0185185185185185};
cell[7] = V{0.0185185185185185}*newRho*(x33 + x37) + V{-0.0185185185185185};
cell[8] = -V{0.0185185185185185}*newRho*(x30 + x31) + V{-0.0185185185185185};
cell[9] = V{0.0185185185185185}*newRho*(x37 + x38) + V{-0.0185185185185185};
cell[10] = -V{0.00462962962962963}*newRho*(x31 + x32) + V{-0.00462962962962963};
cell[11] = -V{0.00462962962962963}*newRho*(x32 + x39) + V{-0.00462962962962963};
cell[12] = -V{0.00462962962962963}*newRho*(x36 + x38) + V{-0.00462962962962963};
cell[13] = V{0.00462962962962963}*newRho*(x31 + x35) + V{-0.00462962962962963};
cell[14] = V{0.0740740740740741}*newRho*x40 + V{-0.0740740740740741};
cell[15] = V{0.0740740740740741}*newRho*x34 + V{-0.0740740740740741};
cell[16] = V{0.0740740740740741}*newRho*x37 + V{-0.0740740740740741};
cell[17] = V{0.0185185185185185}*newRho*x41 + V{-0.0185185185185185};
cell[18] = V{0.0185185185185185}*newRho*x42 + V{-0.0185185185185185};
cell[19] = V{0.0185185185185185}*newRho*(x31 + x40) + V{-0.0185185185185185};
cell[20] = V{0.0185185185185185}*newRho*(x39 + x40) + V{-0.0185185185185185};
cell[21] = V{0.0185185185185185}*newRho*(x31 + x34) + V{-0.0185185185185185};
cell[22] = V{0.0185185185185185}*newRho*(x34 + x39) + V{-0.0185185185185185};
cell[23] = V{0.00462962962962963}*newRho*(x31 + x41) + V{-0.00462962962962963};
cell[24] = V{0.00462962962962963}*newRho*(x39 + x41) + V{-0.00462962962962963};
cell[25] = V{0.00462962962962963}*newRho*(x31 + x42) + V{-0.00462962962962963};
cell[26] = V{0.00462962962962963}*newRho*(x39 + x42) + V{-0.00462962962962963};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x27 = oldU[0]*oldU[0];
auto x28 = V{1.5}*x27;
auto x29 = oldU[1]*oldU[1];
auto x30 = V{1.5}*x29;
auto x31 = oldU[2]*oldU[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = newU[0]*newU[0];
auto x36 = V{1.5}*x35;
auto x37 = newU[1]*newU[1];
auto x38 = V{1.5}*x37;
auto x39 = newU[2]*newU[2];
auto x40 = V{1.5}*x39;
auto x41 = x38 + x40 + V{-1};
auto x42 = x36 + x41;
auto x43 = V{0.0740740740740741}*oldRho;
auto x44 = V{3}*oldU[0];
auto x45 = V{3}*x27;
auto x46 = V{0.0740740740740741}*newRho;
auto x47 = V{3}*newU[0];
auto x48 = V{3}*x35;
auto x49 = V{3}*oldU[1];
auto x50 = V{3}*x29;
auto x51 = x28 + V{-1};
auto x52 = V{3}*newU[1];
auto x53 = V{3}*x37;
auto x54 = x36 + V{-1};
auto x55 = V{3}*oldU[2];
auto x56 = V{3}*x31;
auto x57 = V{3}*newU[2];
auto x58 = V{3}*x39;
auto x59 = V{0.0185185185185185}*oldRho;
auto x60 = oldU[0] + oldU[1];
auto x61 = V{4.5}*(x60*x60);
auto x62 = x34 + x44;
auto x63 = x49 + x62;
auto x64 = V{0.0185185185185185}*newRho;
auto x65 = newU[0] + newU[1];
auto x66 = V{4.5}*(x65*x65);
auto x67 = x42 + x47;
auto x68 = x52 + x67;
auto x69 = oldU[0] - oldU[1];
auto x70 = -V{4.5}*x69*x69;
auto x71 = -x49;
auto x72 = x62 + x71;
auto x73 = newU[0] - newU[1];
auto x74 = -V{4.5}*x73*x73;
auto x75 = -x52;
auto x76 = x67 + x75;
auto x77 = oldU[0] + oldU[2];
auto x78 = V{4.5}*(x77*x77);
auto x79 = newU[0] + newU[2];
auto x80 = V{4.5}*(x79*x79);
auto x81 = -x55;
auto x82 = -oldU[2];
auto x83 = x82 + oldU[0];
auto x84 = -V{4.5}*x83*x83;
auto x85 = -x57;
auto x86 = -newU[2];
auto x87 = x86 + newU[0];
auto x88 = -V{4.5}*x87*x87;
auto x89 = oldU[1] + oldU[2];
auto x90 = V{4.5}*(x89*x89);
auto x91 = x34 + x49;
auto x92 = x55 + x91;
auto x93 = newU[1] + newU[2];
auto x94 = V{4.5}*(x93*x93);
auto x95 = x42 + x52;
auto x96 = x57 + x95;
auto x97 = x82 + oldU[1];
auto x98 = -V{4.5}*x97*x97;
auto x99 = x86 + newU[1];
auto x100 = -V{4.5}*x99*x99;
auto x101 = V{0.00462962962962963}*oldRho;
auto x102 = x60 + oldU[2];
auto x103 = V{4.5}*(x102*x102);
auto x104 = V{0.00462962962962963}*newRho;
auto x105 = x65 + newU[2];
auto x106 = V{4.5}*(x105*x105);
auto x107 = x60 + x82;
auto x108 = V{4.5}*(x107*x107);
auto x109 = x65 + x86;
auto x110 = V{4.5}*(x109*x109);
auto x111 = x69 + oldU[2];
auto x112 = V{4.5}*(x111*x111);
auto x113 = x73 + newU[2];
auto x114 = V{4.5}*(x113*x113);
auto x115 = x93 - newU[0];
auto x116 = -V{4.5}*x115*x115;
auto x117 = x89 - oldU[0];
auto x118 = -V{4.5}*x117*x117;
auto x119 = -x38;
auto x120 = V{1} - x40;
auto x121 = x119 + x120;
auto x122 = x121 + x47;
auto x123 = -x30;
auto x124 = V{1} - x32;
auto x125 = x123 + x124;
auto x126 = x125 + x44;
auto x127 = -x36;
auto x128 = x127 + x52;
auto x129 = -x28;
auto x130 = x129 + x49;
auto x131 = x127 + x57;
auto x132 = x129 + x55;
auto x133 = x122 + x128;
auto x134 = x126 + x130;
auto x135 = -x44;
auto x136 = -x47;
auto x137 = x122 + x131;
auto x138 = x126 + x132;
auto x139 = x34 + x55;
auto x140 = x42 + x57;
cell[0] = -V{0.296296296296296}*newRho*x42 + V{0.296296296296296}*oldRho*x34 + cell[0];
cell[1] = x43*(x33 + x44 - x45) - x46*(x41 + x47 - x48) + cell[1];
cell[2] = x43*(x32 + x49 - x50 + x51) - x46*(x40 + x52 - x53 + x54) + cell[2];
cell[3] = x43*(x30 + x51 + x55 - x56) - x46*(x38 + x54 + x57 - x58) + cell[3];
cell[4] = x59*(-x61 + x63) - x64*(-x66 + x68) + cell[4];
cell[5] = x59*(x70 + x72) - x64*(x74 + x76) + cell[5];
cell[6] = x59*(x55 + x62 - x78) - x64*(x57 + x67 - x80) + cell[6];
cell[7] = x59*(x62 + x81 + x84) - x64*(x67 + x85 + x88) + cell[7];
cell[8] = x59*(-x90 + x92) - x64*(-x94 + x96) + cell[8];
cell[9] = x59*(x81 + x91 + x98) - x64*(x100 + x85 + x95) + cell[9];
cell[10] = x101*(-x103 + x55 + x63) - x104*(-x106 + x57 + x68) + cell[10];
cell[11] = x101*(-x108 + x63 + x81) - x104*(-x110 + x68 + x85) + cell[11];
cell[12] = x101*(-x112 + x55 + x72) - x104*(-x114 + x57 + x76) + cell[12];
cell[13] = -x101*(-x118 - x72 - x81) + x104*(-x116 - x76 - x85) + cell[13];
cell[14] = -x43*(x126 + x45) + x46*(x122 + x48) + cell[14];
cell[15] = -x43*(x124 + x130 + x50) + x46*(x120 + x128 + x53) + cell[15];
cell[16] = -x43*(x123 + x132 + x56 + V{1}) + x46*(x119 + x131 + x58 + V{1}) + cell[16];
cell[17] = -x59*(x134 + x61) + x64*(x133 + x66) + cell[17];
cell[18] = x59*(x135 + x70 + x91) - x64*(x136 + x74 + x95) + cell[18];
cell[19] = -x59*(x138 + x78) + x64*(x137 + x80) + cell[19];
cell[20] = x59*(x135 + x139 + x84) - x64*(x136 + x140 + x88) + cell[20];
cell[21] = -x59*(x125 + x130 + x55 + x90) + x64*(x121 + x128 + x57 + x94) + cell[21];
cell[22] = x59*(x139 + x71 + x98) - x64*(x100 + x140 + x75) + cell[22];
cell[23] = -x101*(x103 + x134 + x55) + x104*(x106 + x133 + x57) + cell[23];
cell[24] = -x101*(x108 + x134 + x81) + x104*(x110 + x133 + x85) + cell[24];
cell[25] = -x101*(x112 + x138 + x71) + x104*(x114 + x137 + x75) + cell[25];
cell[26] = x101*(x118 + x135 + x92) - x104*(x116 + x136 + x96) + cell[26];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x27 = u[0]*u[0];
auto x28 = V{1.5}*x27;
auto x29 = u[1]*u[1];
auto x30 = V{1.5}*x29;
auto x31 = u[2]*u[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = V{0.0740740740740741}*rho;
auto x36 = V{3}*u[0];
auto x37 = V{3}*x27;
auto x38 = V{0.111111111111111}*pi[3];
auto x39 = V{0.111111111111111}*pi[5] + V{0.0740740740740741};
auto x40 = x38 + x39 - V{0.222222222222222}*pi[0];
auto x41 = V{3}*u[1];
auto x42 = V{3}*x29;
auto x43 = x28 + V{-1};
auto x44 = V{0.111111111111111}*pi[0];
auto x45 = x39 + x44 - V{0.222222222222222}*pi[3];
auto x46 = V{3}*u[2];
auto x47 = V{3}*x31;
auto x48 = x38 + x44 - V{0.222222222222222}*pi[5] + V{0.0740740740740741};
auto x49 = V{0.166666666666667}*pi[1];
auto x50 = V{0.0185185185185185}*rho;
auto x51 = u[0] + u[1];
auto x52 = V{4.5}*(x51*x51);
auto x53 = x34 + x36;
auto x54 = x41 + x53;
auto x55 = V{0.0277777777777778}*pi[5];
auto x56 = -V{0.0555555555555556}*pi[3];
auto x57 = V{0.0185185185185185} - V{0.0555555555555556}*pi[0];
auto x58 = x55 + x56 + x57;
auto x59 = u[0] - u[1];
auto x60 = -x59;
auto x61 = -x41;
auto x62 = x53 + x61;
auto x63 = x49 + x58;
auto x64 = V{0.166666666666667}*pi[2];
auto x65 = u[0] + u[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = V{0.0277777777777778}*pi[3];
auto x68 = -V{0.0555555555555556}*pi[5];
auto x69 = x57 + x67 + x68;
auto x70 = -x46;
auto x71 = -u[2];
auto x72 = x71 + u[0];
auto x73 = -x72;
auto x74 = x64 + x69;
auto x75 = V{0.166666666666667}*pi[4];
auto x76 = u[1] + u[2];
auto x77 = V{4.5}*(x76*x76);
auto x78 = x34 + x41;
auto x79 = x46 + x78;
auto x80 = V{0.0277777777777778}*pi[0];
auto x81 = x56 + x68 + x80 + V{0.0185185185185185};
auto x82 = x71 + u[1];
auto x83 = -x82;
auto x84 = x75 + x81;
auto x85 = V{0.00462962962962963}*rho;
auto x86 = x51 + u[2];
auto x87 = V{4.5}*(x86*x86);
auto x88 = V{0.0416666666666667}*pi[2];
auto x89 = V{0.0416666666666667}*pi[4];
auto x90 = x88 + x89;
auto x91 = V{0.0416666666666667}*pi[1];
auto x92 = V{0.0138888888888889}*pi[0] + V{0.0138888888888889}*pi[3] + V{0.0138888888888889}*pi[5] + V{-0.00462962962962963};
auto x93 = x91 + x92;
auto x94 = x90 + x93;
auto x95 = -x91;
auto x96 = x51 + x71;
auto x97 = -x96;
auto x98 = -V{0.0138888888888889}*pi[0] - V{0.0138888888888889}*pi[3] - V{0.0138888888888889}*pi[5] + V{0.00462962962962963};
auto x99 = x59 + u[2];
auto x100 = -x99;
auto x101 = -x88;
auto x102 = x101 + x89;
auto x103 = x91 + x98;
auto x104 = x76 - u[0];
auto x105 = x92 + x95;
auto x106 = -x30;
auto x107 = V{1} - x32;
auto x108 = x106 + x107;
auto x109 = x108 + x36;
auto x110 = -x28;
auto x111 = x110 + x41;
auto x112 = x110 + x46;
auto x113 = V{0.0555555555555556}*pi[3];
auto x114 = x109 + x111;
auto x115 = V{0.0555555555555556}*pi[0] + V{-0.0185185185185185};
auto x116 = -x36;
auto x117 = V{0.0555555555555556}*pi[5];
auto x118 = x109 + x112;
auto x119 = x34 + x46;
auto x120 = -x89;
auto x121 = x120 + x88;
auto x122 = -x104;
cell[0] = -V{0.296296296296296}*rho*x34 - V{0.444444444444444}*pi[0] - V{0.444444444444444}*pi[3] - V{0.444444444444444}*pi[5] + V{-0.296296296296296};
cell[1] = -x35*(x33 + x36 - x37) - x40;
cell[2] = -x35*(x32 + x41 - x42 + x43) - x45;
cell[3] = -x35*(x30 + x43 + x46 - x47) - x48;
cell[4] = x49 - x50*(-x52 + x54) - x58;
cell[5] = -x50*(x62 - V{4.5}*x60*x60) - x63;
cell[6] = -x50*(x46 + x53 - x66) + x64 - x69;
cell[7] = -x50*(x53 + x70 - V{4.5}*x73*x73) - x74;
cell[8] = -x50*(-x77 + x79) + x75 - x81;
cell[9] = -x50*(x70 + x78 - V{4.5}*x83*x83) - x84;
cell[10] = -x85*(x46 + x54 - x87) + x94;
cell[11] = -x85*(x54 + x70 - V{4.5}*x97*x97) - x90 - x95 - x98;
cell[12] = -x102 - x103 - x85*(x46 + x62 - V{4.5}*x100*x100);
cell[13] = x102 + x105 - x85*(x62 + x70 - V{4.5}*x104*x104);
cell[14] = V{0.0740740740740741}*rho*(x109 + x37) - x40;
cell[15] = V{0.0740740740740741}*rho*(x107 + x111 + x42) - x45;
cell[16] = V{0.0740740740740741}*rho*(x106 + x112 + x47 + V{1}) - x48;
cell[17] = x113 + x115 + x49 + x50*(x114 + x52) - x55;
cell[18] = -x50*(x116 + x78 - V{4.5}*x59*x59) - x63;
cell[19] = x115 + x117 + x50*(x118 + x66) + x64 - x67;
cell[20] = -x50*(x116 + x119 - V{4.5}*x72*x72) - x74;
cell[21] = x113 + x117 + x50*(x108 + x111 + x46 + x77) + x75 - x80 + V{-0.0185185185185185};
cell[22] = -x50*(x119 + x61 - V{4.5}*x82*x82) - x84;
cell[23] = x85*(x114 + x46 + x87) + x94;
cell[24] = x101 + x120 + x85*(x114 + x70 + V{4.5}*(x96*x96)) + x93;
cell[25] = x105 + x121 + x85*(x118 + x61 + V{4.5}*(x99*x99));
cell[26] = -x103 - x121 - x85*(x116 + x79 - V{4.5}*x122*x122);

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = V{1}*cell[12];
auto x1 = V{1}*cell[25];
auto x2 = -cell[23];
auto x3 = x2 + cell[10] + cell[12] - cell[19] - cell[25] + cell[6];
auto x4 = -cell[13] - cell[21] + cell[26] + cell[8];
auto x5 = cell[20] + cell[22] + cell[24] + cell[3];
auto x6 = x3 + x4 + x5 - cell[11] - cell[16] - cell[7] - cell[9];
auto x7 = cell[11] - cell[17] - cell[24] + cell[4];
auto x8 = cell[13] + cell[1] + cell[5] + cell[7];
auto x9 = x3 + x7 + x8 - cell[14] - cell[18] - cell[20] - cell[26];
auto x10 = cell[10] + cell[18] + cell[25] + cell[2] + cell[9];
auto x11 = x10 + x5 + x8 + cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8];
auto x12 = V{1} / (x11 + V{1});
auto x13 = x12*(x11 + V{1});
auto x14 = V{0.5}*x13;
auto x15 = x12*x9;
auto x16 = V{1}*x15;
auto x17 = V{1}*cell[13];
auto x18 = V{1}*cell[26];
auto x19 = -V{1}*cell[10];
auto x20 = -V{1}*cell[23];
auto x21 = x17 + x18 + x19 + x20;
auto x22 = V{1}*cell[11];
auto x23 = V{1}*cell[24];
auto x24 = x22 + x23;
auto x25 = x0 + x1 - x14*(x6*force[0] + x9*force[2]) - x16*x6 - x21 - x24 + V{1}*cell[19] - V{1}*cell[20] + V{1}*cell[6] - V{1}*cell[7];
auto x26 = x10 + x2 + x4 + x7 - cell[12] - cell[15] - cell[22] - cell[5];
auto x27 = x26*force[0] + x9*force[1];
auto x28 = x0 + x1;
auto x29 = V{2}*cell[13];
auto x30 = V{2}*cell[26];
auto x31 = V{2}*cell[11];
auto x32 = V{2}*cell[24];
auto x33 = V{1}*x13;
auto x34 = V{2}*x26;
auto x35 = -V{2}*cell[10] + V{2}*cell[12] - V{2}*cell[23] + V{2}*cell[25];
auto x36 = x26*force[2] + x6*force[1];
auto x37 = x12*x6;
auto x38 = V{0.666666666666667}*cell[10];
auto x39 = V{0.666666666666667}*cell[11];
auto x40 = V{0.666666666666667}*cell[12];
auto x41 = V{0.666666666666667}*cell[13];
auto x42 = V{0.666666666666667}*cell[23];
auto x43 = V{0.666666666666667}*cell[24];
auto x44 = V{0.666666666666667}*cell[25];
auto x45 = V{0.666666666666667}*cell[26];
auto x46 = -V{0.333333333333333}*cell[0];
auto x47 = x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45 + x46 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x48 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x49 = -x12*x9*x9 - x13*x9*force[0] + x47 + x48 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x50 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x51 = -x12*x26*x26 - x13*x26*force[1] + x47 + x50 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x52 = -x12*x6*x6 - x13*x6*force[2] + x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x48 + x50 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return (-x14*x27 - x16*x26 - x21 + x22 + x23 - x28 + V{1}*cell[17] - V{1}*cell[18] + V{1}*cell[4] - V{1}*cell[5])*(-x15*x34 - x27*x33 - x29 - x30 + x31 + x32 - x35 + V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[4] - V{2}*cell[5]) + (x29 + x30 - x31 - x32 - x33*x36 - x34*x37 - x35 + V{2}*cell[21] - V{2}*cell[22] + V{2}*cell[8] - V{2}*cell[9])*(-x14*x36 + x17 + x18 - x19 - x20 - x24 - V{1}*x26*x37 - x28 + V{1}*cell[21] - V{1}*cell[22] + V{1}*cell[8] - V{1}*cell[9]) + 2*(x25*x25) + x49*x49 + x51*x51 + x52*x52;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = -cell[11];
auto x1 = -cell[12];
auto x2 = -cell[17];
auto x3 = -cell[24];
auto x4 = x2 + x3 + cell[18];
auto x5 = cell[10] + cell[11] + cell[4];
auto x6 = -cell[23];
auto x7 = -cell[13] - cell[21];
auto x8 = x6 + x7 + cell[26] + cell[8];
auto x9 = x1 + x4 + x5 + x8 - cell[15] - cell[22] + cell[25] + cell[2] - cell[5] + cell[9];
auto x10 = -cell[19];
auto x11 = -cell[25];
auto x12 = x10 + x11 + cell[7];
auto x13 = x6 - cell[26];
auto x14 = cell[12] + cell[6];
auto x15 = cell[13] + cell[1];
auto x16 = x12 + x13 + x14 + x15 + x2 + x3 + x5 - cell[14] - cell[18] - cell[20] + cell[5];
auto x17 = cell[12] + cell[25];
auto x18 = x17 + cell[5];
auto x19 = cell[11] + cell[24];
auto x20 = x19 + cell[20];
auto x21 = cell[22] + cell[9];
auto x22 = cell[10] + cell[3];
auto x23 = V{1} / (x15 + x18 + x20 + x21 + x22 + cell[0] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[21] + cell[23] + cell[26] + cell[2] + cell[4] + cell[6] + cell[7] + cell[8] + V{1});
auto x24 = x16*x23;
auto x25 = -cell[10];
auto x26 = x25 + x6 + cell[13] + cell[26];
auto x27 = x0 + x18 + x24*x9 + x26 + x4 - cell[4];
auto x28 = x0 + x10 + x11 + x14 + x22 + x8 - cell[16] + cell[20] + cell[22] + cell[24] - cell[7] - cell[9];
auto x29 = x1 + x12 + x20 + x24*x28 + x26 - cell[6];
auto x30 = x13 + x17 + x19 + x21 + x23*x28*x9 + x25 + x7 - cell[8];
auto x31 = V{1}*x23;
auto x32 = V{0.666666666666667}*cell[10];
auto x33 = V{0.666666666666667}*cell[11];
auto x34 = V{0.666666666666667}*cell[12];
auto x35 = V{0.666666666666667}*cell[13];
auto x36 = V{0.666666666666667}*cell[23];
auto x37 = V{0.666666666666667}*cell[24];
auto x38 = V{0.666666666666667}*cell[25];
auto x39 = V{0.666666666666667}*cell[26];
auto x40 = -V{0.333333333333333}*cell[0];
auto x41 = x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x42 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x43 = -x31*x16*x16 + x41 + x42 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x44 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x45 = -x31*x9*x9 + x41 + x44 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x46 = -x31*x28*x28 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 + x42 + x44 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return V{2}*(x27*x27) + V{2}*(x29*x29) + V{2}*(x30*x30) + x43*x43 + x45*x45 + x46*x46;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x27 = force[0]*u[0];
auto x28 = force[1]*u[1];
auto x29 = force[2]*u[2];
auto x30 = rho*(V{0.5}*omega + V{-1});
auto x31 = V{6}*u[0];
auto x32 = x31 + V{-3};
auto x33 = V{0.0740740740740741}*force[0];
auto x34 = V{0.222222222222222}*x28;
auto x35 = V{0.222222222222222}*x29;
auto x36 = x34 + x35;
auto x37 = V{6}*u[1];
auto x38 = x37 + V{-3};
auto x39 = V{0.0740740740740741}*force[1];
auto x40 = V{0.222222222222222}*x27;
auto x41 = x35 + x40;
auto x42 = V{6}*u[2];
auto x43 = x42 + V{-3};
auto x44 = V{0.0740740740740741}*force[2];
auto x45 = x34 + x40;
auto x46 = V{9}*u[1];
auto x47 = x32 + x46;
auto x48 = V{0.0185185185185185}*force[0];
auto x49 = V{9}*u[0];
auto x50 = x38 + x49;
auto x51 = V{0.0185185185185185}*force[1];
auto x52 = V{0.0555555555555556}*x29;
auto x53 = -x52;
auto x54 = V{3} - x31;
auto x55 = x46 + x54;
auto x56 = -x49;
auto x57 = x37 + V{3};
auto x58 = x56 + x57;
auto x59 = V{9}*u[2];
auto x60 = x32 + x59;
auto x61 = x43 + x49;
auto x62 = V{0.0185185185185185}*force[2];
auto x63 = V{0.0555555555555556}*x28;
auto x64 = -x63;
auto x65 = x42 + V{3};
auto x66 = x56 + x65;
auto x67 = V{0.0555555555555556}*x27;
auto x68 = -x67;
auto x69 = -x37;
auto x70 = x69 + V{3};
auto x71 = x59 + x70;
auto x72 = -x46;
auto x73 = x65 + x72;
auto x74 = V{0.00462962962962963}*x30;
auto x75 = -x59;
auto x76 = -x42;
auto x77 = x49 + V{-3};
auto x78 = x31 + V{3};
auto x79 = x46 + x78;
auto x80 = x49 + x57;
auto x81 = x49 + x70;
auto x82 = x72 + x78;
auto x83 = x49 + x65;
auto x84 = x76 + V{3};
auto x85 = x49 + x84;
cell[0] = V{0.888888888888889}*x30*(x27 + x28 + x29) + cell[0];
cell[1] = x30*(-x32*x33 + x36) + cell[1];
cell[2] = x30*(-x38*x39 + x41) + cell[2];
cell[3] = x30*(-x43*x44 + x45) + cell[3];
cell[4] = -x30*(x47*x48 + x50*x51 + x53) + cell[4];
cell[5] = -x30*(-x48*x55 - x52 + V{0.0185185185185185}*x58*force[1]) + cell[5];
cell[6] = -x30*(x48*x60 + x61*x62 + x64) + cell[6];
cell[7] = -x30*(-x48*(x54 + x59) - x63 + V{0.0185185185185185}*x66*force[2]) + cell[7];
cell[8] = -x30*(x51*(x38 + x59) + x62*(x43 + x46) + x68) + cell[8];
cell[9] = -x30*(-x51*x71 - x67 + V{0.0185185185185185}*x73*force[2]) + cell[9];
cell[10] = -x74*((x46 + x61)*force[2] + (x47 + x59)*force[0] + (x50 + x59)*force[1]) + cell[10];
cell[11] = -x74*((x47 + x75)*force[0] + (x50 + x75)*force[1] - (x46 + x76 + x77)*force[2]) + cell[11];
cell[12] = -x74*((x60 + x72)*force[0] + (x61 + x72)*force[2] - (x59 + x69 + x77)*force[1]) + cell[12];
cell[13] = -x74*((x46 + x66)*force[2] - (x55 + x59)*force[0] + (x58 + x59)*force[1]) + cell[13];
cell[14] = x30*(-x33*x78 + x36) + cell[14];
cell[15] = x30*(-x39*x57 + x41) + cell[15];
cell[16] = x30*(-x44*x65 + x45) + cell[16];
cell[17] = -x30*(x48*x79 + x51*x80 + x53) + cell[17];
cell[18] = -x30*(-x51*x81 - x52 + V{0.0185185185185185}*x82*force[0]) + cell[18];
cell[19] = -x30*(x48*(x59 + x78) + x62*x83 + x64) + cell[19];
cell[20] = -x30*(-x62*x85 - x63 + V{0.0185185185185185}*(x75 + x78)*force[0]) + cell[20];
cell[21] = -x30*(x51*(x57 + x59) + x62*(x46 + x65) + x68) + cell[21];
cell[22] = -x30*(-x62*(x46 + x84) - x67 + V{0.0185185185185185}*(x57 + x75)*force[1]) + cell[22];
cell[23] = -x74*((x46 + x83)*force[2] + (x59 + x79)*force[0] + (x59 + x80)*force[1]) + cell[23];
cell[24] = -x74*(-(x46 + x85)*force[2] + (x75 + x79)*force[0] + (x75 + x80)*force[1]) + cell[24];
cell[25] = -x74*(-(x49 + x71)*force[1] + (x49 + x73)*force[2] + (x59 + x82)*force[0]) + cell[25];
cell[26] = -x74*(-(x72 + x85)*force[2] - (x75 + x81)*force[1] + (x75 + x82)*force[0]) + cell[26];

}

template <typename CELL, typename RHO, typename NABLARHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addLiangForce(CELL& cell, RHO& rho, NABLARHO& nablarho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x27 = rho*force[0];
auto x28 = nablarho[0]*u[0];
auto x29 = V{0.5}*omega + V{-1};
auto x30 = V{0.222222222222222}*x29;
auto x31 = rho*force[1];
auto x32 = nablarho[1]*u[1];
auto x33 = rho*force[2];
auto x34 = nablarho[2]*u[2];
auto x35 = u[0] + u[1];
auto x36 = x35*nablarho[0];
auto x37 = x35*nablarho[1];
auto x38 = x27 + x31;
auto x39 = V{0.0555555555555556}*x29;
auto x40 = u[0] - u[1];
auto x41 = x40*nablarho[1];
auto x42 = x27 - x31;
auto x43 = u[0] + u[2];
auto x44 = x43*nablarho[0];
auto x45 = x43*nablarho[2];
auto x46 = x27 + x33;
auto x47 = -u[2];
auto x48 = x47 + u[0];
auto x49 = x48*nablarho[2];
auto x50 = -x33;
auto x51 = x27 + x50;
auto x52 = u[1] + u[2];
auto x53 = x52*nablarho[1];
auto x54 = x52*nablarho[2];
auto x55 = x31 + x33;
auto x56 = x47 + u[1];
auto x57 = x56*nablarho[2];
auto x58 = x31 + x50;
auto x59 = x35 + u[2];
auto x60 = x59*nablarho[0];
auto x61 = x59*nablarho[1];
auto x62 = x59*nablarho[2];
auto x63 = x33 + x38;
auto x64 = V{0.0138888888888889}*x29;
auto x65 = x35 + x47;
auto x66 = x65*nablarho[2];
auto x67 = x38 + x50;
auto x68 = x40 + u[2];
auto x69 = x68*nablarho[1];
auto x70 = x33 + x42;
auto x71 = x52 - u[0];
auto x72 = -x71*nablarho[0] + x71*nablarho[1] + x71*nablarho[2];
cell[0] = cell[0];
cell[1] = x30*(x27 - x28) + cell[1];
cell[2] = x30*(x31 - x32) + cell[2];
cell[3] = x30*(x33 - x34) + cell[3];
cell[4] = x39*(-x36 - x37 + x38) + cell[4];
cell[5] = -x39*(x40*nablarho[0] - x41 - x42) + cell[5];
cell[6] = x39*(-x44 - x45 + x46) + cell[6];
cell[7] = -x39*(x48*nablarho[0] - x49 - x51) + cell[7];
cell[8] = x39*(-x53 - x54 + x55) + cell[8];
cell[9] = -x39*(x56*nablarho[1] - x57 - x58) + cell[9];
cell[10] = x64*(-x60 - x61 - x62 + x63) + cell[10];
cell[11] = -x64*(x65*nablarho[0] + x65*nablarho[1] - x66 - x67) + cell[11];
cell[12] = -x64*(x68*nablarho[0] + x68*nablarho[2] - x69 - x70) + cell[12];
cell[13] = -x64*(-x27 + x55 + x72) + cell[13];
cell[14] = -x30*(x27 + x28) + cell[14];
cell[15] = -x30*(x31 + x32) + cell[15];
cell[16] = -x30*(x33 + x34) + cell[16];
cell[17] = -x39*(x36 + x37 + x38) + cell[17];
cell[18] = -x39*(x40*nablarho[0] - x41 + x42) + cell[18];
cell[19] = -x39*(x44 + x45 + x46) + cell[19];
cell[20] = -x39*(x48*nablarho[0] - x49 + x51) + cell[20];
cell[21] = -x39*(x53 + x54 + x55) + cell[21];
cell[22] = -x39*(x56*nablarho[1] - x57 + x58) + cell[22];
cell[23] = -x64*(x60 + x61 + x62 + x63) + cell[23];
cell[24] = -x64*(x65*nablarho[0] + x65*nablarho[1] - x66 + x67) + cell[24];
cell[25] = -x64*(x68*nablarho[0] + x68*nablarho[2] - x69 + x70) + cell[25];
cell[26] = -x64*(x42 + x50 + x72) + cell[26];

}

template <typename CELL, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addAllenCahnForce(CELL& cell, OMEGA& omega, FORCE& force) any_platform
{
auto x27 = V{0.5}*omega + V{-1};
auto x28 = V{0.222222222222222}*x27;
auto x29 = x28*force[0];
auto x30 = x28*force[1];
auto x31 = x28*force[2];
auto x32 = force[0] + force[1];
auto x33 = V{0.0555555555555556}*x27;
auto x34 = x32*x33;
auto x35 = force[0] - force[1];
auto x36 = x33*x35;
auto x37 = x33*(force[0] + force[2]);
auto x38 = -force[2];
auto x39 = x33*(x38 + force[0]);
auto x40 = force[1] + force[2];
auto x41 = x33*x40;
auto x42 = x33*(x38 + force[1]);
auto x43 = V{0.0138888888888889}*x27;
auto x44 = x43*(x32 + force[2]);
auto x45 = x43*(x32 + x38);
auto x46 = x43*(x35 + force[2]);
auto x47 = x43*(x40 - force[0]);
cell[0] = cell[0];
cell[1] = x29 + cell[1];
cell[2] = x30 + cell[2];
cell[3] = x31 + cell[3];
cell[4] = x34 + cell[4];
cell[5] = x36 + cell[5];
cell[6] = x37 + cell[6];
cell[7] = x39 + cell[7];
cell[8] = x41 + cell[8];
cell[9] = x42 + cell[9];
cell[10] = x44 + cell[10];
cell[11] = x45 + cell[11];
cell[12] = x46 + cell[12];
cell[13] = -x47 + cell[13];
cell[14] = -x29 + cell[14];
cell[15] = -x30 + cell[15];
cell[16] = -x31 + cell[16];
cell[17] = -x34 + cell[17];
cell[18] = -x36 + cell[18];
cell[19] = -x37 + cell[19];
cell[20] = -x39 + cell[20];
cell[21] = -x41 + cell[21];
cell[22] = -x42 + cell[22];
cell[23] = -x44 + cell[23];
cell[24] = -x45 + cell[24];
cell[25] = -x46 + cell[25];
cell[26] = x47 + cell[26];

}

template <typename CELL, typename OMEGA, typename SOURCE, typename V=typename CELL::value_t>
static void addAllenCahnSource(CELL& cell, OMEGA& omega, SOURCE& source) any_platform
{
auto x27 = V{0.0740740740740741}*source;
auto x28 = V{0.0185185185185185}*source;
auto x29 = V{0.00462962962962963}*source;
cell[0] = V{0.296296296296296}*source + cell[0];
cell[1] = x27 + cell[1];
cell[2] = x27 + cell[2];
cell[3] = x27 + cell[3];
cell[4] = x28 + cell[4];
cell[5] = x28 + cell[5];
cell[6] = x28 + cell[6];
cell[7] = x28 + cell[7];
cell[8] = x28 + cell[8];
cell[9] = x28 + cell[9];
cell[10] = x29 + cell[10];
cell[11] = x29 + cell[11];
cell[12] = x29 + cell[12];
cell[13] = x29 + cell[13];
cell[14] = x27 + cell[14];
cell[15] = x27 + cell[15];
cell[16] = x27 + cell[16];
cell[17] = x28 + cell[17];
cell[18] = x28 + cell[18];
cell[19] = x28 + cell[19];
cell[20] = x28 + cell[20];
cell[21] = x28 + cell[21];
cell[22] = x28 + cell[22];
cell[23] = x29 + cell[23];
cell[24] = x29 + cell[24];
cell[25] = x29 + cell[25];
cell[26] = x29 + cell[26];

}

};

}

#endif

#endif
