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

#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

namespace olb {

namespace fd {

// -------------------- First derivatives -------------------------------------

/// Second-order central gradient (u_p1 = u(x+1))
template<typename T>
constexpr T centralGradient(T u_p1, T u_m1) any_platform
{
  return (u_p1 - u_m1) / T{2};
}

/// Second-order asymmetric gradient (u_1 = u(x+1))
template<typename T>
constexpr T boundaryGradient(T u_0, T u_1, T u_2) any_platform
{
  return (-T{3}*u_0 + T{4}*u_1 - T{1}*u_2) / T{2};
}

/// Forward second-order first derivative
template<typename T>
constexpr T FSGradient(T u_0, T u_1, T u_2)
{
  return boundaryGradient(u_0, u_1, u_2);
}

/// Backward second-order first derivative
template<typename T>
constexpr T BSGradient(T u_0, T u_1, T u_2)
{
  return T{-1} * boundaryGradient(u_0, u_1, u_2);
}

/// Value at u_0 for which asymmetric gradient is zero (u_1 = u(x+1))
template<typename T>
constexpr T boundaryZeroGradient(T u_1, T u_2)
{
  return T{4}/T{3}*u_1 - T{1}/T{3}*u_2;
}


// -------------------- Second derivatives ------------------------------------
// cf. Abramowitz, Stegun: Handbook of Mathematical Functions, 1964, p. 884

/// Second order central second derivative (u_p1 = u(x+1))
template<typename T>
constexpr T centralSecondDeriv(T u_m1, T u_0, T u_p1) any_platform
{
  return u_m1 - T{2}*u_0 + u_p1;
}

/// Forth order central second derivative (u_p1 = u(x+1))
template<typename T>
constexpr T centralSecondDeriv(T u_m2, T u_m1, T u_0, T u_p1, T u_p2) any_platform
{
  return (- u_m2 + T{16}*u_m1 - T{30}*u_0 + T{16}*u_p1 - u_p2) / T{12};
}


// -------------------- Interpolation -----------------------------------------


/// Linear interpolation (yields u0 at pos=0 and u1 at pos=1)
template<typename T>
constexpr T linearInterpolate(T u_0, T u_1, T pos)
{
  return (T{1}-pos)*u_0 + pos*u_1;
}

}  // namespace fd

}  // namespace olb


#endif
