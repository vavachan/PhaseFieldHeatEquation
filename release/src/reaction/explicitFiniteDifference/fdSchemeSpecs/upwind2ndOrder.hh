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
 * Specializations for the 2nd-order upwind finite-difference differentiation scheme.
 *  -- generic implementation
 */
#ifndef FD_SCHEMES_UPWIND_2ND_HH
#define FD_SCHEMES_UPWIND_2ND_HH

namespace olb {

namespace fd {

template <unsigned D, typename T>
template <typename PARAMETERS>
T AdvectionScheme<D,T,tag::UPWIND_2_ORDER>::
apply(T& f0, T f[], T F[], T u[], PARAMETERS& params)
{
  T fNew = 0.;
  for (unsigned iD=0; iD<D; ++iD) {
    fNew += u[iD] * ( u[iD]==0. ? T(0.)
                      : ( u[iD] >0. ? 3.*f0 - 4.*f[getArrayPos<2>(0,iD)] + f[getArrayPos<2>(1,iD)]
                          : -F[getArrayPos<2>(1,iD)] + 4.*f[getArrayPos<2>(0,iD)] - 3.*f0 ) );
  }
  return 0.5 * fNew;
}

}  // namespace fd

}  // namespace olb

#endif
