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
 * Specializations for the upwind finite-difference differentiation scheme.
 *  -- generic implementation
 */
#ifndef FD_SCHEMES_UPWIND_HH
#define FD_SCHEMES_UPWIND_HH

namespace olb {

namespace fd {

template <unsigned D, typename T>
template <typename PARAMETERS>
T AdvectionScheme<D,T,tag::UPWIND>::
apply(T& f0, T f[], T F[], T u[], PARAMETERS& params)
{
  T fNew = 0.;
  for (unsigned iD=0; iD<D; ++iD) {
    fNew += u[iD] * ( u[iD]==0. ? T(0.)
                      : ( u[iD] >0. ? f0-f[getArrayPos<1>(0,iD)]
                          : F[getArrayPos<1>(0,iD)]-f0 ) );
  }
  return fNew;
}

template <unsigned D, typename T>
constexpr int AdNeumannZeroBoundaryScheme<D,T,tag::UPWIND>::getExtraExtent()
{
  return 0;
}

template <unsigned D, typename T>
template <typename PARAMETERS>
void AdNeumannZeroBoundaryScheme<D,T,tag::UPWIND>::
apply(T fOut[], T& f0, T fIn[], int normal[], T u[], PARAMETERS& params)
{
  for (unsigned iD=0; iD<D; ++iD) {
    fOut[getArrayPos<1>(0,iD)] = ( u[iD]*normal[iD]==0 ? fIn[getArrayPos<1>(0,iD)] : f0 );
  }
}

}  // namespace fd

}  // namespace olb

#endif
