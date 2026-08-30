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
 * Finite-difference models.
 *  -- generic implementation
 */
#ifndef FD_MODEL_HH
#define FD_MODEL_HH

namespace olb {


template<typename T, typename SCHEME_ADV, typename SCHEME_DIFF>
constexpr int FdAdvectionDiffusionModel<T,SCHEME_ADV,SCHEME_DIFF>::extent()
{
  return util::max<int>(SCHEME_ADV::extent(), SCHEME_DIFF::extent());
}

template<typename T, typename SCHEME_ADV, typename SCHEME_DIFF>
template <typename PARAMETERS>
void FdAdvectionDiffusionModel<T,SCHEME_ADV,SCHEME_DIFF>::
apply(T* fNew, T* f0, T f[], T F[], T u[], PARAMETERS& params)
{
  *fNew = *f0
        -  SCHEME_ADV::template apply<PARAMETERS>(*f0, f, F, u, params)
        + SCHEME_DIFF::template apply<PARAMETERS>(*f0, f, F, u, params);
}


}  // namespace olb

#endif
