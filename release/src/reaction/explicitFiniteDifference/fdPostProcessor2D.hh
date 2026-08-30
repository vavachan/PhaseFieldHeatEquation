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
 * 2D postprocessor to locally update the finite-diffusion advection-diffusion external field.
 *  -- generic implementation
 */
#ifndef AD_POST_PROCESSOR_2D_DEV03_HH
#define AD_POST_PROCESSOR_2D_DEV03_HH

namespace olb {


////////  FdBasePostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::FdBasePostProcessor2D()
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
template <typename CELL>
void FdBasePostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::applySourceTerm(T* fNew, CELL& cell)
{
  if constexpr (! std::is_void<SOURCE>::value) {
    *fNew += cell.template getFieldPointer<SOURCE>()[0];
  }
}


////////  FdPostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMS, typename FIELD, typename SOURCE>
FdPostProcessor2D <T,DESCRIPTOR,MODEL,PARAMS,FIELD,SOURCE>::FdPostProcessor2D()
  :  FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>()
{
}

template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMS, typename FIELD, typename SOURCE>
int FdPostProcessor2D <T,DESCRIPTOR,MODEL,PARAMS,FIELD,SOURCE>::getPriority() const
{
  return 1;
}

template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMS, typename FIELD, typename SOURCE>
template <typename CELL, typename PARAMETERS>
void FdPostProcessor2D<T,DESCRIPTOR,MODEL,PARAMS,FIELD,SOURCE>::apply(CELL& cell, PARAMETERS& vars)
{
  std::size_t iT = vars.template get<fd::fdParams::Timestep>();
  T f[MODEL::extent()*DESCRIPTOR::d], F[MODEL::extent()*DESCRIPTOR::d];
  T u[DESCRIPTOR::d];
  T* fNew = fd::accessNew<T,FIELD>(cell, iT);
  T* f0   = fd::accessOld<T,FIELD>(cell, iT);
  for (int iN=1; iN<=MODEL::extent(); ++iN) {
    f[fd::getArrayPos<MODEL::extent()>(iN-1,0)] = *fd::accessOld<T,FIELD>( cell.neighbor({-iN,   0}), iT );
    F[fd::getArrayPos<MODEL::extent()>(iN-1,0)] = *fd::accessOld<T,FIELD>( cell.neighbor({+iN,   0}), iT );
    f[fd::getArrayPos<MODEL::extent()>(iN-1,1)] = *fd::accessOld<T,FIELD>( cell.neighbor({  0, -iN}), iT );
    F[fd::getArrayPos<MODEL::extent()>(iN-1,1)] = *fd::accessOld<T,FIELD>( cell.neighbor({  0, +iN}), iT );
  }
  cell.computeU(u);
  MODEL::template apply<PARAMETERS>(fNew, f0, f, F, u, vars);
  this->applySourceTerm(fNew, cell);
}

}  // namespace olb

#endif
