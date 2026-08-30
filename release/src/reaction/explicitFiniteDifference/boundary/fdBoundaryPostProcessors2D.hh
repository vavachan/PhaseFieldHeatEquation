/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Davide Dapelo
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

#ifndef FD_NO_PENETRATION_BOUNDARIES_2D_DEV03_HH
#define FD_NO_PENETRATION_BOUNDARIES_2D_DEV03_HH

namespace olb {


////////  FdBoundaryPostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMS, typename FIELD, typename SOURCE>
FdBoundaryPostProcessor2D <T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMS,FIELD,SOURCE>::FdBoundaryPostProcessor2D()
  :  FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>()
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::NORMAL_X>() );
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::NORMAL_Y>() );
}

template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMS, typename FIELD, typename SOURCE>
int FdBoundaryPostProcessor2D <T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMS,FIELD,SOURCE>::getPriority() const
{
  return 1;
}

template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMS, typename FIELD, typename SOURCE>
template <typename CELL, typename PARAMETERS>
void FdBoundaryPostProcessor2D<T,DESCRIPTOR,MODEL,SCHEME_BOUND,PARAMS,FIELD,SOURCE>::apply(CELL& cell, PARAMETERS& vars)
{
  int normalX = (int)( cell.template getField<descriptors::NORMAL_X>() );
  int normalY = (int)( cell.template getField<descriptors::NORMAL_Y>() );
  if (normalX!=0 || normalY!=0) { // does absolutely nothing if the cell does not belong to a border
    constexpr int extraExtent = MODEL::extent() + SCHEME_BOUND::getExtraExtent();
    std::size_t iT = vars.template get<fd::fdParams::Timestep>();
    T f[MODEL::extent()*DESCRIPTOR::d], F[MODEL::extent()*DESCRIPTOR::d];
    T fGhost[MODEL::extent()*DESCRIPTOR::d], fNormal[extraExtent*DESCRIPTOR::d];
    T u[DESCRIPTOR::d];
    int normal[] {normalX, normalY};
    std::fill(std::begin(fGhost), std::end(fGhost), 0);
    T* fNew = fd::accessNew<T,FIELD>(cell, iT);
    T* f0   = fd::accessOld<T,FIELD>(cell, iT);
    for (int iN=1; iN<=extraExtent; ++iN) {
      fNormal[fd::getArrayPos<extraExtent>(iN-1,0)] = *fd::accessOld<T,FIELD>( cell.neighbor({normalX*iN,           0}), iT );
      fNormal[fd::getArrayPos<extraExtent>(iN-1,1)] = *fd::accessOld<T,FIELD>( cell.neighbor({          0, normalY*iN}), iT );
    }
    cell.computeU(u);
    SCHEME_BOUND::apply(fGhost, *f0, fNormal, normal, u, vars);
    for (int iN=1; iN<=MODEL::extent(); ++iN) {
      f[fd::getArrayPos<MODEL::extent()>(iN-1,0)] = normalX<=0 ? *fd::accessOld<T,FIELD>( cell.neighbor({-iN,   0}), iT )
                                                               : fGhost[fd::getArrayPos<MODEL::extent()+SCHEME_BOUND::getExtraExtent()>(iN-1,0)];
      F[fd::getArrayPos<MODEL::extent()>(iN-1,0)] = normalX>=0 ? *fd::accessOld<T,FIELD>( cell.neighbor({+iN,   0}), iT )
                                                               : fGhost[fd::getArrayPos<MODEL::extent()+SCHEME_BOUND::getExtraExtent()>(iN-1,0)];
      f[fd::getArrayPos<MODEL::extent()>(iN-1,1)] = normalY<=0 ? *fd::accessOld<T,FIELD>( cell.neighbor({  0, -iN}), iT )
                                                               : fGhost[fd::getArrayPos<MODEL::extent()+SCHEME_BOUND::getExtraExtent()>(iN-1,1)];
      F[fd::getArrayPos<MODEL::extent()>(iN-1,1)] = normalY>=0 ? *fd::accessOld<T,FIELD>( cell.neighbor({  0, +iN}), iT )
                                                               : fGhost[fd::getArrayPos<MODEL::extent()+SCHEME_BOUND::getExtraExtent()>(iN-1,1)];
    }
    MODEL::apply(fNew, f0, f, F, u, vars);
    this->applySourceTerm(fNew, cell);
  }
}


}  // namespace olb

#endif
