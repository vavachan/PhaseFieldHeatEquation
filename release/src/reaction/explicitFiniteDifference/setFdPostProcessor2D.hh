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
 * postprocessor to set finite-diference postprocessors.
 *  -- generic implementation
 */
#ifndef SET_FD_POST_PROCESSOR_2D_DEV03_HH
#define SET_FD_POST_PROCESSOR_2D_DEV03_HH

namespace olb {

template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdPostProcessor2D(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,DESCRIPTOR::d>& superGeometry, int material)
{
  setFdPostProcessor2D<T,DESCRIPTOR,MODEL,PARAMETERS,FIELD,SOURCE>(sLattice, superGeometry.getMaterialIndicator(material));
}


template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdPostProcessor2D(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setFdPostProcessor2D<T,DESCRIPTOR,MODEL,PARAMETERS,FIELD,SOURCE>(sLattice.getBlock(iC), indicator->getBlockIndicatorF(iC));
  }
}


template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMETERS, typename FIELD, typename SOURCE>
void setFdPostProcessor2D(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF<T,DESCRIPTOR::d>& indicator)
{
  block.addPostProcessor(typeid(stage::PostStream), indicator, meta::id<FdPostProcessor2D<T,DESCRIPTOR,MODEL,PARAMETERS,FIELD,SOURCE>>{});
}

}  // namespace olb

#endif
