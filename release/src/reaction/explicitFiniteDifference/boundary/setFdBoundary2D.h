/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, Davide Dapelo
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

//This file contains the Interpolated Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_FD_BOUNDARY_2D_DEV03_H
#define SET_FD_BOUNDARY_2D_DEV03_H

namespace olb {
///Initialising the setFdBoundary2D function on the superLattice domain
///Interpolated Boundaries use the BGKdynamics collision-operator
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdBoundary2D(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material);

///Initialising the setFdBoundary2D function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdBoundary2D(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator);

/// Set interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMETERS, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdBoundary2D(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);

}//namespace olb
#endif
