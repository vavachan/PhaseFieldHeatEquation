/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Luiz Eduardo Czelusniak, Tim Bingert
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

#ifndef OLB_BOUNDARY_PHASE_FIELD_WALL_H_
#define OLB_BOUNDARY_PHASE_FIELD_WALL_H_

#include "dynamics/freeEnergyDynamics.h"
#include "dynamics/dynamics.h"

namespace olb {

namespace boundary {

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = BounceBackBulkDensityADE<T,DESCRIPTOR>
>
struct PhaseFieldWall;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = BounceBackBulkDensityADE<T,DESCRIPTOR>
>
struct PhaseFieldCurvedWall;

template <
concepts::BaseType T,
concepts::LatticeDescriptor DESCRIPTOR,
typename MixinDynamics = BounceBackBulkDensityWellBalanced<T,DESCRIPTOR>
>
struct WellBalancedWall;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = BounceBackBulkDensity<T,DESCRIPTOR>
>
struct FreeEnergyWallMomentum;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = FreeEnergyWallDynamics<T,DESCRIPTOR>
>
struct FreeEnergyWallOrderParameter;

/*template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = NoDynamics<T,DESCRIPTOR>
>
struct SignedDistanceBoundary; //TODO: wait for stage flexibility*/

}

}

#endif
