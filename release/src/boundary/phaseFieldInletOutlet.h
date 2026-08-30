/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

#ifndef BOUNDARY_PHASE_FIELD_INLET_OUTLET_H_
#define BOUNDARY_PHASE_FIELD_INLET_OUTLET_H_

namespace olb {

namespace boundary {

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = MultiPhaseIncompressbileTRTdynamics<T,DESCRIPTOR>
>
struct IncompressibleZouHeVelocity;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = MultiPhaseIncompressbileTRTdynamics<T,DESCRIPTOR>
>
struct IncompressibleZouHePressure;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = MultiPhaseIncompressbileTRTdynamics<T,DESCRIPTOR>
>
struct IncompressibleConvective;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = AllenCahnBGKdynamics<T,DESCRIPTOR>
>
struct PhaseFieldInlet;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = AllenCahnBGKdynamics<T,DESCRIPTOR>
>
struct PhaseFieldConvective;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = RLBdynamics<T, DESCRIPTOR>
>
struct FreeEnergyVelocity;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = RLBdynamics<T, DESCRIPTOR>
>
struct FreeEnergyPressure;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = RLBdynamics<T, DESCRIPTOR>
>
struct FreeEnergyPressureConvective;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR
>
struct FreeEnergyOrderParameter;

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR
>
struct FreeEnergyOrderParameterConvective;

}

}

#endif
