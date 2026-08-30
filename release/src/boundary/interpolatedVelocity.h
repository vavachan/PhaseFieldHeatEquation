/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

#ifndef BOUNDARY_INTERPOLATED_VELOCITY_H_
#define BOUNDARY_INTERPOLATED_VELOCITY_H_

namespace olb {

namespace boundary {

/** Boundary condition details:
 *
 * Given name: Skordos
 *
 * References:
 * finite difference boundary method [Skordos, 1993, 10.1103/PhysRevE.48.4823].
 * finite difference velocity gradient method [Latt et al., 2008, doi: 10.1103/PhysRevE.77.056703]
 *
 * Equation: Navier-Stokes
 *
 * Type: onLattice
 *
 * Condition on first moment (velocity): Dirichlet type
 */
template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = BGKdynamics<T,DESCRIPTOR>
>
struct InterpolatedVelocity;

}

}

#endif
