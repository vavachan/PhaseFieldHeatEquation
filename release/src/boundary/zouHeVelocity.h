/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Stephan Simonis
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

#ifndef BOUNDARY_ZOUHE_VELOCITY
#define BOUNDARY_ZOUHE_VELOCITY

namespace olb {

namespace boundary {

/** Boundary condition details:
 *
 * Given name: Zou-He
 *
 * References:
 * non-equilibrium bounce back method [Zou and He, 1997, doi: 10.1063/1.869307]
 *
 * Equation: Navier-Stokes
 *
 * Type: onLattice
 *
 * Condition on first moment (velocity): Dirichlet type
 *
 * This implementation is lattice independent.
 */
template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = BGKdynamics<T,DESCRIPTOR>
>
struct ZouHeVelocity;

}

}

#endif