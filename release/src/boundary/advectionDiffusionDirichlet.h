/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef OLB_BOUNDARY_ADVECTION_DIFFUSION_DIRICHLET_H_
#define OLB_BOUNDARY_ADVECTION_DIFFUSION_DIRICHLET_H_

namespace olb {

namespace boundary {

/** Boundary condition details:
 *
 * Given name: Allen-Reis, Zou-He
 *
 * References:
 * D2Q5/D3Q7: moment-based BC [Allen and Reis, 2016, doi: 10.1504/PCFD.2016.077296].
 * D2Q9/D3Q19: non-equilibrium bounce back method [Krüger et al., 2017, doi: 10.1007/978-3-319-44649-3]
 *
 * Equation: advection-diffusion
 *
 * Type: onLattice
 *
 * Condition on zeroth moment (scalar): Dirichlet type
 */
template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR,
  typename MixinDynamics = AdvectionDiffusionRLBdynamics<T,DESCRIPTOR>
>
struct AdvectionDiffusionDirichlet;

}

}

#endif