/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 František Prinz, Nicolas Hafen, Mathias J. Krause
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


/**
 *  \file latticeStokesSpheroidDragForce.h
 *  \brief $(Drag Force and Torque in Stokes flow acting on subgrid rotational elipsoids (spheroid))
 *  Assumptions of linear shear flow - Stokes flow regime aproximated in the particle vicinity,
 *  fluid velocity considered only in the spheroid middle point
 *  and then Taylor expansion to the linear term is used. Thus, this force is working only for small particles or
 *  bigger particles in pure laminar flow.
 *  Force based on
 *  Brenner, H. 1963. The Stokes resistance of an arbitrary particle. Chem
 *  Eng Sci, 18: 1–25.
 *  Torque based on
 *  Jeffery, G. B. 1922. The motion of ellipsoidal particles immersed in a
 *  viscous fluid. P Roy Soc A, 102: 161–179.
 *
 *
 */

#ifndef LATTICE_STOKES_SPHEROID_DRAG_FORCE_H
#define LATTICE_STOKES_SPHEROID_DRAG_FORCE_H


namespace olb {


//Forward declaration
namespace particles{
template<typename T, typename PARTICLETYPE> class ParticleSystem;
template<typename T, typename PARTICLETYPE> class Particle;
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class BlockLatticeStokesSpheroidDragForce final : public BlockLatticePhysF<T,DESCRIPTOR> {
private:
  const BlockGeometry<T,DESCRIPTOR::d>& _blockGeometry;
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  particles::ParticleSystem<T,PARTICLETYPE>& _particleSystem;
  PhysR<T,DESCRIPTOR::d> _cellMin;
  PhysR<T,DESCRIPTOR::d> _cellMax;
  Vector<bool,DESCRIPTOR::d> _periodic;
  std::size_t _iP0;
public:
  BlockLatticeStokesSpheroidDragForce( BlockLattice<T,DESCRIPTOR>& blockLattice,
                               const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                               particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
                               const UnitConverter<T,DESCRIPTOR>& converter,
                               PhysR<T,DESCRIPTOR::d> cellMin = PhysR<T,DESCRIPTOR::d> (0.),
                               PhysR<T,DESCRIPTOR::d> cellMax = PhysR<T,DESCRIPTOR::d> (0.),
                               Vector<bool,DESCRIPTOR::d> periodic = Vector<bool,DESCRIPTOR::d> (false),
                               std::size_t iP0=0 );
  void evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP);
  bool operator() (T output[], const int input[]) override;
  static constexpr bool serializeForce = false;
};



}

#endif
