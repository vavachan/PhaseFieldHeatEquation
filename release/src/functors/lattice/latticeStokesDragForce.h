/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef LATTICE_STOKES_DRAG_FORCE_H
#define LATTICE_STOKES_DRAG_FORCE_H

#include <unordered_set>

namespace olb {


//Forward declaration
namespace particles{
template<typename T, typename PARTICLETYPE> class ParticleSystem;
template<typename T, typename PARTICLETYPE> class Particle;
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE, bool serialize=true>
class BlockLatticeStokesDragForce final : public BlockLatticePhysF<T,DESCRIPTOR> {
  using F = std::function<void(particles::Particle<T,PARTICLETYPE>&, const PhysR<T, DESCRIPTOR::d>&, const Vector<T,DESCRIPTOR::d>&, const Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation>&)>;

private:
  const BlockGeometry<T,DESCRIPTOR::d>& _blockGeometry;
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  particles::ParticleSystem<T,PARTICLETYPE>& _particleSystem;
  PhysR<T,DESCRIPTOR::d> _cellMin;
  PhysR<T,DESCRIPTOR::d> _cellMax;
  Vector<bool,DESCRIPTOR::d> _periodic;
  std::size_t _iP0;
  const std::set<int> _ignoredMaterials;
  const F _f;
  //Precalculated constants
  T _delTinv;
  T _C1;
public:
  BlockLatticeStokesDragForce( BlockLattice<T,DESCRIPTOR>& blockLattice,
                               const BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                               particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
                               const UnitConverter<T,DESCRIPTOR>& converter,
                               PhysR<T,DESCRIPTOR::d> cellMin = PhysR<T,DESCRIPTOR::d> (0.),
                               PhysR<T,DESCRIPTOR::d> cellMax = PhysR<T,DESCRIPTOR::d> (0.),
                               Vector<bool,DESCRIPTOR::d> periodic = Vector<bool,DESCRIPTOR::d> (false),
                               std::size_t iP0=0,
                               const std::set<int>& ignoredMaterials = std::set<int>{},
                               const F f = [](auto&, const auto&, const auto&, const auto&){}
                               );
  void evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP);
  bool operator() (T output[], const int input[]) override;
  static constexpr bool serializeForce = serialize; //Return serialized force values via output[]
};



}
#endif
