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
 *  \file latticeStokesSpheroidLiftForce.h
 *  \brief $(Lift Force in Stokes flow acting on subgrid rotational elipsoids by Harper and Chang (1968))
 *  Force is shear-induced
 *  Assumptions of linear shear flow - Stokes flow regime aproximated in the particle vicinity,
 *  fluid velocity considered only in the spheroid middle point
 *  and then Taylor expansion to the linear term is used. Thus, this force is working only for small particles or
 *  bigger particles in pure laminar flow.
 *  Force based on
 *  Harper, E.Y., & Chang, I.D. (1968). Maximum dissipating resulting from lift in a slow viscous flow. Journal of Fluid Mechanics, 33, 209–225.
 *  impolemented from
 *  Jeffery, G. B. 1922. The motion of ellipsoidal particles immersed in a
 *  viscous fluid. P Roy Soc A, 102: 161–179.
 *  Tian, L., Ahmadi, G., Wang, Z., Hopke, P.K., 2012. Transport and deposition of ellipsoidal fibers in low Reynolds number flows. Journal of Aerosol Science 45, 1–18.
    https://doi.org/10.1016/j.jaerosci.2011.09.001
 *  and
 *  Shachar-Berman, L., Ostrovski, Y., Koshiyama, K., Wada, S., Kassinos, S.C., Sznitman, J., 2019. Targeting inhaled fibers to the pulmonary acinus:
 *  Opportunities for augmented delivery from in silico simulations. European Journal of Pharmaceutical Sciences 137, 105003. https://doi.org/10.1016/j.ejps.2019.105003
 *  Supplementary material 1
 *
 *
 */

#ifndef LATTICE_SPHEROID_Lift_FORCE_H
#define LATTICE_SPHEROID_Lift_FORCE_H


namespace olb {


//Forward declaration
namespace particles{
template<typename T, typename PARTICLETYPE> class ParticleSystem;
template<typename T, typename PARTICLETYPE> class Particle;
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class BlockLatticeSpheroidLiftForce final : public BlockLatticePhysF<T,DESCRIPTOR> {
private:
  const BlockGeometry<T,DESCRIPTOR::d>& _blockGeometry;
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  particles::ParticleSystem<T,PARTICLETYPE>& _particleSystem;
  PhysR<T,DESCRIPTOR::d> _cellMin;
  PhysR<T,DESCRIPTOR::d> _cellMax;
  Vector<bool,DESCRIPTOR::d> _periodic;
  std::size_t _iP0;
public:
  BlockLatticeSpheroidLiftForce( BlockLattice<T,DESCRIPTOR>& blockLattice,
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
