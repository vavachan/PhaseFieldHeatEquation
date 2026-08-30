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

#ifndef SUPER_LATTICE_INTERACTION_HH
#define SUPER_LATTICE_INTERACTION_HH

#include "superLatticeInteraction.h"

namespace olb {

namespace particles {

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setSuperParticleField( const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            UnitConverter<T,DESCRIPTOR> const& converter,
                            Particle<T,PARTICLETYPE>& particle,
                            const Vector<bool,DESCRIPTOR::d>& periodicity )
{
  constexpr unsigned D = DESCRIPTOR::d;
  const PhysR<T,D> min = communication::getCuboidMin<T,D>(sGeometry.getCuboidDecomposition());
  const PhysR<T,D> max = communication::getCuboidMax<T,D>(sGeometry.getCuboidDecomposition(), min);


  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    if ( isPeriodic(periodicity) ) {
      setBlockParticleField( sGeometry.getBlockGeometry(iC),
                           sLattice.getBlock(iC), converter, min, max,
                           particle, periodicity);
    }
    else {
      setBlockParticleField( sGeometry.getBlockGeometry(iC),
                           sLattice.getBlock(iC), converter, particle);
    }
  }
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
  typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE, typename F>
void setSuperParticleField( const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                            const PhysR<T,DESCRIPTOR::d>& min,
                            const PhysR<T,DESCRIPTOR::d>& max,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            UnitConverter<T,DESCRIPTOR> const& converter,
                            ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& contactContainer,
                            const size_t iP,
                            Particle<T,PARTICLETYPE>& particle,
                            std::vector<SolidBoundary<T,DESCRIPTOR::d>>& solidBoundaries,
                            F getSetupPeriodicity,
                            int globaliC)
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    if(globaliC == sLattice.getLoadBalancer().glob(iC)
        || !particles::access::providesParallelization<PARTICLETYPE>()) {
      if constexpr ( isPeriodic(getSetupPeriodicity()) ) {
        setBlockParticleField( sGeometry.getBlockGeometry(iC),
                               sLattice.getBlock(iC), converter, min, max,
                               particleSystem, contactContainer, iP,
                               particle, solidBoundaries, getSetupPeriodicity);
      }
      else {
        setBlockParticleField( sGeometry.getBlockGeometry(iC),
            sLattice.getBlock(iC), converter, particleSystem, contactContainer,
            iP, particle, solidBoundaries);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void resetSuperParticleField( SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                              SuperLattice<T, DESCRIPTOR>& sLattice)
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    resetBlockParticleField( sGeometry.getBlockGeometry(iC), sLattice.getBlock(iC));
  }
}

//TODO: HOTFIX ONLY
template<typename T, typename DESCRIPTOR>
void resetContactField( SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                              SuperLattice<T, DESCRIPTOR>& sLattice)
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    resetBlockContactField( sGeometry.getBlockGeometry(iC), sLattice.getBlock(iC));
  }
}

} // namespace particles

} // namespace olb

#endif
