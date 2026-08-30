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


#ifndef SUPER_LATTICE_INTERACTION_H
#define SUPER_LATTICE_INTERACTION_H

#include "core/core.h"
#include "functors/analytical/analyticalBaseF.h"
#include "geometry/superGeometry.h"
#include "particles/particles.h"

// All OpenLB code is contained in this namespace.
namespace olb {

namespace particles {

/// Set particle field with peridic support
//TODO: remove material 1 from hardcoded version
//TODO: also remove additional material 1 check in set setBlockParticleField()
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setSuperParticleField( const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            UnitConverter<T,DESCRIPTOR> const& converter,
                            Particle<T,PARTICLETYPE>& particle,
                            const Vector<bool,DESCRIPTOR::d>& periodicity );


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
  typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE, typename F>
void setSuperParticleField( const SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                            const PhysR<T,DESCRIPTOR::d>& min,
                            const PhysR<T,DESCRIPTOR::d>& max,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            UnitConverter<T,DESCRIPTOR> const& converter,
                            ParticleSystem<T,PARTICLETYPE>& particleSystem,
                            contact::ContactContainer<T,PARTICLECONTACTTYPE,WALLCONTACTTYPE>& particleContacts,
                            size_t iP,
                            Particle<T,PARTICLETYPE>& particle,
                            std::vector<SolidBoundary<T,DESCRIPTOR::d>>& solidBoundaries,
                            F getSetupPeriodicity,
                            int globiC = -1 );


/// Reset particle field
template<typename T, typename DESCRIPTOR>
void resetSuperParticleField( SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                              SuperLattice<T, DESCRIPTOR>& sLattice);

//TODO: HOTFIX ONLY
template<typename T, typename DESCRIPTOR>
void resetContactField( SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                              SuperLattice<T, DESCRIPTOR>& sLattice);

} // namespace particles

} // namespace olb

#endif
