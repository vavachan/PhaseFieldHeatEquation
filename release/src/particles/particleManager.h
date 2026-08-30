/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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


/* The particle manager is intended to provide a semi-generic wrapper for dynamic functions in
 * struct form. Its main purpose is the reduction of loops over all particles. Fitting structs
 * can both be administered by the particle manager or by calling them directly it, if desired.
*/

#include "communication/particleCommunicator.h"

#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H

namespace olb {

namespace particles {

namespace dynamics {

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class ParticleManager{
private:
  XParticleSystem<T,PARTICLETYPE>& _xParticleSystem;
  SuperGeometry<T,DESCRIPTOR::d>& _sGeometry;
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  UnitConverter<T,DESCRIPTOR> const& _converter;
  Vector<T,PARTICLETYPE::d> _externalAcceleration;
  Vector<bool,PARTICLETYPE::d> _periodic;
  communication::ParticleCommunicator _communicator = communication::ParticleCommunicator();

  //Condition for TASKs requiring no particle loop
  template <typename TASK>
  using requires_no_loop = std::integral_constant<bool, !TASK::particleLoop>;

  //Unpack tasks requiring loop
  template<typename taskList, typename ISEQ>
  void unpackTasksLooped(Particle<T,PARTICLETYPE>& particle, T timeStepSize, ISEQ seq, int globiC=0);

public:
  //Constructor
  ParticleManager(
    XParticleSystem<T,PARTICLETYPE>& xParticleSystem,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<T,PARTICLETYPE::d> externalAcceleration = Vector<T,PARTICLETYPE::d>(0.),
    Vector<bool,PARTICLETYPE::d> periodic = Vector<bool,PARTICLETYPE::d>(false) );

  //Execute task and pass time step
  template<typename ...TASKLIST>
  void execute(T timeStepSize);

  //Execute task and use timestep provided by the converter
  template<typename ...TASKLIST>
  void execute();

  // Getter for particle communicator
  const communication::ParticleCommunicator& getParticleCommunicator();

  /// Set external acceleration
  void setExternalAcceleration(const Vector<T,PARTICLETYPE::d>& externalAcceleration);
};


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
