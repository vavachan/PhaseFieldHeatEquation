/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Jan E. Marquardt
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

#ifndef PARTICLE_COMMUNICATOR_H
#define PARTICLE_COMMUNICATOR_H

#include "communication/mpiManager.h"

namespace olb {

namespace particles {

namespace communication {

struct ParticleCommunicator {
  ParticleCommunicator();
  ~ParticleCommunicator();

  ParticleCommunicator(const ParticleCommunicator&)            = delete;
  ParticleCommunicator& operator=(const ParticleCommunicator&) = delete;
  ParticleCommunicator(ParticleCommunicator&&)                 = delete;
  ParticleCommunicator& operator=(ParticleCommunicator&&)      = delete;

#ifdef PARALLEL_MODE_MPI
  MPI_Comm
      particleDistribution; /// Communicator for particle distribution communication
  MPI_Comm surfaceForceComm; /// Communicator for surface force communication
  MPI_Comm
      wallContactDetectionComm; /// Communicator for wall contact detection communication
  MPI_Comm
      particleContactDetectionComm; /// Communicator for particle contact detection communication
  MPI_Comm
      contactTreatmentComm; /// Communicator for contact treatment results communication
  MPI_Comm
      equationsOfMotionComm; /// Communicator for equations of motion results communication
#endif
};

} // namespace communication
} // namespace particles
} // namespace olb

#endif
