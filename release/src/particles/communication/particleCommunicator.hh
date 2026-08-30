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

#ifndef PARTICLE_COMMUNICATOR_HH
#define PARTICLE_COMMUNICATOR_HH

#include "particleCommunicator.h"

namespace olb {
namespace particles {
namespace communication {

ParticleCommunicator::ParticleCommunicator()
{
#ifdef PARALLEL_MODE_MPI
  if (MPI_Comm_dup(MPI_COMM_WORLD, &particleDistribution) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &surfaceForceComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &wallContactDetectionComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &particleContactDetectionComm) !=
      MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &contactTreatmentComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &equationsOfMotionComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
#endif
}

ParticleCommunicator::~ParticleCommunicator()
{
#ifdef PARALLEL_MODE_MPI
  MPI_Comm_free(&particleDistribution);
  MPI_Comm_free(&surfaceForceComm);
  MPI_Comm_free(&wallContactDetectionComm);
  MPI_Comm_free(&particleContactDetectionComm);
  MPI_Comm_free(&contactTreatmentComm);
  MPI_Comm_free(&equationsOfMotionComm);
#endif
}

} // namespace communication
} // namespace particles
} // namespace olb

#endif
