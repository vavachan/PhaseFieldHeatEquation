/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerländer, Nicolas Hafen,
 *  Mathias J. Krause
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

#ifndef MPI_GROUP_HH
#define MPI_GROUP_HH

#include <stdexcept>

namespace olb {

// *INDENT-OFF*

#ifdef PARALLEL_MODE_MPI

MPI_Group_Wrapper::MPI_Group_Wrapper()
{
  if (MPI_Comm_dup(MPI_COMM_WORLD, &_commGroup) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
}

MPI_Group_Wrapper::~MPI_Group_Wrapper()
{
  MPI_Comm_free(&_commGroup);
}

MPI_Comm& MPI_Group_Wrapper::getComm(){
  return _commGroup;
}

#endif

// *INDENT-ON*

}

#endif
