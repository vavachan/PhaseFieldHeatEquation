/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Pavel Eichler, Louis Kronberg, Stephan Simonis
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

#ifndef CUM_DYNAMICS_H
#define CUM_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "collisionCUM.h"

namespace olb {

/**
 * Implementation partially based on:
 * Geier, Martin, et al. "The cumulant lattice Boltzmann equation in three dimensions: Theory and validation."
 * Computers & Mathematics with Applications 70.4 (2015): 507-547.
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using CUMdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::CUM
>;

}
#endif
