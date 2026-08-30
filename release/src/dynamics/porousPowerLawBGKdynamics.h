/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcek, Davide Dapelo
 *                2022 Nando Suntoyo, Adrian Kummerlaender
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

#ifndef POROUS_POWER_LAW_BGK_DYNAMICS_H
#define POROUS_POWER_LAW_BGK_DYNAMICS_H

#include "powerLawBGKdynamics.h"
#include "porousBGKdynamics.h"
#include "porousForcedBGKDynamics.h"

namespace olb {

/// BGK collision using Power Law collision frequency
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticlePowerLawBGKdynamics = dynamics::Tuple<
    T, DESCRIPTOR,
    MOMENTA,
    equilibria::SecondOrder,
    powerlaw::OmegaFromCell<collision::PorousParticle<collision::BGK,false>,false>
    >;

/// BGK collision using Power Law collision frequency with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticlePowerLawForcedBGKdynamics = dynamics::Tuple<
    T, DESCRIPTOR,
    MOMENTA,
    equilibria::SecondOrder,
    powerlaw::OmegaFromCell<collision::BGK,false>,
    forcing::PorousParticleKupershtokh<false>
    >;

/// BGK collision using Power Law (Herschel Bulkley) collision frequency
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticlePowerLawHerschelBulkleyBGKdynamics = dynamics::Tuple<
    T, DESCRIPTOR,
    MOMENTA,
    equilibria::SecondOrder,
    powerlaw::OmegaFromCell<collision::PorousParticle<collision::BGK,false>,true>
    >;

/// BGK collision using Power Law (Herschel Bulkley) collision frequency with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticlePowerLawHerschelBulkleyForcedBGKdynamics = dynamics::Tuple<
    T, DESCRIPTOR,
    MOMENTA,
    equilibria::SecondOrder,
    powerlaw::OmegaFromCell<collision::BGK,true>,
    forcing::PorousParticleKupershtokh<false>
    >;

}

#endif // namespace olb
