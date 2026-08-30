/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2015 Mathias J. Krause
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

/** \file
 * Groups all the 3D include files in the complexDynamics directory.
*/

#include "advectionDiffusionDynamics.h"
#include "advectionDiffusionForces.h"
#include "advectionDiffusionReactionCouplingPostProcessor3D.h"
#include "dynamics.h"
#include "entropicDynamics.h"
#include "entropicLbHelpers.h"
#include "freeEnergyDynamics.h"
#include "freeEnergyCoupling3D.h"
#include "freeEnergyPostProcessor3D.h"
#include "freeSurfacePostProcessor3D.h"
#include "interactionPotential.h"
#include "momenta/aliases.h"
#include "mrtDynamics.h"
#include "navierStokesAdvectionDiffusionCoupling.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor3D.h"
#include "nonNewtonianBGKdynamics.h"
#include "shanChenForcedSingleComponentCoupling.h"
#include "phaseFieldCoupling.h"
#include "powerLawBGKdynamics.h"
#include "porousPowerLawBGKdynamics.h"
//#include "rtlbmDynamics.h"
#include "shanChenDynOmegaForcedPostProcessor3D.h"
#include "shanChenForcedPostProcessor.h"
#include "shanChenForcedSingleComponentPostProcessor3D.h"
#include "smagorinskyBGKdynamics.h"
#include "smagorinskyMRTdynamics.h"
#include "stochasticSGSdynamics.h"
#include "porousForcedBGKDynamics.h"
#include "kbcDynamics.h"
#include "cumulantDynamics.h"
#include "guoZhaoDynamics.h"
#include "thermalCreepBouzidiCoupling.h"
#include "multiphaseBGKdynamics.h"
#include "kEpsilonRANSDynamics.h"
#include "spongeLayerDynamics.h"

#include "dynamics3D.hh"

#include "fsi/dynamics.h"
