/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 the OpenLB project
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
 * Groups all the 3D include files in the boundaryConditions directory.
 */

#include "setBoundary.h"

#include "boundary/postprocessor/advectionDiffusionBoundaryPostProcessor3D.h"
#include "boundaryPostProcessors3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "offBoundaryPostProcessors3D.h"
#include "zouHeVelocity3D.h"
#include "localVelocity3D.h"
#include "phaseFieldWall3D.h"
#include "phaseFieldInletOutlet3D.h"
#include "interpolatedVelocity3D.h"
#include "interpolatedPressure3D.h"
#include "setZeroGradientBoundary3D.h"
#include "zeroDistribution.h"
#include "setBouzidiBoundary.h"
#include "bounceBack.h"
#include "bounceBackVelocity.h"
#include "legacy/defineU3D.h"
#include "slip3D.h"
#include "partialSlip3D.h"
#include "robin3D.h"
#include "helper.h"
#include "vortexMethod.h"
#include "interpolatedConvection3D.h"
#include "externalField.h"
#include "advectionDiffusionDirichlet3D.h"
#include "legacy/setWallFunctionBoundary3D.h"
#include "legacy/wallFunctionBoundaryPostProcessors3D.h"
#include "legacy/setBouzidiVelocityBoundary3D.h"
#include "legacy/setBouzidiZeroVelocityBoundary3D.h"
#include "advectionDiffusionDirichlet3D.h"
#include "zouHePressure.h"
#include "zouHePressure3D.h"
#include "setTurbulentWallModel.h"


#include "boundary3D.hh"
