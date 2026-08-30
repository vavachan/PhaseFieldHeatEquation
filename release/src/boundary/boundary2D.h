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
 * Groups all the 2D include files in the boundaryConditions directory.
 */

#include "setBoundary.h"

#include "boundaryPostProcessors2D.h"
#include "extendedFiniteDifferenceBoundary2D.h"
#include "offBoundaryPostProcessors2D.h"
#include "legacy/defineU2D.h"
#include "zouHeVelocity2D.h"
#include "localVelocity2D.h"
#include "regularizedHeatFlux2D.h"
#include "phaseFieldWall2D.h"
#include "phaseFieldInletOutlet2D.h"
#include "regularizedTemperature2D.h"
#include "interpolatedVelocity2D.h"
#include "interpolatedPressure2D.h"
#include "setSignedDistanceBoundary2D.h"
#include "slip2D.h"
#include "partialSlip2D.h"
#include "helper.h"
#include "interpolatedConvection2D.h"
#include "setBouzidiBoundary.h"
#include "advectionDiffusionDirichlet2D.h"
#include "bounceBack.h"
#include "bounceBackVelocity.h"
#include "zouHePressure2D.h"
#include "setTurbulentWallModel.h"

#include "boundary2D.hh"
