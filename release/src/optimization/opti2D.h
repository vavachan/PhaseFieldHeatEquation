/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
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
 * Groups all the 2D include files for the directory optimization.
 */

#include "optimization/core/latticeData.h"
#include "optimization/core/latticeResults.h"
#include "optimization/core/serialize.h"
#include "optimization/core/optiCase.h"
#include "optimization/core/optiCaseAD.h"
#include "optimization/core/optiCaseDual.h"
#include "optimization/core/optimizer.h"
#include "optimization/core/optimizerBarzilaiBorwein.h"
#include "optimization/core/optimizerLineSearch.h"
#include "optimization/core/optimizerLBFGS.h"
#include "optimization/core/optimizerSteepestDecent.h"
#include "optimization/core/projection.h"
#include "optimization/dynamics/dualDynamics.h"
#include "optimization/functors/dualFunctors3D.h"
#include "optimization/primitives/primitives.h"
#include "optimization/solver/objective.h"
#include "optimization/solver/optiSolverParameters.h"
#include "optimization/solver/serialization.h"
#include "utilities/aDiff.h"
