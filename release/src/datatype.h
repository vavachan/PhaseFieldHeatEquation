/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Julius Jeßberger
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
 * Sets the global type variable FLOATING_POINT_TYPE which is configured in
 * config.mk.
 */

#ifndef DATATYPE_H
#define DATATYPE_H

#ifdef DEFAULT_FLOATING_POINT_TYPE
using ADF_DOUBLE_1 = olb::util::ADf<double,1>;
using FLOATING_POINT_TYPE = DEFAULT_FLOATING_POINT_TYPE;
#else
using FLOATING_POINT_TYPE = double;
#endif

#endif
