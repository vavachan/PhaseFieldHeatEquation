/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Davide Dapelo
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef FD_NO_PENETRATION_BOUNDARY_POST_PROCESSORS_3D_DEV03_H
#define FD_NO_PENETRATION_BOUNDARY_POST_PROCESSORS_3D_DEV03_H

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to treat no-penetration boundaries in the finite-difference advection-diffusion external field.
 */
template<typename T, typename DESCRIPTOR, typename MODEL, typename SCHEME_BOUND, typename PARAMS, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdBoundaryPostProcessor3D final : public FdBasePostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE> {
public:
  using parameters = PARAMS;
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  FdBoundaryPostProcessor3D();
  int getPriority() const;
  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& vars) any_platform;
};

}

#endif
