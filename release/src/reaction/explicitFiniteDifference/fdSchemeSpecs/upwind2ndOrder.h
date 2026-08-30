/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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
 * Specializations for the 2nd-order upwind finite-difference differentiation scheme.
 *  -- header file
 */
#ifndef FD_SCHEMES_UPWIND_2ND_H
#define FD_SCHEMES_UPWIND_2ND_H

namespace olb {

namespace fd {

namespace tag {

struct UPWIND_2_ORDER final : FD_TAG {
  UPWIND_2_ORDER() = delete;
  static constexpr int extent() {return 2;}
};

} // namespace tag

template <unsigned D, typename T>
struct AdvectionScheme<D,T,tag::UPWIND_2_ORDER> final : FdScheme<tag::UPWIND_2_ORDER> {
  AdvectionScheme() = delete;
  template <typename PARAMETERS>
  static T apply(T& f0, T f[], T F[], T u[], PARAMETERS& params);
};

}  // namespace fd

}  // namespace olb

#endif
