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
 * Specializations for the central finite-difference differentiation scheme.
 *  -- header file
 */
#ifndef FD_SCHEMES_CENTRAL_WITH_ANTIDIFFUSIVITY_H
#define FD_SCHEMES_CENTRAL_WITH_ANTIDIFFUSIVITY_H

namespace olb {

namespace fd {

namespace tag {

struct CENTRAL_WITH_ANTIDIFFUSIVITY final : FD_TAG {
  CENTRAL_WITH_ANTIDIFFUSIVITY() = delete;
  static constexpr int extent() {return 1;}
};

} // namespace tag

template <unsigned D, typename T>
struct DiffusionScheme<D,T,tag::CENTRAL_WITH_ANTIDIFFUSIVITY> final : FdScheme<tag::CENTRAL_WITH_ANTIDIFFUSIVITY> {
  DiffusionScheme() = delete;
  template <typename PARAMETERS>
  static T apply(T& f0, T f[], T F[], T u[], PARAMETERS& params);
};
}  // namespace fd

}  // namespace olb

#endif
