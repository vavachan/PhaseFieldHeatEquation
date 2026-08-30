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
 * Finite-difference models.
 *  -- header file
 */
#ifndef FD_MODEL_H
#define FD_MODEL_H

namespace olb {

/*
 * Class for finite-difference advection-diffusion model
 */
template<typename T, typename SCHEME_ADV, typename SCHEME_DIFF>
class FdAdvectionDiffusionModel {
public:
  FdAdvectionDiffusionModel() = delete;
  static constexpr int extent();
  template <typename PARAMETERS>
  static void apply(T* fNew, T* f0, T f[], T F[], T u[], PARAMETERS& params);
};

}  // namespace olb

#endif
