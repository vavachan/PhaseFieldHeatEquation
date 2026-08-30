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
 * 2D postprocessor to locally update the finite-diffusion advection-diffusion external field.
 *  -- header file
 */
#ifndef AD_POST_PROCESSOR_2D_DEV03_H
#define AD_POST_PROCESSOR_2D_DEV03_H

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Virtual parent to all the finite-difference postprocessors.
 * Its raison d'etre is to enforce all the finite-difference postprocessors
 * to possess applySourceTerm().
 */
template<typename T, typename DESCRIPTOR, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdBasePostProcessor2D {
protected:
  FdBasePostProcessor2D();
  template <typename CELL>
  void applySourceTerm(T* fNew, CELL& cell) any_platform;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to update the finite-difference advection-diffusion external field.
 */
template<typename T, typename DESCRIPTOR, typename MODEL, typename PARAMS, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
class FdPostProcessor2D : public FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE> {
public:
  using parameters = PARAMS;
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  FdPostProcessor2D();
  int getPriority() const;
  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& vars) any_platform;
};

}  // namespace olb

#endif
