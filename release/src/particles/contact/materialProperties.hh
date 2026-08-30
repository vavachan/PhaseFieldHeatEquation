/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Jan E. Marquardt, Mathias J. Krause
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

#ifndef MATERIAL_PROPERTIES_HH
#define MATERIAL_PROPERTIES_HH

#include "materialProperties.h"

namespace olb {
namespace particles {
namespace contact {

template <typename T, unsigned N>
void MaterialProperties<T, N>::set(const unsigned material,
                                   const T youngsModulus, const T poissonRatio)
{
  materialProperties[material][0] = youngsModulus;
  materialProperties[material][1] = poissonRatio;
}

template <typename T, unsigned N>
const T
MaterialProperties<T, N>::getYoungsModulus(const unsigned material) const
{
  return materialProperties[material][0];
}

template <typename T, unsigned N>
const T
MaterialProperties<T, N>::getPoissonsRatio(const unsigned material) const
{
  return materialProperties[material][1];
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
