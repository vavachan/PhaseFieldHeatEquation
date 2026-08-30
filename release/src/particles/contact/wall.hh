/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt, Mathias J. Krause
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

#ifndef WALL_HH
#define WALL_HH

#include "wall.h"

namespace olb {

template <typename T, unsigned D>
constexpr SolidBoundary<T, D>::SolidBoundary(
    std::unique_ptr<IndicatorF<T, D>> indPtr, unsigned latticeMaterial,
    unsigned contactMaterial, T enlargementForContact)
    : _indPtr(std::move(indPtr))
    , _latticeMaterial(latticeMaterial)
    , _contactMaterial(contactMaterial)
    , _enlargementForContact(enlargementForContact)
{}

template <typename T, unsigned D>
IndicatorF<T, D>* SolidBoundary<T, D>::getIndicator()
{
  return _indPtr.get();
}

template <typename T, unsigned D>
constexpr unsigned SolidBoundary<T, D>::getLatticeMaterial() const
{
  return _latticeMaterial;
}

template <typename T, unsigned D>
constexpr unsigned SolidBoundary<T, D>::getContactMaterial() const
{
  return _contactMaterial;
}

template <typename T, unsigned D>
constexpr T SolidBoundary<T, D>::getEnlargementForContact() const
{
  return _enlargementForContact;
}

} // namespace olb

#endif
