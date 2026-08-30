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

#ifndef WALL_H
#define WALL_H

#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "utilities/aliases.h"
#include "utilities/functorPtr.h"

namespace olb {

template <typename T, unsigned D>
struct SolidBoundary {
private:
  /// Indicator providing the surface and volume information
  std::unique_ptr<IndicatorF<T, D>> _indPtr;
  /// Material number on the lattice for identification and handling of boundaries
  const unsigned _latticeMaterial;
  /// Material number that corresponds to mechanical properties
  const unsigned _contactMaterial;
  /// Virtual enlargement during contact treatment
  const T _enlargementForContact;

public:
  /// Constructor
  constexpr SolidBoundary() = delete;
  /// Constructor
  constexpr SolidBoundary(std::unique_ptr<IndicatorF<T, D>> indPtr,
                          unsigned latticeMaterial, unsigned contactMaterial,
                          T enlargementForContact = T {0});

  IndicatorF<T, D>*  getIndicator();
  constexpr unsigned getLatticeMaterial() const;
  constexpr unsigned getContactMaterial() const;
  constexpr T        getEnlargementForContact() const;
};


/// Get material numbers of multiple solid boundaries in std::vector
/// - can be used to e.g. limit setBoundaryField() by material number
template<typename T, unsigned D>
std::set<int> getLatticeMaterials(
  const std::vector<SolidBoundary<T,D>>& solidBoundaries )
{
  std::set<int> materials;
  for (const SolidBoundary<T,D>& solidBoundary : solidBoundaries ){
    materials.insert( solidBoundary.getLatticeMaterial() );
  }
  return materials;
}


} // namespace olb

#endif
