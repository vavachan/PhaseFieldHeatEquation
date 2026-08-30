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

#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H

namespace olb {
namespace particles {
namespace contact {

/**Class storing properties that are necessary for the computation of the contact orce
 * N = number of different materials
 * [0]: modulus of elasticity
 * [1]: Poisson's ratio
 */
// TODO: Remove storage overhead (remove duplication of entries by using a one-dimensional array)
template <typename T, unsigned N>
struct MaterialProperties {
private:
  T materialProperties[N][2];

public:
  /// Constructor
  MaterialProperties() = default;

  /// Set material properties
  void set(const unsigned material, const T youngsModulus,
                             const T poissonRatio);

  /// Get modulus of elasticity
  const T getYoungsModulus(const unsigned material) const;
  /// Get Poisson's ratio
  const T getPoissonsRatio(const unsigned material) const;
};
} // namespace contact
} // namespace particles
} // namespace olb
#endif
