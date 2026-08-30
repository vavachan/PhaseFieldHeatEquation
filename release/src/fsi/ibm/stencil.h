/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Timm Krüger, Shota Ito
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

#ifndef STENCIL_H
#define STENCIL_H

namespace olb {

namespace ibm {

template <typename T, unsigned int Width>
struct Stencil
{
  static_assert(Width < 2 || Width > 4, "Invalid stencil size for IBM. Options: {2,3,4}");
};

template <typename T>
struct Stencil<T,2>
{
  static T weight(T& latticeDist) any_platform {
    if (latticeDist < 0) {
      latticeDist *= -1;
    }
    return (1. - latticeDist);
  }
};

template <typename T>
struct Stencil<T,3>
{
  static T weight(T& latticeDist) any_platform {
    if (latticeDist < 0) {
      latticeDist *= -1.;
    }
    if (2 * latticeDist > 1.) {
      return ((5. - 3. * latticeDist - util::sqrt(-2. + 6. * latticeDist - 3. * util::pow(latticeDist, 2.))) / 6.);
    } else {
      return ((1. + util::sqrt(1. - 3. * util::pow(latticeDist, 2.))) / 3.);
    };
  }
};

template <typename T>
struct Stencil<T,4>
{
  static T weight(T& latticeDist) any_platform {
    if (latticeDist < 0) {
      latticeDist *= -1.;
    }
    if (latticeDist >= 1.) {
      return ((5. - 2. * latticeDist - util::sqrt(-7. + 12. * latticeDist - 4. * util::pow(latticeDist, 2.))) / 8.);
    } else {
      return ((3. - 2. * latticeDist + util::sqrt(1. + 4. * latticeDist - 4. * util::pow(latticeDist, 2.))) / 8.);
    };
  }
};

template<int width>
using StencilWidth = std::integral_constant<int,width>;

}

}

#endif
