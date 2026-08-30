/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef OLB_BOUNDARY_PARTIAL_SLIP
#define OLB_BOUNDARY_PARTIAL_SLIP

namespace olb {

namespace descriptors {

struct TUNER : public descriptors::FIELD_BASE<1,0,0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<TUNER>()>(1.);
  }
};
}

namespace boundary {

template <
  concepts::BaseType T,
  concepts::LatticeDescriptor DESCRIPTOR
>
struct PartialSlip;

}

}

#endif
