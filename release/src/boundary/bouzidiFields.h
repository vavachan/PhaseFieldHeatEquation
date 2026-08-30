/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Michael Crocoll, Adrian Kummerlaender, Shota Ito, Anas Selmi
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

#ifndef BOUZIDI_FIELDS_H
#define BOUZIDI_FIELDS_H

namespace olb {

namespace descriptors {

/// Interpolated Bounce Back (Bouzidi) distance field
struct BOUZIDI_DISTANCE : public descriptors::FIELD_BASE<0,0,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<BOUZIDI_DISTANCE>()>(-1);
  }
};

/// Interpolated Bounce Back (Bouzidi) velocity coefficient field
struct BOUZIDI_VELOCITY : public descriptors::FIELD_BASE<0,0,1> { };

/// Interpolated Bounce Back (Bouzidi) for ADE Dirichlet field
struct BOUZIDI_ADE_DIRICHLET : public descriptors::FIELD_BASE<0,0,1> { };

struct NORMAL : public descriptors::FIELD_BASE<0,1,0> { };

struct BOUZIDI_NORMAL : public FIELD_BASE_CUSTOM_SIZE {
  template <unsigned D, unsigned Q >
  static constexpr unsigned size() {
    return Q * D;
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<BOUZIDI_NORMAL>()>{};
  }
};

struct BOUZIDI_SAMPLING_DISTANCE : public FIELD_BASE_CUSTOM_SIZE {
  template <unsigned D, unsigned Q >
  static constexpr unsigned size() {
    return Q * D;
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<BOUZIDI_SAMPLING_DISTANCE>()>{};
  }
};

}

}

#endif
