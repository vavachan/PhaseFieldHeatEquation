/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef FSI_FIELDS_H
#define FSI_FIELDS_H

#include "descriptor/fields.h"

namespace olb {

namespace fields {

namespace fsi {

// General fields shared by all elements

struct ELEMENT_TAG : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct ELEMENT_PIVOT  : public descriptors::FIELD_BASE<0,1> { };
struct REDUCED_ELEMENT_TAG : public ELEMENT_TAG { };

struct ELEMENT_F_CURR : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_F_PREV : public descriptors::FIELD_BASE<0,1> { };

struct ELEMENT_FORCE  : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_TORQUE : public descriptors::ROTATION_FIELD_BASE { };
struct ELEMENT_STRESS : public descriptors::FIELD_BASE<0,1> { };

struct ELEMENTS_COUNT : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };
struct REDUCED_ELEMENTS_COUNT : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };

// Element-specific fields

struct ELEMENT_U_TRANSLATION : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_U_ROTATION : public descriptors::ROTATION_FIELD_BASE { };

struct ELEMENT_LOWER : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_UPPER : public descriptors::FIELD_BASE<0,1> { };
struct ELEMENT_ROTATION  : public descriptors::ROTATION_FIELD_BASE { };

struct ELEMENT_ROTATION_MATRIX : public descriptors::FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return DESCRIPTOR::d*DESCRIPTOR::d;
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};

struct ELEMENT_REFERENCE_DELTA_X : public descriptors::FIELD_BASE<1> { };
struct ELEMENT_REFERENCE_EXTENT : public descriptors::TYPED_FIELD_BASE<std::size_t,0,1> { };
struct ELEMENT_REFERENCE_PROJECTION : public descriptors::TYPED_FIELD_BASE<std::size_t,0,1> { };
struct ELEMENT_REFERENCE_POROSITY : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };
struct ELEMENT_REFERENCE_Y1 : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };

}

}

struct ElementParameters {
  using parameters = meta::list<
    fields::fsi::ELEMENT_TAG,
    fields::fsi::ELEMENT_PIVOT,
    fields::fsi::ELEMENT_LOWER,
    fields::fsi::ELEMENT_ROTATION,
    fields::fsi::ELEMENT_REFERENCE_DELTA_X,
    fields::fsi::ELEMENT_REFERENCE_EXTENT,
    fields::fsi::ELEMENT_REFERENCE_POROSITY
  >;

  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteParametersD<T,DESCRIPTOR,PLATFORM,parameters>;
};

namespace descriptors {

namespace fsi {

/// Descriptor for common FSI element data
template <unsigned D>
struct ELEMENTS : public SPATIAL_DESCRIPTOR<
  D,
  fields::fsi::ELEMENT_TAG,
  fields::fsi::ELEMENT_PIVOT
> { };

/// Descriptor for common reduced properties of FSI coupled elements
template <unsigned D>
struct REDUCED_ELEMENTS : public SPATIAL_DESCRIPTOR<
 D,
 fields::fsi::ELEMENT_TAG,
 fields::fsi::ELEMENT_FORCE,
 fields::fsi::ELEMENT_TORQUE,
 fields::fsi::ELEMENT_STRESS
> { };

/// Descriptor for a FSI reference element lattice
template <unsigned D>
struct REFERENCE_POROSITY_ELEMENT : public SPATIAL_DESCRIPTOR<
  D,
  descriptors::POROSITY
> { };

/// Descriptor for a FSI reference element lattice with wall model
template <unsigned D>
struct REFERENCE_POROSITY_ELEMENT_WITH_WALL_MODEL : public SPATIAL_DESCRIPTOR<
  D,
  descriptors::POROSITY,
  descriptors::Y1
> { };

}

}

}

#endif
