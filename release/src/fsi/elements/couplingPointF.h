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

#ifndef FSI_COUPLING_POINT_F_H
#define FSI_COUPLING_POINT_F_H

#include "fsi/fields.h"

#include "utilities/geometricOperations.h"

namespace olb {

struct CouplingPointF {
  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_PIVOT
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_TAG>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned currElement) any_platform {
    const auto tags = params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>();
    return tags[currElement];
  }

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const Vector<V,DESCRIPTOR::d> currR{pivots, static_cast<unsigned>(iElement)};

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    const V d = sdf::sphere(physR - currR, physDeltaX);
    const V eps = V{0.5}*physDeltaX;
    if (d < eps) {
      return util::max(d / eps, 0);
    } else {
      return V{1};
    }
  }

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    Vector<V,DESCRIPTOR::d> u{};
    return u;
  }
};


}

#endif
