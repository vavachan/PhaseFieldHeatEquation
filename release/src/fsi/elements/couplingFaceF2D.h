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

#ifndef FSI_COUPLING_FACE_F2D_H
#define FSI_COUPLING_FACE_F2D_H

#include "fsi/fields.h"
#include "functors/primitive/sdf.h"
#include "utilities/geometricOperations.h"

namespace olb {

template <typename T>
bool isCouplingFaceInterior(Vector<T,2> p, Vector<T,2> a, Vector<T,2> b, Vector<T,2> c, T w) any_platform
{
  const bool isFrontAB = (a != b) && (((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) >= 0);
  const bool isFrontBC = (c != b) && (((c[0] - b[0]) * (p[1] - b[1]) - (c[1] - b[1]) * (p[0] - b[0])) >= 0);

  return !(isFrontAB || isFrontBC);
}

struct CouplingFaceF2D {
  struct NEIGHBORS : public descriptors::TYPED_FIELD_BASE<int,2> { };

  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_PIVOT,
    fields::fsi::ELEMENT_U_TRANSLATION,
    NEIGHBORS
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_TAG>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    fields::array_of<NEIGHBORS>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned currElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto tags = params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto neighbors = params.template get<fields::array_of<NEIGHBORS>>();

    const int prevElement = neighbors[0][currElement];
    const int nextElement = neighbors[1][currElement];

    const Vector<V,DESCRIPTOR::d> prevR{pivots, static_cast<unsigned>(prevElement)};
    const Vector<V,DESCRIPTOR::d> currR{pivots, static_cast<unsigned>(currElement)};
    const Vector<V,DESCRIPTOR::d> nextR{pivots, static_cast<unsigned>(nextElement)};

    const auto dPrev = norm(physR - prevR);
    const auto dCurr = norm(physR - currR);
    const auto dNext = norm(physR - nextR);

    return (dPrev <  dCurr && dPrev <  dNext) * tags[prevElement]
         + (dCurr <= dPrev && dCurr <= dNext) * tags[currElement]
         + (dNext <  dCurr && dNext <  dPrev) * tags[nextElement];
  }

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto neighbors = params.template get<fields::array_of<NEIGHBORS>>();

    const int prevElement = neighbors[0][iElement];
    const int nextElement = neighbors[1][iElement];

    const Vector<V,DESCRIPTOR::d> prevR{pivots, static_cast<unsigned>(prevElement)};
    const Vector<V,DESCRIPTOR::d> currR{pivots, static_cast<unsigned>(iElement)};
    const Vector<V,DESCRIPTOR::d> nextR{pivots, static_cast<unsigned>(nextElement)};

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    return isCouplingFaceInterior(physR, prevR, currR, nextR, physDeltaX);
  }

  template <concepts::Parameters PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto neighbors = params.template get<fields::array_of<NEIGHBORS>>();

    const int prevElement = neighbors[0][iElement];
    const int nextElement = neighbors[1][iElement];

    const Vector<V,DESCRIPTOR::d> prevR{pivots, static_cast<unsigned>(prevElement)};
    const Vector<V,DESCRIPTOR::d> currR{pivots, static_cast<unsigned>(iElement)};
    const Vector<V,DESCRIPTOR::d> nextR{pivots, static_cast<unsigned>(nextElement)};

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    const V d = util::min(sdf::line(physR, prevR, currR, physDeltaX),
                          sdf::line(physR, currR, nextR, physDeltaX));
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

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto neighbors = params.template get<fields::array_of<NEIGHBORS>>();
    const auto translationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>>();

    const auto prevElement = neighbors[0][iElement];
    const auto nextElement = neighbors[1][iElement];

    const Vector<V,DESCRIPTOR::d> prevR{pivots, static_cast<unsigned>(prevElement)};
    const Vector<V,DESCRIPTOR::d> currR{pivots, static_cast<unsigned>(iElement)};
    const Vector<V,DESCRIPTOR::d> nextR{pivots, static_cast<unsigned>(nextElement)};

    const V distPrev = norm(physR - prevR);
    const V distCurr = norm(physR - currR);
    const V distNext = norm(physR - nextR);

    const Vector<V,DESCRIPTOR::d> prevU{translationUs, static_cast<unsigned>(prevElement)};
    const Vector<V,DESCRIPTOR::d> currU{translationUs, static_cast<unsigned>(iElement)};
    const Vector<V,DESCRIPTOR::d> nextU{translationUs, static_cast<unsigned>(nextElement)};

    const auto totalPrev = distPrev + distCurr;
    const auto totalNext = distNext + distCurr;

    auto u = (distPrev <= distNext) * ((1 - distCurr/totalPrev) * currU + (distCurr/totalPrev) * prevU)
           + (distPrev >  distNext) * ((1 - distCurr/totalNext) * currU + (distCurr/totalNext) * nextU);
    u *= params.template get<fields::converter::PHYS_VELOCITY>();
    return Vector<V,DESCRIPTOR::d>{};
  }
};


}

#endif
