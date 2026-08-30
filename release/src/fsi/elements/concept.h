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

#ifndef FSI_ELEMENTS_CONCEPT_H
#define FSI_ELEMENTS_CONCEPT_H

#include "fsi/fields.h"
#include "core/meta.h"

namespace olb {

namespace concepts {

template <typename PorosityF>
concept BasePorosityElementF = requires() {
  // Declares which per-element data it requires
  requires std::is_base_of_v<meta::list_base, typename PorosityF::data>;
  // Declares which parameters it requires
  requires std::is_base_of_v<meta::list_base, typename PorosityF::parameters>;
};

template <typename PorosityF>
concept PorosityElementF2D = requires(placeholder::Parameters2D parameters,
                                      Vector<placeholder::Parameters2D::value_t,
                                             placeholder::Parameters2D::descriptor_t::d> physR,
                                      unsigned iElement) {
  requires BasePorosityElementF<PorosityF>;

  // Returns the element tag to be set at physR
  /**
   * Usually returns its own tag but may vary for neighborhood-aware
   * connected elements (e.g. coupling faces).
   **/
  { PorosityF().tag(parameters, physR, iElement) } -> std::same_as<int>;

  // Returns if physR is in the non-surface interior
  /**
   * May return false always. Introduced to disable momentum
   * exchange on the non-fluid side of the coupling face.
   **/
  { PorosityF().isInterior(parameters, physR, iElement) } -> std::same_as<bool>;

  // Returns porosity at physR
  { PorosityF().compute(parameters, physR, iElement) } -> std::same_as<typename placeholder::Parameters::value_t>;

};

template <typename PorosityF>
concept PorosityElementF3D = requires(placeholder::Parameters3D parameters,
                                      Vector<placeholder::Parameters3D::value_t,
                                             placeholder::Parameters3D::descriptor_t::d> physR,
                                      unsigned iElement) {
  requires BasePorosityElementF<PorosityF>;

  // Returns the element tag to be set at physR
  /**
   * Usually returns its own tag but may vary for neighborhood-aware
   * connected elements (e.g. coupling faces).
   **/
  { PorosityF().tag(parameters, physR, iElement) } -> std::same_as<int>;

  // Returns if physR is in the non-surface interior
  /**
   * May return false always. Introduced to disable momentum
   * exchange on the non-fluid side of the coupling face.
   **/
  { PorosityF().isInterior(parameters, physR, iElement) } -> std::same_as<bool>;

  // Returns porosity at physR
  { PorosityF().compute(parameters, physR, iElement) } -> std::same_as<typename placeholder::Parameters::value_t>;

};

template <typename PorosityF>
concept PorosityElementF = PorosityElementF2D<PorosityF> || PorosityElementF3D<PorosityF>;

template <typename PorosityF>
concept PorosityElementF2DcomputesU = requires(placeholder::Parameters2D parameters,
                                               Vector<placeholder::Parameters2D::value_t,
                                                      placeholder::Parameters2D::descriptor_t::d> physR,
                                               unsigned iElement) {
  // Returns solid velocity at physR
  { PorosityF().computeU(parameters, physR, iElement) } -> std::same_as<
    Vector<typename placeholder::Parameters::value_t,
           placeholder::Parameters::descriptor_t::d>
  >;
};

template <typename PorosityF>
concept PorosityElementF3DcomputesU = requires(placeholder::Parameters3D parameters,
                                               Vector<placeholder::Parameters3D::value_t,
                                                      placeholder::Parameters3D::descriptor_t::d> physR,
                                               unsigned iElement) {
  // Returns solid velocity at physR
  { PorosityF().computeU(parameters, physR, iElement) } -> std::same_as<
    Vector<typename placeholder::Parameters::value_t,
           placeholder::Parameters::descriptor_t::d>
  >;
};

template <typename PorosityF>
concept PorosityElementFcomputesU = requires() {
  requires PorosityElementF2DcomputesU<PorosityF> || PorosityElementF3DcomputesU<PorosityF>;
};

template <typename PorosityF>
concept PorosityElementF2DcomputesY1 = requires(placeholder::Parameters2D parameters,
                                                Vector<placeholder::Parameters2D::value_t,
                                                       placeholder::Parameters2D::descriptor_t::d> physR,
                                                unsigned iElement) {
  // Returns solid velocity at physR
  { PorosityF().computeY1(parameters, physR, iElement) } -> std::same_as<
    Vector<typename placeholder::Parameters::value_t,
           placeholder::Parameters::descriptor_t::d>
  >;
};

template <typename PorosityF>
concept PorosityElementF3DcomputesY1 = requires(placeholder::Parameters3D parameters,
                                                Vector<placeholder::Parameters3D::value_t,
                                                       placeholder::Parameters3D::descriptor_t::d> physR,
                                                unsigned iElement) {
  // Returns solid velocity at physR
  { PorosityF().computeY1(parameters, physR, iElement) } -> std::same_as<
    Vector<typename placeholder::Parameters::value_t,
           placeholder::Parameters::descriptor_t::d>
  >;
};

template <typename PorosityF>
concept PorosityElementFcomputesY1 = PorosityElementF2DcomputesY1<PorosityF>
                                  || PorosityElementF3DcomputesY1<PorosityF>;

}

}

#endif
