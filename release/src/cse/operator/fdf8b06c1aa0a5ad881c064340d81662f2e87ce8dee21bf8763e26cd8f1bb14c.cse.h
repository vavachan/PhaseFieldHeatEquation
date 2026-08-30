/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-24 Adrian Kummerlaender, Shota Ito
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

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace operators {

template <typename T, typename... FIELDS>
struct CSE_O<StraightFdBoundaryProcessor2D<T, descriptors::D2Q9<FIELDS...>, 1, 1>,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x9 = cell.getDynamics().getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
V v_S0 = cell.computeRho();
auto x11 = v_S0/x9;
auto x12 = V{0.222222222222222}*x11;
auto x15 = V{0.0277777777777778}*x11;
auto x19 = V{0.0555555555555556}*x11;
auto x20 = V{0.111111111111111}*x11;
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_V1 [DESCRIPTOR::d]; cell.neighbor({1,0}).computeU(v_V1);
V v_V2 [DESCRIPTOR::d]; cell.neighbor({-1,0}).computeU(v_V2);
auto x10 = v_V1[0] - v_V2[0];
auto x16 = x10*x15;
V v_V3 [DESCRIPTOR::d]; cell.computeU(v_V3);
V v_V4 [DESCRIPTOR::d]; cell.neighbor({0,-1}).computeU(v_V4);
V v_V5 [DESCRIPTOR::d]; cell.neighbor({0,-2}).computeU(v_V5);
auto x14 = V{0.0833333333333333}*x11*(V{0.5}*v_V1[1] - V{0.5}*v_V2[1] + V{1.5}*v_V3[0] - V{2}*v_V4[0] + V{0.5}*v_V5[0]);
auto x13 = V{3}*v_V3[1] - V{4}*v_V4[1] + V{1}*v_V5[1];
auto x17 = x13*x15;
auto x18 = x14 - x16 - x17;
auto x21 = -x10*x20 + x13*x19;
auto x22 = x14 + x16 + x17;
auto x23 = x10*x19 - x13*x20;
V v_P0 [DESCRIPTOR::q]; cell.getDynamics().computeEquilibrium(cell, v_S0, v_V0, v_P0);
cell[0] = v_P0[0] + x10*x12 + x12*x13;
cell[1] = v_P0[1] + x18;
cell[2] = v_P0[2] + x21;
cell[3] = v_P0[3] - x22;
cell[4] = v_P0[4] + x23;
cell[5] = v_P0[5] + x18;
cell[6] = v_P0[6] + x21;
cell[7] = v_P0[7] - x22;
cell[8] = v_P0[8] + x23;
}
};

}

}
