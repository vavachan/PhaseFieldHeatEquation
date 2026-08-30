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
struct CSE_O<OuterVelocityCornerProcessor3D<T, descriptors::D3Q19<FIELDS...>, -1, 1, 1>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x19 = cell.getDynamics().getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
V v_S0 = cell.neighbor({1,0,0}).computeRho();
V v_S1 = cell.neighbor({0,-1,0}).computeRho();
V v_S2 = cell.neighbor({0,0,-1}).computeRho();
V v_S3 = cell.neighbor({2,0,0}).computeRho();
V v_S4 = cell.neighbor({0,-2,0}).computeRho();
V v_S5 = cell.neighbor({0,0,-2}).computeRho();
auto x21 = (V{0.148148148148148}*v_S0 + V{0.148148148148148}*v_S1 + V{0.148148148148148}*v_S2 - V{0.037037037037037}*v_S3 - V{0.037037037037037}*v_S4 - V{0.037037037037037}*v_S5)/x19;
auto x22 = V{0.5}*x21;
auto x25 = V{0.166666666666667}*x21;
auto x26 = V{0.0833333333333333}*x21;
auto x33 = V{0.25}*x21;
auto x35 = V{0.0833333333333333}*x21;
auto x37 = V{0.0416666666666667}*x21;
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_V1 [DESCRIPTOR::d]; cell.neighbor({1,0,0}).computeU(v_V1);
V v_V2 [DESCRIPTOR::d]; cell.neighbor({2,0,0}).computeU(v_V2);
auto x24 = V{3}*v_V0[0] - V{4}*v_V1[0] + V{1}*v_V2[0];
auto x30 = -x24*x26;
auto x36 = x24*x35;
V v_V3 [DESCRIPTOR::d]; cell.computeU(v_V3);
V v_V4 [DESCRIPTOR::d]; cell.neighbor({0,-1,0}).computeU(v_V4);
V v_V5 [DESCRIPTOR::d]; cell.neighbor({0,-2,0}).computeU(v_V5);
auto x34 = x33*(V{1.5}*v_V0[1] - V{2}*v_V1[1] + V{0.5}*v_V2[1] - V{1.5}*v_V3[0] + V{2}*v_V4[0] - V{0.5}*v_V5[0]);
auto x20 = V{3}*v_V3[1] - V{4}*v_V4[1] + V{1}*v_V5[1];
auto x27 = x20*x26;
auto x38 = x20*x35;
V v_V6 [DESCRIPTOR::d]; cell.computeU(v_V6);
V v_V7 [DESCRIPTOR::d]; cell.neighbor({0,0,-1}).computeU(v_V7);
V v_V8 [DESCRIPTOR::d]; cell.neighbor({0,0,-2}).computeU(v_V8);
auto x42 = x33*(V{1.5}*v_V0[2] - V{2}*v_V1[2] + V{0.5}*v_V2[2] - V{1.5}*v_V6[0] + V{2}*v_V7[0] - V{0.5}*v_V8[0]);
auto x47 = x33*(V{1.5}*v_V3[2] - V{2}*v_V4[2] + V{0.5}*v_V5[2] + V{1.5}*v_V6[1] - V{2}*v_V7[1] + V{0.5}*v_V8[1]);
auto x23 = V{3}*v_V6[2] - V{4}*v_V7[2] + V{1}*v_V8[2];
auto x28 = x23*x26;
auto x29 = x24*x25 + x27 + x28;
auto x31 = -x20*x25 + x28 + x30;
auto x32 = -x23*x25 + x27 + x30;
auto x39 = x23*x37 + x36 - x38;
auto x40 = x34 + x39;
auto x41 = -x34 + x39;
auto x43 = x23*x35;
auto x44 = x20*x37 + x36 - x43;
auto x45 = x42 + x44;
auto x46 = -x42 + x44;
auto x48 = x24*x37 + x38 + x43;
auto x49 = x47 + x48;
auto x50 = -x47 + x48;
V v_V9 [DESCRIPTOR::d]; cell.computeU(v_V9);
V v_S6 = 0.444444444444444*v_S0 + 0.444444444444444*v_S1 + 0.444444444444444*v_S2 - 0.111111111111111*v_S3 - 0.111111111111111*v_S4 - 0.111111111111111*v_S5;
V v_P0 [DESCRIPTOR::q]; cell.getDynamics().computeEquilibrium(cell, v_S6, v_V9, v_P0);
cell[0] = v_P0[0] + x20*x22 + x22*x23 - x22*x24;
cell[1] = v_P0[1] + x29;
cell[2] = v_P0[2] + x31;
cell[3] = v_P0[3] + x32;
cell[4] = v_P0[4] + x40;
cell[5] = v_P0[5] + x41;
cell[6] = v_P0[6] + x45;
cell[7] = v_P0[7] + x46;
cell[8] = v_P0[8] - x49;
cell[9] = v_P0[9] - x50;
cell[10] = v_P0[10] + x29;
cell[11] = v_P0[11] + x31;
cell[12] = v_P0[12] + x32;
cell[13] = v_P0[13] + x40;
cell[14] = v_P0[14] + x41;
cell[15] = v_P0[15] + x45;
cell[16] = v_P0[16] + x46;
cell[17] = v_P0[17] - x49;
cell[18] = v_P0[18] - x50;
}
};

}

}
