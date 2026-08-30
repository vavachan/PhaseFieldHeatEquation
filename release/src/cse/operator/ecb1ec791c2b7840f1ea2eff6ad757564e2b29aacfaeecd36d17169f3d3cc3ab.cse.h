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
struct CSE_O<OuterVelocityEdgeProcessor3D<T, descriptors::D3Q19<FIELDS...>, 1, 1, 1>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x19 = cell.getDynamics().getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
auto x21 = V{1} / (x19);
V v_S0 = cell.neighbor({0,0,-1}).computeRho();
V v_S1 = cell.neighbor({-1,0,0}).computeRho();
V v_S2 = cell.neighbor({0,0,-2}).computeRho();
V v_S3 = cell.neighbor({-2,0,0}).computeRho();
auto x22 = V{0.222222222222222}*v_S0 + V{0.222222222222222}*v_S1 - V{0.0555555555555556}*v_S2 - V{0.0555555555555556}*v_S3;
auto x23 = x21*x22;
auto x24 = V{0.5}*x23;
auto x27 = V{0.0833333333333333}*x23;
auto x30 = V{0.166666666666667}*x23;
auto x35 = V{0.25}*x23;
auto x37 = V{0.0833333333333333}*x23;
auto x41 = V{0.0416666666666667}*x23;
V v_V0 [DESCRIPTOR::d]; cell.neighbor({0,1,0}).computeU(v_V0);
V v_V1 [DESCRIPTOR::d]; cell.neighbor({0,-1,0}).computeU(v_V1);
auto x20 = v_V0[1] - v_V1[1];
auto x28 = x20*x27;
auto x38 = x20*x37;
auto x42 = -x38;
V v_V2 [DESCRIPTOR::d]; cell.computeU(v_V2);
V v_V3 [DESCRIPTOR::d]; cell.neighbor({0,0,-1}).computeU(v_V3);
V v_V4 [DESCRIPTOR::d]; cell.neighbor({0,0,-2}).computeU(v_V4);
auto x50 = x35*(V{0.5}*v_V0[2] - V{0.5}*v_V1[2] + V{1.5}*v_V2[1] - V{2}*v_V3[1] + V{0.5}*v_V4[1]);
auto x25 = V{3}*v_V2[2] - V{4}*v_V3[2] + V{1}*v_V4[2];
auto x29 = x25*x27;
auto x46 = x25*x37;
auto x48 = -x46;
V v_V5 [DESCRIPTOR::d]; cell.computeU(v_V5);
V v_V6 [DESCRIPTOR::d]; cell.neighbor({-1,0,0}).computeU(v_V6);
V v_V7 [DESCRIPTOR::d]; cell.neighbor({-2,0,0}).computeU(v_V7);
auto x26 = V{3}*v_V5[0] - V{4}*v_V6[0] + V{1}*v_V7[0];
auto x31 = -x26*x30 + x28 + x29;
auto x32 = x26*x27;
auto x33 = -x20*x30 + x29 + x32;
auto x34 = -x25*x30 + x28 + x32;
auto x39 = x26*x37;
auto x43 = -x39;
auto x51 = -V{0.0416666666666667}*x21*x22*x26 + x38 + x46 + x50;
auto x52 = x26*x41 + x42 + x48 + x50;
auto x36 = x35*(V{0.5}*v_V0[0] - V{0.5}*v_V1[0] + V{1.5}*v_V5[1] - V{2}*v_V6[1] + V{0.5}*v_V7[1]);
auto x40 = -V{0.0416666666666667}*x21*x22*x25 + x36 + x38 + x39;
auto x44 = x25*x41 + x36 + x42 + x43;
auto x45 = x35*(V{1.5}*v_V2[0] - V{2}*v_V3[0] + V{0.5}*v_V4[0] + V{1.5}*v_V5[2] - V{2}*v_V6[2] + V{0.5}*v_V7[2]);
auto x47 = -V{0.0416666666666667}*x20*x21*x22 + x39 + x45 + x46;
auto x49 = x20*x41 + x43 + x45 + x48;
V v_V8 [DESCRIPTOR::d]; cell.computeU(v_V8);
V v_S4 = 0.666666666666667*v_S0 + 0.666666666666667*v_S1 - 0.166666666666667*v_S2 - 0.166666666666667*v_S3;
V v_P0 [DESCRIPTOR::q]; cell.getDynamics().computeEquilibrium(cell, v_S4, v_V8, v_P0);
cell[0] = v_P0[0] + x20*x24 + x24*x25 + x24*x26;
cell[1] = v_P0[1] + x31;
cell[2] = v_P0[2] + x33;
cell[3] = v_P0[3] + x34;
cell[4] = v_P0[4] - x40;
cell[5] = v_P0[5] + x44;
cell[6] = v_P0[6] - x47;
cell[7] = v_P0[7] + x49;
cell[8] = v_P0[8] - x51;
cell[9] = v_P0[9] + x52;
cell[10] = v_P0[10] + x31;
cell[11] = v_P0[11] + x33;
cell[12] = v_P0[12] + x34;
cell[13] = v_P0[13] - x40;
cell[14] = v_P0[14] + x44;
cell[15] = v_P0[15] - x47;
cell[16] = v_P0[16] + x49;
cell[17] = v_P0[17] - x51;
cell[18] = v_P0[18] + x52;
}
};

}

}
