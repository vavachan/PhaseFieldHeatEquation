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

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[10] + cell[14];
auto x22 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x21 + x22;
auto x24 = x23 + V{1};
auto x25 = V{1} / (x24);
auto x26 = V{0.5}*x25;
auto x27 = cell[13] - cell[4];
auto x28 = cell[15] - cell[6];
auto x29 = x27 + x28;
auto x30 = -cell[1];
auto x31 = cell[16] - cell[7];
auto x32 = x30 + x31;
auto x33 = -cell[5] + x21;
auto x34 = x29 + x32 + x33;
auto x35 = x34*x34;
auto x36 = cell[17] - cell[8];
auto x37 = x27 + x36;
auto x38 = cell[18] - cell[9];
auto x39 = -cell[2];
auto x40 = cell[11] - cell[14] + cell[5] + x39;
auto x41 = x37 + x38 + x40;
auto x42 = x41*x41;
auto x43 = x28 + x36;
auto x44 = -cell[3];
auto x45 = -cell[18] + cell[9];
auto x46 = x44 + x45;
auto x47 = -cell[16] + x22;
auto x48 = x43 + x46 + x47;
auto x49 = x48*x48;
auto x50 = V{1} / ((x24)*(x24));
auto x51 = V{1.5}*x50;
auto x52 = x35*x51;
auto x53 = x42*x51;
auto x54 = x49*x51;
auto x55 = x53 + x54 + V{-1};
auto x56 = x52 + x55;
auto x57 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x58 = V{3}*cell[14];
auto x59 = V{3}*cell[16];
auto x60 = V{3}*cell[5];
auto x61 = V{3}*cell[7];
auto x62 = V{3}*cell[13] - V{3}*cell[4];
auto x63 = V{3}*cell[15] - V{3}*cell[6];
auto x64 = x25*(V{3}*cell[10] - V{3}*cell[1] + x58 + x59 - x60 - x61 + x62 + x63);
auto x65 = V{3}*x50;
auto x66 = x35*x65;
auto x67 = V{0.166666666666667}*x25;
auto x68 = -V{6.93889390390723e-18}*cell[0];
auto x69 = V{0.0833333333333333}*x25;
auto x70 = -V{0.0833333333333333}*cell[12] - V{0.0833333333333333}*cell[3];
auto x71 = V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + x49*x69 + x68 + x70;
auto x72 = -V{0.0833333333333333}*cell[11] - V{0.0833333333333333}*cell[2];
auto x73 = V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + x42*x69 + x72;
auto x74 = x20*(V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] - x35*x67 + x71 + x73) + V{0.0555555555555556};
auto x75 = V{3}*cell[18];
auto x76 = V{3}*cell[9];
auto x77 = V{3}*cell[17] - V{3}*cell[8];
auto x78 = x25*(V{3}*cell[11] - V{3}*cell[2] - x58 + x60 + x62 + x75 - x76 + x77);
auto x79 = x42*x65;
auto x80 = x52 + V{-1};
auto x81 = -V{0.0833333333333333}*cell[10] - V{0.0833333333333333}*cell[1];
auto x82 = V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9] + x35*x69 + x81;
auto x83 = x20*(V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.166666666666667}*cell[7] - x42*x67 + x71 + x82) + V{0.0555555555555556};
auto x84 = x25*(V{3}*cell[12] - V{3}*cell[3] - x59 + x61 + x63 - x75 + x76 + x77);
auto x85 = x49*x65;
auto x86 = x20*(V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x49*x67 + x68 + x73 + x82) + V{0.0555555555555556};
auto x87 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x88 = V{4.5}*x50;
auto x89 = cell[10] + x32;
auto x90 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x38 + x39 + x43 + x89;
auto x91 = x88*(x90*x90);
auto x92 = x56 + x64;
auto x93 = V{0.25}*x25*x34;
auto x94 = x41*x93;
auto x95 = V{0.0416666666666667}*x25;
auto x96 = -V{0.0416666666666667}*cell[0];
auto x97 = V{0.0833333333333333}*x25;
auto x98 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] - x35*x97 + x96;
auto x99 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x42*x97;
auto x100 = x49*x95 + x70 + x98 + x99;
auto x101 = x20*(V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] + x100 - x94) + V{0.0277777777777778};
auto x102 = -x78;
auto x103 = -cell[17] + cell[8];
auto x104 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x103 + x28 + x45 + x89;
auto x105 = -x104;
auto x106 = x20*(-V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x100 + x94) + V{0.0277777777777778};
auto x107 = x30 + x33;
auto x108 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x107 + x37 + x46;
auto x109 = x88*(x108*x108);
auto x110 = x48*x93;
auto x111 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x49*x97;
auto x112 = x111 + x42*x95 + x72 + x98;
auto x113 = x20*(V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x110 + x112) + V{0.0277777777777778};
auto x114 = -x84;
auto x115 = -cell[12] + cell[3] + x27;
auto x116 = V{2}*cell[16] - V{2}*cell[7] + x103 + x107 + x115 + x38;
auto x117 = -x116;
auto x118 = x20*(-V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x110 + x112) + V{0.0277777777777778};
auto x119 = V{2}*cell[17] - V{2}*cell[8] + x29 + x40 + x44 + x47;
auto x120 = x88*(x119*x119);
auto x121 = x56 + x78;
auto x122 = V{0.25}*x25*x41*x48;
auto x123 = x111 + x35*x95 + x81 + x96 + x99;
auto x124 = x20*(V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x122 + x123) + V{0.0277777777777778};
auto x125 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x115 + x31 + x40;
auto x126 = -x125;
auto x127 = x20*(-V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x122 + x123) + V{0.0277777777777778};
auto x128 = x23 + V{1};
auto x129 = -x53;
auto x130 = V{1} - x54;
auto x131 = x129 + x130;
auto x132 = x131 + x64;
auto x133 = -x52;
auto x134 = x133 + x78;
auto x135 = x133 + x84;
auto x136 = -x64;
auto x137 = x56 + x84;
auto x0 = x20*(-V{0.5}*cell[0] + V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] - x26*x35 - x26*x42 - x26*x49) - x56*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = -x57*(x55 + x64 - x66) - x74;
auto x2 = -x57*(x54 + x78 - x79 + x80) - x83;
auto x3 = -x57*(x53 + x80 + x84 - x85) - x86;
auto x4 = -x101 - x87*(x78 - x91 + x92);
auto x5 = -(x106 + x87*(x102 - x88*x105*x105 + x92));
auto x6 = -x113 - x87*(-x109 + x84 + x92);
auto x7 = -(x118 + x87*(x114 - x88*x117*x117 + x92));
auto x8 = -x124 - x87*(-x120 + x121 + x84);
auto x9 = -(x127 + x87*(x114 + x121 - x88*x126*x126));
auto x10 = V{0.0555555555555556}*x128*(x132 + x66) - x74;
auto x11 = V{0.0555555555555556}*x128*(x130 + x134 + x79) - x83;
auto x12 = V{0.0555555555555556}*x128*(x129 + x135 + x85 + V{1}) - x86;
auto x13 = -x101 + V{0.0277777777777778}*x128*(x132 + x134 + x91);
auto x14 = -(x106 + x87*(x121 + x136 - x88*x104*x104));
auto x15 = -x113 + V{0.0277777777777778}*x128*(x109 + x132 + x135);
auto x16 = -(x118 + x87*(x136 + x137 - x88*x116*x116));
auto x17 = -x124 + V{0.0277777777777778}*x128*(x120 + x131 + x134 + x84);
auto x18 = -(x127 + x87*(x102 + x137 - x88*x125*x125));
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x24, V{1}*x50*(x35 + x42 + x49) };
}
};

}

}
