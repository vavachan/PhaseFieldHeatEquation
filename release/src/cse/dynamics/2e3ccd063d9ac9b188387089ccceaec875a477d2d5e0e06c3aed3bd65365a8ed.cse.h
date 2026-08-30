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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x23 = x22 + V{-1};
auto x24 = V{0.5}*cell[13];
auto x25 = V{0.5}*cell[14];
auto x26 = V{0.5}*cell[15];
auto x27 = V{0.5}*cell[16];
auto x28 = V{0.5}*cell[17];
auto x29 = V{0.5}*cell[18];
auto x30 = V{0.5}*cell[4];
auto x31 = V{0.5}*cell[5];
auto x32 = V{0.5}*cell[6];
auto x33 = V{0.5}*cell[7];
auto x34 = V{0.5}*cell[8];
auto x35 = V{0.5}*cell[9];
auto x36 = V{0.5}*cell[0];
auto x37 = V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + x24 + x25 + x26 + x27 + x28 + x29 + x30 + x31 + x32 + x33 + x34 + x35 + x36 + V{0.5};
auto x38 = x19*x19;
auto x39 = x20*x20;
auto x40 = x21*x21;
auto x41 = V{1.5}*x38;
auto x42 = V{1.5}*x39;
auto x43 = V{1.5}*x40;
auto x44 = x42 + x43 + V{-1};
auto x45 = x41 + x44;
auto x46 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x47 = V{3}*x19;
auto x48 = V{3}*x38;
auto x49 = V{0.166666666666667}*cell[10];
auto x50 = V{0.166666666666667}*cell[1];
auto x51 = V{0.166666666666667}*cell[17];
auto x52 = V{0.166666666666667}*cell[18];
auto x53 = V{0.166666666666667}*cell[8];
auto x54 = V{0.166666666666667}*cell[9];
auto x55 = V{0.166666666666667}*cell[11];
auto x56 = V{0.166666666666667}*cell[12];
auto x57 = V{0.166666666666667}*cell[13];
auto x58 = V{0.166666666666667}*cell[14];
auto x59 = V{0.166666666666667}*cell[15];
auto x60 = V{0.166666666666667}*cell[16];
auto x61 = V{0.166666666666667}*cell[2];
auto x62 = V{0.166666666666667}*cell[3];
auto x63 = V{0.166666666666667}*cell[4];
auto x64 = V{0.166666666666667}*cell[5];
auto x65 = V{0.166666666666667}*cell[6];
auto x66 = V{0.166666666666667}*cell[7];
auto x67 = V{0.166666666666667}*cell[0] + x49 + x50 + x51 + x52 + x53 + x54 + x55 + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64 + x65 + x66 + V{0.166666666666667};
auto x68 = V{0.0833333333333333}*cell[13];
auto x69 = V{0.0833333333333333}*cell[14];
auto x70 = V{0.0833333333333333}*cell[4];
auto x71 = V{0.0833333333333333}*cell[5];
auto x72 = -V{6.93889390390723e-18}*cell[0];
auto x73 = V{0.0833333333333333}*cell[0] + V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333333}*cell[2] + V{0.0833333333333333}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9] + V{0.0833333333333333};
auto x74 = V{0.0833333333333333}*cell[12];
auto x75 = V{0.0833333333333333}*cell[3];
auto x76 = -x74 - x75;
auto x77 = x40*x73 + x68 + x69 + x70 + x71 + x72 + x76;
auto x78 = V{0.0833333333333333}*cell[15];
auto x79 = V{0.0833333333333333}*cell[16];
auto x80 = V{0.0833333333333333}*cell[6];
auto x81 = V{0.0833333333333333}*cell[7];
auto x82 = V{0.0833333333333333}*cell[11];
auto x83 = V{0.0833333333333333}*cell[2];
auto x84 = -x82 - x83;
auto x85 = x39*x73 + x78 + x79 + x80 + x81 + x84;
auto x86 = x23*(-x38*x67 + x49 + x50 - x51 - x52 - x53 - x54 + x77 + x85) + V{0.0555555555555556};
auto x87 = V{3}*x20;
auto x88 = V{3}*x39;
auto x89 = x41 + V{-1};
auto x90 = V{0.0833333333333333}*cell[17];
auto x91 = V{0.0833333333333333}*cell[18];
auto x92 = V{0.0833333333333333}*cell[8];
auto x93 = V{0.0833333333333333}*cell[9];
auto x94 = V{0.0833333333333333}*cell[10];
auto x95 = V{0.0833333333333333}*cell[1];
auto x96 = -x94 - x95;
auto x97 = x38*x73 + x90 + x91 + x92 + x93 + x96;
auto x98 = x23*(-x39*x67 + x55 - x59 - x60 + x61 - x65 - x66 + x77 + x97) + V{0.0555555555555556};
auto x99 = V{3}*x21;
auto x100 = V{3}*x40;
auto x101 = x23*(-x40*x67 + x56 - x57 - x58 + x62 - x63 - x64 + x72 + x85 + x97) + V{0.0555555555555556};
auto x102 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x103 = x19 + x20;
auto x104 = V{4.5}*(x103*x103);
auto x105 = x45 + x47;
auto x106 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x107 = x106 + V{1};
auto x108 = V{0.25}*x107*x19;
auto x109 = x108*x20;
auto x110 = V{0.0416666666666667}*cell[0] + V{0.0416666666666667}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0416666666666667}*cell[2] + V{0.0416666666666667}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + V{0.0416666666666667};
auto x111 = -V{0.0416666666666667}*cell[0];
auto x112 = V{0.0833333333333333}*cell[0] + x68 + x69 + x70 + x71 + x74 + x75 + x78 + x79 + x80 + x81 + x82 + x83 + x90 + x91 + x92 + x93 + x94 + x95 + V{0.0833333333333333};
auto x113 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x111 - x112*x38;
auto x114 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x112*x39;
auto x115 = x110*x40 + x113 + x114 + x76;
auto x116 = x23*(V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x109 + x115) + V{0.0277777777777778};
auto x117 = -x87;
auto x118 = x19 - x20;
auto x119 = -x118;
auto x120 = x23*(-V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x109 + x115) + V{0.0277777777777778};
auto x121 = x19 + x21;
auto x122 = V{4.5}*(x121*x121);
auto x123 = x108*x21;
auto x124 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x112*x40;
auto x125 = x110*x39 + x113 + x124 + x84;
auto x126 = x23*(V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x123 + x125) + V{0.0277777777777778};
auto x127 = -x99;
auto x128 = -x21;
auto x129 = x128 + x19;
auto x130 = -x129;
auto x131 = x23*(-V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x123 + x125) + V{0.0277777777777778};
auto x132 = x20 + x21;
auto x133 = V{4.5}*(x132*x132);
auto x134 = x45 + x87;
auto x135 = V{0.25}*x107*x20*x21;
auto x136 = x110*x38 + x111 + x114 + x124 + x96;
auto x137 = x23*(V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x135 + x136) + V{0.0277777777777778};
auto x138 = x128 + x20;
auto x139 = -x138;
auto x140 = x23*(-V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x135 + x136) + V{0.0277777777777778};
auto x141 = x106 + V{1};
auto x142 = -x42;
auto x143 = V{1} - x43;
auto x144 = x142 + x143;
auto x145 = x144 + x47;
auto x146 = -x41;
auto x147 = x146 + x87;
auto x148 = x146 + x99;
auto x149 = -x47;
auto x150 = x45 + x99;
auto x0 = x23*(V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + x24 + x25 + x26 + x27 + x28 + x29 + x30 + x31 + x32 + x33 + x34 + x35 - x36 - x37*x38 - x37*x39 - x37*x40) - x45*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = -x46*(x44 + x47 - x48) - x86;
auto x2 = -x46*(x43 + x87 - x88 + x89) - x98;
auto x3 = -x101 - x46*(-x100 + x42 + x89 + x99);
auto x4 = -x102*(-x104 + x105 + x87) - x116;
auto x5 = -(x102*(x105 + x117 - V{4.5}*x119*x119) + x120);
auto x6 = -x102*(x105 - x122 + x99) - x126;
auto x7 = -(x102*(x105 + x127 - V{4.5}*x130*x130) + x131);
auto x8 = -x102*(-x133 + x134 + x99) - x137;
auto x9 = -(x102*(x127 + x134 - V{4.5}*x139*x139) + x140);
auto x10 = V{0.0555555555555556}*x141*(x145 + x48) - x86;
auto x11 = V{0.0555555555555556}*x141*(x143 + x147 + x88) - x98;
auto x12 = -x101 + V{0.0555555555555556}*x141*(x100 + x142 + x148 + V{1});
auto x13 = -x116 + V{0.0277777777777778}*x141*(x104 + x145 + x147);
auto x14 = -(x102*(x134 + x149 - V{4.5}*x118*x118) + x120);
auto x15 = -x126 + V{0.0277777777777778}*x141*(x122 + x145 + x148);
auto x16 = -(x102*(x149 + x150 - V{4.5}*x129*x129) + x131);
auto x17 = -x137 + V{0.0277777777777778}*x141*(x133 + x144 + x147 + x99);
auto x18 = -(x102*(x117 + x150 - V{4.5}*x138*x138) + x140);
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
return { x107, x38 + x39 + x40 };
}
};

}

}
