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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x23 = x22 + V{-1};
auto x24 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x25 = x24 + V{1};
auto x26 = x24 + V{1};
auto x27 = V{1} / (x26);
auto x28 = V{1}*cell[14];
auto x29 = V{1}*cell[16];
auto x30 = V{1}*cell[5];
auto x31 = V{1}*cell[7];
auto x32 = V{1}*cell[13] - V{1}*cell[4];
auto x33 = V{1}*cell[15] - V{1}*cell[6];
auto x34 = V{0.5}*x19 + x27*(V{1}*cell[10] - V{1}*cell[1] + x28 + x29 - x30 - x31 + x32 + x33);
auto x35 = x34*x34;
auto x36 = V{1.5}*x35;
auto x37 = V{0.5}*x20;
auto x38 = V{1}*cell[18];
auto x39 = V{1}*cell[9];
auto x40 = V{1}*cell[17] - V{1}*cell[8];
auto x41 = x27*(V{1}*cell[11] - V{1}*cell[2] - x28 + x30 + x32 + x38 - x39 + x40);
auto x42 = x37 + x41;
auto x43 = x42*x42;
auto x44 = V{1.5}*x43;
auto x45 = V{1}*cell[12] - V{1}*cell[3] - x29 + x31 + x33 - x38 + x39 + x40;
auto x46 = V{0.5}*x21 + x27*x45;
auto x47 = x46*x46;
auto x48 = V{1.5}*x47;
auto x49 = x36 + x44 + x48 + V{-1};
auto x50 = V{0.5}*x22 + V{-1};
auto x51 = V{1.5}*x19;
auto x52 = V{3}*cell[14];
auto x53 = V{3}*cell[16];
auto x54 = V{3}*cell[5];
auto x55 = V{3}*cell[7];
auto x56 = V{3}*cell[13] - V{3}*cell[4];
auto x57 = V{3}*cell[15] - V{3}*cell[6];
auto x58 = x27*(V{3}*cell[10] - V{3}*cell[1] + x52 + x53 - x54 - x55 + x56 + x57);
auto x59 = x51 + x58;
auto x60 = x19*x59;
auto x61 = V{1.5}*x20;
auto x62 = V{3}*cell[18];
auto x63 = V{3}*cell[9];
auto x64 = V{3}*cell[17] - V{3}*cell[8];
auto x65 = x27*(V{3}*cell[11] - V{3}*cell[2] - x52 + x54 + x56 + x62 - x63 + x64);
auto x66 = x61 + x65;
auto x67 = x20*x66;
auto x68 = V{1.5}*x21;
auto x69 = x27*(V{3}*cell[12] - V{3}*cell[3] - x53 + x55 + x57 - x62 + x63 + x64);
auto x70 = x68 + x69;
auto x71 = x21*x70;
auto x72 = x67 + x71;
auto x73 = V{0.0555555555555556}*x22;
auto x74 = V{4.5}*cell[14];
auto x75 = V{4.5}*cell[16];
auto x76 = V{4.5}*cell[5];
auto x77 = V{4.5}*cell[7];
auto x78 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x79 = V{4.5}*cell[15] - V{4.5}*cell[6];
auto x80 = V{2.25}*x19 + x27*(V{4.5}*cell[10] - V{4.5}*cell[1] + x74 + x75 - x76 - x77 + x78 + x79);
auto x81 = x34*x80;
auto x82 = x49 + x59;
auto x83 = V{3}*x19;
auto x84 = V{6}*cell[14];
auto x85 = V{6}*cell[16];
auto x86 = V{6}*cell[5];
auto x87 = V{6}*cell[7];
auto x88 = V{6}*cell[13] - V{6}*cell[4];
auto x89 = V{6}*cell[15] - V{6}*cell[6];
auto x90 = x27*(V{6}*cell[10] - V{6}*cell[1] + x84 + x85 - x86 - x87 + x88 + x89);
auto x91 = x83 + x90;
auto x92 = x91 + V{-3};
auto x93 = x26*x50;
auto x94 = V{0.0555555555555556}*x93;
auto x95 = V{2.25}*x20;
auto x96 = V{4.5}*cell[18];
auto x97 = V{4.5}*cell[9];
auto x98 = V{4.5}*cell[17] - V{4.5}*cell[8];
auto x99 = x27*(V{4.5}*cell[11] - V{4.5}*cell[2] - x74 + x76 + x78 + x96 - x97 + x98);
auto x100 = x95 + x99;
auto x101 = x100*x42;
auto x102 = x49 + x66;
auto x103 = V{3}*x20;
auto x104 = V{6}*cell[18];
auto x105 = V{6}*cell[9];
auto x106 = V{6}*cell[17] - V{6}*cell[8];
auto x107 = x27*(V{6}*cell[11] - V{6}*cell[2] + x104 - x105 + x106 - x84 + x86 + x88);
auto x108 = x103 + x107;
auto x109 = x108 + V{-3};
auto x110 = x60 + x71;
auto x111 = V{4.5}*cell[12] - V{4.5}*cell[3] - x75 + x77 + x79 - x96 + x97 + x98;
auto x112 = x111*x27 + V{2.25}*x21;
auto x113 = x112*x46;
auto x114 = x49 + x70;
auto x115 = V{3}*x21;
auto x116 = x27*(V{6}*cell[12] - V{6}*cell[3] - x104 + x105 + x106 - x85 + x87 + x89);
auto x117 = x115 + x116;
auto x118 = x117 + V{-3};
auto x119 = x60 + x67;
auto x120 = V{0.0277777777777778}*x22;
auto x121 = (x100 + x80)*(x34 + x42);
auto x122 = V{4.5}*x20;
auto x123 = V{9}*cell[18];
auto x124 = V{9}*cell[5];
auto x125 = V{9}*cell[14];
auto x126 = V{9}*cell[9];
auto x127 = V{9}*cell[13] - V{9}*cell[4];
auto x128 = V{9}*cell[17] - V{9}*cell[8];
auto x129 = x27*(V{9}*cell[11] - V{9}*cell[2] + x123 + x124 - x125 - x126 + x127 + x128);
auto x130 = x122 + x129;
auto x131 = V{4.5}*x19;
auto x132 = V{9}*cell[16];
auto x133 = V{9}*cell[7];
auto x134 = V{9}*cell[15] - V{9}*cell[6];
auto x135 = x27*(V{9}*cell[10] - V{9}*cell[1] - x124 + x125 + x127 + x132 - x133 + x134);
auto x136 = x131 + x135;
auto x137 = -x71;
auto x138 = V{0.0277777777777778}*x93;
auto x139 = x34 - x37 - x41;
auto x140 = x80 - x95 - x99;
auto x141 = -x61 - x65;
auto x142 = x130 + V{3};
auto x143 = -x83 - x90;
auto x144 = x108 + V{3};
auto x145 = -x131 - x135;
auto x146 = (x112 + x80)*(x34 + x46);
auto x147 = V{4.5}*x21;
auto x148 = x27*(V{9}*cell[12] - V{9}*cell[3] - x123 + x126 + x128 - x132 + x133 + x134);
auto x149 = x147 + x148;
auto x150 = -x67;
auto x151 = -V{0.5}*x21 - x27*x45;
auto x152 = x151 + x34;
auto x153 = -x111*x27 - V{2.25}*x21;
auto x154 = x153 + x80;
auto x155 = -x68 - x69;
auto x156 = x149 + V{3};
auto x157 = x117 + V{3};
auto x158 = (x100 + x112)*(x42 + x46);
auto x159 = -x60;
auto x160 = x151 + x42;
auto x161 = x100 + x153;
auto x162 = -x103 - x107;
auto x163 = -x122 - x129;
auto x164 = -x36 - x44 - x48 + V{1};
auto x165 = x164 + x59;
auto x166 = x91 + V{3};
auto x167 = x164 + x66;
auto x168 = -x51 - x58;
auto x169 = x136 + V{3};
auto x170 = -x115 - x116;
auto x171 = -x147 - x148;
auto x0 = -cell[0]*x23 - V{0.333333333333333}*x22*(x25*x49 + V{1}) + V{0.333333333333333}*x26*x50*(x60 + x72);
auto x1 = -cell[1]*x23 - x73*(x25*(-x81 + x82) + V{1}) - x94*(x19*x92 - x72);
auto x2 = -cell[2]*x23 - x73*(x25*(-x101 + x102) + V{1}) - x94*(x109*x20 - x110);
auto x3 = -cell[3]*x23 - x73*(x25*(-x113 + x114) + V{1}) - x94*(x118*x21 - x119);
auto x4 = -cell[4]*x23 - x120*(x25*(-x121 + x66 + x82) + V{1}) - x138*(x137 + x19*(x130 + x92) + x20*(x109 + x136));
auto x5 = -cell[5]*x23 - x120*(x25*(-x139*x140 + x141 + x82) + V{1}) - x138*(-x19*(x142 + x143) + x20*(x144 + x145) - x71);
auto x6 = -cell[6]*x23 - x120*(x25*(-x146 + x70 + x82) + V{1}) - x138*(x150 + x19*(x149 + x92) + x21*(x118 + x136));
auto x7 = -cell[7]*x23 - x120*(x25*(-x152*x154 + x155 + x82) + V{1}) - x138*(-x19*(x143 + x156) + x21*(x145 + x157) - x67);
auto x8 = -cell[8]*x23 - x120*(x25*(x102 - x158 + x70) + V{1}) - x138*(x159 + x20*(x109 + x149) + x21*(x118 + x130));
auto x9 = -cell[9]*x23 - x120*(x25*(x102 + x155 - x160*x161) + V{1}) - x138*(-x20*(x156 + x162) + x21*(x157 + x163) - x60);
auto x10 = -cell[10]*x23 + V{0.0555555555555556}*x22*(x25*(x165 + x81) + V{-1}) - x94*(x166*x19 - x72);
auto x11 = -cell[11]*x23 + V{0.0555555555555556}*x22*(x25*(x101 + x167) + V{-1}) - x94*(-x110 + x144*x20);
auto x12 = -cell[12]*x23 + V{0.0555555555555556}*x22*(x25*(x113 + x164 + x70) + V{-1}) - x94*(-x119 + x157*x21);
auto x13 = -cell[13]*x23 - x138*(x137 + x19*(x130 + x166) + x20*(x136 + x144)) + V{0.0277777777777778}*x22*(x25*(x121 + x165 + x66) + V{-1});
auto x14 = -cell[14]*x23 - x120*(x25*(x102 - x139*x140 + x168) + V{1}) - x138*(x19*(x163 + x166) - x20*(x162 + x169) - x71);
auto x15 = -cell[15]*x23 - x138*(x150 + x19*(x149 + x166) + x21*(x136 + x157)) + V{0.0277777777777778}*x22*(x25*(x146 + x165 + x70) + V{-1});
auto x16 = -cell[16]*x23 - x120*(x25*(x114 - x152*x154 + x168) + V{1}) - x138*(x19*(x166 + x171) - x21*(x169 + x170) - x67);
auto x17 = -cell[17]*x23 - x138*(x159 + x20*(x144 + x149) + x21*(x130 + x157)) + V{0.0277777777777778}*x22*(x25*(x158 + x167 + x70) + V{-1});
auto x18 = -cell[18]*x23 - x120*(x25*(x114 + x141 - x160*x161) + V{1}) - x138*(x20*(x144 + x171) - x21*(x142 + x170) - x60);
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
return { x26, x35 + x43 + x47 };
}
};

}

}
