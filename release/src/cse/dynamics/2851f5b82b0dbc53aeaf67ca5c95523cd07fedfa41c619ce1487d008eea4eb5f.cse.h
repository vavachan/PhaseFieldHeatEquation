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
struct CSE<dynamics::Tuple<T, descriptors::D3Q15<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x17 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x16 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x23 = parameters.template get<multiphase::RHO_LIQUID>();
auto x20 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x15 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x21 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x22 = parameters.template get<multiphase::RHO_VAPOR>();
auto x18 = -x23;
auto x19 = x18 + x22;
auto x24 = V{1} / (x19);
auto x25 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x26 = x25 + V{1};
auto x27 = x20*(x18 + x26);
auto x28 = x24*x27;
auto x29 = x21*(-x22 + x26);
auto x30 = -x29/x19;
auto x31 = x28 + x30 + V{-1};
auto x32 = x27 - x29;
auto x33 = x24*x32;
auto x34 = x25 + V{1};
auto x35 = V{0.5}*x15;
auto x36 = V{1} / (x26);
auto x37 = V{1}*cell[14];
auto x38 = V{1}*cell[7];
auto x39 = V{1}*cell[11];
auto x40 = V{1}*cell[12];
auto x41 = -V{1}*cell[4];
auto x42 = V{1}*cell[5];
auto x43 = x39 + x40 + x41 - x42;
auto x44 = V{1}*cell[13];
auto x45 = V{1}*cell[6];
auto x46 = x44 - x45;
auto x47 = x36*(-V{1}*cell[1] + V{1}*cell[8] + x37 - x38 + x43 + x46);
auto x48 = x35 + x47;
auto x49 = x48*x48;
auto x50 = V{1.5}*x49;
auto x51 = V{0.5}*x16;
auto x52 = -x37 + x38;
auto x53 = x36*(-V{1}*cell[2] + V{1}*cell[9] + x43 - x44 + x45 + x52);
auto x54 = x51 + x53;
auto x55 = x54*x54;
auto x56 = V{1.5}*x55;
auto x57 = V{0.5}*x17;
auto x58 = x36*(V{1}*cell[10] - V{1}*cell[3] + x39 - x40 + x41 + x42 + x46 + x52);
auto x59 = x57 + x58;
auto x60 = x59*x59;
auto x61 = V{1.5}*x60;
auto x62 = x50 + x56 + x61 + V{-1};
auto x63 = V{0.5}*x28 + V{0.5}*x30 + V{-1};
auto x64 = V{1.5}*x15;
auto x65 = V{3}*cell[14];
auto x66 = V{3}*cell[7];
auto x67 = V{3}*cell[11];
auto x68 = V{3}*cell[12];
auto x69 = -V{3}*cell[4];
auto x70 = V{3}*cell[5];
auto x71 = x67 + x68 + x69 - x70;
auto x72 = V{3}*cell[13];
auto x73 = V{3}*cell[6];
auto x74 = x72 - x73;
auto x75 = x36*(-V{3}*cell[1] + V{3}*cell[8] + x65 - x66 + x71 + x74);
auto x76 = x64 + x75;
auto x77 = x15*x76;
auto x78 = V{1.5}*x16;
auto x79 = -x65 + x66;
auto x80 = x36*(-V{3}*cell[2] + V{3}*cell[9] + x71 - x72 + x73 + x79);
auto x81 = x78 + x80;
auto x82 = x16*x81;
auto x83 = V{1.5}*x17;
auto x84 = x36*(V{3}*cell[10] - V{3}*cell[3] + x67 - x68 + x69 + x70 + x74 + x79);
auto x85 = x83 + x84;
auto x86 = x17*x85;
auto x87 = x82 + x86;
auto x88 = V{0.111111111111111}*x33;
auto x89 = V{2.25}*x15;
auto x90 = V{4.5}*cell[14];
auto x91 = V{4.5}*cell[7];
auto x92 = V{4.5}*cell[11];
auto x93 = V{4.5}*cell[12];
auto x94 = -V{4.5}*cell[4];
auto x95 = V{4.5}*cell[5];
auto x96 = x92 + x93 + x94 - x95;
auto x97 = V{4.5}*cell[13];
auto x98 = V{4.5}*cell[6];
auto x99 = x97 - x98;
auto x100 = x36*(-V{4.5}*cell[1] + V{4.5}*cell[8] + x90 - x91 + x96 + x99);
auto x101 = x100 + x89;
auto x102 = x101*x48;
auto x103 = x62 + x76;
auto x104 = V{3}*x15;
auto x105 = V{6}*cell[14];
auto x106 = V{6}*cell[7];
auto x107 = V{6}*cell[11];
auto x108 = V{6}*cell[12];
auto x109 = -V{6}*cell[4];
auto x110 = V{6}*cell[5];
auto x111 = x107 + x108 + x109 - x110;
auto x112 = V{6}*cell[13];
auto x113 = V{6}*cell[6];
auto x114 = x112 - x113;
auto x115 = x36*(-V{6}*cell[1] + V{6}*cell[8] + x105 - x106 + x111 + x114);
auto x116 = x104 + x115;
auto x117 = x116 + V{-3};
auto x118 = x26*x63;
auto x119 = V{0.111111111111111}*x118;
auto x120 = V{2.25}*x16;
auto x121 = -x90 + x91;
auto x122 = x36*(-V{4.5}*cell[2] + V{4.5}*cell[9] + x121 + x96 - x97 + x98);
auto x123 = x120 + x122;
auto x124 = x123*x54;
auto x125 = x62 + x81;
auto x126 = V{3}*x16;
auto x127 = -x105 + x106;
auto x128 = x36*(-V{6}*cell[2] + V{6}*cell[9] + x111 - x112 + x113 + x127);
auto x129 = x126 + x128;
auto x130 = x129 + V{-3};
auto x131 = x77 + x86;
auto x132 = V{2.25}*x17;
auto x133 = x36*(V{4.5}*cell[10] - V{4.5}*cell[3] + x121 + x92 - x93 + x94 + x95 + x99);
auto x134 = x132 + x133;
auto x135 = x134*x59;
auto x136 = V{3}*x17;
auto x137 = x36*(V{6}*cell[10] - V{6}*cell[3] + x107 - x108 + x109 + x110 + x114 + x127);
auto x138 = x136 + x137;
auto x139 = x138 + V{-3};
auto x140 = x77 + x82;
auto x141 = V{0.0138888888888889}*x33;
auto x142 = x48 + x54;
auto x143 = x101 + x123;
auto x144 = (x134 + x143)*(x142 + x59);
auto x145 = x103 + x81;
auto x146 = V{4.5}*x17;
auto x147 = V{9}*cell[5];
auto x148 = V{9}*cell[12];
auto x149 = V{9}*cell[11];
auto x150 = V{9}*cell[7];
auto x151 = V{9}*cell[14];
auto x152 = -V{9}*cell[4];
auto x153 = x149 + x150 - x151 + x152;
auto x154 = V{9}*cell[13];
auto x155 = V{9}*cell[6];
auto x156 = x154 - x155;
auto x157 = x36*(V{9}*cell[10] - V{9}*cell[3] + x147 - x148 + x153 + x156);
auto x158 = x146 + x157;
auto x159 = V{4.5}*x16;
auto x160 = -x147 + x148;
auto x161 = x36*(-V{9}*cell[2] + V{9}*cell[9] + x153 - x154 + x155 + x160);
auto x162 = x159 + x161;
auto x163 = x117 + x162;
auto x164 = V{4.5}*x15;
auto x165 = x36*(-V{9}*cell[1] + V{9}*cell[8] + x149 - x150 + x151 + x152 + x156 + x160);
auto x166 = x164 + x165;
auto x167 = x130 + x158;
auto x168 = x162 + x166;
auto x169 = V{0.0138888888888889}*x118;
auto x170 = x142 - x57 - x58;
auto x171 = -x132 - x133 + x143;
auto x172 = -x83 - x84;
auto x173 = -x146 - x157;
auto x174 = x166 + x173;
auto x175 = -x136 - x137 + x168;
auto x176 = x48 - x51 - x53 + x59;
auto x177 = x101 - x120 - x122 + x134;
auto x178 = -x78 - x80 + x85;
auto x179 = -x159 - x161;
auto x180 = x158 + x179;
auto x181 = x166 + x179;
auto x182 = x158 + x166;
auto x183 = -x126 - x128 + x182;
auto x184 = x129 + V{3};
auto x185 = -x164 - x165;
auto x186 = x138 + V{3};
auto x187 = x162 + x185;
auto x188 = x158 + x162;
auto x189 = -x104 - x115 + x188;
auto x190 = -x35 - x47 + x54 + x59;
auto x191 = -x100 + x123 + x134 - x89;
auto x192 = -x50 - x56 - x61 + V{1};
auto x193 = x192 + x81;
auto x194 = -x64 - x75 + x85;
auto x195 = x116 + V{3};
auto x196 = x192 + x76;
auto x197 = x196 + x81;
auto x0 = -cell[0]*x31 + V{0.222222222222222}*x26*x63*(x77 + x87) - V{0.222222222222222}*x33*(x34*x62 + V{1});
auto x1 = -cell[1]*x31 - x119*(x117*x15 - x87) - x88*(x34*(-x102 + x103) + V{1});
auto x2 = -cell[2]*x31 - x119*(x130*x16 - x131) - x88*(x34*(-x124 + x125) + V{1});
auto x3 = -cell[3]*x31 - x119*(x139*x17 - x140) - x88*(x34*(-x135 + x62 + x85) + V{1});
auto x4 = -cell[4]*x31 - x141*(x34*(-x144 + x145 + x85) + V{1}) - x169*(x15*(x158 + x163) + x16*(x166 + x167) + x17*(x139 + x168));
auto x5 = -cell[5]*x31 - x141*(x34*(x145 - x170*x171 + x172) + V{1}) - x169*(x15*(x163 + x173) + x16*(x130 + x174) - x17*(x175 + V{-3}));
auto x6 = -cell[6]*x31 - x141*(x34*(x103 - x176*x177 + x178) + V{1}) - x169*(x15*(x117 + x180) - x16*(x183 + V{-3}) + x17*(x139 + x181));
auto x7 = -cell[7]*x31 - x169*(-x15*(x189 + V{3}) + x16*(x158 + x184 + x185) + x17*(x186 + x187)) + V{0.0138888888888889}*x24*x32*(x34*(x190*x191 + x193 + x194) + V{-1});
auto x8 = -cell[8]*x31 - x119*(x15*x195 - x87) + V{0.111111111111111}*x24*x32*(x34*(x102 + x196) + V{-1});
auto x9 = -cell[9]*x31 - x119*(-x131 + x16*x184) + V{0.111111111111111}*x24*x32*(x34*(x124 + x193) + V{-1});
auto x10 = -cell[10]*x31 - x119*(-x140 + x17*x186) + V{0.111111111111111}*x24*x32*(x34*(x135 + x192 + x85) + V{-1});
auto x11 = -cell[11]*x31 - x169*(x15*(x188 + x195) + x16*(x182 + x184) + x17*(x168 + x186)) + V{0.0138888888888889}*x24*x32*(x34*(x144 + x197 + x85) + V{-1});
auto x12 = -cell[12]*x31 - x169*(x15*(x162 + x173 + x195) + x16*(x174 + x184) - x17*(x175 + V{3})) + V{0.0138888888888889}*x24*x32*(x34*(x170*x171 + x172 + x197) + V{-1});
auto x13 = -cell[13]*x31 - x169*(x15*(x180 + x195) - x16*(x183 + V{3}) + x17*(x181 + x186)) + V{0.0138888888888889}*x24*x32*(x34*(x176*x177 + x178 + x196) + V{-1});
auto x14 = -cell[14]*x31 - x141*(x34*(x125 - x190*x191 + x194) + V{1}) - x169*(-x15*(x189 + V{-3}) + x16*(x167 + x185) + x17*(x139 + x187));
cell.template getFieldPointer<descriptors::OMEGA>()[0] = x33;
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
return { x26, x49 + x55 + x60 };
}
};

}

}
