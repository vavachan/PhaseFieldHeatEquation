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
struct CSE<SpongeLayerDynamics<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<descriptors::UX>(0);
auto x23 = cell.template getFieldComponent<descriptors::UZ>(0);
auto x22 = cell.template getFieldComponent<descriptors::UY>(0);
auto x24 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<descriptors::DAMPING>(0);
auto x20 = cell.template getFieldComponent<descriptors::DENSITY>(0);
auto x25 = cell[10] + cell[14];
auto x26 = cell[12] + cell[7];
auto x27 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x25 + x26;
auto x28 = x27 + V{1};
auto x29 = V{1} / ((x28)*(x28));
auto x30 = V{1.5}*x29;
auto x31 = -cell[4];
auto x32 = cell[13] + x31;
auto x33 = -cell[6];
auto x34 = cell[15] + x33;
auto x35 = x32 + x34;
auto x36 = -cell[1];
auto x37 = -cell[7];
auto x38 = cell[16] + x37;
auto x39 = x36 + x38;
auto x40 = -cell[5];
auto x41 = x25 + x40;
auto x42 = x35 + x39 + x41;
auto x43 = x42*x42;
auto x44 = x30*x43;
auto x45 = -cell[8];
auto x46 = cell[17] + x45;
auto x47 = x32 + x46;
auto x48 = -cell[9];
auto x49 = cell[18] + x48;
auto x50 = -cell[2];
auto x51 = -cell[14];
auto x52 = cell[11] + cell[5] + x50 + x51;
auto x53 = x47 + x49 + x52;
auto x54 = x53*x53;
auto x55 = x30*x54;
auto x56 = x34 + x46;
auto x57 = -cell[3];
auto x58 = -cell[18];
auto x59 = cell[9] + x58;
auto x60 = x57 + x59;
auto x61 = -cell[16];
auto x62 = x26 + x61;
auto x63 = x56 + x60 + x62;
auto x64 = x63*x63;
auto x65 = x30*x64;
auto x66 = x55 + x65 + V{-1};
auto x67 = x44 + x66;
auto x68 = x21*x21;
auto x69 = V{1.5}*x68;
auto x70 = x22*x22;
auto x71 = V{1.5}*x70;
auto x72 = x23*x23;
auto x73 = V{1.5}*x72;
auto x74 = x71 + x73 + V{-1};
auto x75 = x69 + x74;
auto x76 = x27 + V{1};
auto x77 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x78 = V{1} / (x28);
auto x79 = V{3}*cell[14];
auto x80 = V{3}*cell[16];
auto x81 = V{3}*cell[5];
auto x82 = V{3}*cell[7];
auto x83 = V{3}*cell[13] - V{3}*cell[4];
auto x84 = V{3}*cell[15] - V{3}*cell[6];
auto x85 = x78*(V{3}*cell[10] - V{3}*cell[1] + x79 + x80 - x81 - x82 + x83 + x84);
auto x86 = V{3}*x29;
auto x87 = x43*x86;
auto x88 = x66 + x85 - x87;
auto x89 = V{0.0555555555555556}*x19;
auto x90 = V{3}*x21;
auto x91 = V{3}*x68;
auto x92 = V{3}*cell[18];
auto x93 = V{3}*cell[9];
auto x94 = V{3}*cell[17] - V{3}*cell[8];
auto x95 = x78*(V{3}*cell[11] - V{3}*cell[2] - x79 + x81 + x83 + x92 - x93 + x94);
auto x96 = x54*x86;
auto x97 = x44 + V{-1};
auto x98 = x65 + x95 - x96 + x97;
auto x99 = V{3}*x22;
auto x100 = V{3}*x70;
auto x101 = x69 + V{-1};
auto x102 = x78*(V{3}*cell[12] - V{3}*cell[3] - x80 + x82 + x84 - x92 + x93 + x94);
auto x103 = x64*x86;
auto x104 = x102 - x103 + x55 + x97;
auto x105 = V{3}*x23;
auto x106 = V{3}*x72;
auto x107 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x108 = V{4.5}*x29;
auto x109 = cell[10] + x39;
auto x110 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x109 + x49 + x50 + x56;
auto x111 = x108*(x110*x110);
auto x112 = x67 + x85;
auto x113 = -x111 + x112 + x95;
auto x114 = V{0.0277777777777778}*x19;
auto x115 = x21 + x22;
auto x116 = V{4.5}*(x115*x115);
auto x117 = x75 + x90;
auto x118 = -x95;
auto x119 = -cell[17] + cell[8];
auto x120 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x109 + x119 + x34 + x59;
auto x121 = -x120;
auto x122 = -x108*x121*x121 + x112 + x118;
auto x123 = -x99;
auto x124 = x21 - x22;
auto x125 = -x124;
auto x126 = x36 + x41;
auto x127 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x126 + x47 + x60;
auto x128 = x108*(x127*x127);
auto x129 = x102 + x112 - x128;
auto x130 = x21 + x23;
auto x131 = V{4.5}*(x130*x130);
auto x132 = -x102;
auto x133 = -cell[12] + cell[3] + x32;
auto x134 = V{2}*cell[16] - V{2}*cell[7] + x119 + x126 + x133 + x49;
auto x135 = -x134;
auto x136 = -x108*x135*x135 + x112 + x132;
auto x137 = -x105;
auto x138 = -x23;
auto x139 = x138 + x21;
auto x140 = -x139;
auto x141 = V{2}*cell[17] - V{2}*cell[8] + x35 + x52 + x57 + x62;
auto x142 = x108*(x141*x141);
auto x143 = x67 + x95;
auto x144 = x102 - x142 + x143;
auto x145 = x22 + x23;
auto x146 = V{4.5}*(x145*x145);
auto x147 = x75 + x99;
auto x148 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x133 + x38 + x52;
auto x149 = -x148;
auto x150 = -x108*x149*x149 + x132 + x143;
auto x151 = x138 + x22;
auto x152 = -x151;
auto x153 = -x55;
auto x154 = V{1} - x65;
auto x155 = x153 + x154;
auto x156 = x155 + x85;
auto x157 = x156 + x87;
auto x158 = -x71;
auto x159 = V{1} - x73;
auto x160 = x158 + x159;
auto x161 = x160 + x90;
auto x162 = -x44;
auto x163 = x162 + x95;
auto x164 = x154 + x163 + x96;
auto x165 = -x69;
auto x166 = x165 + x99;
auto x167 = x102 + x162;
auto x168 = x103 + x153 + x167 + V{1};
auto x169 = x105 + x165;
auto x170 = x111 + x156 + x163;
auto x171 = -x85;
auto x172 = -x108*x120*x120 + x143 + x171;
auto x173 = -x90;
auto x174 = x128 + x156 + x167;
auto x175 = x102 + x67;
auto x176 = -x108*x134*x134 + x171 + x175;
auto x177 = x105 + x75;
auto x178 = x102 + x142 + x155 + x163;
auto x179 = -x108*x148*x148 + x118 + x175;
auto x0 = cell[0] - V{0.333333333333333}*x19*(x20*x75 - x67*x76) - x24*(cell[0] + x67*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = -x24*(cell[1] + x77*x88 + V{0.0555555555555556}) - x36 - x89*(x20*(x74 + x90 - x91) - x76*x88);
auto x2 = -x24*(cell[2] + x77*x98 + V{0.0555555555555556}) - x50 - x89*(x20*(-x100 + x101 + x73 + x99) - x76*x98);
auto x3 = -x24*(cell[3] + x104*x77 + V{0.0555555555555556}) - x57 - x89*(-x104*x76 + x20*(x101 + x105 - x106 + x71));
auto x4 = -x114*(-x113*x76 + x20*(-x116 + x117 + x99)) - x24*(cell[4] + x107*x113 + V{0.0277777777777778}) - x31;
auto x5 = -(x114*(-x122*x76 + x20*(x117 + x123 - V{4.5}*x125*x125)) + x24*(cell[5] + x107*x122 + V{0.0277777777777778}) + x40);
auto x6 = -x114*(-x129*x76 + x20*(x105 + x117 - x131)) - x24*(cell[6] + x107*x129 + V{0.0277777777777778}) - x33;
auto x7 = -(x114*(-x136*x76 + x20*(x117 + x137 - V{4.5}*x140*x140)) + x24*(cell[7] + x107*x136 + V{0.0277777777777778}) + x37);
auto x8 = -x114*(-x144*x76 + x20*(x105 - x146 + x147)) - x24*(cell[8] + x107*x144 + V{0.0277777777777778}) - x45;
auto x9 = -(x114*(-x150*x76 + x20*(x137 + x147 - V{4.5}*x152*x152)) + x24*(cell[9] + x107*x150 + V{0.0277777777777778}) + x48);
auto x10 = cell[10] - x24*(cell[10] - x157*x77 + V{0.0555555555555556}) + x89*(-x157*x76 + x20*(x161 + x91));
auto x11 = cell[11] - x24*(cell[11] - x164*x77 + V{0.0555555555555556}) + x89*(-x164*x76 + x20*(x100 + x159 + x166));
auto x12 = cell[12] - x24*(cell[12] - x168*x77 + V{0.0555555555555556}) + x89*(-x168*x76 + x20*(x106 + x158 + x169 + V{1}));
auto x13 = cell[13] + x114*(-x170*x76 + x20*(x116 + x161 + x166)) - x24*(cell[13] - x107*x170 + V{0.0277777777777778});
auto x14 = -(x114*(-x172*x76 + x20*(x147 + x173 - V{4.5}*x124*x124)) + x24*(cell[14] + x107*x172 + V{0.0277777777777778}) + x51);
auto x15 = cell[15] + x114*(-x174*x76 + x20*(x131 + x161 + x169)) - x24*(cell[15] - x107*x174 + V{0.0277777777777778});
auto x16 = -(x114*(-x176*x76 + x20*(x173 + x177 - V{4.5}*x139*x139)) + x24*(cell[16] + x107*x176 + V{0.0277777777777778}) + x61);
auto x17 = cell[17] + x114*(-x178*x76 + x20*(x105 + x146 + x160 + x166)) - x24*(cell[17] - x107*x178 + V{0.0277777777777778});
auto x18 = -(x114*(-x179*x76 + x20*(x123 + x177 - V{4.5}*x151*x151)) + x24*(cell[18] + x107*x179 + V{0.0277777777777778}) + x58);
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
return { x28, V{1}*x29*(x43 + x54 + x64) };
}
};

}

}
