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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x19 = x22 + V{-1};
auto x20 = cell[10] + cell[14];
auto x21 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x20 + x21;
auto x24 = x23 + V{1};
auto x25 = x23 + V{1};
auto x26 = V{1} / ((x25)*(x25));
auto x27 = V{1.5}*x26;
auto x28 = cell[13] - cell[4];
auto x29 = cell[15] - cell[6];
auto x30 = x28 + x29;
auto x31 = -cell[1];
auto x32 = cell[16] - cell[7];
auto x33 = x31 + x32;
auto x34 = -cell[5] + x20;
auto x35 = x30 + x33 + x34;
auto x36 = x35*x35;
auto x37 = x27*x36;
auto x38 = cell[17] - cell[8];
auto x39 = x28 + x38;
auto x40 = cell[18] - cell[9];
auto x41 = -cell[2];
auto x42 = cell[11] - cell[14] + cell[5] + x41;
auto x43 = x39 + x40 + x42;
auto x44 = x43*x43;
auto x45 = x27*x44;
auto x46 = x29 + x38;
auto x47 = -cell[3];
auto x48 = -cell[18] + cell[9];
auto x49 = x47 + x48;
auto x50 = -cell[16] + x21;
auto x51 = x46 + x49 + x50;
auto x52 = x51*x51;
auto x53 = x27*x52;
auto x54 = x45 + x53 + V{-1};
auto x55 = x37 + x54;
auto x56 = cell[0]*x19 + V{0.333333333333333}*x22*(x24*x55 + V{1});
auto x57 = cell[1]*x19;
auto x58 = V{0.0555555555555556}*x22;
auto x59 = V{1} / (x25);
auto x60 = V{3}*cell[14];
auto x61 = V{3}*cell[16];
auto x62 = V{3}*cell[5];
auto x63 = V{3}*cell[7];
auto x64 = V{3}*cell[13] - V{3}*cell[4];
auto x65 = V{3}*cell[15] - V{3}*cell[6];
auto x66 = x59*(V{3}*cell[10] - V{3}*cell[1] + x60 + x61 - x62 - x63 + x64 + x65);
auto x67 = V{3}*x26;
auto x68 = x36*x67;
auto x69 = x58*(x24*(x54 + x66 - x68) + V{1});
auto x70 = x57 + x69;
auto x71 = cell[2]*x19;
auto x72 = V{3}*cell[18];
auto x73 = V{3}*cell[9];
auto x74 = V{3}*cell[17] - V{3}*cell[8];
auto x75 = x59*(V{3}*cell[11] - V{3}*cell[2] - x60 + x62 + x64 + x72 - x73 + x74);
auto x76 = x44*x67;
auto x77 = x37 + V{-1};
auto x78 = x58*(x24*(x53 + x75 - x76 + x77) + V{1});
auto x79 = x71 + x78;
auto x80 = cell[3]*x19;
auto x81 = x59*(V{3}*cell[12] - V{3}*cell[3] - x61 + x63 + x65 - x72 + x73 + x74);
auto x82 = x52*x67;
auto x83 = x58*(x24*(x45 + x77 + x81 - x82) + V{1});
auto x84 = x80 + x83;
auto x85 = cell[4]*x19;
auto x86 = V{0.0277777777777778}*x22;
auto x87 = V{4.5}*x26;
auto x88 = cell[10] + x33;
auto x89 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x40 + x41 + x46 + x88;
auto x90 = x87*(x89*x89);
auto x91 = x55 + x66;
auto x92 = x86*(x24*(x75 - x90 + x91) + V{1});
auto x93 = x85 + x92;
auto x94 = cell[5]*x19;
auto x95 = -x75;
auto x96 = -cell[17] + cell[8];
auto x97 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x29 + x48 + x88 + x96;
auto x98 = -x97;
auto x99 = x86*(x24*(-x87*x98*x98 + x91 + x95) + V{1});
auto x100 = x94 + x99;
auto x101 = cell[6]*x19;
auto x102 = x31 + x34;
auto x103 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x102 + x39 + x49;
auto x104 = x87*(x103*x103);
auto x105 = x86*(x24*(-x104 + x81 + x91) + V{1});
auto x106 = x101 + x105;
auto x107 = cell[7]*x19;
auto x108 = -x81;
auto x109 = -cell[12] + cell[3] + x28;
auto x110 = V{2}*cell[16] - V{2}*cell[7] + x102 + x109 + x40 + x96;
auto x111 = -x110;
auto x112 = x86*(x24*(x108 - x87*x111*x111 + x91) + V{1});
auto x113 = x107 + x112;
auto x114 = cell[8]*x19;
auto x115 = V{2}*cell[17] - V{2}*cell[8] + x30 + x42 + x47 + x50;
auto x116 = x87*(x115*x115);
auto x117 = x55 + x75;
auto x118 = x86*(x24*(-x116 + x117 + x81) + V{1});
auto x119 = x114 + x118;
auto x120 = cell[9]*x19;
auto x121 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x109 + x32 + x42;
auto x122 = -x121;
auto x123 = x86*(x24*(x108 + x117 - x87*x122*x122) + V{1});
auto x124 = x120 + x123;
auto x125 = cell[10]*x19;
auto x126 = -x45;
auto x127 = V{1} - x53;
auto x128 = x126 + x127;
auto x129 = x128 + x66;
auto x130 = x58*(x24*(x129 + x68) + V{-1});
auto x131 = cell[11]*x19;
auto x132 = -x37;
auto x133 = x132 + x75;
auto x134 = x58*(x24*(x127 + x133 + x76) + V{-1});
auto x135 = cell[12]*x19;
auto x136 = x132 + x81;
auto x137 = x58*(x24*(x126 + x136 + x82 + V{1}) + V{-1});
auto x138 = cell[13]*x19;
auto x139 = x86*(x24*(x129 + x133 + x90) + V{-1});
auto x140 = cell[14]*x19;
auto x141 = -x66;
auto x142 = x86*(x24*(x117 + x141 - x87*x97*x97) + V{1});
auto x143 = x140 + x142;
auto x144 = cell[15]*x19;
auto x145 = x86*(x24*(x104 + x129 + x136) + V{-1});
auto x146 = cell[16]*x19;
auto x147 = x55 + x81;
auto x148 = x86*(x24*(x141 + x147 - x87*x110*x110) + V{1});
auto x149 = x146 + x148;
auto x150 = cell[17]*x19;
auto x151 = x86*(x24*(x116 + x128 + x133 + x81) + V{-1});
auto x152 = cell[18]*x19;
auto x153 = x86*(x24*(x147 - x87*x121*x121 + x95) + V{1});
auto x154 = x152 + x153;
auto x155 = V{1} / (-x100 - x106 - x113 - x119 - x124 - x125 + x130 - x131 + x134 - x135 + x137 - x138 + x139 - x143 - x144 + x145 - x149 - x150 + x151 - x154 - x56 - x70 - x79 - x84 - x93 + V{1});
auto x156 = V{1}*x94;
auto x157 = V{1}*x107;
auto x158 = V{1}*x140;
auto x159 = V{1}*x146;
auto x160 = -V{1}*cell[13]*x19 + x139 + V{1}*x85 + x92;
auto x161 = -V{1}*cell[15]*x19 + V{1}*x101 + x105 + x145;
auto x162 = V{1}*x120;
auto x163 = V{1}*x152;
auto x164 = -V{1}*cell[17]*x19 + V{1}*x114 + x118 + x151;
cell[0] = -x56;
cell[1] = -x70;
cell[2] = -x79;
cell[3] = -x84;
cell[4] = -x93;
cell[5] = -x100;
cell[6] = -x106;
cell[7] = -x113;
cell[8] = -x119;
cell[9] = -x124;
cell[10] = -x125 + x130;
cell[11] = -x131 + x134;
cell[12] = -x135 + x137;
cell[13] = -x138 + x139;
cell[14] = -x143;
cell[15] = -x144 + x145;
cell[16] = -x149;
cell[17] = -x150 + x151;
cell[18] = -x154;
cell.template getFieldPointer<descriptors::VELOCITY>()[0] = -x155*(V{1}*cell[10]*x19 - x112 - x130 + x142 + x148 - x156 - x157 + x158 + x159 - x160 - x161 - V{1}*x57 - x69 - x99);
cell.template getFieldPointer<descriptors::VELOCITY>()[1] = -x155*(V{1}*cell[11]*x19 - x123 - x134 - x142 + x153 + x156 - x158 - x160 - x162 + x163 - x164 - V{1}*x71 - x78 + x99);
cell.template getFieldPointer<descriptors::VELOCITY>()[2] = -x155*(V{1}*cell[12]*x19 + x112 + x123 - x137 - x148 - x153 + x157 - x159 - x161 + x162 - x163 - x164 - V{1}*x80 - x83);
return { x25, V{1}*x26*(x36 + x44 + x52) };
}
};

}

}
