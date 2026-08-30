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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x20 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell[15] + cell[17];
auto x22 = cell[12] + x21;
auto x23 = cell[11] + cell[18];
auto x24 = cell[10] + cell[14] + cell[16];
auto x25 = cell[2] + cell[8] + cell[9];
auto x26 = cell[13] + cell[3];
auto x27 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + cell[7] + x22 + x23 + x24 + x25 + x26 + V{1};
auto x28 = V{1} / (x27);
auto x29 = V{0.333333333333333}*cell[13];
auto x30 = V{0.333333333333333}*cell[14];
auto x31 = V{0.333333333333333}*cell[4];
auto x32 = V{0.333333333333333}*cell[5];
auto x33 = V{1}*x28;
auto x34 = -cell[18];
auto x35 = -cell[3];
auto x36 = -cell[8];
auto x37 = cell[9] + x36;
auto x38 = x34 + x35 + x37;
auto x39 = -cell[6];
auto x40 = cell[7] + x39;
auto x41 = -cell[16] + x40;
auto x42 = x22 + x38 + x41;
auto x43 = x42*x42;
auto x44 = V{0.333333333333333}*cell[0];
auto x45 = V{0.333333333333333}*cell[10];
auto x46 = V{0.333333333333333}*cell[1];
auto x47 = -V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x44 + x45 + x46;
auto x48 = V{0.333333333333333}*cell[11];
auto x49 = V{0.333333333333333}*cell[2];
auto x50 = -V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7] + x48 + x49;
auto x51 = -V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x29 + x30 + x31 + x32 + x33*x43 + x47 + x50;
auto x52 = V{0.333333333333333}*cell[15];
auto x53 = V{0.333333333333333}*cell[16];
auto x54 = V{0.333333333333333}*cell[6];
auto x55 = V{0.333333333333333}*cell[7];
auto x56 = cell[13] + cell[17];
auto x57 = -cell[2];
auto x58 = -cell[9];
auto x59 = x23 + x36 + x57 + x58;
auto x60 = -cell[4];
auto x61 = cell[5] + x60;
auto x62 = -cell[14] + x61;
auto x63 = x56 + x59 + x62;
auto x64 = x63*x63;
auto x65 = V{0.333333333333333}*cell[12];
auto x66 = V{0.333333333333333}*cell[3];
auto x67 = -V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5] + x65 + x66;
auto x68 = -V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x33*x64 + x47 + x52 + x53 + x54 + x55 + x67;
auto x69 = V{0.333333333333333}*cell[17];
auto x70 = V{0.333333333333333}*cell[18];
auto x71 = V{0.333333333333333}*cell[8];
auto x72 = V{0.333333333333333}*cell[9];
auto x73 = cell[13] + cell[15];
auto x74 = -cell[1];
auto x75 = -cell[7];
auto x76 = x39 + x74 + x75;
auto x77 = -cell[5] + x60;
auto x78 = x24 + x73 + x76 + x77;
auto x79 = x78*x78;
auto x80 = -V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x33*x79 + x44 + x50 + x67 + x69 + x70 + x71 + x72;
auto x81 = x28*x78;
auto x82 = -cell[13] + cell[14] + x61 + x63*x81;
auto x83 = -cell[15] + cell[16];
auto x84 = x40 + x42*x81 + x83;
auto x85 = -cell[17];
auto x86 = cell[18] + x85;
auto x87 = x28*x42*x63 + x37 + x86;
auto x88 = V{1} / (V{3.00000046417339}*util::sqrt(x28*(x20*x20)*util::sqrt(V{0.5}*(x51*x51) + V{0.5}*(x68*x68) + V{0.5}*(x80*x80) + x82*x82 + x84*x84 + x87*x87) + V{0.0277777691819762}/((x19)*(x19))) + V{0.5}/x19);
auto x89 = V{1} / ((x27)*(x27));
auto x90 = V{1.5}*x89;
auto x91 = x79*x90;
auto x92 = x64*x90;
auto x93 = x43*x90;
auto x94 = x92 + x93 + V{-1};
auto x95 = x91 + x94;
auto x96 = V{1} - x88;
auto x97 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x98 = V{3}*cell[14];
auto x99 = V{3}*cell[16];
auto x100 = V{3}*cell[5];
auto x101 = V{3}*cell[7];
auto x102 = V{3}*cell[13] - V{3}*cell[4];
auto x103 = V{3}*cell[15] - V{3}*cell[6];
auto x104 = x28*(V{3}*cell[10] - V{3}*cell[1] - x100 - x101 + x102 + x103 + x98 + x99);
auto x105 = V{3}*x89;
auto x106 = x105*x79;
auto x107 = V{3}*cell[18];
auto x108 = V{3}*cell[9];
auto x109 = V{3}*cell[17] - V{3}*cell[8];
auto x110 = x28*(V{3}*cell[11] - V{3}*cell[2] + x100 + x102 + x107 - x108 + x109 - x98);
auto x111 = x105*x64;
auto x112 = x91 + V{-1};
auto x113 = x28*(V{3}*cell[12] - V{3}*cell[3] + x101 + x103 - x107 + x108 + x109 - x99);
auto x114 = x105*x43;
auto x115 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x116 = V{4.5}*x89;
auto x117 = cell[10] + cell[16] + x76;
auto x118 = V{2}*cell[13] - V{2}*cell[4] + x117 + x21 + x59;
auto x119 = x116*(x118*x118);
auto x120 = x104 + x95;
auto x121 = -x110;
auto x122 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x117 + x25 + x34 + x85;
auto x123 = -x122;
auto x124 = cell[10] + cell[14] + x74 + x77;
auto x125 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x124 + x38 + x56;
auto x126 = x116*(x125*x125);
auto x127 = -x113;
auto x128 = -cell[12] + x26;
auto x129 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x124 + x128 + x58 + x86;
auto x130 = -x129;
auto x131 = cell[11] + x57 + x62;
auto x132 = cell[12] + V{2}*cell[17] - V{2}*cell[8] + x131 + x35 + x41 + x73;
auto x133 = x116*(x132*x132);
auto x134 = x110 + x95;
auto x135 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x128 + x131 + x75 + x83;
auto x136 = -x135;
auto x137 = -x92;
auto x138 = V{1} - x93;
auto x139 = x137 + x138;
auto x140 = x104 + x139;
auto x141 = -x91;
auto x142 = x110 + x141;
auto x143 = x113 + x141;
auto x144 = -x104;
auto x145 = x113 + x95;
auto x0 = V{1}*cell[0]*x96 - x88*(x95*(x29 + x30 + x31 + x32 + x44 + x45 + x46 + x48 + x49 + x52 + x53 + x54 + x55 + x65 + x66 + x69 + x70 + x71 + x72 + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = V{1}*cell[1]*x96 - x88*(x97*(x104 - x106 + x94) + V{0.0555555555555556});
auto x2 = V{1}*cell[2]*x96 - x88*(x97*(x110 - x111 + x112 + x93) + V{0.0555555555555556});
auto x3 = V{1}*cell[3]*x96 - x88*(x97*(x112 + x113 - x114 + x92) + V{0.0555555555555556});
auto x4 = V{1}*cell[4]*x96 - x88*(x115*(x110 - x119 + x120) + V{0.0277777777777778});
auto x5 = V{1}*cell[5]*x96 - x88*(x115*(-x116*x123*x123 + x120 + x121) + V{0.0277777777777778});
auto x6 = V{1}*cell[6]*x96 - x88*(x115*(x113 + x120 - x126) + V{0.0277777777777778});
auto x7 = V{1}*cell[7]*x96 - x88*(x115*(-x116*x130*x130 + x120 + x127) + V{0.0277777777777778});
auto x8 = V{1}*cell[8]*x96 - x88*(x115*(x113 - x133 + x134) + V{0.0277777777777778});
auto x9 = V{1}*cell[9]*x96 - x88*(x115*(-x116*x136*x136 + x127 + x134) + V{0.0277777777777778});
auto x10 = V{1}*cell[10]*x96 + x88*(x97*(x106 + x140) + V{-0.0555555555555556});
auto x11 = V{1}*cell[11]*x96 + x88*(x97*(x111 + x138 + x142) + V{-0.0555555555555556});
auto x12 = V{1}*cell[12]*x96 + x88*(x97*(x114 + x137 + x143 + V{1}) + V{-0.0555555555555556});
auto x13 = V{1}*cell[13]*x96 + x88*(x115*(x119 + x140 + x142) + V{-0.0277777777777778});
auto x14 = V{1}*cell[14]*x96 - x88*(x115*(-x116*x122*x122 + x134 + x144) + V{0.0277777777777778});
auto x15 = V{1}*cell[15]*x96 + x88*(x115*(x126 + x140 + x143) + V{-0.0277777777777778});
auto x16 = V{1}*cell[16]*x96 - x88*(x115*(-x116*x129*x129 + x144 + x145) + V{0.0277777777777778});
auto x17 = V{1}*cell[17]*x96 + x88*(x115*(x113 + x133 + x139 + x142) + V{-0.0277777777777778});
auto x18 = V{1}*cell[18]*x96 - x88*(x115*(-x116*x135*x135 + x121 + x145) + V{0.0277777777777778});
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
return { x27, V{1}*x89*(x43 + x64 + x79) };
}
};

}

}
