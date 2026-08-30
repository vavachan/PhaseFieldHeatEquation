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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x20 = parameters.template get<statistics::AVERAGE_RHO>();
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x21 = x19 + V{-1};
auto x22 = cell[10] + cell[14];
auto x23 = cell[12] + cell[7];
auto x24 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x22 + x23;
auto x25 = x24 + V{1};
auto x26 = V{1} / ((x25)*(x25));
auto x27 = V{1.5}*x26;
auto x28 = cell[13] - cell[4];
auto x29 = cell[15] - cell[6];
auto x30 = x28 + x29;
auto x31 = -cell[1];
auto x32 = cell[16] - cell[7];
auto x33 = x31 + x32;
auto x34 = -cell[5] + x22;
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
auto x50 = -cell[16] + x23;
auto x51 = x46 + x49 + x50;
auto x52 = x51*x51;
auto x53 = x27*x52;
auto x54 = x45 + x53 + V{-1};
auto x55 = x37 + x54;
auto x56 = V{1} / (x25);
auto x57 = -x56*(x20 + V{-1}) + V{1};
auto x58 = x24 + V{1};
auto x59 = x57*x58;
auto x60 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x61 = V{3}*cell[14];
auto x62 = V{3}*cell[16];
auto x63 = V{3}*cell[5];
auto x64 = V{3}*cell[7];
auto x65 = V{3}*cell[13] - V{3}*cell[4];
auto x66 = V{3}*cell[15] - V{3}*cell[6];
auto x67 = x56*(V{3}*cell[10] - V{3}*cell[1] + x61 + x62 - x63 - x64 + x65 + x66);
auto x68 = V{3}*x26;
auto x69 = x36*x68;
auto x70 = x54 + x67 - x69;
auto x71 = V{0.0555555555555556}*x59;
auto x72 = V{3}*cell[18];
auto x73 = V{3}*cell[9];
auto x74 = V{3}*cell[17] - V{3}*cell[8];
auto x75 = x56*(V{3}*cell[11] - V{3}*cell[2] - x61 + x63 + x65 + x72 - x73 + x74);
auto x76 = x44*x68;
auto x77 = x37 + V{-1};
auto x78 = x53 + x75 - x76 + x77;
auto x79 = x56*(V{3}*cell[12] - V{3}*cell[3] - x62 + x64 + x66 - x72 + x73 + x74);
auto x80 = x52*x68;
auto x81 = x45 + x77 + x79 - x80;
auto x82 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x83 = V{4.5}*x26;
auto x84 = cell[10] + x33;
auto x85 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x40 + x41 + x46 + x84;
auto x86 = x83*(x85*x85);
auto x87 = x55 + x67;
auto x88 = x75 - x86 + x87;
auto x89 = V{0.0277777777777778}*x59;
auto x90 = -x75;
auto x91 = -cell[17] + cell[8];
auto x92 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x29 + x48 + x84 + x91;
auto x93 = -x92;
auto x94 = -x83*x93*x93 + x87 + x90;
auto x95 = x31 + x34;
auto x96 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x39 + x49 + x95;
auto x97 = x83*(x96*x96);
auto x98 = x79 + x87 - x97;
auto x99 = -x79;
auto x100 = -cell[12] + cell[3] + x28;
auto x101 = V{2}*cell[16] - V{2}*cell[7] + x100 + x40 + x91 + x95;
auto x102 = -x101;
auto x103 = -x83*x102*x102 + x87 + x99;
auto x104 = V{2}*cell[17] - V{2}*cell[8] + x30 + x42 + x47 + x50;
auto x105 = x83*(x104*x104);
auto x106 = x55 + x75;
auto x107 = -x105 + x106 + x79;
auto x108 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x100 + x32 + x42;
auto x109 = -x108;
auto x110 = x106 - x83*x109*x109 + x99;
auto x111 = -x45;
auto x112 = V{1} - x53;
auto x113 = x111 + x112;
auto x114 = x113 + x67;
auto x115 = x114 + x69;
auto x116 = -x37;
auto x117 = x116 + x75;
auto x118 = x112 + x117 + x76;
auto x119 = x116 + x79;
auto x120 = x111 + x119 + x80 + V{1};
auto x121 = x114 + x117 + x86;
auto x122 = -x67;
auto x123 = x106 + x122 - x83*x92*x92;
auto x124 = x114 + x119 + x97;
auto x125 = x55 + x79;
auto x126 = x122 + x125 - x83*x101*x101;
auto x127 = x105 + x113 + x117 + x79;
auto x128 = x125 - x83*x108*x108 + x90;
auto x0 = -x21*(cell[0] + x55*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{0.333333333333333}) - V{0.333333333333333}*x55*x59 + V{-0.333333333333333};
auto x1 = -x21*(cell[1] + x60*x70 + V{0.0555555555555556}) - x70*x71 + V{-0.0555555555555556};
auto x2 = -x21*(cell[2] + x60*x78 + V{0.0555555555555556}) - x71*x78 + V{-0.0555555555555556};
auto x3 = -x21*(cell[3] + x60*x81 + V{0.0555555555555556}) - x71*x81 + V{-0.0555555555555556};
auto x4 = -x21*(cell[4] + x82*x88 + V{0.0277777777777778}) - x88*x89 + V{-0.0277777777777778};
auto x5 = -x21*(cell[5] + x82*x94 + V{0.0277777777777778}) - x89*x94 + V{-0.0277777777777778};
auto x6 = -x21*(cell[6] + x82*x98 + V{0.0277777777777778}) - x89*x98 + V{-0.0277777777777778};
auto x7 = -x103*x89 - x21*(cell[7] + x103*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x8 = -x107*x89 - x21*(cell[8] + x107*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x9 = -x110*x89 - x21*(cell[9] + x110*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x10 = V{0.0555555555555556}*x115*x57*x58 - x21*(cell[10] - x115*x60 + V{0.0555555555555556}) + V{-0.0555555555555556};
auto x11 = V{0.0555555555555556}*x118*x57*x58 - x21*(cell[11] - x118*x60 + V{0.0555555555555556}) + V{-0.0555555555555556};
auto x12 = V{0.0555555555555556}*x120*x57*x58 - x21*(cell[12] - x120*x60 + V{0.0555555555555556}) + V{-0.0555555555555556};
auto x13 = V{0.0277777777777778}*x121*x57*x58 - x21*(cell[13] - x121*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x14 = -x123*x89 - x21*(cell[14] + x123*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x15 = V{0.0277777777777778}*x124*x57*x58 - x21*(cell[15] - x124*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x16 = -x126*x89 - x21*(cell[16] + x126*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x17 = V{0.0277777777777778}*x127*x57*x58 - x21*(cell[17] - x127*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
auto x18 = -x128*x89 - x21*(cell[18] + x128*x82 + V{0.0277777777777778}) + V{-0.0277777777777778};
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
return { -x20 + x24 + V{2}, V{1}*x26*(x36 + x44 + x52) };
}
};

}

}
