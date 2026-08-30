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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[10] + cell[14];
auto x22 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x21 + x22;
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
auto x34 = -cell[5] + x21;
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
auto x50 = -cell[16] + x22;
auto x51 = x46 + x49 + x50;
auto x52 = x51*x51;
auto x53 = x27*x52;
auto x54 = x45 + x53 + V{-1};
auto x55 = x37 + x54;
auto x56 = V{0.0555555555555556}*x19;
auto x57 = V{1} / (x25);
auto x58 = V{3}*cell[14];
auto x59 = V{3}*cell[16];
auto x60 = V{3}*cell[5];
auto x61 = V{3}*cell[7];
auto x62 = V{3}*cell[13] - V{3}*cell[4];
auto x63 = V{3}*cell[15] - V{3}*cell[6];
auto x64 = x57*(V{3}*cell[10] - V{3}*cell[1] + x58 + x59 - x60 - x61 + x62 + x63);
auto x65 = V{3}*x26;
auto x66 = x36*x65;
auto x67 = V{3}*cell[18];
auto x68 = V{3}*cell[9];
auto x69 = V{3}*cell[17] - V{3}*cell[8];
auto x70 = x57*(V{3}*cell[11] - V{3}*cell[2] - x58 + x60 + x62 + x67 - x68 + x69);
auto x71 = x44*x65;
auto x72 = x37 + V{-1};
auto x73 = x57*(V{3}*cell[12] - V{3}*cell[3] - x59 + x61 + x63 - x67 + x68 + x69);
auto x74 = x52*x65;
auto x75 = V{0.0277777777777778}*x19;
auto x76 = V{4.5}*x26;
auto x77 = cell[10] + x33;
auto x78 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x40 + x41 + x46 + x77;
auto x79 = x76*(x78*x78);
auto x80 = x55 + x64;
auto x81 = -x70;
auto x82 = -cell[17] + cell[8];
auto x83 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x29 + x48 + x77 + x82;
auto x84 = -x83;
auto x85 = x31 + x34;
auto x86 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x39 + x49 + x85;
auto x87 = x76*(x86*x86);
auto x88 = -x73;
auto x89 = -cell[12] + cell[3] + x28;
auto x90 = V{2}*cell[16] - V{2}*cell[7] + x40 + x82 + x85 + x89;
auto x91 = -x90;
auto x92 = V{2}*cell[17] - V{2}*cell[8] + x30 + x42 + x47 + x50;
auto x93 = x76*(x92*x92);
auto x94 = x55 + x70;
auto x95 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x32 + x42 + x89;
auto x96 = -x95;
auto x97 = -x45;
auto x98 = V{1} - x53;
auto x99 = x97 + x98;
auto x100 = x64 + x99;
auto x101 = -x37;
auto x102 = x101 + x70;
auto x103 = x101 + x73;
auto x104 = -x64;
auto x105 = x55 + x73;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x24*x55 + V{1});
auto x1 = -cell[1]*x20 - x56*(x24*(x54 + x64 - x66) + V{1});
auto x2 = -cell[2]*x20 - x56*(x24*(x53 + x70 - x71 + x72) + V{1});
auto x3 = -cell[3]*x20 - x56*(x24*(x45 + x72 + x73 - x74) + V{1});
auto x4 = -cell[4]*x20 - x75*(x24*(x70 - x79 + x80) + V{1});
auto x5 = -(cell[5]*x20 + x75*(x24*(-x76*x84*x84 + x80 + x81) + V{1}));
auto x6 = -cell[6]*x20 - x75*(x24*(x73 + x80 - x87) + V{1});
auto x7 = -(cell[7]*x20 + x75*(x24*(-x76*x91*x91 + x80 + x88) + V{1}));
auto x8 = -cell[8]*x20 - x75*(x24*(x73 - x93 + x94) + V{1});
auto x9 = -(cell[9]*x20 + x75*(x24*(-x76*x96*x96 + x88 + x94) + V{1}));
auto x10 = -cell[10]*x20 + x56*(x24*(x100 + x66) + V{-1});
auto x11 = -cell[11]*x20 + x56*(x24*(x102 + x71 + x98) + V{-1});
auto x12 = -cell[12]*x20 + x56*(x24*(x103 + x74 + x97 + V{1}) + V{-1});
auto x13 = -cell[13]*x20 + x75*(x24*(x100 + x102 + x79) + V{-1});
auto x14 = -(cell[14]*x20 + x75*(x24*(x104 - x76*x83*x83 + x94) + V{1}));
auto x15 = -cell[15]*x20 + x75*(x24*(x100 + x103 + x87) + V{-1});
auto x16 = -(cell[16]*x20 + x75*(x24*(x104 + x105 - x76*x90*x90) + V{1}));
auto x17 = -cell[17]*x20 + x75*(x24*(x102 + x73 + x93 + x99) + V{-1});
auto x18 = -(cell[18]*x20 + x75*(x24*(x105 - x76*x95*x95 + x81) + V{1}));
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
return { x25, V{1}*x26*(x36 + x44 + x52) };
}
};

}

}
