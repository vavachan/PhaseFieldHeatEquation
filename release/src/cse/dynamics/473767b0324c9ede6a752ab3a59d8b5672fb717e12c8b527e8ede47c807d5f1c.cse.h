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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x20 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x21 = x20 + V{-1};
auto x22 = cell[10] + cell[14];
auto x23 = cell[12] + cell[7];
auto x24 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x22 + x23;
auto x25 = x24 + V{1};
auto x26 = x24 + V{1};
auto x27 = V{1} / ((x26)*(x26));
auto x28 = V{1.5}*x27;
auto x29 = x19*x19;
auto x30 = cell[13] - cell[4];
auto x31 = cell[15] - cell[6];
auto x32 = x30 + x31;
auto x33 = -cell[1];
auto x34 = cell[16] - cell[7];
auto x35 = x33 + x34;
auto x36 = -cell[5] + x22;
auto x37 = x32 + x35 + x36;
auto x38 = x37*x37;
auto x39 = x29*x38;
auto x40 = x28*x39;
auto x41 = cell[17] - cell[8];
auto x42 = x30 + x41;
auto x43 = cell[18] - cell[9];
auto x44 = -cell[2];
auto x45 = cell[11] - cell[14] + cell[5] + x44;
auto x46 = x42 + x43 + x45;
auto x47 = x46*x46;
auto x48 = x29*x47;
auto x49 = x28*x48;
auto x50 = x31 + x41;
auto x51 = -cell[3];
auto x52 = -cell[18] + cell[9];
auto x53 = x51 + x52;
auto x54 = -cell[16] + x23;
auto x55 = x50 + x53 + x54;
auto x56 = x55*x55;
auto x57 = x29*x56;
auto x58 = x28*x57;
auto x59 = x49 + x58 + V{-1};
auto x60 = x40 + x59;
auto x61 = V{0.0555555555555556}*x20;
auto x62 = V{3}*x19/x26;
auto x63 = x37*x62;
auto x64 = V{3}*x27;
auto x65 = x39*x64;
auto x66 = x46*x62;
auto x67 = x48*x64;
auto x68 = x40 + V{-1};
auto x69 = x55*x62;
auto x70 = x57*x64;
auto x71 = V{0.0277777777777778}*x20;
auto x72 = V{4.5}*x27;
auto x73 = cell[10] + x35;
auto x74 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x43 + x44 + x50 + x73;
auto x75 = x29*x72*(x74*x74);
auto x76 = x60 + x63;
auto x77 = -x66;
auto x78 = -cell[17] + cell[8];
auto x79 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x31 + x52 + x73 + x78;
auto x80 = -x79;
auto x81 = x33 + x36;
auto x82 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x42 + x53 + x81;
auto x83 = x29*x72*(x82*x82);
auto x84 = -x69;
auto x85 = -cell[12] + cell[3] + x30;
auto x86 = V{2}*cell[16] - V{2}*cell[7] + x43 + x78 + x81 + x85;
auto x87 = -x86;
auto x88 = V{2}*cell[17] - V{2}*cell[8] + x32 + x45 + x51 + x54;
auto x89 = x29*x72*(x88*x88);
auto x90 = x60 + x66;
auto x91 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x34 + x45 + x85;
auto x92 = -x91;
auto x93 = -x49;
auto x94 = V{1} - x58;
auto x95 = x93 + x94;
auto x96 = x63 + x95;
auto x97 = -x40;
auto x98 = x66 + x97;
auto x99 = x69 + x97;
auto x100 = -x63;
auto x101 = x60 + x69;
auto x0 = -cell[0]*x21 - V{0.333333333333333}*x20*(x25*x60 + V{1});
auto x1 = -cell[1]*x21 - x61*(x25*(x59 + x63 - x65) + V{1});
auto x2 = -cell[2]*x21 - x61*(x25*(x58 + x66 - x67 + x68) + V{1});
auto x3 = -cell[3]*x21 - x61*(x25*(x49 + x68 + x69 - x70) + V{1});
auto x4 = -cell[4]*x21 - x71*(x25*(x66 - x75 + x76) + V{1});
auto x5 = -(cell[5]*x21 + x71*(x25*(-x29*x72*x80*x80 + x76 + x77) + V{1}));
auto x6 = -cell[6]*x21 - x71*(x25*(x69 + x76 - x83) + V{1});
auto x7 = -(cell[7]*x21 + x71*(x25*(-x29*x72*x87*x87 + x76 + x84) + V{1}));
auto x8 = -cell[8]*x21 - x71*(x25*(x69 - x89 + x90) + V{1});
auto x9 = -(cell[9]*x21 + x71*(x25*(-x29*x72*x92*x92 + x84 + x90) + V{1}));
auto x10 = -cell[10]*x21 + x61*(x25*(x65 + x96) + V{-1});
auto x11 = -cell[11]*x21 + x61*(x25*(x67 + x94 + x98) + V{-1});
auto x12 = -cell[12]*x21 + x61*(x25*(x70 + x93 + x99 + V{1}) + V{-1});
auto x13 = -cell[13]*x21 + x71*(x25*(x75 + x96 + x98) + V{-1});
auto x14 = -(cell[14]*x21 + x71*(x25*(x100 - x29*x72*x79*x79 + x90) + V{1}));
auto x15 = -cell[15]*x21 + x71*(x25*(x83 + x96 + x99) + V{-1});
auto x16 = -(cell[16]*x21 + x71*(x25*(x100 + x101 - x29*x72*x86*x86) + V{1}));
auto x17 = -cell[17]*x21 + x71*(x25*(x69 + x89 + x95 + x98) + V{-1});
auto x18 = -(cell[18]*x21 + x71*(x25*(x101 - x29*x72*x91*x91 + x77) + V{1}));
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
return { x26, V{1}*x27*x29*(x38 + x47 + x56) };
}
};

}

}
