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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::fsi::HLBM>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x9 = cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x11 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x12 = parameters.template get<descriptors::OMEGA>();
auto x13 = x12 + V{-1};
auto x14 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x15 = x14 + V{1};
auto x16 = x14 + V{1};
auto x17 = V{1} / ((x16)*(x16));
auto x18 = V{1.5}*x17;
auto x19 = cell[1] - cell[5];
auto x20 = -cell[4] + cell[8];
auto x21 = -cell[3] + cell[7] + x19 + x20;
auto x22 = x21*x21;
auto x23 = x18*x22;
auto x24 = cell[2] - cell[6];
auto x25 = cell[3] - cell[7] + x19 + x24;
auto x26 = -x25;
auto x27 = x18*(x26*x26) + V{-1};
auto x28 = x23 + x27;
auto x29 = V{1} / (x16);
auto x30 = V{1}*cell[3];
auto x31 = V{1}*cell[1] - V{1}*cell[5];
auto x32 = x29*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x30 + x31);
auto x33 = x9 + V{-1};
auto x34 = -x33;
auto x35 = -x32;
auto x36 = x11 + x35;
auto x37 = x34*x36;
auto x38 = x32 + x37;
auto x39 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x30 + x31;
auto x40 = -x29*x39;
auto x41 = x34*(x10 - x40);
auto x42 = x40 + x41;
auto x43 = V{-1} + V{1.5}*(x38*x38) + V{1.5}*(x42*x42);
auto x44 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x45 = V{3}*cell[3];
auto x46 = V{3}*cell[7];
auto x47 = V{3}*cell[1] - V{3}*cell[5];
auto x48 = V{3}*cell[2] - V{3}*cell[6] + x45 - x46 + x47;
auto x49 = x29*x48;
auto x50 = V{4.5}*x17;
auto x51 = V{2}*cell[1] - V{2}*cell[5] + x20 + x24;
auto x52 = -x51;
auto x53 = V{1} - x23;
auto x54 = x29*(-V{3}*cell[4] + V{3}*cell[8] - x45 + x46 + x47);
auto x55 = x25*x25;
auto x56 = x18*x55;
auto x57 = x54 - x56;
auto x58 = x53 + x57;
auto x59 = x49 + x50*(x52*x52) + x58;
auto x60 = x21*x29;
auto x61 = -x60;
auto x62 = x26*x29 + x41;
auto x63 = x34*x36 - x61 - x62;
auto x64 = -x54;
auto x65 = V{3}*x37;
auto x66 = -x29*x48;
auto x67 = V{3}*x41 + x43 + x66;
auto x68 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x69 = V{3}*x17;
auto x70 = x53 + x55*x69;
auto x71 = x49 + x70;
auto x72 = V{0.0277777777777778}*x12;
auto x73 = -V{2}*cell[3] - cell[4] + V{2}*cell[7] + cell[8] - x24;
auto x74 = x50*(x73*x73);
auto x75 = x28 + x54 + x66 - x74;
auto x76 = x37 + x60;
auto x77 = x62 + x76;
auto x78 = x54 + x65;
auto x79 = V{0.111111111111111}*x12;
auto x80 = x22*x69;
auto x81 = x27 + x54 - x80;
auto x82 = x49 + V{-1};
auto x83 = x54 + x82;
auto x84 = -(-V{4.5}*x17*x51*x51 + x23 + x56 + x83);
auto x85 = x33*x36;
auto x86 = V{3}*x85;
auto x87 = -x86;
auto x88 = x29*x39;
auto x89 = x33*(x10 + x88);
auto x90 = x25*x29 + x89;
auto x91 = -x60 + x85 - x90;
auto x92 = x88 + x89;
auto x93 = V{1.5}*(x92*x92);
auto x94 = -x35 - x85;
auto x95 = V{1.5}*(x94*x94);
auto x96 = V{3}*x89 + x93 + x95;
auto x97 = -x49;
auto x98 = x70 + x97;
auto x99 = x82 + x96;
auto x100 = x58 + x74 + x97;
auto x101 = x61 + x85;
auto x102 = -x101 - x90;
auto x103 = x57 + x80 + V{1};
auto x104 = -x101;
auto x0 = -cell[0]*x13 - V{0.444444444444444}*x12*(x15*x28 + V{1}) + V{0.444444444444444}*x15*x28 - x43*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444});
auto x1 = -(cell[1]*x13 - V{0.0277777777777778}*x12*(x15*x59 + V{-1}) + x44*x59 + x44*(x64 - x65 + x67 - V{4.5}*x63*x63));
auto x2 = -(cell[2]*x13 - V{0.111111111111111}*x12*(x15*x71 + V{-1}) + x68*x71 + x68*(x67 - V{4.5}*x62*x62));
auto x3 = -(cell[3]*x13 - V{0.0277777777777778}*x15*x75 + x44*(x67 + x78 - V{4.5}*x77*x77) + x72*(x15*x75 + V{1}));
auto x4 = -(cell[4]*x13 - V{0.111111111111111}*x15*x81 + x68*(x43 + x78 - V{4.5}*x76*x76) + x79*(x15*x81 + V{1}));
auto x5 = -cell[5]*x13 - x44*x84 - x44*(x83 + x87 + x96 - V{4.5}*x91*x91) + x72*(x15*x84 + V{-1});
auto x6 = -cell[6]*x13 - x68*x98 - x68*(x99 - V{4.5}*x90*x90) + x79*(x15*x98 + V{-1});
auto x7 = -cell[7]*x13 - x100*x44 - x44*(x64 + x86 + x99 - V{4.5}*x102*x102) + x72*(x100*x15 + V{-1});
auto x8 = -cell[8]*x13 - x103*x68 + x68*(x54 + x87 - x93 - x95 + V{1} + V{4.5}*(x104*x104)) + x79*(x103*x15 + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x16, V{1}*x17*(x22 + x55) };
}
};

}

}
