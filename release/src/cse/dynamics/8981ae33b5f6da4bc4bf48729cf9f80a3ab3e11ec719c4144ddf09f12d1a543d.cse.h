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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x25 = parameters.template get<descriptors::OMEGA>();
auto x23 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x24 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x26 = x25 + V{-1};
auto x27 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x28 = x27 + V{1};
auto x29 = x22*x22;
auto x30 = V{1.5}*x29;
auto x31 = x23*x23;
auto x32 = V{1.5}*x31;
auto x33 = x24*x24;
auto x34 = V{1.5}*x33;
auto x35 = x32 + x34 + V{-1};
auto x36 = x30 + x35;
auto x37 = V{0.5}*x25 + V{-1};
auto x38 = x19*x22;
auto x39 = x20*x23;
auto x40 = x21*x24;
auto x41 = x27 + V{1};
auto x42 = V{0.0555555555555556}*x25;
auto x43 = V{3}*x22;
auto x44 = V{3}*x29;
auto x45 = V{6}*x22;
auto x46 = x45 + V{-3};
auto x47 = V{0.0555555555555556}*x19;
auto x48 = V{0.166666666666667}*x39;
auto x49 = V{0.166666666666667}*x40;
auto x50 = x48 + x49;
auto x51 = V{3}*x23;
auto x52 = V{3}*x31;
auto x53 = x30 + V{-1};
auto x54 = V{6}*x23;
auto x55 = x54 + V{-3};
auto x56 = V{0.0555555555555556}*x20;
auto x57 = V{0.166666666666667}*x38;
auto x58 = x49 + x57;
auto x59 = V{3}*x24;
auto x60 = V{3}*x33;
auto x61 = V{6}*x24;
auto x62 = x61 + V{-3};
auto x63 = V{0.0555555555555556}*x21;
auto x64 = x48 + x57;
auto x65 = V{0.0277777777777778}*x25;
auto x66 = x22 + x23;
auto x67 = V{4.5}*(x66*x66);
auto x68 = x36 + x43;
auto x69 = V{9}*x23;
auto x70 = V{0.0277777777777778}*x19;
auto x71 = V{9}*x22;
auto x72 = V{0.0277777777777778}*x20;
auto x73 = V{0.0833333333333333}*x40;
auto x74 = -x73;
auto x75 = x37*x41;
auto x76 = -x51;
auto x77 = x22 - x23;
auto x78 = -x77;
auto x79 = V{3} - x45;
auto x80 = -x71;
auto x81 = x54 + V{3};
auto x82 = x22 + x24;
auto x83 = V{4.5}*(x82*x82);
auto x84 = V{9}*x24;
auto x85 = V{0.0277777777777778}*x21;
auto x86 = V{0.0833333333333333}*x39;
auto x87 = -x86;
auto x88 = -x59;
auto x89 = -x24;
auto x90 = x22 + x89;
auto x91 = -x90;
auto x92 = x61 + V{3};
auto x93 = x23 + x24;
auto x94 = V{4.5}*(x93*x93);
auto x95 = x36 + x51;
auto x96 = V{0.0833333333333333}*x38;
auto x97 = -x96;
auto x98 = x23 + x89;
auto x99 = -x98;
auto x100 = V{3} - x54;
auto x101 = -x69;
auto x102 = -x32;
auto x103 = V{1} - x34;
auto x104 = x102 + x103;
auto x105 = x104 + x43;
auto x106 = x45 + V{3};
auto x107 = -x30;
auto x108 = x107 + x51;
auto x109 = x107 + x59;
auto x110 = -x43;
auto x111 = x36 + x59;
auto x112 = V{3} - x61;
auto x113 = -x84;
auto x0 = -cell[0]*x26 - V{0.333333333333333}*x25*(x28*x36 + V{1}) + V{1}*x37*x41*(x38 + x39 + x40);
auto x1 = -cell[1]*x26 + x37*x41*(-x46*x47 + x50) - x42*(x28*(x35 + x43 - x44) + V{1});
auto x2 = -cell[2]*x26 + x37*x41*(-x55*x56 + x58) - x42*(x28*(x34 + x51 - x52 + x53) + V{1});
auto x3 = -cell[3]*x26 + x37*x41*(-x62*x63 + x64) - x42*(x28*(x32 + x53 + x59 - x60) + V{1});
auto x4 = -cell[4]*x26 - x65*(x28*(x51 - x67 + x68) + V{1}) - x75*(x70*(x46 + x69) + x72*(x55 + x71) + x74);
auto x5 = -(cell[5]*x26 + x65*(x28*(x68 + x76 - V{4.5}*x78*x78) + V{1}) + x75*(V{0.0277777777777778}*x20*(x80 + x81) - x70*(x69 + x79) - x73));
auto x6 = -cell[6]*x26 - x65*(x28*(x59 + x68 - x83) + V{1}) - x75*(x70*(x46 + x84) + x85*(x62 + x71) + x87);
auto x7 = -(cell[7]*x26 + x65*(x28*(x68 + x88 - V{4.5}*x91*x91) + V{1}) + x75*(V{0.0277777777777778}*x21*(x80 + x92) - x70*(x79 + x84) - x86));
auto x8 = -cell[8]*x26 - x65*(x28*(x59 - x94 + x95) + V{1}) - x75*(x72*(x55 + x84) + x85*(x62 + x69) + x97);
auto x9 = -(cell[9]*x26 + x65*(x28*(x88 + x95 - V{4.5}*x99*x99) + V{1}) + x75*(V{0.0277777777777778}*x21*(x101 + x92) - x72*(x100 + x84) - x96));
auto x10 = -cell[10]*x26 + x42*(x28*(x105 + x44) + V{-1}) + x75*(-x106*x47 + x50);
auto x11 = -cell[11]*x26 + x42*(x28*(x103 + x108 + x52) + V{-1}) + x75*(-x56*x81 + x58);
auto x12 = -cell[12]*x26 + x42*(x28*(x102 + x109 + x60 + V{1}) + V{-1}) + x75*(-x63*x92 + x64);
auto x13 = -cell[13]*x26 + V{0.0277777777777778}*x25*(x28*(x105 + x108 + x67) + V{-1}) - x75*(x70*(x106 + x69) + x72*(x71 + x81) + x74);
auto x14 = -(cell[14]*x26 + x65*(x28*(x110 + x95 - V{4.5}*x77*x77) + V{1}) + x75*(V{0.0277777777777778}*x19*(x101 + x106) - x72*(x100 + x71) - x73));
auto x15 = -cell[15]*x26 + V{0.0277777777777778}*x25*(x28*(x105 + x109 + x83) + V{-1}) - x75*(x70*(x106 + x84) + x85*(x71 + x92) + x87);
auto x16 = -(cell[16]*x26 + x65*(x28*(x110 + x111 - V{4.5}*x90*x90) + V{1}) + x75*(V{0.0277777777777778}*x19*(x106 + x113) - x85*(x112 + x71) - x86));
auto x17 = -cell[17]*x26 + V{0.0277777777777778}*x25*(x28*(x104 + x108 + x59 + x94) + V{-1}) - x75*(x72*(x81 + x84) + x85*(x69 + x92) + x97);
auto x18 = -(cell[18]*x26 + x65*(x28*(x111 + x76 - V{4.5}*x98*x98) + V{1}) + x75*(V{0.0277777777777778}*x20*(x113 + x81) - x85*(x112 + x69) - x96));
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
return { x41, x29 + x31 + x33 };
}
};

}

}
