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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::LocalSmagorinskyEffectiveOmega<collision::HRR> >, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x16 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x19 = cell.template getFieldComponent<descriptors::WMVELOCITY>(1);
auto x18 = cell.template getFieldComponent<descriptors::WMVELOCITY>(0);
auto x15 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x11 = cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x17 = cell.template getFieldComponent<descriptors::WMPOROSITY>(0);
auto x20 = parameters.template get<descriptors::OMEGA>();
auto x14 = cell.template getFieldComponent<descriptors::TENSOR>(2);
auto x22 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x12 = cell.template getFieldComponent<descriptors::TENSOR>(0);
auto x13 = cell.template getFieldComponent<descriptors::TENSOR>(1);
auto x9 = cell.template getFieldComponent<collision::HYBRID>(0);
auto x10 = V{0.5}/x20;
auto x21 = V{0.0277777691819762}/((x20)*(x20));
auto x23 = cell[1] + cell[7] + cell[8];
auto x24 = cell[0] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + x23;
auto x25 = x24 + V{1};
auto x26 = V{1} / (x25);
auto x27 = x22*x22;
auto x28 = -cell[3];
auto x29 = x11 + V{-1};
auto x30 = -x29;
auto x31 = cell[1] - cell[7];
auto x32 = -cell[5];
auto x33 = cell[3] + x32;
auto x34 = cell[2] - cell[6] + x31 + x33;
auto x35 = V{1}*x11*x26;
auto x36 = x15*x30 - x34*x35;
auto x37 = x17*x36;
auto x38 = x17 + V{-1};
auto x39 = -x38;
auto x40 = x18*x39;
auto x41 = x37 + x40;
auto x42 = x35*(-cell[4] + x23 + x28 + x32);
auto x43 = x16*x30 + x42;
auto x44 = x17*x43;
auto x45 = x19*x39;
auto x46 = x44 + x45;
auto x47 = x25*x41*x46;
auto x48 = cell[5] + x28 + x31 + x47;
auto x49 = V{1}*cell[1];
auto x50 = V{1}*cell[5];
auto x51 = V{1}*cell[3];
auto x52 = V{1}*cell[7];
auto x53 = x47 + x49 + x50 - x51 - x52;
auto x54 = x15*x29 + x34*x35;
auto x55 = x17*x54;
auto x56 = x18*x38;
auto x57 = x55 + x56;
auto x58 = x57*x57;
auto x59 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x60 = V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8] + x59;
auto x61 = -x25*x58 + x60;
auto x62 = -x16*x29 + x42;
auto x63 = x17*x62;
auto x64 = x19*x38;
auto x65 = x63 - x64;
auto x66 = x65*x65;
auto x67 = -V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8] + x59;
auto x68 = -x25*x66 + x67;
auto x69 = V{0.5}*(x61*x61) + V{0.5}*(x68*x68);
auto x70 = x10 + V{3.00000046417339}*util::sqrt(x21 + x26*x27*util::sqrt(x48*x53 + x69));
auto x71 = V{1} - 1/x70;
auto x72 = V{0.666666666666667}*x9;
auto x73 = x12*x25;
auto x74 = x9 + V{-1};
auto x75 = V{0.444444444444444}*x74;
auto x76 = x14*x25;
auto x77 = x41*x41;
auto x78 = x46*x46;
auto x79 = V{1.5}*x77 + V{1.5}*x78 + V{-1};
auto x80 = -x25*x78 + x67;
auto x81 = x80*x9;
auto x82 = -x74;
auto x83 = V{0.666666666666667}*x82;
auto x84 = -x25*x77 + x60;
auto x85 = -x53;
auto x86 = x10 + V{3.00000046417339}*util::sqrt(x21 + x26*x27*util::sqrt(-x48*x85 + V{0.5}*(x80*x80) + V{0.5}*(x84*x84)));
auto x87 = x76*x83*x86;
auto x88 = x81 - x87;
auto x89 = x88*(V{0.166667}*x17*x36 + V{0.166667}*x18*x39);
auto x90 = x85*x9;
auto x91 = x13*x25;
auto x92 = x83*x86*x91;
auto x93 = x90 - x92;
auto x94 = x93*(V{0.333334}*x17*x43 + V{0.333334}*x19*x39);
auto x95 = x84*x9;
auto x96 = x73*x83*x86;
auto x97 = x95 - x96;
auto x98 = x97*(V{0.166667}*x17*x43 + V{0.166667}*x19*x39);
auto x99 = x93*(V{0.333334}*x17*x36 + V{0.333334}*x18*x39);
auto x100 = V{0.25}*x90;
auto x101 = V{0.166666666666667}*x91;
auto x102 = x101*x82*x86;
auto x103 = V{0.0833333333333333}*x9;
auto x104 = V{0.0555555555555556}*x82;
auto x105 = -x103*x80 - x103*x84 + x104*x73*x86 + x104*x76*x86;
auto x106 = x100 - x102 + x105;
auto x107 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x108 = x17*x43 + x19*x39 - x41;
auto x109 = V{3}*x44;
auto x110 = V{3}*x45;
auto x111 = -x109 - x110;
auto x112 = V{3}*x37 + V{3}*x40 + x79;
auto x113 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x114 = -x81 + x87;
auto x115 = -x90 + x92;
auto x116 = V{0.333333333333333}*x9;
auto x117 = -x95 + x96;
auto x118 = V{0.166666666666667}*x9;
auto x119 = V{0.222222222222222}*x74;
auto x120 = V{0.111111111111111}*x74;
auto x121 = x41 + x46;
auto x122 = x24 + V{1};
auto x123 = -x57 - x65;
auto x124 = V{3}*x63;
auto x125 = V{3}*x64;
auto x126 = x124 - x125;
auto x127 = V{1.5}*x58;
auto x128 = V{3}*x55 + V{3}*x56 + V{1.5}*x66 + V{-1};
auto x129 = x127 + x128;
auto x130 = V{0.222222222222222}*x82;
auto x131 = V{0.111111111111111}*x82;
auto x132 = x25*x57*x65;
auto x133 = -x132 + x49 + x50 - x51 - x52;
auto x134 = x10 + V{3.00000046417339}*util::sqrt(x21 + x26*x27*util::sqrt(x133*(cell[1] - cell[7] - x132 - x33) + x69));
auto x135 = V{1} / (x134);
auto x136 = V{0.666666666666667}*x74;
auto x137 = -x136*x70*x91 + x53*x9;
auto x138 = V{0.0555555555555556}*x74;
auto x139 = -x57 + x63 - x64;
auto x0 = -V{1}*x71*(x61*x72 + x68*x72 + x70*x73*x75 + x70*x75*x76) - x79*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{-0.444444444444444};
auto x1 = -(x107*(x111 + x112 - x77*(V{6.000012}*x17*x43 + V{6.000012}*x19*x39) + x78*(V{6.000012}*x17*x36 + V{6.000012}*x18*x39) - V{4.5}*x108*x108) + V{1}*x71*(x106 + x89 + x94 - x98 - x99) + V{0.0277777777777778});
auto x2 = -x113*(V{3}*x17*x36 + V{3}*x18*x39 - x77*(V{2.999997}*x17*x43 + V{2.999997}*x19*x39) - V{3}*x77 - x78*(V{6.000003}*x17*x36 + V{6.000003}*x18*x39) + V{1.5}*x78 + V{-1}) + V{1}*x71*(x114*(V{0.666667}*x17*x54 + V{0.666667}*x18*x38) + x115*(V{0.666666}*x17*x54 + V{0.666666}*x18*x38) - x115*(V{1.333334}*x17*x62 - V{1.333334}*x19*x38) + x116*x61 - x117*(V{0.333333}*x17*x62 - V{0.333333}*x19*x38) - x118*x68 + x119*x70*x73 - x120*x70*x76) + V{-0.111111111111111};
auto x3 = -(x107*(x109 + x110 + x112 + x77*(V{18}*x17*x43 + V{18}*x19*x39) + x78*(V{18}*x17*x36 + V{18}*x18*x39) - V{4.5}*x121*x121) + V{1}*x71*(-x100 + x102 + x105 + x88*(V{0.5}*x17*x36 + V{0.5}*x18*x39) + x93*(V{1}*x17*x36 + V{1}*x18*x39) + x93*(V{1}*x17*x43 + V{1}*x19*x39) + x97*(V{0.5}*x17*x43 + V{0.5}*x19*x39)) + V{0.0277777777777778});
auto x4 = -x113*(-x111 - x77*(V{6.000003}*x17*x43 + V{6.000003}*x19*x39) + V{1.5}*x77 - x78*(V{2.999997}*x17*x36 + V{2.999997}*x18*x39) - V{3}*x78 + V{-1}) + V{1}*x71*(x114*(V{0.333333}*x17*x54 + V{0.333333}*x18*x38) + x115*(V{1.333334}*x17*x54 + V{1.333334}*x18*x38) - x115*(V{0.666666}*x17*x62 - V{0.666666}*x19*x38) + x116*x68 - x117*(V{0.666667}*x17*x62 - V{0.666667}*x19*x38) - x118*x61 + x119*x70*x76 - x120*x70*x73) + V{-0.111111111111111};
auto x5 = -(V{0.0277777777777778}*x122*(x126 + x129 + x58*(V{6.000012}*x17*x62 - V{6.000012}*x19*x38) + x66*(V{6.000012}*x17*x54 + V{6.000012}*x18*x38) - V{4.5}*x123*x123) + V{1}*x71*(x106 - x89 - x94 + x98 + x99) + V{0.0277777777777778});
auto x6 = V{0.111111111111111}*x122*(-x128 + V{6.000003}*x57*x66 - x58*(V{2.999997}*x17*x62 - V{2.999997}*x19*x38) + V{3}*x58) - V{1}*x71*(-x116*x84 + x118*x80 + x130*x73*x86 - x131*x76*x86 + x88*(V{0.666667}*x17*x36 + V{0.666667}*x18*x39) + x93*(V{0.666666}*x17*x36 + V{0.666666}*x18*x39) + x93*(V{1.333334}*x17*x43 + V{1.333334}*x19*x39) + x97*(V{0.333333}*x17*x43 + V{0.333333}*x19*x39)) + V{-0.111111111111111};
auto x7 = -x107*(-x124 + x125 + x129 - V{18}*x58*x65 + x66*(V{18}*x17*x54 + V{18}*x18*x38) - V{4.5}*x139*x139) + V{1}*(V{1} - x135)*(x101*x134*x74 + x103*x61 + x103*x68 - V{0.25}*x133*x9 + x134*x138*x73 + x134*x138*x76 + x137*(V{1}*x17*x54 + V{1}*x18*x38) - x137*(V{1}*x17*x62 - V{1}*x19*x38) - (V{0.5}*x17*x54 + V{0.5}*x18*x38)*(x136*x70*x76 + x68*x9) + (V{0.5}*x17*x62 - V{0.5}*x19*x38)*(x136*x70*x73 + x61*x9)) + V{-0.0277777777777778};
auto x8 = V{0.111111111111111}*x122*(x126 - x127 - x58*(V{6.000003}*x17*x62 - V{6.000003}*x19*x38) + x66*(V{2.999997}*x17*x54 + V{2.999997}*x18*x38) + V{3}*x66 + V{1}) - V{1}*x71*(-x116*x80 + x118*x84 + x130*x76*x86 - x131*x73*x86 + x88*(V{0.333333}*x17*x36 + V{0.333333}*x18*x39) + x93*(V{1.333334}*x17*x36 + V{1.333334}*x18*x39) + x93*(V{0.666666}*x17*x43 + V{0.666666}*x19*x39) + x97*(V{0.666667}*x17*x43 + V{0.666667}*x19*x39)) + V{-0.111111111111111};
cell.template getFieldPointer<descriptors::OMEGA>()[0] = V{1}*x135;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x25, x58 + x66 };
}
};

}

}
