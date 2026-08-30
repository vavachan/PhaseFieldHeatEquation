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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x23 = x22 + V{-1};
auto x24 = V{1} / (x21 + V{1});
auto x25 = V{2}*cell[15] + V{2}*cell[17] + V{2}*cell[7] + V{2}*cell[9];
auto x26 = cell[0] + cell[10] + cell[11] + V{2}*cell[12] + cell[13] + cell[14] + cell[1] + cell[2] + cell[4] + cell[5] + x25 + V{1};
auto x27 = x24*x26;
auto x28 = V{0.0277777777777778}*x27;
auto x29 = x19*x19;
auto x30 = V{3}*x19;
auto x31 = x20*x20;
auto x32 = V{1.5}*x31;
auto x33 = -x32;
auto x34 = x21*x21;
auto x35 = V{1.5}*x34;
auto x36 = V{1} - x35;
auto x37 = x33 + x36;
auto x38 = x30 + x37;
auto x39 = V{3}*x29 + x38;
auto x40 = V{3}*x20;
auto x41 = V{1.5}*x29;
auto x42 = -x41;
auto x43 = x40 + x42;
auto x44 = V{3}*x31 + x36 + x43;
auto x45 = x32 + x35 + V{-1};
auto x46 = -V{3}*x29 + x30 + x45;
auto x47 = -x46;
auto x48 = x41 + V{-1};
auto x49 = -V{3}*x31 + x35 + x40 + x48;
auto x50 = -x49;
auto x51 = x19 + x20;
auto x52 = x51*x51;
auto x53 = x38 + x43 + V{4.5}*x52;
auto x54 = -x30;
auto x55 = x19 - x20;
auto x56 = -V{4.5}*x55*x55;
auto x57 = x41 + x45;
auto x58 = x40 + x57;
auto x59 = x54 + x56 + x58;
auto x60 = -x59;
auto x61 = -x55;
auto x62 = -x40;
auto x63 = x30 + x57;
auto x64 = x62 + x63;
auto x65 = x64 - V{4.5}*x61*x61;
auto x66 = -x65;
auto x67 = x30 - V{4.5}*x52 + x58;
auto x68 = -x67;
auto x69 = V{0.0555555555555556}*x27;
auto x70 = V{3}*x21;
auto x71 = x42 + x70;
auto x72 = x33 + V{3}*x34 + x71 + V{1};
auto x73 = x19 + x21;
auto x74 = V{4.5}*(x73*x73);
auto x75 = x38 + x71 + x74;
auto x76 = x20 + x21;
auto x77 = V{4.5}*(x76*x76);
auto x78 = x37 + x43 + x70 + x77;
auto x79 = -x21;
auto x80 = x19 + x79;
auto x81 = -x80;
auto x82 = -x70;
auto x83 = x63 + x82;
auto x84 = x83 - V{4.5}*x81*x81;
auto x85 = -x84;
auto x86 = x20 + x79;
auto x87 = -x86;
auto x88 = x58 + x82;
auto x89 = x88 - V{4.5}*x87*x87;
auto x90 = -x89;
auto x91 = V{0.333333333333333}*x27;
auto x92 = -x57;
auto x93 = x91*x92;
auto x94 = V{0.0277777777777778}*x27;
auto x95 = V{0.0555555555555555}*x27;
auto x96 = V{0.0555555555555555}*x27;
auto x97 = V{0.0277777777777778}*x27;
auto x98 = V{0.0555555555555555}*x27;
auto x99 = V{0.0277777777777778}*x27;
auto x100 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{4.44089209850063e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{6.66133814775094e-16}*cell[15] + V{4.44089209850063e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{6.66133814775094e-16}*cell[7] + V{0.111111111111111}*x27*x72 + x39*x69 + x44*x69 + x53*x99 + x75*x98 + x78*x95 + V{2.22044604925031e-16};
auto x101 = x100 + x47*x96 + x50*x95 + x60*x99 + x66*x94 + x68*x97 + x85*x98 + x90*x95 + x93;
auto x102 = V{0.0462962962962963}*x27;
auto x103 = V{0.00925925925925926}*x27;
auto x104 = V{0.166666666666667}*cell[15];
auto x105 = V{0.166666666666667}*cell[7];
auto x106 = V{0.0833333333333333}*cell[11];
auto x107 = V{0.0833333333333333}*cell[2];
auto x108 = V{0.00462962962962963}*x27;
auto x109 = x108*x44;
auto x110 = V{0.00925925925925926}*x27;
auto x111 = V{0.00462962962962963}*x27;
auto x112 = x108*x49;
auto x113 = x111*x75;
auto x114 = V{0.00231481481481482}*x27;
auto x115 = -V{0.166666666666667}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x110*x72 - x114*x53 + x114*x59 + x114*x65 + x114*x67 + V{-0.0555555555555555};
auto x116 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] + V{0.166666666666667}*cell[1] - V{0.333333333333333}*cell[9] + x104 + x105 - x106 - x107 + x109 + x110*x78 - x110*x89 + x111*x84 - x112 - x113 + x115;
auto x117 = V{0.0555555555555556}*x22;
auto x118 = V{0.166666666666667}*cell[17];
auto x119 = V{0.166666666666667}*cell[9];
auto x120 = V{0.0833333333333333}*cell[10];
auto x121 = V{0.0833333333333333}*cell[1];
auto x122 = x108*x39;
auto x123 = x108*x46;
auto x124 = x111*x78;
auto x125 = V{0.166666666666667}*cell[11] - V{0.333333333333333}*cell[15] + V{0.166666666666667}*cell[2] - V{0.333333333333333}*cell[7] + x110*x75 - x110*x84 + x111*x89 + x115 + x118 + x119 - x120 - x121 + x122 - x123 - x124;
auto x126 = -x32 + V{3}*x34 - x48 - x70;
auto x127 = -V{0.333333333333333}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x104 - x105 + x106 + x107 - x108*x53 - x109 + x113 - x118 - x119 + x120 + x121 - x122 + x124 + V{0.0555555555555555};
auto x128 = V{0.0115740740740741}*x27;
auto x129 = V{0.0162037037037037}*x27;
auto x130 = V{0.00231481481481481}*x27;
auto x131 = V{0.00462962962962963}*x27;
auto x132 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[9] - x114*x78 + x114*x89 - x131*x39 + x131*x46 + V{0.0138888888888889};
auto x133 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[7] - x114*x75 + x114*x84 - x131*x44 + x131*x49;
auto x134 = -V{0.0833333333333333}*cell[12] + x108*x72 + x132 + x133;
auto x135 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x130*x59 - x130*x65 + x134;
auto x136 = V{0.0277777777777778}*x22;
auto x137 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x130*x53 - x130*x67 + x134;
auto x138 = x63 + x70 - x74;
auto x139 = V{0.0231481481481481}*x27;
auto x140 = V{0.00462962962962963}*x27;
auto x141 = V{0.00231481481481481}*x27;
auto x142 = V{0.00115740740740741}*x27;
auto x143 = V{0.166666666666667}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x103*x72 - x142*x53 + x142*x59 + x142*x65 + x142*x67;
auto x144 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x132 + x141*x44 - x141*x49 + x143;
auto x145 = V{0.833333333333333}*cell[15] - V{0.166666666666667}*cell[7] - x140*x84 + x144;
auto x146 = -V{0.166666666666667}*cell[15] + V{0.833333333333333}*cell[7] + x140*x75 + x144;
auto x147 = x58 + x70 - x77;
auto x148 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x133 + x141*x39 - x141*x46 + x143 + V{0.0138888888888889};
auto x149 = V{0.833333333333333}*cell[17] - V{0.166666666666667}*cell[9] - x140*x89 + x148;
auto x150 = -V{0.166666666666667}*cell[17] + V{0.833333333333333}*cell[9] + x140*x78 + x148;
auto x151 = -V{4.5}*x80*x80;
auto x152 = x57 + x70;
auto x153 = x151 + x152 + x54;
auto x154 = -V{4.5}*x86*x86;
auto x155 = x152 + x154 + x62;
auto x0 = -V{0.333333333333333}*x22*(-x101*x24*x92 + V{1}) + x23*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{1}*cell[4] + V{1}*cell[5] + x25 - x28*x39 - x28*x44 - x28*x47 - x28*x50 - x28*x53 - x28*x60 - x28*x66 - x28*x68 - x69*x72 - x69*x75 - x69*x78 - x69*x85 - x69*x90 - x93 + V{0.833333333333333});
auto x1 = -x117*(-x101*x24*x47 + V{1}) - x23*(-x102*x46 - x103*x39 + x116);
auto x2 = -x117*(-x101*x24*x50 + V{1}) - x23*(-x102*x49 - x103*x44 + x125);
auto x3 = -x117*(-x101*x126*x24 + V{1}) + x23*(-x108*x47 - x108*x50 - x108*x60 - x108*x66 - x108*x68 + x111*x85 + x111*x90 - x126*x69 + x127 + V{0.0185185185185185}*x27*x72);
auto x4 = -x136*(-x101*x24*x68 + V{1}) - x23*(-x128*x53 - x129*x67 + x135);
auto x5 = -x136*(-x101*x24*x66 + V{1}) - x23*(x128*x59 - x129*x65 + x137);
auto x6 = -x136*(x101*x138*x24 + V{1}) - x23*(-x138*x28 - x139*x75 + x145);
auto x7 = -x136*(-x101*x24*x85 + V{1}) - x23*(-x108*x84 + x146);
auto x8 = -x136*(x101*x147*x24 + V{1}) - x23*(-x139*x78 - x147*x28 + x149);
auto x9 = -x136*(-x101*x24*x90 + V{1}) - x23*(-x108*x89 + x150);
auto x10 = -x117*(-x101*x24*x39 + V{1}) - x23*(x102*x39 + x103*x46 + x116);
auto x11 = -x117*(-x101*x24*x44 + V{1}) - x23*(x102*x44 + x103*x49 + x125);
auto x12 = -x117*(-x101*x24*x72 + V{1}) - x23*(-x108*x59 - x108*x65 - x108*x67 - x112 - x123 - x127 + V{0.037037037037037}*x24*x26*x72 + V{0.00462962962962963}*x24*x26*x84 + V{0.00462962962962963}*x24*x26*x89);
auto x13 = -x136*(-x101*x24*x53 + V{1}) - x23*(x128*x67 + x129*x53 + x135);
auto x14 = -x136*(-x101*x24*x60 + V{1}) - x23*(x128*x65 - x129*x59 + x137);
auto x15 = -x136*(-x101*x24*x75 + V{1}) - x23*(x108*x75 + x145);
auto x16 = -x136*(x101*x153*x24 + V{1}) - x23*(x139*x84 + x146 - x153*x28);
auto x17 = -x136*(-x101*x24*x78 + V{1}) - x23*(x108*x78 + x149);
auto x18 = -x136*(x101*x155*x24 + V{1}) - x23*(x139*x89 + x150 - x155*x28);
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
return { V{1}*x24*(x100 - x46*x96 - x49*x95 - x57*x91 - x59*x99 - x67*x97 - x94*(x56 + x64) - x95*(x154 + x88) - x98*(x151 + x83)), x29 + x31 + x34 };
}
};

}

}
