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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x23 = x22 + V{-1};
auto x24 = V{1} / (x19 + V{1});
auto x25 = V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16];
auto x26 = cell[0] + V{2}*cell[10] + cell[11] + cell[12] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + x25 + V{1};
auto x27 = x24*x26;
auto x28 = V{0.0277777777777778}*x27;
auto x29 = x20*x20;
auto x30 = V{3}*x20;
auto x31 = x21*x21;
auto x32 = V{1.5}*x31;
auto x33 = -x32;
auto x34 = x19*x19;
auto x35 = V{1.5}*x34;
auto x36 = V{1} - x35;
auto x37 = x33 + x36;
auto x38 = x30 + x37;
auto x39 = V{3}*x29 + x38;
auto x40 = V{3}*x21;
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
auto x51 = x20 + x21;
auto x52 = x51*x51;
auto x53 = x38 + x43 + V{4.5}*x52;
auto x54 = -x30;
auto x55 = -x21;
auto x56 = x20 + x55;
auto x57 = -V{4.5}*x56*x56;
auto x58 = x41 + x45;
auto x59 = x40 + x58;
auto x60 = x54 + x57 + x59;
auto x61 = -x60;
auto x62 = -x56;
auto x63 = -x40;
auto x64 = x30 + x58;
auto x65 = x63 + x64;
auto x66 = x65 - V{4.5}*x62*x62;
auto x67 = -x66;
auto x68 = x30 - V{4.5}*x52 + x59;
auto x69 = -x68;
auto x70 = V{0.0555555555555556}*x27;
auto x71 = V{3}*x34;
auto x72 = V{3}*x19;
auto x73 = x42 + x72;
auto x74 = x33 + x71 + x73 + V{1};
auto x75 = x19 + x20;
auto x76 = V{4.5}*(x75*x75);
auto x77 = x38 + x73 + x76;
auto x78 = x19 + x21;
auto x79 = V{4.5}*(x78*x78);
auto x80 = x37 + x43 + x72 + x79;
auto x81 = -x72;
auto x82 = x19 - x20;
auto x83 = x64 + x81 - V{4.5}*x82*x82;
auto x84 = -x83;
auto x85 = x19 + x55;
auto x86 = x59 + x81 - V{4.5}*x85*x85;
auto x87 = -x86;
auto x88 = V{0.333333333333333}*x27;
auto x89 = -x58;
auto x90 = x88*x89;
auto x91 = V{0.0277777777777778}*x27;
auto x92 = V{0.0555555555555555}*x27;
auto x93 = V{0.0555555555555555}*x27;
auto x94 = V{0.0277777777777778}*x27;
auto x95 = V{0.0555555555555555}*x27;
auto x96 = V{0.0277777777777778}*x27;
auto x97 = V{4.44089209850063e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{6.66133814775094e-16}*cell[14] + V{8.88178419700125e-16}*cell[15] + V{4.44089209850063e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{0.111111111111111}*x27*x74 + x39*x70 + x44*x70 + x53*x96 + x77*x95 + x80*x92 + V{2.22044604925031e-16};
auto x98 = x47*x93 + x50*x92 + x61*x96 + x67*x91 + x69*x94 + x84*x95 + x87*x92 + x90 + x97;
auto x99 = x32 + x48 - x71 + x72;
auto x100 = V{0.0833333333333333}*cell[11];
auto x101 = V{0.0833333333333333}*cell[12];
auto x102 = V{0.0833333333333333}*cell[2];
auto x103 = V{0.0833333333333333}*cell[3];
auto x104 = V{0.00462962962962963}*x27;
auto x105 = x104*x49;
auto x106 = x104*x46;
auto x107 = V{0.00462962962962963}*x27;
auto x108 = x107*x77;
auto x109 = x107*x80;
auto x110 = -V{0.333333333333333}*cell[10] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x100 + x101 + x102 + x103 + x104*x60 + x104*x66 + x104*x68 + x105 + x106 + x108 + x109 - V{0.00462962962962963}*x24*x26*x39 - V{0.00462962962962963}*x24*x26*x44 - V{0.00462962962962963}*x24*x26*x53 - V{0.00462962962962963}*x24*x26*x83 - V{0.00462962962962963}*x24*x26*x86 + V{0.0555555555555555};
auto x111 = V{0.0555555555555556}*x22;
auto x112 = V{0.0462962962962963}*x27;
auto x113 = V{0.00925925925925926}*x27;
auto x114 = V{0.00925925925925926}*x27;
auto x115 = V{0.00231481481481482}*x27;
auto x116 = -V{0.166666666666667}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x114*x74 - x115*x53 + x115*x60 + x115*x66 + x115*x68 + V{-0.0555555555555555};
auto x117 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2] - x101 - x103 + x104*x44 - x105 + x107*x83 - x108 + x114*x80 - x114*x86 + x116;
auto x118 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[3] - x100 - x102 + x104*x39 - x106 + x107*x86 - x109 + x114*x77 - x114*x83 + x116;
auto x119 = x64 + x72 - x76;
auto x120 = V{0.0231481481481481}*x27;
auto x121 = V{0.00462962962962963}*x27;
auto x122 = V{0.00231481481481481}*x27;
auto x123 = V{0.00462962962962963}*x27;
auto x124 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] - x115*x80 + x115*x86 - x123*x39 + x123*x46 + V{0.0138888888888889};
auto x125 = V{0.00115740740740741}*x27;
auto x126 = V{0.166666666666667}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] - x113*x74 - x125*x53 + x125*x60 + x125*x66 + x125*x68;
auto x127 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x122*x44 - x122*x49 + x124 + x126;
auto x128 = V{0.833333333333333}*cell[13] - V{0.166666666666667}*cell[14] - x121*x83 + x127;
auto x129 = V{0.0277777777777778}*x22;
auto x130 = -x82;
auto x131 = x58 + x72;
auto x132 = x131 + x54 - V{4.5}*x130*x130;
auto x133 = -V{0.166666666666667}*cell[13] + V{0.833333333333333}*cell[14] + x121*x77 + x127;
auto x134 = x59 + x72 - x79;
auto x135 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - x115*x77 + x115*x83 - x123*x44 + x123*x49;
auto x136 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x122*x39 - x122*x46 + x126 + x135 + V{0.0138888888888889};
auto x137 = V{0.833333333333333}*cell[15] - V{0.166666666666667}*cell[16] - x121*x86 + x136;
auto x138 = -x85;
auto x139 = x131 + x63 - V{4.5}*x138*x138;
auto x140 = -V{0.166666666666667}*cell[15] + V{0.833333333333333}*cell[16] + x121*x80 + x136;
auto x141 = V{0.0115740740740741}*x27;
auto x142 = V{0.0162037037037037}*x27;
auto x143 = V{0.00231481481481481}*x27;
auto x144 = -V{0.0833333333333333}*cell[10] + x104*x74 + x124 + x135;
auto x145 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x143*x60 - x143*x66 + x144;
auto x146 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x143*x53 - x143*x68 + x144;
auto x0 = -V{0.333333333333333}*x22*(-x24*x89*x98 + V{1}) + x23*(V{1}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[8] + V{1}*cell[9] + x25 - x28*x39 - x28*x44 - x28*x47 - x28*x50 - x28*x53 - x28*x61 - x28*x67 - x28*x69 - x70*x74 - x70*x77 - x70*x80 - x70*x84 - x70*x87 - x90 + V{0.833333333333333});
auto x1 = -x111*(x24*x98*x99 + V{1}) - x23*(-x110 - V{0.0185185185185185}*x27*x74 - x70*x99);
auto x2 = -x111*(-x24*x47*x98 + V{1}) - x23*(-x112*x46 - x113*x39 + x117);
auto x3 = -x111*(-x24*x50*x98 + V{1}) - x23*(-x112*x49 - x113*x44 + x118);
auto x4 = -x129*(x119*x24*x98 + V{1}) - x23*(-x119*x28 - x120*x77 + x128);
auto x5 = -x129*(x132*x24*x98 + V{1}) - x23*(x120*x83 - x132*x28 + x133);
auto x6 = -x129*(x134*x24*x98 + V{1}) - x23*(-x120*x80 - x134*x28 + x137);
auto x7 = -x129*(x139*x24*x98 + V{1}) - x23*(x120*x86 - x139*x28 + x140);
auto x8 = -x129*(-x24*x69*x98 + V{1}) - x23*(-x141*x53 - x142*x68 + x145);
auto x9 = -x129*(-x24*x67*x98 + V{1}) - x23*(x141*x60 - x142*x66 + x146);
auto x10 = -x111*(-x24*x74*x98 + V{1}) - x23*(-x110 + V{0.037037037037037}*x24*x26*x74);
auto x11 = -x111*(-x24*x39*x98 + V{1}) - x23*(x112*x39 + x113*x46 + x117);
auto x12 = -x111*(-x24*x44*x98 + V{1}) - x23*(x112*x44 + x113*x49 + x118);
auto x13 = -x129*(-x24*x77*x98 + V{1}) - x23*(x104*x77 + x128);
auto x14 = -x129*(-x24*x84*x98 + V{1}) - x23*(-x104*x83 + x133);
auto x15 = -x129*(-x24*x80*x98 + V{1}) - x23*(x104*x80 + x137);
auto x16 = -x129*(-x24*x87*x98 + V{1}) - x23*(-x104*x86 + x140);
auto x17 = -x129*(-x24*x53*x98 + V{1}) - x23*(x141*x68 + x142*x53 + x145);
auto x18 = -x129*(-x24*x61*x98 + V{1}) - x23*(x141*x66 - x142*x60 + x146);
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
return { V{1}*x24*(-x46*x93 - x49*x92 - x58*x88 - x60*x96 - x68*x94 - x83*x95 - x86*x92 - x91*(x57 + x65) + x97), x29 + x31 + x34 };
}
};

}

}
