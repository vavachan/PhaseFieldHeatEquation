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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x23 = x22 + V{-1};
auto x24 = x21 + V{-1};
auto x25 = -1/x24;
auto x26 = V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[6] + V{2}*cell[8];
auto x27 = cell[0] + cell[10] + cell[11] + cell[13] + cell[14] + cell[1] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + x26 + V{1};
auto x28 = x25*x27;
auto x29 = V{0.0277777777777778}*x28;
auto x30 = x19*x19;
auto x31 = V{3}*x19;
auto x32 = x20*x20;
auto x33 = V{1.5}*x32;
auto x34 = -x33;
auto x35 = x21*x21;
auto x36 = V{1.5}*x35;
auto x37 = V{1} - x36;
auto x38 = x34 + x37;
auto x39 = x31 + x38;
auto x40 = V{3}*x30 + x39;
auto x41 = V{3}*x20;
auto x42 = V{1.5}*x30;
auto x43 = -x42;
auto x44 = x41 + x43;
auto x45 = V{3}*x32 + x37 + x44;
auto x46 = x33 + x36 + V{-1};
auto x47 = -V{3}*x30 + x31 + x46;
auto x48 = -x47;
auto x49 = x42 + V{-1};
auto x50 = -V{3}*x32 + x36 + x41 + x49;
auto x51 = -x50;
auto x52 = x19 + x20;
auto x53 = x52*x52;
auto x54 = x39 + x44 + V{4.5}*x53;
auto x55 = -x31;
auto x56 = x19 - x20;
auto x57 = -V{4.5}*x56*x56;
auto x58 = x42 + x46;
auto x59 = x41 + x58;
auto x60 = x55 + x57 + x59;
auto x61 = -x60;
auto x62 = -x56;
auto x63 = -x41 + x58;
auto x64 = x31 + x63;
auto x65 = x64 - V{4.5}*x62*x62;
auto x66 = -x65;
auto x67 = x31 - V{4.5}*x53 + x59;
auto x68 = -x67;
auto x69 = V{0.0555555555555556}*x28;
auto x70 = V{3}*x21;
auto x71 = x33 - V{3}*x35 + x49 + x70;
auto x72 = -x71;
auto x73 = -x21;
auto x74 = x19 + x73;
auto x75 = x58 + x70;
auto x76 = x55 + x75 - V{4.5}*x74*x74;
auto x77 = -x76;
auto x78 = x20 + x73;
auto x79 = x63 + x70 - V{4.5}*x78*x78;
auto x80 = -x79;
auto x81 = x19 + x21;
auto x82 = x81*x81;
auto x83 = x31 + x75 - V{4.5}*x82;
auto x84 = -x83;
auto x85 = x20 + x21;
auto x86 = x85*x85;
auto x87 = x59 + x70 - V{4.5}*x86;
auto x88 = -x87;
auto x89 = -x58;
auto x90 = V{0.333333333333333}*x28*x89;
auto x91 = V{0.0555555555555555}*x28;
auto x92 = V{0.0555555555555555}*x28;
auto x93 = V{0.0277777777777778}*x28;
auto x94 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{6.66133814775094e-16}*cell[16] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{4.44089209850063e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{2.22044604925031e-16};
auto x95 = V{0.0555555555555555}*x28*x48 + V{0.0277777777777778}*x28*x66 + V{0.0277777777777778}*x28*x68 + V{0.111111111111111}*x28*x72 + x40*x69 + x45*x69 + x51*x91 + x54*x93 + x61*x93 + x77*x92 + x80*x91 + x84*x92 + x88*x91 + x90 + x94;
auto x96 = V{0.166666666666667}*cell[16];
auto x97 = V{0.166666666666667}*cell[6];
auto x98 = V{0.0833333333333333}*cell[11];
auto x99 = V{0.0833333333333333}*cell[2];
auto x100 = V{0.333333333333333}*cell[18];
auto x101 = V{0.333333333333333}*cell[8];
auto x102 = V{1} / (x24);
auto x103 = x102*x27;
auto x104 = V{0.00462962962962963}*x103;
auto x105 = V{0.00925925925925926}*x103;
auto x106 = V{0.0462962962962963}*x103;
auto x107 = V{0.00925925925925926}*x103;
auto x108 = x104*x45;
auto x109 = V{0.00462962962962963}*x103;
auto x110 = x109*x83;
auto x111 = x109*x76;
auto x112 = V{0.166666666666667}*cell[3];
auto x113 = V{0.00231481481481482}*x103;
auto x114 = x113*x67;
auto x115 = x113*x65;
auto x116 = x113*x60;
auto x117 = V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x105*x71 - x112 + x113*x54 - x114 - x115 - x116 + V{-0.0555555555555555};
auto x118 = V{0.0555555555555556}*x22;
auto x119 = V{0.166666666666667}*cell[18];
auto x120 = V{0.166666666666667}*cell[8];
auto x121 = V{0.0833333333333333}*cell[10];
auto x122 = V{0.0833333333333333}*cell[1];
auto x123 = V{0.333333333333333}*cell[16];
auto x124 = V{0.333333333333333}*cell[6];
auto x125 = x104*x40;
auto x126 = x109*x87;
auto x127 = x109*x79;
auto x128 = -x96;
auto x129 = -x97;
auto x130 = -V{0.00462962962962963}*x102*x27*x50 + x108 + x110 + x111 + x128 + x129 + x98 + x99 + V{0.0555555555555555};
auto x131 = -x119;
auto x132 = -x120;
auto x133 = -V{0.00462962962962963}*x102*x27*x47 + x121 + x122 + x125 + x126 + x127 + x131 + x132;
auto x134 = V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.333333333333333}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5];
auto x135 = V{0.0115740740740741}*x103;
auto x136 = V{0.0162037037037037}*x103;
auto x137 = V{0.00231481481481481}*x103;
auto x138 = V{0.00462962962962963}*x103;
auto x139 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] - x113*x79 - x113*x87 + x138*x40 - x138*x47 + V{0.0138888888888889};
auto x140 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] - x113*x76 - x113*x83 + x138*x45 - x138*x50;
auto x141 = -V{0.0833333333333333}*cell[3] + x104*x71 + x139 + x140;
auto x142 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x137*x60 + x137*x65 + x141;
auto x143 = V{0.0277777777777778}*x22;
auto x144 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x137*x54 + x137*x67 + x141;
auto x145 = V{0.00462962962962963}*x103;
auto x146 = V{0.00231481481481481}*x103;
auto x147 = V{0.00115740740740741}*x103;
auto x148 = V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.166666666666667}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x107*x71 + x147*x54 - x147*x60 - x147*x65 - x147*x67;
auto x149 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x139 - x146*x45 + x146*x50 + x148;
auto x150 = -V{0.166666666666667}*cell[16] + V{0.833333333333333}*cell[6] + x145*x76 + x149;
auto x151 = V{0.0277777777777778}*x103;
auto x152 = -x70;
auto x153 = -x74;
auto x154 = x152 + x31 + x58 - V{4.5}*x153*x153;
auto x155 = V{0.0231481481481481}*x103;
auto x156 = V{0.833333333333333}*cell[16] - V{0.166666666666667}*cell[6] + x145*x83 + x149;
auto x157 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x140 - x146*x40 + x146*x47 + x148 + V{0.0138888888888889};
auto x158 = -V{0.166666666666667}*cell[18] + V{0.833333333333333}*cell[8] + x145*x79 + x157;
auto x159 = -x78;
auto x160 = x152 + x59 - V{4.5}*x159*x159;
auto x161 = V{0.833333333333333}*cell[18] - V{0.166666666666667}*cell[8] + x145*x87 + x157;
auto x162 = -V{0.0833333333333334}*cell[13] - V{0.0833333333333334}*cell[14] - V{0.0833333333333334}*cell[4] - V{0.0833333333333334}*cell[5] - V{0.00231481481481482}*x102*x27*x54 - V{0.00925925925925926}*x102*x27*x71 + x112 + x114 + x115 + x116;
auto x163 = V{0.00462962962962963}*x28;
auto x164 = V{0.00462962962962963}*x28;
auto x165 = x43 + x70;
auto x166 = x165 + x34 + V{3}*x35 + V{1};
auto x167 = x165 + x39 + V{4.5}*x82;
auto x168 = x38 + x44 + x70 + V{4.5}*x86;
auto x169 = V{0.0555555555555555}*x103;
auto x170 = V{0.0555555555555555}*x103;
auto x171 = V{0.0277777777777778}*x103;
auto x172 = V{0.0555555555555556}*x103;
auto x0 = -V{0.333333333333333}*x22*(-x25*x89*x95 + V{1}) + x23*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[13] + V{1}*cell[14] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{1}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + x26 - x29*x40 - x29*x45 - x29*x48 - x29*x51 - x29*x54 - x29*x61 - x29*x66 - x29*x68 - x69*x72 - x69*x77 - x69*x80 - x69*x84 - x69*x88 - x90 + V{0.833333333333333});
auto x1 = -x118*(-x25*x48*x95 + V{1}) - x23*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x100 - x101 + x104*x50 + x105*x79 + x105*x87 + x106*x47 + x107*x40 - x108 - x110 - x111 + x117 + x96 + x97 - x98 - x99);
auto x2 = -x118*(-x25*x51*x95 + V{1}) - x23*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x104*x47 + x105*x76 + x105*x83 + x106*x50 + x107*x45 + x117 + x119 + x120 - x121 - x122 - x123 - x124 - x125 - x126 - x127);
auto x3 = -x118*(-x25*x72*x95 + V{1}) - x23*(V{0.00462962962962963}*x102*x27*x60 + V{0.00462962962962963}*x102*x27*x65 + V{0.00462962962962963}*x102*x27*x67 + V{0.037037037037037}*x102*x27*x71 - x104*x54 - x130 - x133 - x134);
auto x4 = -x143*(-x25*x68*x95 + V{1}) - x23*(x135*x54 + x136*x67 + x142);
auto x5 = -x143*(-x25*x66*x95 + V{1}) - x23*(-x135*x60 + x136*x65 + x144);
auto x6 = -x143*(-x25*x84*x95 + V{1}) - x23*(x104*x83 + x150);
auto x7 = -x143*(x154*x25*x95 + V{1}) - x23*(x151*x154 - x155*x76 + x156);
auto x8 = -x143*(-x25*x88*x95 + V{1}) - x23*(x104*x87 + x158);
auto x9 = -x143*(x160*x25*x95 + V{1}) - x23*(x151*x160 - x155*x79 + x161);
auto x10 = -x118*(-x25*x40*x95 + V{1}) - x23*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x100 - x101 + V{0.00925925925925926}*x102*x27*x79 + V{0.00925925925925926}*x102*x27*x87 - x106*x40 - x107*x47 - x130 - x162);
auto x11 = -x118*(-x25*x45*x95 + V{1}) - x23*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.00925925925925926}*x102*x27*x76 + V{0.00925925925925926}*x102*x27*x83 - x106*x45 - x107*x50 - x123 - x124 - x133 - x162 + V{-0.0555555555555555});
auto x12 = -x118*(-x166*x25*x95 + V{1}) + x23*(x121 + x122 + x128 + x129 + x131 + x132 + x134 + x163*x77 + x163*x80 + x163*x84 + x163*x88 - x164*x40 - x164*x45 - x164*x48 - x164*x51 - x164*x54 - x164*x61 - x164*x66 - x164*x68 - x166*x69 + V{0.0185185185185185}*x28*x72 + x98 + x99 + V{0.0555555555555555});
auto x13 = -x143*(-x25*x54*x95 + V{1}) - x23*(-x135*x67 - x136*x54 + x142);
auto x14 = -x143*(-x25*x61*x95 + V{1}) - x23*(-x135*x65 + x136*x60 + x144);
auto x15 = -x143*(-x167*x25*x95 + V{1}) - x23*(x150 - x151*x167 - x155*x83);
auto x16 = -x143*(-x25*x77*x95 + V{1}) - x23*(x104*x76 + x156);
auto x17 = -x143*(-x168*x25*x95 + V{1}) - x23*(-x151*x168 - x155*x87 + x158);
auto x18 = -x143*(-x25*x80*x95 + V{1}) - x23*(x104*x79 + x161);
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
return { -V{1}*x102*(V{0.0555555555555555}*x103*x47 + V{0.333333333333333}*x103*x58 + V{0.0277777777777778}*x103*x67 + V{0.111111111111111}*x103*x71 + V{0.0277777777777778}*x103*(x57 + x64) + x169*x50 + x169*x79 + x169*x87 + x170*x76 + x170*x83 - x171*x54 + x171*x60 - x172*x40 - x172*x45 + x94), x30 + x32 + x35 };
}
};

}

}
