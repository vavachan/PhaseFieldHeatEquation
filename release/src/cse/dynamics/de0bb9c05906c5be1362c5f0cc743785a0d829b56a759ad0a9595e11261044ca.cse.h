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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x23 = x22 + V{-1};
auto x24 = x20 + V{-1};
auto x25 = -1/x24;
auto x26 = V{2}*cell[14] + V{2}*cell[4] + V{2}*cell[8] + V{2}*cell[9];
auto x27 = cell[0] + cell[10] + cell[12] + cell[15] + cell[16] + cell[1] + V{2}*cell[2] + cell[3] + cell[6] + cell[7] + x26 + V{1};
auto x28 = x25*x27;
auto x29 = V{0.0277777777777778}*x28;
auto x30 = x19*x19;
auto x31 = V{3}*x19;
auto x32 = x21*x21;
auto x33 = V{1.5}*x32;
auto x34 = -x33;
auto x35 = x20*x20;
auto x36 = V{1.5}*x35;
auto x37 = V{1} - x36;
auto x38 = x34 + x37;
auto x39 = x31 + x38;
auto x40 = V{3}*x30 + x39;
auto x41 = V{3}*x21;
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
auto x52 = x19 + x21;
auto x53 = x52*x52;
auto x54 = x39 + x44 + V{4.5}*x53;
auto x55 = -x31;
auto x56 = -x21;
auto x57 = x19 + x56;
auto x58 = -V{4.5}*x57*x57;
auto x59 = x42 + x46;
auto x60 = x41 + x59;
auto x61 = x55 + x58 + x60;
auto x62 = -x61;
auto x63 = -x57;
auto x64 = -x41 + x59;
auto x65 = x31 + x64;
auto x66 = x65 - V{4.5}*x63*x63;
auto x67 = -x66;
auto x68 = x31 - V{4.5}*x53 + x60;
auto x69 = -x68;
auto x70 = V{0.0555555555555556}*x28;
auto x71 = V{3}*x20;
auto x72 = x33 - V{3}*x35 + x49 + x71;
auto x73 = -x72;
auto x74 = x19 - x20;
auto x75 = x59 + x71;
auto x76 = x55 + x75 - V{4.5}*x74*x74;
auto x77 = -x76;
auto x78 = x20 + x56;
auto x79 = -x78;
auto x80 = x64 + x71;
auto x81 = x80 - V{4.5}*x79*x79;
auto x82 = -x81;
auto x83 = x19 + x20;
auto x84 = x83*x83;
auto x85 = x31 + x75 - V{4.5}*x84;
auto x86 = -x85;
auto x87 = x20 + x21;
auto x88 = x87*x87;
auto x89 = x60 + x71 - V{4.5}*x88;
auto x90 = -x89;
auto x91 = -x59;
auto x92 = V{0.333333333333333}*x28*x91;
auto x93 = V{0.0555555555555555}*x28;
auto x94 = V{0.0555555555555555}*x28;
auto x95 = V{0.0277777777777778}*x28;
auto x96 = V{1.66533453693773e-16}*cell[10] + V{9.71445146547012e-17}*cell[12] + V{6.66133814775094e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{1.66533453693773e-16}*cell[1] + V{4.44089209850063e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{2.22044604925031e-16}*cell[4] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{4.44089209850063e-16}*cell[9] + V{2.22044604925031e-16};
auto x97 = V{0.0555555555555555}*x28*x48 + V{0.0277777777777778}*x28*x67 + V{0.0277777777777778}*x28*x69 + V{0.111111111111111}*x28*x73 + x40*x70 + x45*x70 + x51*x93 + x54*x95 + x62*x95 + x77*x94 + x82*x93 + x86*x94 + x90*x93 + x92 + x96;
auto x98 = V{0.166666666666667}*cell[14];
auto x99 = V{0.166666666666667}*cell[4];
auto x100 = V{0.0833333333333333}*cell[12];
auto x101 = V{0.0833333333333333}*cell[3];
auto x102 = V{0.333333333333333}*cell[8];
auto x103 = V{0.333333333333333}*cell[9];
auto x104 = V{1} / (x24);
auto x105 = x104*x27;
auto x106 = V{0.00462962962962963}*x105;
auto x107 = V{0.00925925925925926}*x105;
auto x108 = V{0.0462962962962963}*x105;
auto x109 = V{0.00925925925925926}*x105;
auto x110 = x106*x45;
auto x111 = V{0.00462962962962963}*x105;
auto x112 = x111*x85;
auto x113 = x111*x76;
auto x114 = V{0.166666666666667}*cell[2];
auto x115 = V{0.00231481481481482}*x105;
auto x116 = x115*x68;
auto x117 = x115*x66;
auto x118 = x115*x61;
auto x119 = V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x107*x72 - x114 + x115*x54 - x116 - x117 - x118 + V{-0.0555555555555555};
auto x120 = V{0.0555555555555556}*x22;
auto x121 = -x98;
auto x122 = -x99;
auto x123 = x100 + x101 - V{0.00462962962962963}*x104*x27*x50 + x110 + x112 + x113 + x121 + x122 + V{0.0555555555555555};
auto x124 = V{0.0833333333333333}*cell[10];
auto x125 = V{0.0833333333333333}*cell[1];
auto x126 = V{0.166666666666667}*cell[8];
auto x127 = -x126;
auto x128 = V{0.166666666666667}*cell[9];
auto x129 = -x128;
auto x130 = x106*x40;
auto x131 = x111*x89;
auto x132 = x111*x81;
auto x133 = -V{0.00462962962962963}*x104*x27*x47 + x124 + x125 + x127 + x129 + x130 + x131 + x132;
auto x134 = V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.333333333333333}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7];
auto x135 = V{0.333333333333333}*cell[14];
auto x136 = V{0.333333333333333}*cell[4];
auto x137 = V{0.00462962962962963}*x105;
auto x138 = V{0.00231481481481481}*x105;
auto x139 = V{0.00462962962962963}*x105;
auto x140 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x115*x81 - x115*x89 + x139*x40 - x139*x47 + V{0.0138888888888889};
auto x141 = V{0.00115740740740741}*x105;
auto x142 = V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.166666666666667}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x109*x72 + x141*x54 - x141*x61 - x141*x66 - x141*x68;
auto x143 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x138*x45 + x138*x50 + x140 + x142;
auto x144 = -V{0.166666666666667}*cell[14] + V{0.833333333333333}*cell[4] + x137*x76 + x143;
auto x145 = V{0.0277777777777778}*x22;
auto x146 = V{0.0277777777777778}*x105;
auto x147 = -x71;
auto x148 = -x74;
auto x149 = x147 + x31 + x59 - V{4.5}*x148*x148;
auto x150 = V{0.0231481481481481}*x105;
auto x151 = V{0.833333333333333}*cell[14] - V{0.166666666666667}*cell[4] + x137*x85 + x143;
auto x152 = V{0.0115740740740741}*x105;
auto x153 = V{0.0162037037037037}*x105;
auto x154 = V{0.00231481481481481}*x105;
auto x155 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] - x115*x76 - x115*x85 + x139*x45 - x139*x50;
auto x156 = -V{0.0833333333333333}*cell[2] + x106*x72 + x140 + x155;
auto x157 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x154*x61 + x154*x66 + x156;
auto x158 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x154*x54 + x154*x68 + x156;
auto x159 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x138*x40 + x138*x47 + x142 + x155 + V{0.0138888888888889};
auto x160 = V{0.833333333333333}*cell[8] - V{0.166666666666667}*cell[9] + x137*x81 + x159;
auto x161 = -V{0.166666666666667}*cell[8] + V{0.833333333333333}*cell[9] + x137*x89 + x159;
auto x162 = -V{0.0833333333333334}*cell[15] - V{0.0833333333333334}*cell[16] - V{0.0833333333333334}*cell[6] - V{0.0833333333333334}*cell[7] - V{0.00231481481481482}*x104*x27*x54 - V{0.00925925925925926}*x104*x27*x72 + x114 + x116 + x117 + x118;
auto x163 = V{0.00462962962962963}*x28;
auto x164 = V{0.00462962962962963}*x28;
auto x165 = x43 + x71;
auto x166 = x165 + x34 + V{3}*x35 + V{1};
auto x167 = x165 + x39 + V{4.5}*x84;
auto x168 = x38 + x44 + x71 + V{4.5}*x88;
auto x169 = -V{4.5}*x78*x78;
auto x170 = x147 + x169 + x60;
auto x171 = V{0.0555555555555555}*x105;
auto x172 = V{0.0555555555555555}*x105;
auto x173 = V{0.0277777777777778}*x105;
auto x174 = V{0.0555555555555556}*x105;
auto x0 = -V{0.333333333333333}*x22*(-x25*x91*x97 + V{1}) + x23*(V{0.5}*cell[10] + V{0.5}*cell[12] + V{1}*cell[15] + V{1}*cell[16] + V{0.5}*cell[1] + V{1}*cell[2] + V{0.5}*cell[3] + V{1}*cell[6] + V{1}*cell[7] + x26 - x29*x40 - x29*x45 - x29*x48 - x29*x51 - x29*x54 - x29*x62 - x29*x67 - x29*x69 - x70*x73 - x70*x77 - x70*x82 - x70*x86 - x70*x90 - x92 + V{0.833333333333333});
auto x1 = -x120*(-x25*x48*x97 + V{1}) - x23*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x100 - x101 - x102 - x103 + x106*x50 + x107*x81 + x107*x89 + x108*x47 + x109*x40 - x110 - x112 - x113 + x119 + x98 + x99);
auto x2 = -x120*(-x25*x73*x97 + V{1}) - x23*(V{0.00462962962962963}*x104*x27*x61 + V{0.00462962962962963}*x104*x27*x66 + V{0.00462962962962963}*x104*x27*x68 + V{0.037037037037037}*x104*x27*x72 - x106*x54 - x123 - x133 - x134);
auto x3 = -x120*(-x25*x51*x97 + V{1}) - x23*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x106*x47 + x107*x76 + x107*x85 + x108*x50 + x109*x45 + x119 - x124 - x125 + x126 + x128 - x130 - x131 - x132 - x135 - x136);
auto x4 = -x145*(-x25*x86*x97 + V{1}) - x23*(x106*x85 + x144);
auto x5 = -x145*(x149*x25*x97 + V{1}) - x23*(x146*x149 - x150*x76 + x151);
auto x6 = -x145*(-x25*x69*x97 + V{1}) - x23*(x152*x54 + x153*x68 + x157);
auto x7 = -x145*(-x25*x67*x97 + V{1}) - x23*(-x152*x61 + x153*x66 + x158);
auto x8 = -x145*(-x25*x90*x97 + V{1}) - x23*(x106*x89 + x160);
auto x9 = -x145*(-x25*x82*x97 + V{1}) - x23*(x106*x81 + x161);
auto x10 = -x120*(-x25*x40*x97 + V{1}) - x23*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x102 - x103 + V{0.00925925925925926}*x104*x27*x81 + V{0.00925925925925926}*x104*x27*x89 - x108*x40 - x109*x47 - x123 - x162);
auto x11 = -x120*(-x166*x25*x97 + V{1}) + x23*(x100 + x101 + x121 + x122 + x124 + x125 + x127 + x129 + x134 + x163*x77 + x163*x82 + x163*x86 + x163*x90 - x164*x40 - x164*x45 - x164*x48 - x164*x51 - x164*x54 - x164*x62 - x164*x67 - x164*x69 - x166*x70 + V{0.0185185185185185}*x28*x73 + V{0.0555555555555555});
auto x12 = -x120*(-x25*x45*x97 + V{1}) - x23*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.00925925925925926}*x104*x27*x76 + V{0.00925925925925926}*x104*x27*x85 - x108*x45 - x109*x50 - x133 - x135 - x136 - x162 + V{-0.0555555555555555});
auto x13 = -x145*(-x167*x25*x97 + V{1}) - x23*(x144 - x146*x167 - x150*x85);
auto x14 = -x145*(-x25*x77*x97 + V{1}) - x23*(x106*x76 + x151);
auto x15 = -x145*(-x25*x54*x97 + V{1}) - x23*(-x152*x68 - x153*x54 + x157);
auto x16 = -x145*(-x25*x62*x97 + V{1}) - x23*(-x152*x66 + x153*x61 + x158);
auto x17 = -x145*(-x168*x25*x97 + V{1}) - x23*(-x146*x168 - x150*x89 + x160);
auto x18 = -x145*(x170*x25*x97 + V{1}) - x23*(x146*x170 - x150*x81 + x161);
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
return { -V{1}*x104*(V{0.0555555555555555}*x105*x47 + V{0.333333333333333}*x105*x59 + V{0.0277777777777778}*x105*x68 + V{0.111111111111111}*x105*x72 + V{0.0277777777777778}*x105*(x58 + x65) + x171*x50 + x171*x89 + x171*(x169 + x80) + x172*x76 + x172*x85 - x173*x54 + x173*x61 - x174*x40 - x174*x45 + x96), x30 + x32 + x35 };
}
};

}

}
