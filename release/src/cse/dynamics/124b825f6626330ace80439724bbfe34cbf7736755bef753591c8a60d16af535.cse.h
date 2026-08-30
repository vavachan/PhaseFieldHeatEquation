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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x23 = x22 + V{-1};
auto x24 = x19 + V{-1};
auto x25 = -1/x24;
auto x26 = V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7];
auto x27 = cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + cell[3] + cell[8] + cell[9] + x26 + V{1};
auto x28 = x25*x27;
auto x29 = V{0.0277777777777778}*x28;
auto x30 = x20*x20;
auto x31 = V{3}*x20;
auto x32 = x21*x21;
auto x33 = V{1.5}*x32;
auto x34 = -x33;
auto x35 = x19*x19;
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
auto x52 = x20 + x21;
auto x53 = x52*x52;
auto x54 = x39 + x44 + V{4.5}*x53;
auto x55 = -x31;
auto x56 = -x21;
auto x57 = x20 + x56;
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
auto x71 = V{3}*x19;
auto x72 = x33 - V{3}*x35 + x49 + x71;
auto x73 = -x72;
auto x74 = x19 - x20;
auto x75 = -x74;
auto x76 = x59 + x71;
auto x77 = x55 + x76;
auto x78 = x77 - V{4.5}*x75*x75;
auto x79 = -x78;
auto x80 = x19 + x56;
auto x81 = -x80;
auto x82 = x64 + x71;
auto x83 = x82 - V{4.5}*x81*x81;
auto x84 = -x83;
auto x85 = x19 + x20;
auto x86 = x85*x85;
auto x87 = x31 + x76 - V{4.5}*x86;
auto x88 = -x87;
auto x89 = x19 + x21;
auto x90 = x89*x89;
auto x91 = x60 + x71 - V{4.5}*x90;
auto x92 = -x91;
auto x93 = -x59;
auto x94 = V{0.333333333333333}*x28*x93;
auto x95 = V{0.0555555555555555}*x28;
auto x96 = V{0.0555555555555555}*x28;
auto x97 = V{0.0277777777777778}*x28;
auto x98 = V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{4.44089209850063e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{6.66133814775094e-16}*cell[4] + V{6.66133814775094e-16}*cell[5] + V{8.88178419700125e-16}*cell[6] + V{4.44089209850063e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{2.22044604925031e-16};
auto x99 = V{0.0555555555555555}*x28*x48 + V{0.0277777777777778}*x28*x67 + V{0.0277777777777778}*x28*x69 + V{0.111111111111111}*x28*x73 + x40*x70 + x45*x70 + x51*x95 + x54*x97 + x62*x97 + x79*x96 + x84*x95 + x88*x96 + x92*x95 + x94 + x98;
auto x100 = V{1} / (x24);
auto x101 = x100*x27;
auto x102 = V{0.00462962962962963}*x101;
auto x103 = V{0.0833333333333333}*cell[12];
auto x104 = V{0.0833333333333333}*cell[3];
auto x105 = V{0.166666666666667}*cell[4];
auto x106 = -x105;
auto x107 = V{0.166666666666667}*cell[5];
auto x108 = -x107;
auto x109 = x102*x45;
auto x110 = V{0.00462962962962963}*x101;
auto x111 = x110*x87;
auto x112 = x110*x78;
auto x113 = -V{0.00462962962962963}*x100*x27*x50 + x103 + x104 + x106 + x108 + x109 + x111 + x112 + V{0.0555555555555555};
auto x114 = V{0.0833333333333333}*cell[11];
auto x115 = V{0.0833333333333333}*cell[2];
auto x116 = V{0.166666666666667}*cell[6];
auto x117 = -x116;
auto x118 = V{0.166666666666667}*cell[7];
auto x119 = -x118;
auto x120 = x102*x40;
auto x121 = x110*x91;
auto x122 = x110*x83;
auto x123 = -V{0.00462962962962963}*x100*x27*x47 + x114 + x115 + x117 + x119 + x120 + x121 + x122;
auto x124 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.333333333333333}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x125 = V{0.0555555555555556}*x22;
auto x126 = V{0.333333333333333}*cell[6];
auto x127 = V{0.333333333333333}*cell[7];
auto x128 = V{0.00925925925925926}*x101;
auto x129 = V{0.0462962962962963}*x101;
auto x130 = V{0.00925925925925926}*x101;
auto x131 = V{0.166666666666667}*cell[1];
auto x132 = V{0.00231481481481482}*x101;
auto x133 = x132*x68;
auto x134 = x132*x66;
auto x135 = x132*x61;
auto x136 = V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x128*x72 - x131 + x132*x54 - x133 - x134 - x135 + V{-0.0555555555555555};
auto x137 = V{0.333333333333333}*cell[4];
auto x138 = V{0.333333333333333}*cell[5];
auto x139 = V{0.00462962962962963}*x101;
auto x140 = V{0.00231481481481481}*x101;
auto x141 = V{0.00462962962962963}*x101;
auto x142 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] - x132*x83 - x132*x91 + x141*x40 - x141*x47 + V{0.0138888888888889};
auto x143 = V{0.00115740740740741}*x101;
auto x144 = V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.166666666666667}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] - x130*x72 + x143*x54 - x143*x61 - x143*x66 - x143*x68;
auto x145 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x140*x45 + x140*x50 + x142 + x144;
auto x146 = V{0.833333333333333}*cell[4] - V{0.166666666666667}*cell[5] + x139*x78 + x145;
auto x147 = V{0.0277777777777778}*x22;
auto x148 = -V{0.166666666666667}*cell[4] + V{0.833333333333333}*cell[5] + x139*x87 + x145;
auto x149 = V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] - x132*x78 - x132*x87 + x141*x45 - x141*x50;
auto x150 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x140*x40 + x140*x47 + x144 + x149 + V{0.0138888888888889};
auto x151 = V{0.833333333333333}*cell[6] - V{0.166666666666667}*cell[7] + x139*x83 + x150;
auto x152 = -V{0.166666666666667}*cell[6] + V{0.833333333333333}*cell[7] + x139*x91 + x150;
auto x153 = V{0.0115740740740741}*x101;
auto x154 = V{0.0162037037037037}*x101;
auto x155 = V{0.00231481481481481}*x101;
auto x156 = -V{0.0833333333333333}*cell[1] + x102*x72 + x142 + x149;
auto x157 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x155*x61 + x155*x66 + x156;
auto x158 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x155*x54 + x155*x68 + x156;
auto x159 = V{0.00462962962962963}*x28;
auto x160 = V{0.00462962962962963}*x28;
auto x161 = x43 + x71;
auto x162 = x161 + x34 + V{3}*x35 + V{1};
auto x163 = -V{0.0833333333333334}*cell[17] - V{0.0833333333333334}*cell[18] - V{0.0833333333333334}*cell[8] - V{0.0833333333333334}*cell[9] - V{0.00231481481481482}*x100*x27*x54 - V{0.00925925925925926}*x100*x27*x72 + x131 + x133 + x134 + x135;
auto x164 = V{0.0277777777777778}*x101;
auto x165 = x161 + x39 + V{4.5}*x86;
auto x166 = V{0.0231481481481481}*x101;
auto x167 = -x71;
auto x168 = -V{4.5}*x74*x74;
auto x169 = x167 + x168 + x31 + x59;
auto x170 = x38 + x44 + x71 + V{4.5}*x90;
auto x171 = -V{4.5}*x80*x80;
auto x172 = x167 + x171 + x60;
auto x173 = V{0.0555555555555555}*x101;
auto x174 = V{0.0555555555555555}*x101;
auto x175 = V{0.0277777777777778}*x101;
auto x176 = V{0.0555555555555556}*x101;
auto x0 = -V{0.333333333333333}*x22*(-x25*x93*x99 + V{1}) + x23*(V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[17] + V{1}*cell[18] + V{1}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[8] + V{1}*cell[9] + x26 - x29*x40 - x29*x45 - x29*x48 - x29*x51 - x29*x54 - x29*x62 - x29*x67 - x29*x69 - x70*x73 - x70*x79 - x70*x84 - x70*x88 - x70*x92 - x94 + V{0.833333333333333});
auto x1 = -x125*(-x25*x73*x99 + V{1}) - x23*(V{0.00462962962962963}*x100*x27*x61 + V{0.00462962962962963}*x100*x27*x66 + V{0.00462962962962963}*x100*x27*x68 + V{0.037037037037037}*x100*x27*x72 - x102*x54 - x113 - x123 - x124);
auto x2 = -x125*(-x25*x48*x99 + V{1}) - x23*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x102*x50 - x103 - x104 + x105 + x107 - x109 - x111 - x112 - x126 - x127 + x128*x83 + x128*x91 + x129*x47 + x130*x40 + x136);
auto x3 = -x125*(-x25*x51*x99 + V{1}) - x23*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x102*x47 - x114 - x115 + x116 + x118 - x120 - x121 - x122 + x128*x78 + x128*x87 + x129*x50 + x130*x45 + x136 - x137 - x138);
auto x4 = -x147*(-x25*x88*x99 + V{1}) - x23*(x102*x87 + x146);
auto x5 = -x147*(-x25*x79*x99 + V{1}) - x23*(x102*x78 + x148);
auto x6 = -x147*(-x25*x92*x99 + V{1}) - x23*(x102*x91 + x151);
auto x7 = -x147*(-x25*x84*x99 + V{1}) - x23*(x102*x83 + x152);
auto x8 = -x147*(-x25*x69*x99 + V{1}) - x23*(x153*x54 + x154*x68 + x157);
auto x9 = -x147*(-x25*x67*x99 + V{1}) - x23*(-x153*x61 + x154*x66 + x158);
auto x10 = -x125*(-x162*x25*x99 + V{1}) + x23*(x103 + x104 + x106 + x108 + x114 + x115 + x117 + x119 + x124 + x159*x79 + x159*x84 + x159*x88 + x159*x92 - x160*x40 - x160*x45 - x160*x48 - x160*x51 - x160*x54 - x160*x62 - x160*x67 - x160*x69 - x162*x70 + V{0.0185185185185185}*x28*x73 + V{0.0555555555555555});
auto x11 = -x125*(-x25*x40*x99 + V{1}) - x23*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.00925925925925926}*x100*x27*x83 + V{0.00925925925925926}*x100*x27*x91 - x113 - x126 - x127 - x129*x40 - x130*x47 - x163);
auto x12 = -x125*(-x25*x45*x99 + V{1}) - x23*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.00925925925925926}*x100*x27*x78 + V{0.00925925925925926}*x100*x27*x87 - x123 - x129*x45 - x130*x50 - x137 - x138 - x163 + V{-0.0555555555555555});
auto x13 = -x147*(-x165*x25*x99 + V{1}) - x23*(x146 - x164*x165 - x166*x87);
auto x14 = -x147*(x169*x25*x99 + V{1}) - x23*(x148 + x164*x169 - x166*x78);
auto x15 = -x147*(-x170*x25*x99 + V{1}) - x23*(x151 - x164*x170 - x166*x91);
auto x16 = -x147*(x172*x25*x99 + V{1}) - x23*(x152 + x164*x172 - x166*x83);
auto x17 = -x147*(-x25*x54*x99 + V{1}) - x23*(-x153*x68 - x154*x54 + x157);
auto x18 = -x147*(-x25*x62*x99 + V{1}) - x23*(-x153*x66 + x154*x61 + x158);
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
return { -V{1}*x100*(V{0.0555555555555555}*x101*x47 + V{0.333333333333333}*x101*x59 + V{0.0277777777777778}*x101*x68 + V{0.111111111111111}*x101*x72 + V{0.0277777777777778}*x101*(x58 + x65) + x173*x50 + x173*x91 + x173*(x171 + x82) + x174*x87 + x174*(x168 + x77) - x175*x54 + x175*x61 - x176*x40 - x176*x45 + x98), x30 + x32 + x35 };
}
};

}

}
