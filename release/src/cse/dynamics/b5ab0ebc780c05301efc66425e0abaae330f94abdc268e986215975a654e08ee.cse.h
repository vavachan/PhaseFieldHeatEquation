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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x23 = x22 + V{-1};
auto x24 = V{0.5}*cell[10];
auto x25 = V{0.5}*cell[11];
auto x26 = V{2}*cell[13];
auto x27 = V{1} / (x19 + V{1});
auto x28 = V{0.25}*cell[0] + V{0.25}*cell[12] + V{0.5}*cell[13] + V{0.25}*cell[3] + V{0.25};
auto x29 = V{1} / (x20 + V{1});
auto x30 = cell[0] + cell[12] + cell[3] + x26 + V{1};
auto x31 = V{2}*cell[10] + cell[11] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] + x30;
auto x32 = cell[10] + V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + V{2}*cell[5] + cell[6] + cell[7] + x30;
auto x33 = x21*x21;
auto x34 = V{1.5}*x33;
auto x35 = x19*x19;
auto x36 = V{1.5}*x35;
auto x37 = x20*x20;
auto x38 = V{1.5}*x37;
auto x39 = x36 + x38 + V{-1};
auto x40 = x34 + x39;
auto x41 = x40*(V{0.166666666666667}*x27*x31 + V{0.166666666666667}*x29*x32);
auto x42 = V{0.25}*x27*x31 + V{0.25}*x29*x32;
auto x43 = V{0.0138888888888889}*x27*x31 + V{0.0138888888888889}*x29*x32;
auto x44 = x19 + x20;
auto x45 = V{4.5}*(x44*x44);
auto x46 = V{3}*x20;
auto x47 = -x34;
auto x48 = V{1} - x36;
auto x49 = x47 + x48;
auto x50 = x46 + x49;
auto x51 = -x38;
auto x52 = V{3}*x19;
auto x53 = x51 + x52;
auto x54 = x45 + x50 + x53;
auto x55 = x40 + x46;
auto x56 = -x45 + x52 + x55;
auto x57 = -x21;
auto x58 = x20 + x57;
auto x59 = -x58;
auto x60 = V{3}*x21;
auto x61 = -x60;
auto x62 = x55 + x61;
auto x63 = x62 - V{4.5}*x59*x59;
auto x64 = x20 + x21;
auto x65 = V{4.5}*(x64*x64);
auto x66 = x55 + x60 - x65;
auto x67 = -x43*x66;
auto x68 = x51 + x60;
auto x69 = x50 + x65 + x68;
auto x70 = x43*x69;
auto x71 = -x46;
auto x72 = -V{4.5}*x58*x58;
auto x73 = x40 + x60;
auto x74 = x71 + x72 + x73;
auto x75 = -x43*x74;
auto x76 = V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32;
auto x77 = V{3}*x37;
auto x78 = x50 + x77;
auto x79 = x76*x78;
auto x80 = x19 + x21;
auto x81 = V{4.5}*(x80*x80);
auto x82 = x49 + x52 + x68 + x81;
auto x83 = V{0.0555555555555556}*x27*x31 + V{0.0555555555555556}*x29*x32;
auto x84 = V{3}*x35;
auto x85 = x47 + x53 + x84 + V{1};
auto x86 = V{1.11022302462516e-16}*x27*x31 + V{1.11022302462516e-16}*x29*x32;
auto x87 = V{8.32667268468867e-17}*x27*x31 + V{8.32667268468867e-17}*x29*x32;
auto x88 = x34 + V{-1};
auto x89 = x36 + x46 - x77 + x88;
auto x90 = x76*x89;
auto x91 = -x52;
auto x92 = x19 - x20;
auto x93 = -V{4.5}*x92*x92;
auto x94 = x55 + x91 + x93;
auto x95 = x19 + x57;
auto x96 = -V{4.5}*x95*x95;
auto x97 = x73 + x91 + x96;
auto x98 = V{3}*x33;
auto x99 = x39 + x60 - x98;
auto x100 = -x76*x99;
auto x101 = x48 + x68 + x98;
auto x102 = x101*x76;
auto x103 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17};
auto x104 = V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[5] + x100 + x102 - x27*(V{1.11022302462516e-16}*cell[10] + V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x103) - x29*(V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x103) - x33*(V{5.55111512312578e-17}*x27*x31 + V{5.55111512312578e-17}*x29*x32) - x41 + x54*(V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32) - x56*(V{4.62592926927149e-18}*x27*x31 + V{4.62592926927149e-18}*x29*x32) + V{-2.22044604925031e-16};
auto x105 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x104 - x35*x86 - x37*x87 + x67 + x70 + x75 + x76*x82 - x76*x94 - x76*x97 + x79 + x83*x85 - x90;
auto x106 = -x95;
auto x107 = x40 + x52;
auto x108 = x107 + x61;
auto x109 = x108 - V{4.5}*x106*x106;
auto x110 = -x92;
auto x111 = x107 + x71;
auto x112 = x111 - V{4.5}*x110*x110;
auto x113 = x52 + x73 - x81;
auto x114 = -x113*x43;
auto x115 = x43*x82;
auto x116 = -x43*x97;
auto x117 = x76*x85;
auto x118 = x38 + x52 - x84 + x88;
auto x119 = x118*x76;
auto x120 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x104 + x114 + x115 + x116 + x117 - x119 - x35*x87 - x37*x86 + x69*x76 - x74*x76 + x78*x83;
auto x121 = x27*(x105 - x43*x63) + x29*(-x109*x43 - x112*x76 + x120);
auto x122 = V{0.0833333333333333}*cell[11];
auto x123 = V{0.0833333333333333}*cell[2];
auto x124 = V{0.0833333333333334}*cell[15];
auto x125 = V{0.0833333333333334}*cell[16];
auto x126 = V{0.0833333333333334}*cell[6];
auto x127 = V{0.0833333333333334}*cell[7];
auto x128 = V{0.0833333333333333}*x27*x31 + V{0.0833333333333333}*x29*x32;
auto x129 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x130 = x129*x37;
auto x131 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[3] + V{3.46944695195361e-18};
auto x132 = x27*(V{6.93889390390723e-18}*cell[10] + V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x131);
auto x133 = x29*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x131);
auto x134 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x135 = V{0.0833333333333333}*cell[12] - V{0.166666666666667}*cell[13] - V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - V{0.0833333333333334}*cell[5] - x129*x33 + x132 + x133 + x134*x54 + x134*x56 + V{0.0555555555555555};
auto x136 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x122 + x123 - x124 - x125 - x126 - x127 + x128*x35 - x130 + x135;
auto x137 = V{0.0833333333333333}*cell[10];
auto x138 = V{0.0833333333333333}*cell[1];
auto x139 = V{0.0833333333333334}*cell[17];
auto x140 = V{0.0833333333333334}*cell[18];
auto x141 = V{0.0833333333333334}*cell[8];
auto x142 = V{0.0833333333333334}*cell[9];
auto x143 = x129*x35;
auto x144 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x128*x37 + x135 + x137 + x138 - x139 - x140 - x141 - x142 - x143;
auto x145 = V{0.00231481481481481}*x27*x31 + V{0.00231481481481481}*x29*x32;
auto x146 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[5] - x122 - x123 + x124 + x125 + x126 + x127 - x128*x33 + x130 - x132 - x133 - x137 - x138 + x139 + x140 + x141 + x142 + x143 + x145*x54 + x145*x56 + V{-0.0555555555555555};
auto x147 = V{0.00578703703703704}*x27*x31 + V{0.00578703703703704}*x29*x32;
auto x148 = x27*x31 + x29*x32;
auto x149 = V{0.125}*x148*x19;
auto x150 = x149*x20;
auto x151 = V{0.0208333333333333}*x27*x31 + V{0.0208333333333333}*x29*x32;
auto x152 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[3] + V{0.0208333333333333};
auto x153 = -x27*(V{0.0416666666666667}*cell[10] + V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x152);
auto x154 = -x29*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x152);
auto x155 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x156 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x153 + x154 - x155*x35 + V{0.0138888888888889};
auto x157 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x155*x37;
auto x158 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x151*x33 + x156 + x157;
auto x159 = V{0.833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[5] - x150 + x158;
auto x160 = x111 + x93;
auto x161 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x162 = -V{0.166666666666667}*cell[13] + V{0.416666666666667}*cell[14] + V{0.416666666666667}*cell[5] + x150 + x158 + x161*x54 + x161*x56;
auto x163 = x149*x21;
auto x164 = V{0.000578703703703704}*x27*x31 + V{0.000578703703703704}*x29*x32;
auto x165 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[5] - x155*x33 - x164*x54 - x164*x56;
auto x166 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x151*x37 + x156 + x165;
auto x167 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x163 + x166;
auto x168 = -x43*(x108 + x96);
auto x169 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x163 + x166;
auto x170 = V{0.125}*x148*x20*x21;
auto x171 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x151*x35 + x153 + x154 + x157 + x165 + V{0.0138888888888889};
auto x172 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x170 + x171;
auto x173 = -x43*(x62 + x72);
auto x174 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x170 + x171;
auto x175 = x27*(x105 + x173) + x29*(x120 - x160*x76 + x168);
auto x0 = -x22*(V{0.166666666666667}*x121*x40 + V{0.333333333333333}) + x23*(V{0.5}*cell[12] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x24 + x25 + x26 - x27*(V{0.25}*cell[11] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[2] + V{0.25}*cell[8] + V{0.25}*cell[9] + x24 + x28) - x29*(V{0.25}*cell[10] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.5}*cell[5] + V{0.25}*cell[6] + V{0.25}*cell[7] + x25 + x28) - x33*x42 - x35*x42 - x37*x42 + x41 - x43*x54 - x43*x56 + V{0.833333333333333});
auto x1 = -x22*(V{0.0277777777777778}*x118*x121 + V{0.0555555555555556}) + x23*(x119 + x136);
auto x2 = -x22*(V{0.0277777777777778}*x121*x89 + V{0.0555555555555556}) + x23*(x144 + x90);
auto x3 = -x22*(V{0.0277777777777778}*x121*x99 + V{0.0555555555555556}) - x23*(x100 + x146);
auto x4 = -x22*(V{0.0138888888888889}*x121*x56 + V{0.0277777777777778}) - x23*(-x147*x54 + x159 - x56*(V{0.0196759259259259}*x27*x31 + V{0.0196759259259259}*x29*x32));
auto x5 = -x22*(V{0.0138888888888889}*x112*x121 + V{0.0277777777777778}) - x23*(-x160*x43 + x162);
auto x6 = -x22*(V{0.0138888888888889}*x113*x121 + V{0.0277777777777778}) - x23*(x114 + x167);
auto x7 = -x22*(V{0.0138888888888889}*x109*x121 + V{0.0277777777777778}) - x23*(x168 + x169);
auto x8 = -x22*(V{0.0138888888888889}*x121*x66 + V{0.0277777777777778}) - x23*(x172 + x67);
auto x9 = -x22*(V{0.0138888888888889}*x121*x63 + V{0.0277777777777778}) - x23*(x173 + x174);
auto x10 = x22*(V{0.0277777777777778}*x175*x85 + V{-0.0555555555555556}) + x23*(-x117 + x136);
auto x11 = x22*(V{0.0277777777777778}*x175*x78 + V{-0.0555555555555556}) + x23*(x144 - x79);
auto x12 = x22*(V{0.0277777777777778}*x101*x175 + V{-0.0555555555555556}) - x23*(x102 + x146);
auto x13 = x22*(V{0.0138888888888889}*x175*x54 + V{-0.0277777777777778}) - x23*(-x147*x56 + x159 + x54*(V{0.00810185185185185}*x27*x31 + V{0.00810185185185185}*x29*x32));
auto x14 = -x22*(V{0.0138888888888889}*x121*x94 + V{0.0277777777777778}) - x23*(x162 - x43*x94);
auto x15 = x22*(V{0.0138888888888889}*x175*x82 + V{-0.0277777777777778}) - x23*(x115 + x167);
auto x16 = -x22*(V{0.0138888888888889}*x121*x97 + V{0.0277777777777778}) - x23*(x116 + x169);
auto x17 = x22*(V{0.0138888888888889}*x175*x69 + V{-0.0277777777777778}) - x23*(x172 + x70);
auto x18 = -x22*(V{0.0138888888888889}*x121*x74 + V{0.0277777777777778}) - x23*(x174 + x75);
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
return { V{0.5}*x175, x33 + x35 + x37 };
}
};

}

}
