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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x23 = x22 + V{-1};
auto x24 = V{0.5}*cell[11];
auto x25 = V{0.5}*cell[12];
auto x26 = V{2}*cell[17];
auto x27 = V{1} / (x20 + V{1});
auto x28 = V{0.25}*cell[0] + V{0.25}*cell[10] + V{0.5}*cell[17] + V{0.25}*cell[1] + V{0.25};
auto x29 = V{1} / (x21 + V{1});
auto x30 = cell[0] + cell[10] + cell[1] + x26 + V{1};
auto x31 = V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[18] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + x30;
auto x32 = cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x30;
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
auto x44 = V{3}*x20;
auto x45 = x20 + x21;
auto x46 = V{4.5}*(x45*x45);
auto x47 = -x34;
auto x48 = V{1} - x38;
auto x49 = x47 + x48;
auto x50 = V{3}*x21;
auto x51 = -x36;
auto x52 = x50 + x51;
auto x53 = x44 + x46 + x49 + x52;
auto x54 = x40 + x50;
auto x55 = x44 - x46 + x54;
auto x56 = -x21;
auto x57 = x19 + x56;
auto x58 = -x57;
auto x59 = -x50;
auto x60 = V{3}*x19;
auto x61 = x40 + x60;
auto x62 = x59 + x61;
auto x63 = x62 - V{4.5}*x58*x58;
auto x64 = V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32;
auto x65 = x19 - x20;
auto x66 = -x65;
auto x67 = -x44;
auto x68 = x61 + x67;
auto x69 = x68 - V{4.5}*x66*x66;
auto x70 = V{1.11022302462516e-16}*cell[12];
auto x71 = x19 + x21;
auto x72 = V{4.5}*(x71*x71);
auto x73 = x54 + x60 - x72;
auto x74 = -x43*x73;
auto x75 = x49 + x60;
auto x76 = x52 + x72 + x75;
auto x77 = x43*x76;
auto x78 = -x60;
auto x79 = -V{4.5}*x57*x57;
auto x80 = x54 + x78 + x79;
auto x81 = -x43*x80;
auto x82 = V{3}*x33;
auto x83 = x48 + x52 + x82;
auto x84 = x64*x83;
auto x85 = x19 + x20;
auto x86 = V{4.5}*(x85*x85);
auto x87 = x44 + x51;
auto x88 = x75 + x86 + x87;
auto x89 = V{0.0555555555555556}*x27*x31 + V{0.0555555555555556}*x29*x32;
auto x90 = V{3}*x37;
auto x91 = x47 + x87 + x90 + V{1};
auto x92 = V{5.55111512312578e-17}*x27*x31 + V{5.55111512312578e-17}*x29*x32;
auto x93 = V{1.11022302462516e-16}*x27*x31 + V{1.11022302462516e-16}*x29*x32;
auto x94 = x39 + x50 - x82;
auto x95 = x64*x94;
auto x96 = x20 + x56;
auto x97 = -V{4.5}*x96*x96;
auto x98 = x54 + x67 + x97;
auto x99 = V{3}*x35;
auto x100 = x34 + V{-1};
auto x101 = x100 + x38 + x60 - x99;
auto x102 = -x101*x64;
auto x103 = x75 + x99;
auto x104 = x103*x64;
auto x105 = V{1.11022302462516e-16}*cell[11];
auto x106 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17};
auto x107 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] + x102 + x104 - x27*(V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x105 + x106) - x29*(V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[7] + V{1.11022302462516e-16}*cell[9] + x106 + x70) - x35*(V{8.32667268468867e-17}*x27*x31 + V{8.32667268468867e-17}*x29*x32) - x41 + V{-2.22044604925031e-16};
auto x108 = V{2.22044604925031e-16}*cell[11] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{8.88178419700125e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x107 - x33*x92 - x37*x93 + x53*(V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32) - x55*(V{6.16790569236198e-18}*x27*x31 + V{6.16790569236198e-18}*x29*x32) + x64*x88 - x64*x98 + x70 + x74 + x77 + x81 + x84 + x89*x91 - x95;
auto x109 = -x96;
auto x110 = x40 + x44;
auto x111 = x110 + x59;
auto x112 = x111 - V{4.5}*x109*x109;
auto x113 = x44 + x61 - x86;
auto x114 = -x113*x43;
auto x115 = x43*x88;
auto x116 = -V{4.5}*x65*x65;
auto x117 = x110 + x116 + x78;
auto x118 = -x117*x43;
auto x119 = x64*x91;
auto x120 = x100 + x36 + x44 - x90;
auto x121 = x120*x64;
auto x122 = V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x105 + x107 + x114 + x115 + x118 + x119 - x121 - x33*x93 - x37*x92 + x53*(V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32) - x55*(V{3.08395284618099e-18}*x27*x31 + V{3.08395284618099e-18}*x29*x32) + x64*x76 + x83*x89;
auto x123 = x27*(x108 - x43*x63 - x64*x69) + x29*(-x112*x64 + x122 - x43*x69 - x63*x64);
auto x124 = V{0.0833333333333334}*cell[13];
auto x125 = V{0.0833333333333334}*cell[14];
auto x126 = V{0.0833333333333334}*cell[15];
auto x127 = V{0.0833333333333334}*cell[16];
auto x128 = V{0.0833333333333334}*cell[4];
auto x129 = V{0.0833333333333334}*cell[5];
auto x130 = V{0.0833333333333334}*cell[6];
auto x131 = V{0.0833333333333334}*cell[7];
auto x132 = V{0.0833333333333333}*cell[11];
auto x133 = V{0.0833333333333333}*cell[12];
auto x134 = V{0.0833333333333333}*cell[2];
auto x135 = V{0.0833333333333333}*cell[3];
auto x136 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18};
auto x137 = x27*(V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x136);
auto x138 = x29*(V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x136);
auto x139 = V{0.00231481481481481}*x27*x31 + V{0.00231481481481481}*x29*x32;
auto x140 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x141 = x140*x37;
auto x142 = x140*x33;
auto x143 = V{0.0833333333333333}*x27*x31 + V{0.0833333333333333}*x29*x32;
auto x144 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[9] + x124 + x125 + x126 + x127 + x128 + x129 + x130 + x131 - x132 - x133 - x134 - x135 - x137 - x138 + x139*x53 + x139*x55 + x141 + x142 - x143*x35 + V{-0.0555555555555555};
auto x145 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x146 = V{0.0833333333333333}*cell[10] - V{0.166666666666667}*cell[17] - V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - V{0.0833333333333334}*cell[9] + x137 + x138 - x140*x35 + x145*x53 + x145*x55 + V{0.0555555555555555};
auto x147 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x124 - x125 - x128 - x129 + x133 + x135 - x142 + x143*x37 + x146;
auto x148 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x126 - x127 - x130 - x131 + x132 + x134 - x141 + x143*x33 + x146;
auto x149 = x27*x31 + x29*x32;
auto x150 = V{0.125}*x149*x19;
auto x151 = x150*x20;
auto x152 = V{0.0208333333333333}*x27*x31 + V{0.0208333333333333}*x29*x32;
auto x153 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333};
auto x154 = -x27*(V{0.0416666666666667}*cell[11] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x153);
auto x155 = -x29*(V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[12] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[7] + V{0.0416666666666667}*cell[9] + x153);
auto x156 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x157 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x154 + x155 - x156*x37 + V{0.0138888888888889};
auto x158 = V{0.000578703703703704}*x27*x31 + V{0.000578703703703704}*x29*x32;
auto x159 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[9] - x156*x35 - x158*x53 - x158*x55;
auto x160 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x152*x33 + x157 + x159;
auto x161 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x151 + x160;
auto x162 = x116 + x68;
auto x163 = -x162*x43;
auto x164 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x151 + x160;
auto x165 = x150*x21;
auto x166 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x156*x33;
auto x167 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x152*x37 + x154 + x155 + x159 + x166 + V{0.0138888888888889};
auto x168 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x165 + x167;
auto x169 = x62 + x79;
auto x170 = -x169*x43;
auto x171 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x165 + x167;
auto x172 = V{0.00578703703703704}*x27*x31 + V{0.00578703703703704}*x29*x32;
auto x173 = V{0.125}*x149*x20*x21;
auto x174 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x152*x35 + x157 + x166;
auto x175 = V{0.833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[9] - x173 + x174;
auto x176 = x111 + x97;
auto x177 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x178 = -V{0.166666666666667}*cell[17] + V{0.416666666666667}*cell[18] + V{0.416666666666667}*cell[9] + x173 + x174 + x177*x53 + x177*x55;
auto x179 = x27*(x108 - x162*x64 + x170) + x29*(x122 + x163 - x169*x64 - x176*x64);
auto x0 = -x22*(V{0.166666666666667}*x123*x40 + V{0.333333333333333}) + x23*(V{0.5}*cell[10] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[9] + x24 + x25 + x26 - x27*(V{0.25}*cell[12] + V{0.5}*cell[13] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[18] + V{0.25}*cell[3] + V{0.5}*cell[5] + V{0.25}*cell[6] + V{0.25}*cell[7] + x24 + x28) - x29*(V{0.25}*cell[11] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[15] + V{0.25}*cell[2] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[7] + V{0.5}*cell[9] + x25 + x28) - x33*x42 - x35*x42 - x37*x42 + x41 - x43*x53 - x43*x55 + V{0.833333333333333});
auto x1 = -x22*(V{0.0277777777777778}*x101*x123 + V{0.0555555555555556}) - x23*(x102 + x144);
auto x2 = -x22*(V{0.0277777777777778}*x120*x123 + V{0.0555555555555556}) + x23*(x121 + x147);
auto x3 = -x22*(V{0.0277777777777778}*x123*x94 + V{0.0555555555555556}) + x23*(x148 + x95);
auto x4 = -x22*(V{0.0138888888888889}*x113*x123 + V{0.0277777777777778}) - x23*(x114 + x161);
auto x5 = -x22*(V{0.0138888888888889}*x123*x69 + V{0.0277777777777778}) - x23*(x163 + x164);
auto x6 = -x22*(V{0.0138888888888889}*x123*x73 + V{0.0277777777777778}) - x23*(x168 + x74);
auto x7 = -x22*(V{0.0138888888888889}*x123*x63 + V{0.0277777777777778}) - x23*(x170 + x171);
auto x8 = -x22*(V{0.0138888888888889}*x123*x55 + V{0.0277777777777778}) - x23*(-x172*x53 + x175 - x55*(V{0.0196759259259259}*x27*x31 + V{0.0196759259259259}*x29*x32));
auto x9 = -x22*(V{0.0138888888888889}*x112*x123 + V{0.0277777777777778}) - x23*(-x176*x43 + x178);
auto x10 = x22*(V{0.0277777777777778}*x103*x179 + V{-0.0555555555555556}) - x23*(x104 + x144);
auto x11 = x22*(V{0.0277777777777778}*x179*x91 + V{-0.0555555555555556}) + x23*(-x119 + x147);
auto x12 = x22*(V{0.0277777777777778}*x179*x83 + V{-0.0555555555555556}) + x23*(x148 - x84);
auto x13 = x22*(V{0.0138888888888889}*x179*x88 + V{-0.0277777777777778}) - x23*(x115 + x161);
auto x14 = -x22*(V{0.0138888888888889}*x117*x123 + V{0.0277777777777778}) - x23*(x118 + x164);
auto x15 = x22*(V{0.0138888888888889}*x179*x76 + V{-0.0277777777777778}) - x23*(x168 + x77);
auto x16 = -x22*(V{0.0138888888888889}*x123*x80 + V{0.0277777777777778}) - x23*(x171 + x81);
auto x17 = x22*(V{0.0138888888888889}*x179*x53 + V{-0.0277777777777778}) - x23*(-x172*x55 + x175 + x53*(V{0.00810185185185185}*x27*x31 + V{0.00810185185185185}*x29*x32));
auto x18 = -x22*(V{0.0138888888888889}*x123*x98 + V{0.0277777777777778}) - x23*(x178 - x43*x98);
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
return { V{0.5}*x179, x33 + x35 + x37 };
}
};

}

}
