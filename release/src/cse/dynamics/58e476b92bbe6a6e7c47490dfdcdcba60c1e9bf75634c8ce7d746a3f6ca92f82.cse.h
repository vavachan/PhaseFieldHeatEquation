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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x23 = x22 + V{-1};
auto x24 = V{0.5}*cell[11];
auto x25 = V{0.5}*cell[1];
auto x26 = V{2}*cell[5];
auto x27 = V{1} / (x19 + V{-1});
auto x28 = V{0.25}*cell[0] + V{0.25}*cell[12] + V{0.25}*cell[3] + V{0.5}*cell[5] + V{0.25};
auto x29 = V{1} / (x20 + V{1});
auto x30 = cell[0] + cell[12] + cell[3] + x26 + V{1};
auto x31 = cell[10] + V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[6] + cell[7] + x30;
auto x32 = cell[11] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + V{2}*cell[4] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x30;
auto x33 = -V{0.0138888888888889}*x27*x32 + V{0.0138888888888889}*x29*x31;
auto x34 = x19 - x20;
auto x35 = -x34;
auto x36 = V{3}*x20;
auto x37 = -x36;
auto x38 = V{3}*x19;
auto x39 = x21*x21;
auto x40 = V{1.5}*x39;
auto x41 = x19*x19;
auto x42 = V{1.5}*x41;
auto x43 = x20*x20;
auto x44 = V{1.5}*x43;
auto x45 = x42 + x44 + V{-1};
auto x46 = x40 + x45;
auto x47 = x38 + x46;
auto x48 = x37 + x47;
auto x49 = x48 - V{4.5}*x35*x35;
auto x50 = x46*(-V{0.166666666666667}*x27*x32 + V{0.166666666666667}*x29*x31);
auto x51 = -V{0.25}*x27*x32 + V{0.25}*x29*x31;
auto x52 = -V{4.5}*x34*x34;
auto x53 = -x38 + x46;
auto x54 = x36 + x52 + x53;
auto x55 = -x21;
auto x56 = x19 + x55;
auto x57 = -x56;
auto x58 = V{3}*x21;
auto x59 = -x58;
auto x60 = x47 + x59;
auto x61 = x60 - V{4.5}*x57*x57;
auto x62 = -V{0.0277777777777778}*x27*x32 + V{0.0277777777777778}*x29*x31;
auto x63 = -V{0.0277777777777778}*x27*x32 + V{0.0277777777777778}*x29*x31;
auto x64 = V{3}*x39;
auto x65 = x45 + x58 - x64;
auto x66 = -x63*x65;
auto x67 = V{1} - x44;
auto x68 = -x42;
auto x69 = x58 + x68;
auto x70 = x64 + x67 + x69;
auto x71 = x63*x70;
auto x72 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17};
auto x73 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{6.66133814775094e-16}*cell[5] + x27*(V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[6] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x72) - x29*(V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x72) - x39*(-V{5.55111512312578e-17}*x27*x32 + V{5.55111512312578e-17}*x29*x31) - x50 - x54*(-V{4.62592926927149e-18}*x27*x32 + V{4.62592926927149e-18}*x29*x31) + x66 + x71 + V{-2.22044604925031e-16};
auto x74 = -x49*x62 + x73;
auto x75 = V{3}*x41;
auto x76 = x40 + V{-1};
auto x77 = x38 + x44 - x75 + x76;
auto x78 = -x63*x77;
auto x79 = -x40;
auto x80 = x67 + x79;
auto x81 = x38 + x80;
auto x82 = x75 + x81;
auto x83 = x63*x82;
auto x84 = x19 + x21;
auto x85 = V{4.5}*(x84*x84);
auto x86 = x47 + x58 - x85;
auto x87 = -x33*x86;
auto x88 = x69 + x81 + x85;
auto x89 = x33*x88;
auto x90 = -V{4.5}*x56*x56;
auto x91 = x53 + x58 + x90;
auto x92 = -x33*x91;
auto x93 = x19 + x20;
auto x94 = V{4.5}*(x93*x93);
auto x95 = x36 + x68;
auto x96 = x81 + x94 + x95;
auto x97 = x20 + x21;
auto x98 = V{4.5}*(x97*x97);
auto x99 = x36 + x69 + x80 + x98;
auto x100 = -V{0.0555555555555556}*x27*x32 + V{0.0555555555555556}*x29*x31;
auto x101 = V{3}*x43;
auto x102 = x101 + x79 + x95 + V{1};
auto x103 = -V{1.11022302462516e-16}*x27*x32 + V{1.11022302462516e-16}*x29*x31;
auto x104 = -V{8.32667268468867e-17}*x27*x32 + V{8.32667268468867e-17}*x29*x31;
auto x105 = x20 + x55;
auto x106 = -V{4.5}*x105*x105;
auto x107 = x46 + x58;
auto x108 = x106 + x107 + x37;
auto x109 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x100*x102 - x103*x43 - x104*x41 - x108*x63 + x63*x96 + x63*x99 + x78 + x83 + x87 + x89 + x92;
auto x110 = -x105;
auto x111 = x36 + x46 + x59;
auto x112 = x111 - V{4.5}*x110*x110;
auto x113 = x102*x63;
auto x114 = x107 + x36 - x98;
auto x115 = -x114*x33;
auto x116 = x33*x99;
auto x117 = -x108*x33;
auto x118 = -x101 + x36 + x42 + x76;
auto x119 = x118*x63;
auto x120 = x36 + x47 - x94;
auto x121 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x100*x77 - x103*x41 - x104*x43 + x113 + x115 + x116 + x117 - x119 - x120*x63 - x63*x86;
auto x122 = -x27*(-x112*x33 + x121 - x61*x63 + x74) + x29*(x109 - x33*x61 + x74);
auto x123 = -V{0.00115740740740741}*x27*x32 + V{0.00115740740740741}*x29*x31;
auto x124 = x48 + x52;
auto x125 = -V{0.0833333333333333}*x27*x32 + V{0.0833333333333333}*x29*x31;
auto x126 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18};
auto x127 = x27*(V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x126);
auto x128 = x29*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x126);
auto x129 = -x128;
auto x130 = -V{0.0416666666666667}*x27*x32 + V{0.0416666666666667}*x29*x31;
auto x131 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x127 + x129 + x130*x43 + V{-0.0555555555555555};
auto x132 = V{0.0833333333333334}*cell[13];
auto x133 = V{0.0833333333333334}*cell[4];
auto x134 = V{0.166666666666667}*cell[5];
auto x135 = V{0.0833333333333333}*cell[12];
auto x136 = V{0.0833333333333333}*cell[3];
auto x137 = x130*x39;
auto x138 = x123*x54;
auto x139 = x132 + x133 + x134 - x135 - x136 + x137 - x138;
auto x140 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] + x123*x124 - x125*x41 + x131 + x139;
auto x141 = V{0.0833333333333333}*cell[10];
auto x142 = V{0.0833333333333333}*cell[1];
auto x143 = V{0.166666666666667}*cell[11];
auto x144 = V{0.166666666666667}*cell[2];
auto x145 = V{0.0833333333333334}*cell[17];
auto x146 = V{0.0833333333333334}*cell[18];
auto x147 = V{0.0833333333333334}*cell[8];
auto x148 = V{0.0833333333333334}*cell[9];
auto x149 = x123*x49;
auto x150 = x130*x41;
auto x151 = -V{0.00231481481481481}*x27*x32 + V{0.00231481481481481}*x29*x31;
auto x152 = -x141 - x142 + x145 + x146 + x147 + x148 + x150;
auto x153 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.333333333333333}*cell[5] - x124*x151 - x125*x39 + x131 + x151*x54 + x152;
auto x154 = -V{0.00115740740740741}*x27*x32 + V{0.00115740740740741}*x29*x31;
auto x155 = -x27*x32 + x29*x31;
auto x156 = V{0.125}*x155*x19;
auto x157 = x156*x20;
auto x158 = -V{0.0208333333333333}*x27*x32 + V{0.0208333333333333}*x29*x31;
auto x159 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[12] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333};
auto x160 = x27*(V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x159);
auto x161 = -x29*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x159);
auto x162 = -V{0.0416666666666667}*x27*x32 + V{0.0416666666666667}*x29*x31;
auto x163 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x160 + x161 - x162*x41 + V{0.0138888888888889};
auto x164 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x162*x43;
auto x165 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x158*x39 + x163 + x164;
auto x166 = V{0.416666666666667}*cell[13] + V{0.416666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x124*x154 + x154*x54 - x157 + x165;
auto x167 = -V{0.00578703703703704}*x27*x32 + V{0.00578703703703704}*x29*x31;
auto x168 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[4] + V{0.833333333333333}*cell[5] + x157 + x165;
auto x169 = x156*x21;
auto x170 = -V{0.000578703703703704}*x27*x32 + V{0.000578703703703704}*x29*x31;
auto x171 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0833333333333334}*cell[5] + x124*x170 - x162*x39 - x170*x54;
auto x172 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x158*x43 + x163 + x171;
auto x173 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x169 + x172;
auto x174 = x60 + x90;
auto x175 = -x174*x33;
auto x176 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x169 + x172;
auto x177 = V{0.125}*x155*x20*x21;
auto x178 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x158*x41 + x160 + x161 + x164 + x171 + V{0.0138888888888889};
auto x179 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x177 + x178;
auto x180 = -x33*(x106 + x111);
auto x181 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x177 + x178;
auto x182 = -x124*x62 + x73;
auto x183 = -x27*(x121 - x174*x63 + x180 + x182) + x29*(x109 + x175 + x182);
auto x0 = -x22*(V{0.166666666666667}*x122*x46 + V{0.333333333333333}) + x23*(V{0.5}*cell[10] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x24 + x25 + x26 + x27*(V{0.25}*cell[11] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[2] + V{0.5}*cell[4] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.25}*cell[8] + V{0.25}*cell[9] + x25 + x28) - x29*(V{0.25}*cell[10] + V{0.5}*cell[13] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.25}*cell[6] + V{0.25}*cell[7] + x24 + x28) + x33*x49 - x33*x54 - x39*x51 - x41*x51 - x43*x51 + x50 + V{0.833333333333333});
auto x1 = -x22*(V{0.0277777777777778}*x122*x77 + V{0.0555555555555556}) - x23*(x140 + x78);
auto x2 = -x22*(V{0.0277777777777778}*x118*x122 + V{0.0555555555555556}) + x23*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x119 + x125*x43 - x127 + x128 - x132 - x133 - x134 + x135 + x136 - x137 + x138 + x141 + x142 - x143 - x144 - x145 - x146 - x147 - x148 - x149 - x150 + V{0.0555555555555555});
auto x3 = -x22*(V{0.0277777777777778}*x122*x65 + V{0.0555555555555556}) - x23*(x153 + x66);
auto x4 = -x22*(V{0.0138888888888889}*x120*x122 + V{0.0277777777777778}) - x23*(-x120*x33 + x166);
auto x5 = -x22*(V{0.0138888888888889}*x122*x49 + V{0.0277777777777778}) - x23*(-x124*(-V{0.00810185185185185}*x27*x32 + V{0.00810185185185185}*x29*x31) - x167*x54 + x168);
auto x6 = -x22*(V{0.0138888888888889}*x122*x86 + V{0.0277777777777778}) - x23*(x173 + x87);
auto x7 = -x22*(V{0.0138888888888889}*x122*x61 + V{0.0277777777777778}) - x23*(x175 + x176);
auto x8 = -x22*(V{0.0138888888888889}*x114*x122 + V{0.0277777777777778}) - x23*(x115 + x179);
auto x9 = -x22*(V{0.0138888888888889}*x112*x122 + V{0.0277777777777778}) - x23*(x180 + x181);
auto x10 = x22*(V{0.0277777777777778}*x183*x82 + V{-0.0555555555555556}) - x23*(x140 + x83);
auto x11 = x22*(V{0.0277777777777778}*x102*x183 + V{-0.0555555555555556}) + x23*(V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x113 - x127 - x129 - x139 - x143 - x144 - x149 - x152 + V{0.0833333333333333}*x155*x43 + V{0.0555555555555555});
auto x12 = x22*(V{0.0277777777777778}*x183*x70 + V{-0.0555555555555556}) - x23*(x153 + x71);
auto x13 = x22*(V{0.0138888888888889}*x183*x96 + V{-0.0277777777777778}) - x23*(x166 + x33*x96);
auto x14 = -x22*(V{0.0138888888888889}*x122*x54 + V{0.0277777777777778}) - x23*(x124*x167 + x168 - x54*(-V{0.0196759259259259}*x27*x32 + V{0.0196759259259259}*x29*x31));
auto x15 = x22*(V{0.0138888888888889}*x183*x88 + V{-0.0277777777777778}) - x23*(x173 + x89);
auto x16 = -x22*(V{0.0138888888888889}*x122*x91 + V{0.0277777777777778}) - x23*(x176 + x92);
auto x17 = x22*(V{0.0138888888888889}*x183*x99 + V{-0.0277777777777778}) - x23*(x116 + x179);
auto x18 = -x22*(V{0.0138888888888889}*x108*x122 + V{0.0277777777777778}) - x23*(x117 + x181);
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
return { V{0.5}*x183, x39 + x41 + x43 };
}
};

}

}
