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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x23 = x22 + V{-1};
auto x24 = V{0.5}*cell[2];
auto x25 = V{0.5}*cell[3];
auto x26 = V{2}*cell[8];
auto x27 = V{1} / (x20 + V{-1});
auto x28 = V{0.25}*cell[0] + V{0.25}*cell[10] + V{0.25}*cell[1] + V{0.5}*cell[8] + V{0.25};
auto x29 = V{1} / (x21 + V{-1});
auto x30 = cell[0] + cell[10] + cell[1] + x26 + V{1};
auto x31 = cell[12] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[9] + x30;
auto x32 = cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + x30;
auto x33 = V{0.25}*x27*x31 + V{0.25}*x29*x32;
auto x34 = x19*x19;
auto x35 = x20*x20;
auto x36 = x21*x21;
auto x37 = V{0.0138888888888889}*x27*x31 + V{0.0138888888888889}*x29*x32;
auto x38 = V{3}*x20;
auto x39 = x20 + x21;
auto x40 = V{4.5}*(x39*x39);
auto x41 = V{1.5}*x36;
auto x42 = -x41;
auto x43 = V{1.5}*x35;
auto x44 = V{1} - x43;
auto x45 = x42 + x44;
auto x46 = V{3}*x21;
auto x47 = V{1.5}*x34;
auto x48 = -x47;
auto x49 = x46 + x48;
auto x50 = x38 + x40 + x45 + x49;
auto x51 = x43 + x47 + V{-1};
auto x52 = x41 + x51;
auto x53 = x46 + x52;
auto x54 = x38 - x40 + x53;
auto x55 = x52*(V{0.166666666666667}*x27*x31 + V{0.166666666666667}*x29*x32);
auto x56 = V{1.11022302462516e-16}*cell[4];
auto x57 = V{3}*x19;
auto x58 = x19 + x21;
auto x59 = V{4.5}*(x58*x58);
auto x60 = x53 + x57 - x59;
auto x61 = x37*x60;
auto x62 = x45 + x57;
auto x63 = x49 + x59 + x62;
auto x64 = -x37*x63;
auto x65 = -x57;
auto x66 = -x21;
auto x67 = x19 + x66;
auto x68 = -V{4.5}*x67*x67;
auto x69 = x53 + x65 + x68;
auto x70 = x37*x69;
auto x71 = V{1.11022302462516e-16}*x27*x31 + V{1.11022302462516e-16}*x29*x32;
auto x72 = -x67;
auto x73 = -x46;
auto x74 = x52 + x57;
auto x75 = x73 + x74;
auto x76 = x75 - V{4.5}*x72*x72;
auto x77 = V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32;
auto x78 = V{3}*x36;
auto x79 = x46 + x51 - x78;
auto x80 = x77*x79;
auto x81 = x19 + x20;
auto x82 = V{4.5}*(x81*x81);
auto x83 = x38 + x74 - x82;
auto x84 = x19 - x20;
auto x85 = -V{4.5}*x84*x84;
auto x86 = x38 + x52;
auto x87 = x65 + x85 + x86;
auto x88 = x20 + x66;
auto x89 = -x88;
auto x90 = x73 + x86;
auto x91 = x90 - V{4.5}*x89*x89;
auto x92 = V{0.0555555555555556}*x27*x31 + V{0.0555555555555556}*x29*x32;
auto x93 = V{3}*x35;
auto x94 = x41 + V{-1};
auto x95 = x38 + x47 - x93 + x94;
auto x96 = x44 + x49 + x78;
auto x97 = x77*x96;
auto x98 = V{3}*x34;
auto x99 = x43 + x57 + x94 - x98;
auto x100 = x77*x99;
auto x101 = x62 + x98;
auto x102 = -x101*x77;
auto x103 = V{1.11022302462516e-16}*cell[2];
auto x104 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[8] + V{5.55111512312578e-17};
auto x105 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x100 + x102 + x27*(V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[9] + x103 + x104 + x56) + x29*(V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[6] + x104) + x34*(V{8.32667268468867e-17}*x27*x31 + V{8.32667268468867e-17}*x29*x32) + x55 + V{-2.22044604925031e-16};
auto x106 = x37*x83;
auto x107 = x38 + x48;
auto x108 = x107 + x62 + x82;
auto x109 = -x108*x37;
auto x110 = x37*x87;
auto x111 = -x84;
auto x112 = -x38;
auto x113 = x112 + x74;
auto x114 = x113 - V{4.5}*x111*x111;
auto x115 = x77*x95;
auto x116 = -V{4.5}*x88*x88;
auto x117 = x112 + x116 + x53;
auto x118 = x107 + x42 + x93 + V{1};
auto x119 = x118*x77;
auto x120 = x27*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x105 + x35*x71 + x36*(V{4.85722573273506e-17}*x27*x31 + V{4.85722573273506e-17}*x29*x32) + x37*x76 - x50*(V{3.08395284618099e-18}*x27*x31 + V{3.08395284618099e-18}*x29*x32) + x54*(V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32) + x56 + x61 + x64 + x70 + x77*x83 + x77*x87 + x77*x91 + x80 + x92*x95 - x97) + x29*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x103 + x105 + x106 + x109 + x110 + x114*x37 + x115 + x117*x77 - x119 + x35*(V{5.55111512312578e-17}*x27*x31 + V{5.55111512312578e-17}*x29*x32) + x36*x71 + x54*(V{0.0277777777777778}*x27*x31 + V{0.0277777777777778}*x29*x32) + x60*x77 + x69*x77 + x79*x92);
auto x121 = V{0.0833333333333334}*cell[13];
auto x122 = V{0.0833333333333334}*cell[14];
auto x123 = V{0.0833333333333334}*cell[15];
auto x124 = V{0.0833333333333334}*cell[16];
auto x125 = V{0.0833333333333334}*cell[4];
auto x126 = V{0.0833333333333334}*cell[5];
auto x127 = V{0.0833333333333334}*cell[6];
auto x128 = V{0.0833333333333334}*cell[7];
auto x129 = V{0.0833333333333333}*cell[11];
auto x130 = V{0.0833333333333333}*cell[12];
auto x131 = V{0.0833333333333333}*cell[2];
auto x132 = V{0.0833333333333333}*cell[3];
auto x133 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{3.46944695195361e-18};
auto x134 = x27*(V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x133);
auto x135 = x29*(V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + x133);
auto x136 = V{0.00231481481481481}*x27*x31 + V{0.00231481481481481}*x29*x32;
auto x137 = V{0.0833333333333333}*x27*x31 + V{0.0833333333333333}*x29*x32;
auto x138 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x139 = x138*x35;
auto x140 = x138*x36;
auto x141 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.166666666666667}*cell[9] + x121 + x122 + x123 + x124 + x125 + x126 + x127 + x128 - x129 - x130 - x131 - x132 + x134 + x135 + x136*x50 + x136*x54 + x137*x34 - x139 - x140 + V{-0.0555555555555555};
auto x142 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x143 = V{0.0833333333333333}*cell[10] - V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - V{0.166666666666667}*cell[8] - V{0.0833333333333334}*cell[9] - x134 - x135 + x138*x34 + x142*x50 + x142*x54 + V{0.0555555555555555};
auto x144 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x121 - x122 - x125 - x126 + x130 + x132 - x137*x35 + x140 + x143;
auto x145 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x123 - x124 - x127 - x128 + x129 + x131 - x137*x36 + x139 + x143;
auto x146 = x27*x31 + x29*x32;
auto x147 = V{0.125}*x146*x19;
auto x148 = x147*x20;
auto x149 = V{0.0208333333333333}*x27*x31 + V{0.0208333333333333}*x29*x32;
auto x150 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[10] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0208333333333333};
auto x151 = x27*(V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[14] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[2] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + V{0.0416666666666667}*cell[9] + x150);
auto x152 = x29*(V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[3] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[6] + x150);
auto x153 = V{0.0416666666666667}*x27*x31 + V{0.0416666666666667}*x29*x32;
auto x154 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x151 + x152 + x153*x35 + V{0.0138888888888889};
auto x155 = V{0.000578703703703704}*x27*x31 + V{0.000578703703703704}*x29*x32;
auto x156 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0416666666666667}*cell[9] + x153*x34 - x155*x50 - x155*x54;
auto x157 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x149*x36 + x154 + x156;
auto x158 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x148 + x157;
auto x159 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x148 + x157;
auto x160 = x147*x21;
auto x161 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x153*x36;
auto x162 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x149*x35 + x151 + x152 + x156 + x161 + V{0.0138888888888889};
auto x163 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x160 + x162;
auto x164 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x160 + x162;
auto x165 = V{0.00578703703703704}*x27*x31 + V{0.00578703703703704}*x29*x32;
auto x166 = V{0.125}*x146*x20*x21;
auto x167 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x149*x34 + x154 + x161;
auto x168 = -V{0.0833333333333333}*cell[18] + V{0.833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] + x166 + x167;
auto x169 = V{0.00115740740740741}*x27*x31 + V{0.00115740740740741}*x29*x32;
auto x170 = V{0.416666666666667}*cell[18] - V{0.166666666666667}*cell[8] + V{0.416666666666667}*cell[9] - x166 + x167 + x169*x50 + x169*x54;
auto x0 = -x22*(-V{0.166666666666667}*x120*x52 + V{0.333333333333333}) + x23*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[18] + V{0.5}*cell[1] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[9] + x24 + x25 + x26 + x27*(V{0.25}*cell[12] + V{0.5}*cell[14] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.25}*cell[3] + V{0.5}*cell[4] + V{0.25}*cell[6] + V{0.25}*cell[7] + V{0.5}*cell[9] + x24 + x28) + x29*(V{0.25}*cell[11] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[16] + V{0.5}*cell[18] + V{0.25}*cell[2] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[6] + x25 + x28) + x33*x34 + x33*x35 + x33*x36 - x37*x50 - x37*x54 - x55 + V{0.833333333333333});
auto x1 = -x22*(-V{0.0277777777777778}*x120*x99 + V{0.0555555555555556}) - x23*(x100 + x141);
auto x2 = -x22*(-V{0.0277777777777778}*x120*x95 + V{0.0555555555555556}) + x23*(-x115 + x144);
auto x3 = -x22*(-V{0.0277777777777778}*x120*x79 + V{0.0555555555555556}) + x23*(x145 - x80);
auto x4 = -x22*(-V{0.0138888888888889}*x120*x83 + V{0.0277777777777778}) - x23*(x106 + x158);
auto x5 = -x22*(-V{0.0138888888888889}*x114*x120 + V{0.0277777777777778}) - x23*(x159 + x37*(x113 + x85));
auto x6 = -x22*(-V{0.0138888888888889}*x120*x60 + V{0.0277777777777778}) - x23*(x163 + x61);
auto x7 = -x22*(-V{0.0138888888888889}*x120*x76 + V{0.0277777777777778}) - x23*(x164 + x37*(x68 + x75));
auto x8 = -x22*(-V{0.0138888888888889}*x120*x54 + V{0.0277777777777778}) - x23*(-x165*x50 + x168 + x54*(V{0.00810185185185185}*x27*x31 + V{0.00810185185185185}*x29*x32));
auto x9 = -x22*(-V{0.0138888888888889}*x120*x91 + V{0.0277777777777778}) - x23*(x170 + x37*(x116 + x90));
auto x10 = -x22*(V{0.0277777777777778}*x101*x120 + V{0.0555555555555556}) - x23*(x102 + x141);
auto x11 = -x22*(V{0.0277777777777778}*x118*x120 + V{0.0555555555555556}) + x23*(x119 + x144);
auto x12 = -x22*(V{0.0277777777777778}*x120*x96 + V{0.0555555555555556}) + x23*(x145 + x97);
auto x13 = -x22*(V{0.0138888888888889}*x108*x120 + V{0.0277777777777778}) - x23*(x109 + x158);
auto x14 = -x22*(-V{0.0138888888888889}*x120*x87 + V{0.0277777777777778}) - x23*(x110 + x159);
auto x15 = -x22*(V{0.0138888888888889}*x120*x63 + V{0.0277777777777778}) - x23*(x163 + x64);
auto x16 = -x22*(-V{0.0138888888888889}*x120*x69 + V{0.0277777777777778}) - x23*(x164 + x70);
auto x17 = -x22*(V{0.0138888888888889}*x120*x50 + V{0.0277777777777778}) - x23*(-x165*x54 + x168 - x50*(V{0.0196759259259259}*x27*x31 + V{0.0196759259259259}*x29*x32));
auto x18 = -x22*(-V{0.0138888888888889}*x117*x120 + V{0.0277777777777778}) - x23*(x117*x37 + x170);
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
return { -V{0.5}*x120, x34 + x35 + x36 };
}
};

}

}
