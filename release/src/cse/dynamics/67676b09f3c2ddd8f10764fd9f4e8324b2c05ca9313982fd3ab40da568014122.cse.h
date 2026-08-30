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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x25 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(2);
auto x24 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x22 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x26 = parameters.template get<descriptors::OMEGA>();
auto x19 = x26 + V{-1};
auto x20 = V{0.0277777777777778}*x22;
auto x21 = V{1} / (x22);
auto x23 = V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16];
auto x27 = cell[0] + V{2}*cell[10] + cell[11] + cell[12] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + x23 + V{1};
auto x28 = -x21*x27 + V{1};
auto x29 = -x28;
auto x30 = x29*x29;
auto x31 = V{1.5}*x30;
auto x32 = x31 + V{-1};
auto x33 = x24*x24;
auto x34 = V{1.5}*x33;
auto x35 = V{3}*x25;
auto x36 = x25*x25;
auto x37 = V{3}*x36;
auto x38 = x34 + x35 - x37;
auto x39 = x32 + x38;
auto x40 = V{1.5}*x36;
auto x41 = V{3}*x24;
auto x42 = V{3}*x33;
auto x43 = x40 + x41 - x42;
auto x44 = x32 + x43;
auto x45 = x34 + x40;
auto x46 = x45 + V{-1};
auto x47 = x31 + x46;
auto x48 = x41 + x47;
auto x49 = x24 + x25;
auto x50 = V{4.5}*(x49*x49);
auto x51 = x35 - x50;
auto x52 = x48 + x51;
auto x53 = -x35;
auto x54 = x24 - x25;
auto x55 = -x54;
auto x56 = x53 - V{4.5}*x55*x55;
auto x57 = x48 + x56;
auto x58 = -x41;
auto x59 = -V{4.5}*x54*x54;
auto x60 = x35 + x58 + x59;
auto x61 = x47 + x60;
auto x62 = V{0.0555555555555556}*x22;
auto x63 = x21*(V{3}*cell[0] + V{6}*cell[10] + V{3}*cell[11] + V{3}*cell[12] + V{6}*cell[13] + V{6}*cell[14] + V{6}*cell[15] + V{6}*cell[16] + V{3}*cell[17] + V{3}*cell[18] + V{3}*cell[2] + V{3}*cell[3] + V{3}*cell[8] + V{3}*cell[9] + V{3});
auto x64 = x45 + V{2};
auto x65 = -x63 + x64;
auto x66 = -V{3}*x30 + x65;
auto x67 = x24 + x28;
auto x68 = -x67;
auto x69 = x31 + x65;
auto x70 = x41 + x69 - V{4.5}*x68*x68;
auto x71 = x25 + x28;
auto x72 = -x71;
auto x73 = x35 + x69 - V{4.5}*x72*x72;
auto x74 = x21*x27 + V{-1};
auto x75 = x24 + x74;
auto x76 = -V{4.5}*x75*x75;
auto x77 = x58 + x76;
auto x78 = x69 + x77;
auto x79 = x25 + x74;
auto x80 = -V{4.5}*x79*x79;
auto x81 = x53 + x80;
auto x82 = x69 + x81;
auto x83 = V{0.333333333333333}*x22;
auto x84 = x47*x83;
auto x85 = x28*x28;
auto x86 = V{1.5}*x85;
auto x87 = V{1} - x86;
auto x88 = -x40 + x41;
auto x89 = x42 + x88;
auto x90 = x87 + x89;
auto x91 = -x34 + x35;
auto x92 = x37 + x91;
auto x93 = x87 + x92;
auto x94 = x50 + x88 + x91;
auto x95 = x87 + x94;
auto x96 = V{1}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[8] + V{1}*cell[9] + x23 + V{0.833333333333333};
auto x97 = V{0.0277777777777778}*x22;
auto x98 = -x57;
auto x99 = V{0.0555555555555555}*x22;
auto x100 = -x39;
auto x101 = -x82;
auto x102 = -x73;
auto x103 = V{1} - x31;
auto x104 = x103 + x89;
auto x105 = x103 + x92;
auto x106 = V{0.0555555555555555}*x22;
auto x107 = -x44;
auto x108 = V{0.111111111111111}*x22;
auto x109 = -x66;
auto x110 = -x47;
auto x111 = x110*x83;
auto x112 = V{0.0277777777777778}*x22;
auto x113 = -x52;
auto x114 = V{0.0555555555555555}*x22;
auto x115 = -x78;
auto x116 = -x70;
auto x117 = V{0.0277777777777778}*x22;
auto x118 = x103 + x94;
auto x119 = -x61;
auto x120 = V{4.44089209850063e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{6.66133814775094e-16}*cell[14] + V{8.88178419700125e-16}*cell[15] + V{4.44089209850063e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{2.22044604925031e-16};
auto x121 = x21*(x100*x99 + x101*x99 + x102*x99 + x104*x62 + x105*x62 + x106*x107 + x108*x109 + x111 + x112*x113 + x114*x115 + x114*x116 + x117*x118 + x117*x119 + x120 + x97*x98);
auto x122 = x121 + V{-1};
auto x123 = x122*x122;
auto x124 = V{1.5}*x123;
auto x125 = x124 + x46;
auto x126 = V{0.333333333333333}*x26*(x125*x22 + V{1});
auto x127 = x45 + V{-4};
auto x128 = x127 + x63;
auto x129 = V{0.0833333333333333}*cell[11];
auto x130 = V{0.0833333333333333}*cell[12];
auto x131 = V{0.0833333333333333}*cell[2];
auto x132 = V{0.0833333333333333}*cell[3];
auto x133 = V{0.166666666666667}*cell[13];
auto x134 = V{0.166666666666667}*cell[14];
auto x135 = V{0.166666666666667}*cell[15];
auto x136 = V{0.166666666666667}*cell[16];
auto x137 = V{0.00462962962962963}*x22;
auto x138 = x137*x39;
auto x139 = x137*x44;
auto x140 = x137*x90;
auto x141 = x137*x93;
auto x142 = V{0.00462962962962963}*x22;
auto x143 = x142*x70;
auto x144 = x142*x73;
auto x145 = x142*x78;
auto x146 = x142*x82;
auto x147 = -V{0.333333333333333}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x129 + x130 + x131 + x132 - x133 - x134 - x135 - x136 + x137*x52 + x137*x57 + x137*x61 - x137*x95 + x138 + x139 - x140 - x141 - x143 - x144 - x145 - x146 + V{0.0555555555555555};
auto x148 = V{0.0555555555555556}*x26;
auto x149 = V{1} - x121;
auto x150 = x149*x149;
auto x151 = V{0.0833333333333333}*x22;
auto x152 = V{0.166666666666667}*x22;
auto x153 = V{0.166666666666667}*x22;
auto x154 = V{0.166666666666667}*x22;
auto x155 = x21*(V{1.33226762955019e-15}*cell[10] + V{4.9960036108132e-16}*cell[11] + V{3.33066907387547e-16}*cell[12] + V{1.99840144432528e-15}*cell[13] + V{1.99840144432528e-15}*cell[14] + V{2.66453525910038e-15}*cell[15] + V{1.33226762955019e-15}*cell[16] + V{1.41553435639707e-15}*cell[17] + V{1.58206781009085e-15}*cell[18] + V{4.9960036108132e-16}*cell[2] + V{3.33066907387547e-16}*cell[3] + V{1.41553435639707e-15}*cell[8] + V{1.58206781009085e-15}*cell[9] + x100*x152 + x101*x152 + x102*x152 + x104*x154 + x105*x154 + V{0.166666666666667}*x107*x22 + V{0.333333333333333}*x109*x22 + V{1}*x110*x22 + V{0.0833333333333333}*x113*x22 + x115*x153 + x116*x153 + x118*x151 + x119*x151 + V{0.0833333333333333}*x22*x98 + V{6.66133814775094e-16});
auto x156 = x124 + x127 + x155;
auto x157 = V{0.0462962962962963}*x22;
auto x158 = x86 + V{-1};
auto x159 = x158 + x43;
auto x160 = V{0.00925925925925926}*x22;
auto x161 = V{0.166666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2] - x130 - x132 + x133 + x134 + x141;
auto x162 = -x160*x90 + x161;
auto x163 = -V{4.5}*x67*x67;
auto x164 = x65 + x86;
auto x165 = x163 + x164 + x41;
auto x166 = x164 + x77;
auto x167 = x158 + x38;
auto x168 = V{0.00925925925925926}*x22;
auto x169 = -V{4.5}*x71*x71;
auto x170 = x164 + x169 + x35;
auto x171 = x164 + x81;
auto x172 = V{0.00231481481481482}*x22;
auto x173 = x46 + x86;
auto x174 = x41 + x51;
auto x175 = x173 + x174;
auto x176 = x173 + x41 + x53 + x59;
auto x177 = x173 + x60;
auto x178 = x65 - V{3}*x85;
auto x179 = -V{0.166666666666667}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x172*x95 + V{-0.0555555555555555};
auto x180 = -x168*x178 + x172*x175 + x172*x176 + x172*x177 + x179;
auto x181 = -x137*x167 + x142*x165 + x142*x166 - x168*x170 - x168*x171 + x180;
auto x182 = x124 + V{-1};
auto x183 = x148*(x22*(x182 + x43) + V{1});
auto x184 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[3] - x129 - x131 + x135 + x136 + x140;
auto x185 = -x160*x93 + x184;
auto x186 = -x137*x159 + x142*x170 + x142*x171 - x165*x168 - x166*x168 + x180;
auto x187 = x148*(x22*(x182 + x38) + V{1});
auto x188 = V{0.0231481481481481}*x22;
auto x189 = x128 + x86;
auto x190 = V{0.00462962962962963}*x22;
auto x191 = V{0.00231481481481481}*x22;
auto x192 = V{0.00115740740740741}*x22;
auto x193 = V{0.166666666666667}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] - x192*x95;
auto x194 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x191*x93 + x193;
auto x195 = V{0.833333333333333}*cell[13] - V{0.166666666666667}*cell[14] + x194;
auto x196 = V{0.00462962962962963}*x22;
auto x197 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] - x196*x90 + V{0.0138888888888889};
auto x198 = x159*x196 + x170*x172 + x171*x172 + x197;
auto x199 = x160*x178 + x175*x192 + x176*x192 + x177*x192;
auto x200 = -x167*x191 + x198 + x199;
auto x201 = -x165*x190 + x195 + x200;
auto x202 = V{0.0277777777777778}*x26;
auto x203 = x122 + x24;
auto x204 = -x203;
auto x205 = -V{0.166666666666667}*cell[13] + V{0.833333333333333}*cell[14] + x194;
auto x206 = -x166*x190 + x200 + x205;
auto x207 = x149 + x24;
auto x208 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x191*x90 + x193 + V{0.0138888888888889};
auto x209 = V{0.833333333333333}*cell[15] - V{0.166666666666667}*cell[16] + x208;
auto x210 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - x196*x93;
auto x211 = x165*x172 + x166*x172 + x167*x196 + x210;
auto x212 = -x159*x191 + x199 + x211;
auto x213 = -x170*x190 + x209 + x212;
auto x214 = x122 + x25;
auto x215 = -x214;
auto x216 = -V{0.166666666666667}*cell[15] + V{0.833333333333333}*cell[16] + x208;
auto x217 = -x171*x190 + x212 + x216;
auto x218 = x149 + x25;
auto x219 = V{0.0162037037037037}*x22;
auto x220 = V{0.0115740740740741}*x22;
auto x221 = -V{0.0833333333333333}*cell[10];
auto x222 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x221;
auto x223 = -x220*x95 + x222;
auto x224 = V{0.00231481481481481}*x22;
auto x225 = -x137*x178 + x198 + x211;
auto x226 = -x176*x224 - x177*x224 + x225;
auto x227 = x202*(x22*(x125 + x174) + V{1});
auto x228 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x221 + x224*x95;
auto x229 = -x175*x224 + x225 + x228;
auto x230 = x202*(x22*(x125 + x41 + x56) + V{1});
auto x231 = x147 + V{0.037037037037037}*x22*x66;
auto x232 = -x155 + x64;
auto x233 = x22*(-V{3}*x123 + x232) + V{1};
auto x234 = x157*x90 + x161;
auto x235 = -x21*(-x106*x44 - x108*x66 - x112*x52 - x114*x70 - x114*x78 - x117*x61 + x117*x95 + x120 - x39*x99 - x57*x97 + x62*x90 + x62*x93 - x73*x99 - x82*x99 - x84) + V{1};
auto x236 = x235*x235;
auto x237 = V{1} - V{1.5}*x236;
auto x238 = x157*x93 + x184;
auto x239 = x124 + x232;
auto x240 = x22*(x239 + x58 - V{4.5}*x203*x203) + V{1};
auto x241 = -x207;
auto x242 = x22*(x239 + x41 - V{4.5}*x241*x241) + V{1};
auto x243 = x22*(x239 + x53 - V{4.5}*x214*x214) + V{1};
auto x244 = -x218;
auto x245 = x22*(x239 + x35 - V{4.5}*x244*x244) + V{1};
auto x246 = x219*x95 + x222;
auto x247 = x202*(x22*(x125 + x60) + V{1});
auto x248 = x172*x73 + x172*x82 + x196*x44 + x197;
auto x249 = x172*x70 + x172*x78 + x196*x39 + x210;
auto x250 = -x137*x66 + x248 + x249;
auto x251 = -x224*x57 - x224*x61 + x250;
auto x252 = -x224*x52 + x228 + x250;
auto x253 = -x168*x66 + x172*x52 + x172*x57 + x172*x61 + x179;
auto x254 = -x138 + x143 + x145 - x168*x73 - x168*x82 + x253;
auto x255 = -x139 + x144 + x146 - x168*x70 - x168*x78 + x253;
auto x256 = V{2}*x26 + V{-2};
auto x257 = x160*x66 + x192*x52 + x192*x57 + x192*x61;
auto x258 = -x191*x39 + x248 + x257;
auto x259 = -x191*x44 + x249 + x257;
auto x260 = V{1} - V{1.5}*x150;
cell[0] = -x126 + x19*(x20*x39 + x20*x44 + x20*x52 + x20*x57 + x20*x61 - x20*x90 - x20*x93 - x20*x95 + x62*x66 + x62*x70 + x62*x73 + x62*x78 + x62*x82 + x84 + x96);
cell[1] = -x148*(x22*(-V{4.5}*x150 + x156) + V{1}) + x19*(x147 - V{0.0185185185185185}*x22*x66 + x62*(x128 + x31 - V{4.5}*x85));
cell[2] = -x183 - x19*(-x157*x159 + x162 + x181);
cell[3] = -x187 - x19*(-x157*x167 + x185 + x186);
cell[4] = -(x19*(x166*x188 - x20*(x189 + x41 + x76) + x201) + x202*(x22*(x156 + x41 - V{4.5}*x204*x204) + V{1}));
cell[5] = -(x19*(x165*x188 - x20*(x163 + x189 + x58) + x206) + x202*(x22*(x156 + x58 - V{4.5}*x207*x207) + V{1}));
cell[6] = -(x19*(x171*x188 - x20*(x189 + x35 + x80) + x213) + x202*(x22*(x156 + x35 - V{4.5}*x215*x215) + V{1}));
cell[7] = -(x19*(x170*x188 - x20*(x169 + x189 + x53) + x217) + x202*(x22*(x156 + x53 - V{4.5}*x218*x218) + V{1}));
cell[8] = -x19*(-x175*x219 + x223 + x226) - x227;
cell[9] = -x19*(-x176*x219 + x177*x220 + x229) - x230;
cell[10] = -x148*x233 + x19*x231;
cell[11] = x148*(x22*(x237 + x89) + V{-1}) - x19*(x159*x160 + x181 + x234);
cell[12] = x148*(x22*(x237 + x92) + V{-1}) - x19*(x160*x167 + x186 + x238);
cell[13] = -x19*(-x137*x166 + x201) - x202*x240;
cell[14] = -x19*(-x137*x165 + x206) - x202*x242;
cell[15] = -x19*(-x137*x171 + x213) - x202*x243;
cell[16] = -x19*(-x137*x170 + x217) - x202*x245;
cell[17] = -x19*(x175*x220 + x226 + x246) + x202*(x22*(x237 + x94) + V{-1});
cell[18] = -x19*(x176*x220 - x177*x219 + x229) - x247;
cell.template getFieldPointer<descriptors::VELOCITY>()[0] = V{1}*x21*(-x126 - x148*x240 - x148*x242 - x148*x243 - x148*x245 - x183 - x187 - x19*(-x157*x39 + x185 + x255) - x19*(-x157*x44 + x162 + x254) - x19*(x160*x39 + x238 + x255) - x19*(x160*x44 + x234 + x254) - x19*(-x219*x52 + x223 + x251) - x19*(-x219*x57 + x220*x61 + x252) - x19*(-x219*x61 + x220*x57 + x252) - x19*(x220*x52 + x246 + x251) + x19*(-x100*x20 - x101*x62 - x102*x62 - x104*x20 - x105*x20 - x107*x20 - x109*x62 - x111 - x113*x20 - x115*x62 - x116*x62 - x118*x20 - x119*x20 - x20*x98 + x96) - x227 - x230 + x231*x256 - V{0.111111111111111}*x233*x26 - x247 - x256*(-x137*x70 - x190*x78 + x205 + x258) - x256*(-x137*x73 - x190*x82 + x216 + x259) - x256*(-x137*x78 - x190*x70 + x195 + x258) - x256*(-x137*x82 - x190*x73 + x209 + x259) + V{0.0555555555555556}*x26*(x22*(x260 + x89) + V{-1}) + V{0.0555555555555556}*x26*(x22*(x260 + x92) + V{-1}) + V{0.0277777777777778}*x26*(x22*(x260 + x94) + V{-1}) + V{1}) + V{-1};
cell.template getFieldPointer<descriptors::VELOCITY>()[1] = x24;
cell.template getFieldPointer<descriptors::VELOCITY>()[2] = x25;
return { x22, V{1}*x236 + x33 + x36 };
}
};

}

}
