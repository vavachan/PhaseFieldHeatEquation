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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x30 = parameters.template get<descriptors::OMEGA>();
auto x29 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x28 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x31 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x27 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x32 = cell[13] + cell[1] + cell[5] + cell[7];
auto x33 = cell[18] + cell[25] + cell[2] + cell[9];
auto x34 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x35 = cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + x32 + x33 + x34;
auto x36 = x35 + V{1};
auto x37 = V{1} / (x36);
auto x38 = -cell[23];
auto x39 = cell[10] + cell[11] - cell[17] - cell[24] + cell[4] + x38;
auto x40 = cell[12] - cell[19] - cell[25] + cell[6];
auto x41 = -cell[14] - cell[18] - cell[20] - cell[26] + x32 + x39 + x40;
auto x42 = x37*(x35 + V{1});
auto x43 = V{0.666666666666667}*cell[10];
auto x44 = V{0.666666666666667}*cell[11];
auto x45 = V{0.666666666666667}*cell[12];
auto x46 = V{0.666666666666667}*cell[13];
auto x47 = V{0.666666666666667}*cell[23];
auto x48 = V{0.666666666666667}*cell[24];
auto x49 = V{0.666666666666667}*cell[25];
auto x50 = V{0.666666666666667}*cell[26];
auto x51 = -V{0.333333333333333}*cell[0];
auto x52 = -V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + x43 + x44 + x45 + x46 + x47 + x48 + x49 + x50 + x51;
auto x53 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x54 = V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9] - x27*x41*x42 - x37*x41*x41 + x52 + x53;
auto x55 = -cell[13] - cell[21] + cell[26] + cell[8];
auto x56 = -cell[12] - cell[15] - cell[22] - cell[5] + x33 + x39 + x55;
auto x57 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x58 = V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7] - x28*x42*x56 - x37*x56*x56 + x52 + x57;
auto x59 = -cell[11] - cell[16] - cell[7] - cell[9] + x34 + x38 + x40 + x55;
auto x60 = V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5] - x29*x42*x59 - x37*x59*x59 + x43 + x44 + x45 + x46 + x47 + x48 + x49 + x50 + x51 + x53 + x57;
auto x61 = V{1}*cell[18];
auto x62 = V{1}*cell[11];
auto x63 = -x62;
auto x64 = V{1}*cell[4];
auto x65 = x42*(x27*x56 + x28*x41);
auto x66 = x37*x56;
auto x67 = x41*x66;
auto x68 = V{1}*cell[12];
auto x69 = V{1}*cell[25];
auto x70 = V{1}*cell[10];
auto x71 = -x70;
auto x72 = -V{1}*cell[23];
auto x73 = x68 + x69 + x71 + x72;
auto x74 = V{1}*cell[13];
auto x75 = V{1}*cell[26];
auto x76 = x74 + x75;
auto x77 = V{1}*cell[5];
auto x78 = -V{1}*cell[17];
auto x79 = V{1}*cell[24];
auto x80 = -x79;
auto x81 = x77 + x78 + x80;
auto x82 = V{2}*cell[13];
auto x83 = V{2}*cell[26];
auto x84 = V{2}*cell[11];
auto x85 = V{2}*cell[24];
auto x86 = -V{2}*cell[10] + V{2}*cell[12] - V{2}*cell[23] + V{2}*cell[25];
auto x87 = V{1}*cell[22];
auto x88 = -x75;
auto x89 = V{1}*cell[8];
auto x90 = x42*(x28*x59 + x29*x56);
auto x91 = x59*x66;
auto x92 = x62 + x79;
auto x93 = V{1}*cell[9];
auto x94 = -x74;
auto x95 = -V{1}*cell[21];
auto x96 = x93 + x94 + x95;
auto x97 = V{1}*cell[20];
auto x98 = V{1}*cell[6];
auto x99 = -x68 + x72;
auto x100 = V{1}*cell[7];
auto x101 = -V{1}*cell[19];
auto x102 = -x69;
auto x103 = x100 + x101 + x102;
auto x104 = -x103 - V{1}*x37*x41*x59 - V{0.5}*x42*(x27*x59 + x29*x41) - x71 - x76 - x92 - x97 + x98 - x99;
auto x105 = V{1} / (V{3.00000046417339}*util::sqrt(x37*(x31*x31)*util::sqrt(V{0.5}*(-x61 - x63 + x64 - V{0.5}*x65 - V{1}*x67 - x73 - x76 - x81)*(V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[4] - V{2}*cell[5] - V{1}*x65 - V{2}*x67 - x82 - x83 + x84 + x85 - x86) + V{0.5}*(-x73 - x87 - x88 + x89 - V{0.5}*x90 - V{1}*x91 - x92 - x96)*(V{2}*cell[21] - V{2}*cell[22] + V{2}*cell[8] - V{2}*cell[9] + x82 + x83 - x84 - x85 - x86 - V{1}*x90 - V{2}*x91) + x104*x104 + V{0.5}*(x54*x54) + V{0.5}*(x58*x58) + V{0.5}*(x60*x60)) + V{0.0277777691819762}/((x30)*(x30))) + V{0.5}/x30);
auto x106 = V{0.5}*x27;
auto x107 = x62 + x64 + x70;
auto x108 = x68 + x72 + x98;
auto x109 = -V{1}*cell[14] + V{1}*cell[1] + x103 + x107 + x108 - x61 + x74 + x81 + x88 - x97;
auto x110 = -x109*x37;
auto x111 = x106 + x110;
auto x112 = V{0.5}*x28;
auto x113 = x75 + x89;
auto x114 = -V{1}*cell[15] + V{1}*cell[2] + x107 + x113 + x61 + x69 - x77 + x78 + x80 - x87 + x96 + x99;
auto x115 = -x114*x37;
auto x116 = x112 + x115;
auto x117 = V{0.5}*x29;
auto x118 = -V{1}*cell[16] + V{1}*cell[3] - x100 + x101 + x102 + x108 + x113 + x63 + x70 + x79 + x87 - x93 + x94 + x95 + x97;
auto x119 = -x118*x37;
auto x120 = x117 + x119;
auto x121 = V{-1} + V{1.5}*(x111*x111) + V{1.5}*(x116*x116) + V{1.5}*(x120*x120);
auto x122 = V{1} - x105;
auto x123 = V{1.5}*x27;
auto x124 = V{3}*cell[13];
auto x125 = V{3}*cell[5];
auto x126 = V{3}*cell[7];
auto x127 = V{3}*cell[18];
auto x128 = V{3}*cell[20];
auto x129 = V{3}*cell[26];
auto x130 = V{3}*cell[10];
auto x131 = V{3}*cell[11];
auto x132 = -V{3}*cell[23];
auto x133 = V{3}*cell[24];
auto x134 = -V{3}*cell[17] + V{3}*cell[4] + x130 + x131 + x132 - x133;
auto x135 = V{3}*cell[12];
auto x136 = V{3}*cell[25];
auto x137 = -V{3}*cell[19] + V{3}*cell[6] + x135 - x136;
auto x138 = -V{3}*cell[14] + V{3}*cell[1] + x124 + x125 + x126 - x127 - x128 - x129 + x134 + x137;
auto x139 = -x138*x37;
auto x140 = x123 + x139;
auto x141 = x140*x27;
auto x142 = V{1.5}*x28;
auto x143 = V{3}*cell[9];
auto x144 = V{3}*cell[22];
auto x145 = -V{3}*cell[21] + V{3}*cell[8] - x124 + x129;
auto x146 = -V{3}*cell[15] + V{3}*cell[2] - x125 + x127 + x134 - x135 + x136 + x143 - x144 + x145;
auto x147 = -x146*x37;
auto x148 = x142 + x147;
auto x149 = x148*x28;
auto x150 = V{1.5}*x29;
auto x151 = -V{3}*cell[16] + V{3}*cell[3] - x126 + x128 + x130 - x131 + x132 + x133 + x137 - x143 + x144 + x145;
auto x152 = -x151*x37;
auto x153 = x150 + x152;
auto x154 = x153*x29;
auto x155 = x149 + x154;
auto x156 = V{1} - V{0.5}*x105;
auto x157 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x158 = V{2.25}*x27;
auto x159 = V{4.5}*cell[13];
auto x160 = V{4.5}*cell[5];
auto x161 = V{4.5}*cell[7];
auto x162 = V{4.5}*cell[18];
auto x163 = V{4.5}*cell[20];
auto x164 = V{4.5}*cell[26];
auto x165 = V{4.5}*cell[10];
auto x166 = V{4.5}*cell[11];
auto x167 = -V{4.5}*cell[23];
auto x168 = V{4.5}*cell[24];
auto x169 = -V{4.5}*cell[17] + V{4.5}*cell[4] + x165 + x166 + x167 - x168;
auto x170 = V{4.5}*cell[12];
auto x171 = V{4.5}*cell[25];
auto x172 = -V{4.5}*cell[19] + V{4.5}*cell[6] + x170 - x171;
auto x173 = -V{4.5}*cell[14] + V{4.5}*cell[1] + x159 + x160 + x161 - x162 - x163 - x164 + x169 + x172;
auto x174 = -x173*x37;
auto x175 = x158 + x174;
auto x176 = x121 + x140;
auto x177 = V{6}*cell[13];
auto x178 = V{6}*cell[5];
auto x179 = V{6}*cell[7];
auto x180 = V{6}*cell[18];
auto x181 = V{6}*cell[20];
auto x182 = V{6}*cell[26];
auto x183 = V{6}*cell[10];
auto x184 = V{6}*cell[11];
auto x185 = -V{6}*cell[23];
auto x186 = V{6}*cell[24];
auto x187 = -V{6}*cell[17] + V{6}*cell[4] + x183 + x184 + x185 - x186;
auto x188 = V{6}*cell[12];
auto x189 = V{6}*cell[25];
auto x190 = -V{6}*cell[19] + V{6}*cell[6] + x188 - x189;
auto x191 = -V{6}*cell[14] + V{6}*cell[1] + x177 + x178 + x179 - x180 - x181 - x182 + x187 + x190;
auto x192 = x191*x37;
auto x193 = V{3}*x27;
auto x194 = V{3} - x193;
auto x195 = x192 + x194;
auto x196 = V{0.0740740740740741}*x36;
auto x197 = V{2.25}*x28;
auto x198 = V{4.5}*cell[9];
auto x199 = V{4.5}*cell[22];
auto x200 = -V{4.5}*cell[21] + V{4.5}*cell[8] - x159 + x164;
auto x201 = -V{4.5}*cell[15] + V{4.5}*cell[2] - x160 + x162 + x169 - x170 + x171 + x198 - x199 + x200;
auto x202 = -x201*x37;
auto x203 = x197 + x202;
auto x204 = x121 + x148;
auto x205 = V{6}*cell[9];
auto x206 = V{6}*cell[22];
auto x207 = -V{6}*cell[21] + V{6}*cell[8] - x177 + x182;
auto x208 = -V{6}*cell[15] + V{6}*cell[2] - x178 + x180 + x187 - x188 + x189 + x205 - x206 + x207;
auto x209 = x208*x37;
auto x210 = V{3}*x28;
auto x211 = V{3} - x210;
auto x212 = x209 + x211;
auto x213 = x141 + x154;
auto x214 = V{2.25}*x29;
auto x215 = -V{4.5}*cell[16] + V{4.5}*cell[3] - x161 + x163 + x165 - x166 + x167 + x168 + x172 - x198 + x199 + x200;
auto x216 = -x215*x37;
auto x217 = x214 + x216;
auto x218 = x121 + x153;
auto x219 = -V{6}*cell[16] + V{6}*cell[3] - x179 + x181 + x183 - x184 + x185 + x186 + x190 - x205 + x206 + x207;
auto x220 = x219*x37;
auto x221 = V{3}*x29;
auto x222 = V{3} - x221;
auto x223 = x220 + x222;
auto x224 = x141 + x149;
auto x225 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x226 = x111 + x116;
auto x227 = x175 + x203;
auto x228 = x148 + x176;
auto x229 = V{4.5}*x28;
auto x230 = V{9}*cell[18];
auto x231 = V{9}*cell[25];
auto x232 = V{9}*cell[9];
auto x233 = V{9}*cell[12];
auto x234 = V{9}*cell[22];
auto x235 = V{9}*cell[5];
auto x236 = V{9}*cell[10];
auto x237 = V{9}*cell[11];
auto x238 = -V{9}*cell[23];
auto x239 = V{9}*cell[24];
auto x240 = -V{9}*cell[17] + V{9}*cell[4] + x236 + x237 + x238 - x239;
auto x241 = V{9}*cell[26];
auto x242 = V{9}*cell[13];
auto x243 = -V{9}*cell[21] + V{9}*cell[8] + x241 - x242;
auto x244 = -V{9}*cell[15] + V{9}*cell[2] + x230 + x231 + x232 - x233 - x234 - x235 + x240 + x243;
auto x245 = x244*x37;
auto x246 = -x229 + x245;
auto x247 = x195 + x246;
auto x248 = V{4.5}*x27;
auto x249 = V{9}*cell[7];
auto x250 = V{9}*cell[20];
auto x251 = -V{9}*cell[19] + V{9}*cell[6] - x231 + x233;
auto x252 = -V{9}*cell[14] + V{9}*cell[1] - x230 + x235 + x240 - x241 + x242 + x249 - x250 + x251;
auto x253 = x252*x37;
auto x254 = -x248 + x253;
auto x255 = x212 + x254;
auto x256 = x151*x37;
auto x257 = x150 - x256;
auto x258 = x257*x29;
auto x259 = V{0.0185185185185185}*x36;
auto x260 = -x112;
auto x261 = x111 - x115 + x260;
auto x262 = -x197;
auto x263 = x175 - x202 + x262;
auto x264 = -x142;
auto x265 = -x147 + x264;
auto x266 = x176 + x265;
auto x267 = x191*x37 + x194;
auto x268 = x229 - x244*x37;
auto x269 = -x209 + x210 + V{3};
auto x270 = x254 + x269;
auto x271 = V{4.5}*x29;
auto x272 = -V{9}*cell[16] + V{9}*cell[3] - x232 + x234 + x236 - x237 + x238 + x239 + x243 - x249 + x250 + x251;
auto x273 = x272*x37;
auto x274 = -x271 + x273;
auto x275 = x195 + x274;
auto x276 = x146*x37;
auto x277 = x142 - x276;
auto x278 = x277*x28;
auto x279 = x223 + x254;
auto x280 = -x117;
auto x281 = -x119 + x280;
auto x282 = x111 + x281;
auto x283 = -x214;
auto x284 = -x216 + x283;
auto x285 = x175 + x284;
auto x286 = -x150;
auto x287 = -x152 + x286;
auto x288 = x271 - x272*x37;
auto x289 = -x220 + x221 + V{3};
auto x290 = x254 + x289;
auto x291 = x116 + x120;
auto x292 = x203 + x217;
auto x293 = x153 + x204;
auto x294 = x138*x37;
auto x295 = x123 - x294;
auto x296 = x27*x295;
auto x297 = x212 + x274;
auto x298 = x223 + x246;
auto x299 = x116 + x281;
auto x300 = x203 + x284;
auto x301 = x208*x37 + x211;
auto x302 = x246 + x289;
auto x303 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x304 = V{0.00462962962962963}*x36;
auto x305 = x271 - x273;
auto x306 = x229 - x245;
auto x307 = x106 - x109*x37;
auto x308 = x118*x37;
auto x309 = x280 + x308;
auto x310 = x114*x37;
auto x311 = x260 + x310;
auto x312 = x158 - x173*x37;
auto x313 = x215*x37;
auto x314 = x283 + x313;
auto x315 = x201*x37;
auto x316 = x262 + x315;
auto x317 = x307*x307;
auto x318 = x112 - x310;
auto x319 = x318*x318;
auto x320 = x117 - x308;
auto x321 = x320*x320;
auto x322 = V{1.5}*x317 + V{1.5}*x319 + V{1.5}*x321 + V{-1};
auto x323 = x264 + x276;
auto x324 = x322 + x323;
auto x325 = x256 + x286;
auto x326 = x324 + x325;
auto x327 = -x123;
auto x328 = x294 + x322 + x327;
auto x329 = -x192 + x193 + V{3};
auto x330 = x197 - x315;
auto x331 = x214 - x313;
auto x332 = x307 + x318;
auto x333 = x312 + x330;
auto x334 = x323 + x328;
auto x335 = x306 + x329;
auto x336 = x248 - x253;
auto x337 = x269 + x336;
auto x338 = -x139 + x327;
auto x339 = x248 - x252*x37;
auto x340 = x246 + x329;
auto x341 = x307 + x320;
auto x342 = x312 + x331;
auto x343 = x325 + x328;
auto x344 = x289 + x336;
auto x345 = x219*x37 + x222;
auto x0 = V{1}*cell[0]*x122 - x105*(x121*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + V{0.296296296296296}) - V{0.296296296296296}*x156*x36*(x141 + x155);
auto x1 = V{1}*cell[1]*x122 - x105*(x157*(-x111*x175 + x176) + V{0.0740740740740741}) - x156*x196*(x155 + x195*x27);
auto x2 = V{1}*cell[2]*x122 - x105*(x157*(-x116*x203 + x204) + V{0.0740740740740741}) - x156*x196*(x212*x28 + x213);
auto x3 = V{1}*cell[3]*x122 - x105*(x157*(-x120*x217 + x218) + V{0.0740740740740741}) - x156*x196*(x223*x29 + x224);
auto x4 = V{1}*cell[4]*x122 - x105*(x225*(-x226*x227 + x228) + V{0.0185185185185185}) - x156*x259*(x247*x27 + x255*x28 + x258);
auto x5 = V{1}*cell[5]*x122 - x105*(x225*(-x261*x263 + x266) + V{0.0185185185185185}) - x156*x259*(x154 + x27*(x267 + x268) - x270*x28);
auto x6 = V{1}*cell[6]*x122 - x105*(x225*(x153 + x176 - (x111 + x120)*(x175 + x217)) + V{0.0185185185185185}) - x156*x259*(x27*x275 + x278 + x279*x29);
auto x7 = V{1}*cell[7]*x122 - x105*(x225*(x176 - x282*x285 + x287) + V{0.0185185185185185}) - x156*x259*(x149 + x27*(x267 + x288) - x29*x290);
auto x8 = V{1}*cell[8]*x122 - x105*(x225*(-x291*x292 + x293) + V{0.0185185185185185}) - x156*x259*(x28*x297 + x29*x298 + x296);
auto x9 = V{1}*cell[9]*x122 - x105*(x225*(x204 + x287 - x299*x300) + V{0.0185185185185185}) - x156*x259*(x141 + x28*(x288 + x301) - x29*x302);
auto x10 = V{1}*cell[10]*x122 - x105*(x303*(x153 + x228 - (x120 + x226)*(x217 + x227)) + V{0.00462962962962963}) - x156*x304*(x27*(x247 + x274) + x28*(x255 + x274) + x29*(x246 + x279));
auto x11 = V{1}*cell[11]*x122 - x105*(x303*(x228 + x287 - (-x226 - x281)*(-x227 - x284)) + V{0.00462962962962963}) - x156*x304*(x27*(x247 + x305) + x28*(x255 + x305) + x29*(-x246 - x290));
auto x12 = V{1}*cell[12]*x122 - x105*(x303*(x153 + x266 - (-x120 - x261)*(-x217 - x263)) + V{0.00462962962962963}) - x156*x304*(x27*(x275 + x306) + x28*(-x270 - x274) + x29*(x279 + x306));
auto x13 = -x105*(x303*(x295 + x326 - (-x307 - x309 - x311)*(-x312 - x314 - x316)) + V{0.00462962962962963}) + x122*x74 + x156*x304*(-x27*(x195 + x305 + x306) + x28*(x270 + x305) + x29*(x290 + x306));
auto x14 = V{1}*cell[14]*x122 - x105*(x157*(-x307*x312 + x328) + V{0.0740740740740741}) - x156*x196*(x155 - x27*x329);
auto x15 = V{1}*cell[15]*x122 - x105*(x157*(-x318*x330 + x324) + V{0.0740740740740741}) - x156*x196*(x213 - x269*x28);
auto x16 = V{1}*cell[16]*x122 - x105*(x157*(-x320*x331 + x322 + x325) + V{0.0740740740740741}) - x156*x196*(x224 - x289*x29);
auto x17 = V{1}*cell[17]*x122 - x105*(x225*(-x332*x333 + x334) + V{0.0185185185185185}) + x156*x259*(-x258 + x27*x335 + x28*x337);
auto x18 = V{1}*cell[18]*x122 - x105*(x225*(x204 - x261*x263 + x338) + V{0.0185185185185185}) - x156*x259*(x154 - x27*x340 + x28*(x301 + x339));
auto x19 = V{1}*cell[19]*x122 - x105*(x225*(-x341*x342 + x343) + V{0.0185185185185185}) + x156*x259*(x27*(x305 + x329) - x278 + x29*x344);
auto x20 = V{1}*cell[20]*x122 - x105*(x225*(x218 - x282*x285 + x338) + V{0.0185185185185185}) - x156*x259*(x149 - x27*(x274 + x329) + x29*(x339 + x345));
auto x21 = V{1}*cell[21]*x122 - x105*(x225*(x326 - (x318 + x320)*(x330 + x331)) + V{0.0185185185185185}) + x156*x259*(x28*(x269 + x305) + x29*(x289 + x306) - x296);
auto x22 = V{1}*cell[22]*x122 - x105*(x225*(x218 + x265 - x299*x300) + V{0.0185185185185185}) - x156*x259*(x141 - x28*(x269 + x274) + x29*(x268 + x345));
auto x23 = V{1}*cell[23]*x122 - x105*(x303*(x325 + x334 - (x320 + x332)*(x331 + x333)) + V{0.00462962962962963}) + x156*x304*(x27*(x305 + x335) + x28*(x305 + x337) + x29*(x306 + x344));
auto x24 = -x105*(x303*(x257 + x334 - (x309 + x332)*(x314 + x333)) + V{0.00462962962962963}) + x122*x79 + x156*x304*(x27*(x274 + x335) + x28*(x274 + x337) - x29*(x223 + x306 + x336));
auto x25 = -x105*(x303*(x277 + x343 - (x311 + x341)*(x316 + x342)) + V{0.00462962962962963}) + x122*x69 + x156*x304*(x27*(x305 + x340) - x28*(x212 + x305 + x336) + x29*(x302 + x336));
auto x26 = V{1}*cell[26]*x122 - x105*(x303*(x293 + x338 - (x106 + x110 - x291)*(x158 + x174 - x292)) + V{0.00462962962962963}) - x156*x304*(x27*(-x274 - x340) + x28*(x297 + x336) + x29*(x298 + x336));
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { x36, x317 + x319 + x321 };
}
};

}

}
