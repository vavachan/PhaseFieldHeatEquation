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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::BGK> >, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = cell.template getFieldComponent<collision::LES::SMAGORINSKY>(0);
auto x28 = parameters.template get<descriptors::OMEGA>();
auto x29 = cell[18] + cell[25] + cell[26];
auto x30 = cell[10] + cell[4];
auto x31 = cell[12] + x30;
auto x32 = cell[22] + cell[6];
auto x33 = cell[15] + cell[21];
auto x34 = cell[16] + cell[9];
auto x35 = cell[19] + cell[7];
auto x36 = cell[0] + cell[11] + cell[13] + cell[14] + cell[17] + cell[1] + cell[20] + cell[23] + cell[24] + cell[2] + cell[3] + cell[5] + cell[8] + x29 + x31 + x32 + x33 + x34 + x35 + V{1};
auto x37 = V{1} / (x36);
auto x38 = V{1}*x37;
auto x39 = -cell[26];
auto x40 = -cell[25];
auto x41 = -cell[18];
auto x42 = -cell[19];
auto x43 = cell[7] + x42;
auto x44 = -cell[23];
auto x45 = cell[11] + x44;
auto x46 = x43 + x45;
auto x47 = -cell[24];
auto x48 = -cell[17];
auto x49 = cell[5] + x48;
auto x50 = cell[13] + x47 + x49;
auto x51 = -cell[20];
auto x52 = -cell[14] + cell[1];
auto x53 = cell[6] + x51 + x52;
auto x54 = x31 + x39 + x40 + x41 + x46 + x50 + x53;
auto x55 = x54*x54;
auto x56 = V{0.666666666666667}*cell[10];
auto x57 = V{0.666666666666667}*cell[11];
auto x58 = V{0.666666666666667}*cell[12];
auto x59 = V{0.666666666666667}*cell[13];
auto x60 = V{0.666666666666667}*cell[23];
auto x61 = V{0.666666666666667}*cell[24];
auto x62 = V{0.666666666666667}*cell[25];
auto x63 = V{0.666666666666667}*cell[26];
auto x64 = -V{0.333333333333333}*cell[0];
auto x65 = -V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64;
auto x66 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x67 = V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9] - x38*x55 + x65 + x66;
auto x68 = -cell[12];
auto x69 = -cell[13];
auto x70 = cell[8] + x69;
auto x71 = -cell[15] + cell[2];
auto x72 = x47 + x71;
auto x73 = -cell[5] + x48;
auto x74 = -cell[22];
auto x75 = -cell[21];
auto x76 = cell[9] + x75;
auto x77 = x74 + x76;
auto x78 = x29 + x30 + x45 + x68 + x70 + x72 + x73 + x77;
auto x79 = x78*x78;
auto x80 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x81 = V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7] - x38*x79 + x65 + x80;
auto x82 = cell[26] + x40;
auto x83 = cell[20] + x82;
auto x84 = -cell[11];
auto x85 = cell[12] + x44 + x84;
auto x86 = -cell[16] + cell[3];
auto x87 = cell[10] + cell[24] + x86;
auto x88 = -cell[9];
auto x89 = x75 + x88;
auto x90 = -cell[7] + x42;
auto x91 = x32 + x70 + x83 + x85 + x87 + x89 + x90;
auto x92 = x91*x91;
auto x93 = V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5] - x38*x92 + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64 + x66 + x80;
auto x94 = x37*x91;
auto x95 = -cell[8];
auto x96 = cell[22] + x95;
auto x97 = cell[11] + cell[12] + x69;
auto x98 = cell[25] + x39 + x44;
auto x99 = -cell[10];
auto x100 = cell[24] + x99;
auto x101 = x100 + x76 + x78*x94 + x96 + x97 + x98;
auto x102 = -cell[6];
auto x103 = cell[13] + x68;
auto x104 = x100 + x102 + x103 + x46 + x54*x94 + x83;
auto x105 = -cell[4] + x29 + x37*x54*x78 + x50 + x85 + x99;
auto x106 = V{1} / (V{3.00000046417339}*util::sqrt(x37*(x27*x27)*util::sqrt(x101*x101 + x104*x104 + x105*x105 + V{0.5}*(x67*x67) + V{0.5}*(x81*x81) + V{0.5}*(x93*x93)) + V{0.0277777691819762}/((x28)*(x28))) + V{0.5}/x28);
auto x107 = V{1} / ((x36)*(x36));
auto x108 = V{1.5}*x107;
auto x109 = -x54;
auto x110 = x109*x109;
auto x111 = x108*x110;
auto x112 = -x91;
auto x113 = x112*x112;
auto x114 = x108*x113;
auto x115 = -x78;
auto x116 = x115*x115;
auto x117 = x108*x116;
auto x118 = x114 + x117 + V{-1};
auto x119 = x111 + x118;
auto x120 = V{1} - x106;
auto x121 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x122 = V{3}*cell[13];
auto x123 = V{3}*cell[5];
auto x124 = V{3}*cell[7];
auto x125 = V{3}*cell[18];
auto x126 = V{3}*cell[20];
auto x127 = V{3}*cell[26];
auto x128 = V{3}*cell[10];
auto x129 = V{3}*cell[11];
auto x130 = -V{3}*cell[23];
auto x131 = V{3}*cell[24];
auto x132 = -V{3}*cell[17] + V{3}*cell[4] + x128 + x129 + x130 - x131;
auto x133 = V{3}*cell[12];
auto x134 = V{3}*cell[25];
auto x135 = -V{3}*cell[19] + V{3}*cell[6] + x133 - x134;
auto x136 = -V{3}*cell[14] + V{3}*cell[1] + x122 + x123 + x124 - x125 - x126 - x127 + x132 + x135;
auto x137 = -x136*x37;
auto x138 = V{3}*x107;
auto x139 = V{3}*cell[9];
auto x140 = V{3}*cell[22];
auto x141 = -V{3}*cell[21] + V{3}*cell[8] - x122 + x127;
auto x142 = -V{3}*cell[15] + V{3}*cell[2] - x123 + x125 + x132 - x133 + x134 + x139 - x140 + x141;
auto x143 = -x142*x37;
auto x144 = x111 + V{-1};
auto x145 = -V{3}*cell[16] + V{3}*cell[3] - x124 + x126 + x128 - x129 + x130 + x131 + x135 - x139 + x140 + x141;
auto x146 = -x145*x37;
auto x147 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x148 = V{4.5}*x107;
auto x149 = x43 + x53;
auto x150 = -V{2}*cell[17] + V{2}*cell[4];
auto x151 = -V{2}*cell[23];
auto x152 = V{2}*cell[10];
auto x153 = cell[8] + x151 + x152;
auto x154 = V{2}*cell[11] - V{2}*cell[24];
auto x155 = x154 + x71;
auto x156 = x149 + x150 + x153 + x155 + x77;
auto x157 = -x156;
auto x158 = x119 + x137;
auto x159 = x143 + x158;
auto x160 = V{2}*cell[25];
auto x161 = V{2}*cell[12];
auto x162 = -x160 + x161;
auto x163 = V{2}*cell[26];
auto x164 = V{2}*cell[13];
auto x165 = -x163 + x164;
auto x166 = V{2}*cell[18];
auto x167 = V{2}*cell[5];
auto x168 = -cell[2] - x166 + x167;
auto x169 = x149 + x162 + x165 + x168 + x33 + x88 + x96;
auto x170 = -x143;
auto x171 = x158 + x170;
auto x172 = cell[4] + x86;
auto x173 = -V{2}*cell[19] + V{2}*cell[6] + x52;
auto x174 = x41 + x49;
auto x175 = cell[22] + x153 + x162 + x172 + x173 + x174 + x89;
auto x176 = -x175;
auto x177 = -x146;
auto x178 = -cell[3];
auto x179 = cell[4] + x178;
auto x180 = V{2}*cell[20];
auto x181 = V{2}*cell[7];
auto x182 = -x180 + x181 + x52;
auto x183 = cell[21] + x154 + x165 + x174 + x179 + x182 + x34 + x74 + x95;
auto x184 = cell[18] + x73;
auto x185 = -V{2}*cell[21] + V{2}*cell[8];
auto x186 = x185 + x71;
auto x187 = cell[20] + cell[6] + x151 + x152 + x163 - x164 + x172 + x184 + x186 + x90;
auto x188 = -x187;
auto x189 = x119 + x143;
auto x190 = x146 + x189;
auto x191 = V{2}*cell[22];
auto x192 = V{2}*cell[9];
auto x193 = cell[16] - x191 + x192;
auto x194 = x102 + x155 + x160 - x161 + x179 + x184 + x193 + x35 + x51;
auto x195 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x196 = V{3}*cell[10] - V{3}*cell[23] + x150 + x173 + x185 + x72 + x82 + x86 + x97;
auto x197 = -x196;
auto x198 = cell[10] + V{3}*cell[11] - V{3}*cell[24] + x103 + x150 + x178 + x182 + x193 + x71 + x98;
auto x199 = x44 + x84 + x87;
auto x200 = V{3}*cell[12] + cell[13] + cell[15] - V{3}*cell[25] + x168 + x173 + x191 - x192 + x199 + x39;
auto x201 = x136*x37;
auto x202 = -V{3}*cell[13] + cell[14] - cell[1] + cell[25] + V{3}*cell[26] + x166 - x167 + x180 - x181 + x186 + x199 + x68;
auto x203 = -x202;
auto x204 = x145*x37;
auto x205 = x108*x79;
auto x206 = x108*x92 + V{-1};
auto x207 = x205 + x206;
auto x208 = x142*x37;
auto x209 = x108*x55;
auto x210 = x208 + x209;
auto x211 = x204 + x207 + x210;
auto x212 = x201 + x207;
auto x213 = x204 + x209;
auto x214 = x210 + x212;
auto x215 = -x137;
auto x216 = -x169;
auto x217 = x212 + x213;
auto x218 = -x183;
auto x219 = x119 + x146;
auto x220 = -x194;
auto x221 = -x198;
auto x222 = -x200;
auto x0 = V{1}*cell[0]*x120 - x106*(x119*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + V{0.296296296296296});
auto x1 = V{1}*cell[1]*x120 - x106*(x121*(-x110*x138 + x118 + x137) + V{0.0740740740740741});
auto x2 = V{1}*cell[2]*x120 - x106*(x121*(x114 - x116*x138 + x143 + x144) + V{0.0740740740740741});
auto x3 = V{1}*cell[3]*x120 - x106*(x121*(-x113*x138 + x117 + x144 + x146) + V{0.0740740740740741});
auto x4 = V{1}*cell[4]*x120 - x106*(x147*(-x148*x157*x157 + x159) + V{0.0185185185185185});
auto x5 = V{1}*cell[5]*x120 - x106*(x147*(-x148*x169*x169 + x171) + V{0.0185185185185185});
auto x6 = V{1}*cell[6]*x120 - x106*(x147*(x146 - x148*x176*x176 + x158) + V{0.0185185185185185});
auto x7 = V{1}*cell[7]*x120 - x106*(x147*(-x148*x183*x183 + x158 + x177) + V{0.0185185185185185});
auto x8 = V{1}*cell[8]*x120 - x106*(x147*(-x148*x188*x188 + x190) + V{0.0185185185185185});
auto x9 = V{1}*cell[9]*x120 - x106*(x147*(-x148*x194*x194 + x177 + x189) + V{0.0185185185185185});
auto x10 = V{1}*cell[10]*x120 - x106*(x195*(x146 - x148*x197*x197 + x159) + V{0.00462962962962963});
auto x11 = V{1}*cell[11]*x120 - x106*(x195*(-x148*x198*x198 + x159 + x177) + V{0.00462962962962963});
auto x12 = V{1}*cell[12]*x120 - x106*(x195*(x146 - x148*x200*x200 + x171) + V{0.00462962962962963});
auto x13 = V{1}*cell[13]*x120 - x106*(x195*(-x148*x203*x203 - x201 + x211) + V{0.00462962962962963});
auto x14 = V{1}*cell[14]*x120 - x106*(x121*(-x138*x55 + x212) + V{0.0740740740740741});
auto x15 = V{1}*cell[15]*x120 - x106*(x121*(-x138*x79 + x206 + x210) + V{0.0740740740740741});
auto x16 = V{1}*cell[16]*x120 - x106*(x121*(-x138*x92 + x205 + x213 + V{-1}) + V{0.0740740740740741});
auto x17 = V{1}*cell[17]*x120 - x106*(x147*(-x148*x156*x156 + x214) + V{0.0185185185185185});
auto x18 = V{1}*cell[18]*x120 - x106*(x147*(-x148*x216*x216 + x189 + x215) + V{0.0185185185185185});
auto x19 = V{1}*cell[19]*x120 - x106*(x147*(-x148*x175*x175 + x217) + V{0.0185185185185185});
auto x20 = V{1}*cell[20]*x120 - x106*(x147*(-x148*x218*x218 + x215 + x219) + V{0.0185185185185185});
auto x21 = V{1}*cell[21]*x120 - x106*(x147*(-x148*x187*x187 + x211) + V{0.0185185185185185});
auto x22 = V{1}*cell[22]*x120 - x106*(x147*(-x148*x220*x220 + x170 + x219) + V{0.0185185185185185});
auto x23 = V{1}*cell[23]*x120 - x106*(x195*(-x148*x196*x196 + x204 + x214) + V{0.00462962962962963});
auto x24 = V{1}*cell[24]*x120 - x106*(x195*(-x148*x221*x221 - x204 + x214) + V{0.00462962962962963});
auto x25 = V{1}*cell[25]*x120 - x106*(x195*(-x148*x222*x222 - x208 + x217) + V{0.00462962962962963});
auto x26 = V{1}*cell[26]*x120 - x106*(x195*(-x148*x202*x202 + x190 + x215) + V{0.00462962962962963});
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
return { x36, V{1}*x107*(x55 + x79 + x92) };
}
};

}

}
