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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x23 = x22 + V{-1};
auto x24 = V{1} / (x20 + V{1});
auto x25 = V{2}*cell[13] + V{2}*cell[17] + V{2}*cell[18] + V{2}*cell[5];
auto x26 = cell[0] + cell[10] + V{2}*cell[11] + cell[12] + cell[15] + cell[16] + cell[1] + cell[3] + cell[6] + cell[7] + x25 + V{1};
auto x27 = x24*x26;
auto x28 = V{0.0277777777777778}*x27;
auto x29 = x19*x19;
auto x30 = V{3}*x19;
auto x31 = x21*x21;
auto x32 = V{1.5}*x31;
auto x33 = -x32;
auto x34 = x20*x20;
auto x35 = V{1.5}*x34;
auto x36 = V{1} - x35;
auto x37 = x33 + x36;
auto x38 = x30 + x37;
auto x39 = V{3}*x29 + x38;
auto x40 = V{3}*x21;
auto x41 = V{1.5}*x29;
auto x42 = -x41;
auto x43 = x40 + x42;
auto x44 = V{3}*x31 + x36 + x43;
auto x45 = x32 + x35 + V{-1};
auto x46 = -V{3}*x29 + x30 + x45;
auto x47 = -x46;
auto x48 = x41 + V{-1};
auto x49 = -V{3}*x31 + x35 + x40 + x48;
auto x50 = -x49;
auto x51 = x19 + x21;
auto x52 = x51*x51;
auto x53 = x38 + x43 + V{4.5}*x52;
auto x54 = -x30;
auto x55 = -x21;
auto x56 = x19 + x55;
auto x57 = -V{4.5}*x56*x56;
auto x58 = x41 + x45;
auto x59 = x40 + x58;
auto x60 = x54 + x57 + x59;
auto x61 = -x60;
auto x62 = -x56;
auto x63 = -x40;
auto x64 = x30 + x58;
auto x65 = x63 + x64;
auto x66 = x65 - V{4.5}*x62*x62;
auto x67 = -x66;
auto x68 = x30 - V{4.5}*x52 + x59;
auto x69 = -x68;
auto x70 = V{0.0555555555555556}*x27;
auto x71 = V{3}*x20;
auto x72 = x42 + x71;
auto x73 = x33 + V{3}*x34 + x72 + V{1};
auto x74 = x19 + x20;
auto x75 = V{4.5}*(x74*x74);
auto x76 = x38 + x72 + x75;
auto x77 = x20 + x21;
auto x78 = V{4.5}*(x77*x77);
auto x79 = x37 + x43 + x71 + x78;
auto x80 = x19 - x20;
auto x81 = -x80;
auto x82 = -x71;
auto x83 = x64 + x82;
auto x84 = x83 - V{4.5}*x81*x81;
auto x85 = -x84;
auto x86 = x20 + x55;
auto x87 = x59 + x82 - V{4.5}*x86*x86;
auto x88 = -x87;
auto x89 = V{0.333333333333333}*x27;
auto x90 = -x58;
auto x91 = x89*x90;
auto x92 = V{0.0277777777777778}*x27;
auto x93 = V{0.0555555555555555}*x27;
auto x94 = V{0.0555555555555555}*x27;
auto x95 = V{0.0277777777777778}*x27;
auto x96 = V{0.0555555555555555}*x27;
auto x97 = V{0.0277777777777778}*x27;
auto x98 = V{1.66533453693773e-16}*cell[10] + V{4.44089209850063e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{8.88178419700125e-16}*cell[17] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[3] + V{6.66133814775094e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{0.111111111111111}*x27*x73 + x39*x70 + x44*x70 + x53*x97 + x76*x96 + x79*x93 + V{2.22044604925031e-16};
auto x99 = x47*x94 + x50*x93 + x61*x97 + x67*x92 + x69*x95 + x85*x96 + x88*x93 + x91 + x98;
auto x100 = V{0.0462962962962963}*x27;
auto x101 = V{0.00925925925925926}*x27;
auto x102 = V{0.166666666666667}*cell[13];
auto x103 = V{0.166666666666667}*cell[5];
auto x104 = V{0.0833333333333333}*cell[12];
auto x105 = V{0.0833333333333333}*cell[3];
auto x106 = V{0.00462962962962963}*x27;
auto x107 = x106*x44;
auto x108 = V{0.00925925925925926}*x27;
auto x109 = V{0.00462962962962963}*x27;
auto x110 = x106*x49;
auto x111 = x109*x76;
auto x112 = V{0.00231481481481482}*x27;
auto x113 = -V{0.166666666666667}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x108*x73 - x112*x53 + x112*x60 + x112*x66 + x112*x68 + V{-0.0555555555555555};
auto x114 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1] + x102 + x103 - x104 - x105 + x107 + x108*x79 - x108*x87 + x109*x84 - x110 - x111 + x113;
auto x115 = V{0.0555555555555556}*x22;
auto x116 = -x32 + V{3}*x34 - x48 - x71;
auto x117 = V{0.0833333333333333}*cell[10];
auto x118 = V{0.0833333333333333}*cell[1];
auto x119 = V{0.166666666666667}*cell[17];
auto x120 = V{0.166666666666667}*cell[18];
auto x121 = x109*x79;
auto x122 = x106*x39;
auto x123 = -V{0.333333333333333}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x102 - x103 + x104 + x105 - x106*x53 - x107 + x111 + x117 + x118 - x119 - x120 + x121 - x122 + V{0.0555555555555555};
auto x124 = x106*x46;
auto x125 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] + V{0.166666666666667}*cell[3] - V{0.333333333333333}*cell[5] + x108*x76 - x108*x84 + x109*x87 + x113 - x117 - x118 + x119 + x120 - x121 + x122 - x124;
auto x126 = x64 + x71 - x75;
auto x127 = V{0.0231481481481481}*x27;
auto x128 = V{0.00462962962962963}*x27;
auto x129 = V{0.00231481481481481}*x27;
auto x130 = V{0.00462962962962963}*x27;
auto x131 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - x112*x79 + x112*x87 - x130*x39 + x130*x46 + V{0.0138888888888889};
auto x132 = V{0.00115740740740741}*x27;
auto x133 = V{0.166666666666667}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x101*x73 - x132*x53 + x132*x60 + x132*x66 + x132*x68;
auto x134 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x129*x44 - x129*x49 + x131 + x133;
auto x135 = V{0.833333333333333}*cell[13] - V{0.166666666666667}*cell[5] - x128*x84 + x134;
auto x136 = V{0.0277777777777778}*x22;
auto x137 = -V{0.166666666666667}*cell[13] + V{0.833333333333333}*cell[5] + x128*x76 + x134;
auto x138 = V{0.0115740740740741}*x27;
auto x139 = V{0.0162037037037037}*x27;
auto x140 = V{0.00231481481481481}*x27;
auto x141 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[5] - x112*x76 + x112*x84 - x130*x44 + x130*x49;
auto x142 = -V{0.0833333333333333}*cell[11] + x106*x73 + x131 + x141;
auto x143 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x140*x60 - x140*x66 + x142;
auto x144 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x140*x53 - x140*x68 + x142;
auto x145 = x59 + x71 - x78;
auto x146 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x129*x39 - x129*x46 + x133 + x141 + V{0.0138888888888889};
auto x147 = V{0.833333333333333}*cell[17] - V{0.166666666666667}*cell[18] - x128*x87 + x146;
auto x148 = -x86;
auto x149 = x58 + x71;
auto x150 = x149 + x63 - V{4.5}*x148*x148;
auto x151 = -V{0.166666666666667}*cell[17] + V{0.833333333333333}*cell[18] + x128*x79 + x146;
auto x152 = -V{4.5}*x80*x80;
auto x153 = x149 + x152 + x54;
auto x0 = -V{0.333333333333333}*x22*(-x24*x90*x99 + V{1}) + x23*(V{0.5}*cell[10] + V{1}*cell[11] + V{0.5}*cell[12] + V{1}*cell[15] + V{1}*cell[16] + V{0.5}*cell[1] + V{0.5}*cell[3] + V{1}*cell[6] + V{1}*cell[7] + x25 - x28*x39 - x28*x44 - x28*x47 - x28*x50 - x28*x53 - x28*x61 - x28*x67 - x28*x69 - x70*x73 - x70*x76 - x70*x79 - x70*x85 - x70*x88 - x91 + V{0.833333333333333});
auto x1 = -x115*(-x24*x47*x99 + V{1}) - x23*(-x100*x46 - x101*x39 + x114);
auto x2 = -x115*(-x116*x24*x99 + V{1}) + x23*(-x106*x47 - x106*x50 - x106*x61 - x106*x67 - x106*x69 + x109*x85 + x109*x88 - x116*x70 + x123 + V{0.0185185185185185}*x27*x73);
auto x3 = -x115*(-x24*x50*x99 + V{1}) - x23*(-x100*x49 - x101*x44 + x125);
auto x4 = -x136*(x126*x24*x99 + V{1}) - x23*(-x126*x28 - x127*x76 + x135);
auto x5 = -x136*(-x24*x85*x99 + V{1}) - x23*(-x106*x84 + x137);
auto x6 = -x136*(-x24*x69*x99 + V{1}) - x23*(-x138*x53 - x139*x68 + x143);
auto x7 = -x136*(-x24*x67*x99 + V{1}) - x23*(x138*x60 - x139*x66 + x144);
auto x8 = -x136*(x145*x24*x99 + V{1}) - x23*(-x127*x79 - x145*x28 + x147);
auto x9 = -x136*(x150*x24*x99 + V{1}) - x23*(x127*x87 - x150*x28 + x151);
auto x10 = -x115*(-x24*x39*x99 + V{1}) - x23*(x100*x39 + x101*x46 + x114);
auto x11 = -x115*(-x24*x73*x99 + V{1}) - x23*(-x106*x60 - x106*x66 - x106*x68 - x110 - x123 - x124 + V{0.037037037037037}*x24*x26*x73 + V{0.00462962962962963}*x24*x26*x84 + V{0.00462962962962963}*x24*x26*x87);
auto x12 = -x115*(-x24*x44*x99 + V{1}) - x23*(x100*x44 + x101*x49 + x125);
auto x13 = -x136*(-x24*x76*x99 + V{1}) - x23*(x106*x76 + x135);
auto x14 = -x136*(x153*x24*x99 + V{1}) - x23*(x127*x84 + x137 - x153*x28);
auto x15 = -x136*(-x24*x53*x99 + V{1}) - x23*(x138*x68 + x139*x53 + x143);
auto x16 = -x136*(-x24*x61*x99 + V{1}) - x23*(x138*x66 - x139*x60 + x144);
auto x17 = -x136*(-x24*x79*x99 + V{1}) - x23*(x106*x79 + x147);
auto x18 = -x136*(-x24*x88*x99 + V{1}) - x23*(-x106*x87 + x151);
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
return { V{1}*x24*(-x46*x94 - x49*x93 - x58*x89 - x60*x97 - x68*x95 - x87*x93 - x92*(x57 + x65) - x96*(x152 + x83) + x98), x29 + x31 + x34 };
}
};

}

}
