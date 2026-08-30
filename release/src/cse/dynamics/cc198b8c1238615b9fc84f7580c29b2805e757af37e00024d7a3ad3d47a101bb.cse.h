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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x12 = parameters.template get<descriptors::OMEGA>();
auto x11 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x9 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x10 = x12 + V{-1};
auto x13 = V{0.0740740740740741}*x9;
auto x14 = V{3}*x11;
auto x15 = V{1} / (x9);
auto x16 = V{0.111111111111111}*x9;
auto x17 = x15*(cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1});
auto x18 = V{1} - x17;
auto x19 = -x18;
auto x20 = x19*x19;
auto x21 = V{1.5}*x20;
auto x22 = -x21;
auto x23 = x11*x11;
auto x24 = V{3}*x23 + V{1};
auto x25 = x22 + x24;
auto x26 = x14 + x25;
auto x27 = -x14;
auto x28 = x25 + x27;
auto x29 = V{0.444444444444444}*x9;
auto x30 = V{1.5}*x23;
auto x31 = x30 + V{-1};
auto x32 = x21 + x31;
auto x33 = -x32;
auto x34 = V{0.222222222222222}*x9;
auto x35 = V{3}*x20;
auto x36 = x15*(V{3}*cell[0] + V{3}*cell[4] + V{6}*cell[5] + V{6}*cell[6] + V{6}*cell[7] + V{3}*cell[8] + V{3});
auto x37 = -x30 + V{-2};
auto x38 = x36 + x37;
auto x39 = x35 + x38;
auto x40 = V{0.0555555555555555}*x9;
auto x41 = x11 + x17 + V{-1};
auto x42 = V{4.5}*(x41*x41);
auto x43 = x14 + x22 + x38 + x42;
auto x44 = -x11 - x18;
auto x45 = x30 + V{2};
auto x46 = -x36 + x45;
auto x47 = x21 + x46;
auto x48 = x14 + x47 - V{4.5}*x44*x44;
auto x49 = -x48;
auto x50 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x51 = x15*(x16*x26 + x16*x28 + x29*x33 + x34*x39 + x40*x43 + x40*x49 + x50);
auto x52 = V{1} - x51;
auto x53 = x11 + x52;
auto x54 = -x53;
auto x55 = x51 + V{-1};
auto x56 = x55*x55;
auto x57 = V{1.5}*x56;
auto x58 = V{0.166666666666667}*x9;
auto x59 = V{0.333333333333333}*x9;
auto x60 = x15*(V{4.9960036108132e-16}*cell[4] + V{3.99680288865056e-15}*cell[5] + V{1.99840144432528e-15}*cell[6] + V{2.66453525910038e-15}*cell[7] + V{4.9960036108132e-16}*cell[8] + x26*x59 + x28*x59 + V{1.33333333333333}*x33*x9 + V{0.666666666666667}*x39*x9 + x43*x58 + x49*x58 + V{3.33066907387547e-16});
auto x61 = x45 - x60;
auto x62 = x57 + x61;
auto x63 = x14 + x62 - V{4.5}*x54*x54;
auto x64 = x11 + x55;
auto x65 = V{4.5}*(x64*x64);
auto x66 = x27 + x62 - x65;
auto x67 = V{3}*x56;
auto x68 = x61 - x67;
auto x69 = V{1.15648231731787e-17}*x9;
auto x70 = x27 - x42 + x47;
auto x71 = -x35 + x46;
auto x72 = x24 - V{1.5}*x18*x18;
auto x73 = x14 + x72;
auto x74 = x27 + x72;
auto x75 = x52*x52;
auto x76 = x24 - V{1.5}*x75;
auto x77 = x14 + x76;
auto x78 = x27 + x76;
auto x79 = V{0.0277777777777778}*x9;
auto x80 = x30 + x57 + x60 + V{-4};
auto x81 = V{1.01192202765314e-18}*x9;
auto x82 = V{1.87928376564154e-18}*x9;
auto x83 = V{0.0231481481481481}*x9;
auto x84 = -x63;
auto x85 = V{0.00462962962962963}*x9;
auto x86 = -x57;
auto x87 = x37 + x60;
auto x88 = x14 + x65 + x86 + x87;
auto x89 = V{1.73472347597681e-18}*x9;
auto x90 = V{0.00925925925925926}*x9;
auto x91 = x24 + x86;
auto x92 = x14 + x91;
auto x93 = x27 + x91;
auto x94 = x67 + x87;
auto x95 = -V{0.0833333333333334}*cell[4] - V{0.166666666666667}*cell[6] - V{0.0833333333333334}*cell[8] + x26*x89 + x28*x89 + V{4.62592926927149e-18}*x39*x9 + V{0.0185185185185185}*x9*x94 + x90*x92 + x90*x93 + V{-0.0555555555555556};
auto x96 = -x10*(-V{0.833333333333333}*cell[5] + V{0.166666666666667}*cell[7] + x43*x82 + x49*x81 + x83*x84 - x85*x88 + x95) + V{0.0277777777777778};
auto x97 = V{0.333333333333333}*cell[7];
auto x98 = V{0.333333333333334}*cell[5];
auto x99 = V{3.85185988877447e-34}*x9;
auto x100 = V{0.00925925925925926}*x9;
auto x101 = V{1.44560289664734e-18}*x9;
auto x102 = V{0.0185185185185185}*x9;
auto x103 = -x10*(V{0.166666666666667}*cell[4] - V{0.666666666666667}*cell[6] + V{0.166666666666667}*cell[8] + x100*x84 + x100*x88 + x101*x43 + x101*x49 - x102*x92 - x102*x93 + x26*x99 + x28*x99 + V{1.15648231731787e-17}*x39*x9 + V{0.0740740740740741}*x9*x94 - x97 - x98 + V{-0.0555555555555556}) + V{0.111111111111111};
auto x104 = -x64;
auto x105 = -x10*(V{0.166666666666666}*cell[5] - V{0.833333333333333}*cell[7] + x43*x81 + x49*x82 + x83*x88 - x84*x85 + x95) + V{0.0277777777777778};
auto x106 = -x15*(x16*x73 + x16*x74 - x29*x32 - x34*x71 - x40*x48 - x40*x70 + x50) + V{1};
auto x107 = x106*x106;
auto x108 = -V{1.5}*x107 + x24;
auto x109 = x10*(-V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[8] - x100*x63 - x100*x66 - x101*x48 - x101*x70 + V{0.037037037037037}*x68*x9 + V{2.31296463463574e-18}*x71*x9 + V{3.46944695195361e-18}*x73*x9 + V{3.46944695195361e-18}*x74*x9 + V{0.037037037037037}*x77*x9 + V{0.037037037037037}*x78*x9 - x97 - x98 + V{-0.0555555555555556}) + V{-0.111111111111111};
auto x0 = -x10*(-V{0.666666666666667}*cell[4] - V{2.66666666666667}*cell[5] - V{1.33333333333333}*cell[6] - V{2.66666666666667}*cell[7] - V{0.666666666666667}*cell[8] - x13*x63 - x13*x66 - x48*x69 - V{0.148148148148148}*x68*x9 - x69*x70 - V{3.70074341541719e-17}*x71*x9 + V{1.38777878078145e-17}*x73*x9 + V{1.38777878078145e-17}*x74*x9 + V{0.0740740740740741}*x77*x9 + V{0.0740740740740741}*x78*x9 + V{-0.444444444444444}) - x29*(x31 + x57) + V{-0.444444444444444};
auto x1 = -(x79*(x27 + x80 - V{4.5}*x53*x53) + x96);
auto x2 = -x103 - x16*(-V{4.5}*x75 + x80);
auto x3 = -(x105 + x79*(x14 + x80 - V{4.5}*x104*x104));
auto x4 = x109 + x16*(x108 + x27);
auto x5 = -x63*x79 - x96;
auto x6 = -x103 - x16*x68;
auto x7 = -x105 - x66*x79;
auto x8 = x109 + x16*(x108 + x14);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x9, V{1}*x107 + x23 };
}
};

}

}
