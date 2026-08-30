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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x13 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x14 = parameters.template get<descriptors::OMEGA>();
auto x9 = x14 + V{-1};
auto x10 = V{0.0740740740740741}*x11;
auto x12 = V{3}*x13;
auto x15 = V{1} / (x11);
auto x16 = cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1};
auto x17 = -x15*x16 + V{1};
auto x18 = x13 + x17;
auto x19 = -x18;
auto x20 = -x17;
auto x21 = x20*x20;
auto x22 = V{1.5}*x21;
auto x23 = x15*(V{3}*cell[0] + V{3}*cell[4] + V{6}*cell[5] + V{6}*cell[6] + V{6}*cell[7] + V{3}*cell[8] + V{3});
auto x24 = x13*x13;
auto x25 = V{1.5}*x24;
auto x26 = x25 + V{2};
auto x27 = -x23 + x26;
auto x28 = x22 + x27;
auto x29 = x12 + x28 - V{4.5}*x19*x19;
auto x30 = -x12;
auto x31 = x13 + x15*x16 + V{-1};
auto x32 = V{4.5}*(x31*x31);
auto x33 = -x32;
auto x34 = x30 + x33;
auto x35 = x28 + x34;
auto x36 = V{0.148148148148148}*x11;
auto x37 = V{3}*x21;
auto x38 = x27 - x37;
auto x39 = V{0.444444444444444}*x11;
auto x40 = x25 + V{-1};
auto x41 = x22 + x40;
auto x42 = x39*x41;
auto x43 = x17*x17;
auto x44 = V{1.5}*x43;
auto x45 = V{3}*x24 + V{1};
auto x46 = -x44 + x45;
auto x47 = x12 + x46;
auto x48 = x10*x47;
auto x49 = x30 + x46;
auto x50 = x10*x49;
auto x51 = V{0.666666666666667}*cell[4] + V{2.66666666666667}*cell[5] + V{1.33333333333333}*cell[6] + V{2.66666666666667}*cell[7] + V{0.666666666666667}*cell[8] + V{0.888888888888889};
auto x52 = V{0.111111111111111}*x11;
auto x53 = -x22;
auto x54 = x45 + x53;
auto x55 = x12 + x54;
auto x56 = x30 + x54;
auto x57 = -x41;
auto x58 = x39*x57;
auto x59 = V{0.222222222222222}*x11;
auto x60 = x23 - x25 + V{-2};
auto x61 = x37 + x60;
auto x62 = V{0.0555555555555555}*x11;
auto x63 = x12 + x32 + x53 + x60;
auto x64 = -x29;
auto x65 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x66 = x15*(x52*x55 + x52*x56 + x58 + x59*x61 + x62*x63 + x62*x64 + x65);
auto x67 = x66 + V{-1};
auto x68 = x67*x67;
auto x69 = V{1.5}*x68;
auto x70 = V{0.444444444444444}*x14*(x11*(x40 + x69) + V{1});
auto x71 = V{0.0231481481481481}*x11;
auto x72 = x44 - V{4.5}*x18*x18;
auto x73 = x12 + x27 + x72;
auto x74 = V{0.0277777777777778}*x11;
auto x75 = x25 + V{-4};
auto x76 = x23 + x75;
auto x77 = V{0.00462962962962963}*x11;
auto x78 = x27 + x34 + x44;
auto x79 = V{0.0185185185185185}*x11;
auto x80 = x27 - V{3}*x43;
auto x81 = V{0.00925925925925926}*x11;
auto x82 = V{0.0833333333333333}*cell[4] + V{0.166666666666667}*cell[6] + V{0.0833333333333333}*cell[8] - x47*x81 - x49*x81 + V{0.0277777777777778};
auto x83 = x79*x80 + x82;
auto x84 = V{0.833333333333333}*cell[5] - V{0.166666666666667}*cell[7];
auto x85 = -x77*x78 + x83 + x84;
auto x86 = V{0.0277777777777778}*x14;
auto x87 = V{1} - x66;
auto x88 = x13 + x87;
auto x89 = V{0.166666666666667}*x11;
auto x90 = V{0.333333333333333}*x11;
auto x91 = x15*(V{4.9960036108132e-16}*cell[4] + V{3.99680288865056e-15}*cell[5] + V{1.99840144432528e-15}*cell[6] + V{2.66453525910038e-15}*cell[7] + V{4.9960036108132e-16}*cell[8] + V{1.33333333333333}*x11*x57 + V{0.666666666666667}*x11*x61 + x55*x90 + x56*x90 + x63*x89 + x64*x89 + V{3.33066907387547e-16});
auto x92 = x69 + x75 + x91;
auto x93 = V{0.0185185185185185}*x11;
auto x94 = -V{0.166666666666667}*cell[4] + V{0.666666666666667}*cell[6] - V{0.166666666666667}*cell[8] + x47*x93 + x49*x93;
auto x95 = V{0.00925925925925926}*x11;
auto x96 = V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[7] + V{-0.0555555555555555};
auto x97 = x29*x95 + x35*x95 + x96;
auto x98 = x94 + x97;
auto x99 = V{0.111111111111111}*x14;
auto x100 = x87*x87;
auto x101 = -V{0.166666666666667}*cell[5] + V{0.833333333333333}*cell[7];
auto x102 = x101 - x73*x77 + x83;
auto x103 = x13 + x67;
auto x104 = -x103;
auto x105 = V{0.037037037037037}*x11;
auto x106 = V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[8];
auto x107 = -x105*x47 + x106 + x50;
auto x108 = V{0.037037037037037}*x11;
auto x109 = x73*x95 + x78*x95 + x96;
auto x110 = -x108*x80 + x109;
auto x111 = -x15*(-x29*x62 - x35*x62 - x38*x59 - x42 + x47*x52 + x49*x52 + x65) + V{1};
auto x112 = x111*x111;
auto x113 = -V{1.5}*x112 + x45;
auto x114 = V{0.00462962962962963}*x11;
auto x115 = -x88;
auto x116 = x26 - x91;
auto x117 = x116 + x69;
auto x118 = x11*(x117 + x12 - V{4.5}*x115*x115) + V{1};
auto x119 = V{0.037037037037037}*x11;
auto x120 = x11*(x116 - V{3}*x68) + V{1};
auto x121 = x11*(x117 + x30 - V{4.5}*x103*x103) + V{1};
auto x122 = -x105*x49 + x106 + x48;
auto x123 = -x108*x38 + x97;
auto x124 = V{2}*x14 + V{-2};
auto x125 = x38*x79 + x82;
auto x126 = V{0.0555555555555556}*x14;
auto x127 = -V{1.5}*x100 + x45;
cell[0] = -x70 + x9*(x10*x29 + x10*x35 + x36*x38 + x42 - x48 - x50 + x51);
cell[1] = -(x86*(x11*(x30 + x92 - V{4.5}*x88*x88) + V{1}) + x9*(x71*x73 - x74*(x30 + x72 + x76) + x85));
cell[2] = x9*(-V{0.0740740740740741}*x11*x38 + V{0.111111111111111}*x11*(x22 - V{4.5}*x43 + x76) - x98) - x99*(x11*(-V{4.5}*x100 + x92) + V{1});
cell[3] = -(x86*(x11*(x12 + x92 - V{4.5}*x104*x104) + V{1}) + x9*(x102 + x71*x78 - x74*(x12 + x33 + x44 + x76)));
cell[4] = -x9*(x107 + x110) + x99*(x11*(x113 + x30) + V{-1});
cell[5] = -x118*x86 - x9*(-x114*x73 + x85);
cell[6] = -x120*x99 - x9*(x109 - x119*x80 + x94);
cell[7] = -x121*x86 - x9*(x102 - x114*x78);
cell[8] = -x9*(x110 + x122) + x99*(x11*(x113 + x12) + V{-1});
cell.template getFieldPointer<descriptors::VELOCITY>()[0] = V{1}*x15*(-x118*x126 - V{0.222222222222222}*x120*x14 - x121*x126 - x124*(-x119*x38 + x98) - x124*(x101 - x114*x35 + x125 - x29*x77) - x124*(-x114*x29 + x125 - x35*x77 + x84) + V{0.111111111111111}*x14*(x11*(x12 + x127) + V{-1}) + V{0.111111111111111}*x14*(x11*(x127 + x30) + V{-1}) - x70 - x9*(x107 + x123) - x9*(x122 + x123) + x9*(-x10*x55 - x10*x56 - x10*x63 - x10*x64 - x36*x61 + x51 - x58) + V{1}) + V{-1};
cell.template getFieldPointer<descriptors::VELOCITY>()[1] = x13;
return { x11, V{1}*x112 + x24 };
}
};

}

}
