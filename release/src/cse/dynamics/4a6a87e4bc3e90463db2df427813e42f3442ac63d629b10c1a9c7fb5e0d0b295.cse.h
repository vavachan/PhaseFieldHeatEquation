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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x13 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x14 = parameters.template get<descriptors::OMEGA>();
auto x9 = x14 + V{-1};
auto x10 = V{0.0740740740740741}*x11;
auto x12 = V{1} / (x11);
auto x15 = x12*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1});
auto x16 = V{1} - x15;
auto x17 = x13 + x16;
auto x18 = -x17;
auto x19 = V{3}*x13;
auto x20 = x12*(V{3}*cell[0] + V{6}*cell[1] + V{6}*cell[2] + V{6}*cell[3] + V{3}*cell[4] + V{3}*cell[8] + V{3});
auto x21 = V{2} - x20;
auto x22 = x13*x13;
auto x23 = V{1.5}*x22;
auto x24 = x16*x16;
auto x25 = V{1.5}*x24;
auto x26 = x23 + x25;
auto x27 = x21 + x26;
auto x28 = x19 + x27;
auto x29 = x28 - V{4.5}*x18*x18;
auto x30 = -x19;
auto x31 = x13 + x15 + V{-1};
auto x32 = V{4.5}*(x31*x31);
auto x33 = -x32;
auto x34 = x27 + x30 + x33;
auto x35 = -x16;
auto x36 = x27 - V{4.5}*x35*x35;
auto x37 = V{0.444444444444444}*x11;
auto x38 = x23 + V{-1};
auto x39 = x25 + x38;
auto x40 = x37*x39;
auto x41 = V{3}*x22 + V{1};
auto x42 = -x25;
auto x43 = x19 + x42;
auto x44 = x41 + x43;
auto x45 = x10*x44;
auto x46 = x30 + x41;
auto x47 = x42 + x46;
auto x48 = x10*x47;
auto x49 = V{2.66666666666667}*cell[1] + V{1.33333333333333}*cell[2] + V{2.66666666666667}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[8] + x10*x29 + x10*x34 + V{0.148148148148148}*x11*x36 + x40 - x45 - x48 + V{0.888888888888889};
auto x50 = -x39;
auto x51 = V{0.222222222222222}*x11;
auto x52 = -x36;
auto x53 = V{0.0555555555555555}*x11;
auto x54 = x20 - x23 + x32 + x43 + V{-2};
auto x55 = -x29;
auto x56 = V{0.111111111111111}*x11;
auto x57 = V{1.33226762955019e-15}*cell[1] + V{6.66133814775094e-16}*cell[2] + V{8.88178419700125e-16}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[8] + x44*x56 + x47*x56 + V{1.11022302462516e-16};
auto x58 = x12*(x37*x50 + x51*x52 + x53*x54 + x53*x55 + x57);
auto x59 = V{1} - x58;
auto x60 = x59*x59;
auto x61 = V{1.5}*x60;
auto x62 = V{0.444444444444444}*x14*(x11*(x38 + x61) + V{1});
auto x63 = V{0.00462962962962963}*x11;
auto x64 = V{0.00462962962962963}*x11;
auto x65 = -V{4.5}*x17*x17;
auto x66 = x28 + x65;
auto x67 = x23 - V{3}*x24;
auto x68 = x21 + x67;
auto x69 = V{0.00925925925925926}*x11;
auto x70 = V{0.166666666666667}*cell[2] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[8] + V{0.0185185185185185}*x11*x68 - x44*x69 - x47*x69 + V{0.0277777777777778};
auto x71 = V{0.833333333333333}*cell[1] - V{0.166666666666667}*cell[3] - x64*x66 + x70;
auto x72 = -x34*x63 + x71;
auto x73 = V{0.0277777777777778}*x14;
auto x74 = x58 + V{-1};
auto x75 = x13 + x74;
auto x76 = V{0.166666666666667}*x11;
auto x77 = V{0.333333333333333}*x11;
auto x78 = x12*(V{3.99680288865056e-15}*cell[1] + V{1.99840144432528e-15}*cell[2] + V{2.66453525910038e-15}*cell[3] + V{4.9960036108132e-16}*cell[4] + V{4.9960036108132e-16}*cell[8] + V{1.33333333333333}*x11*x50 + V{0.666666666666667}*x11*x52 + x44*x77 + x47*x77 + x54*x76 + x55*x76 + V{3.33066907387547e-16});
auto x79 = x23 + x61 - x78 + V{2};
auto x80 = x11*(x30 + x79 - V{4.5}*x75*x75) + V{1};
auto x81 = V{0.0185185185185185}*x11;
auto x82 = V{0.00925925925925926}*x11;
auto x83 = V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[3] + x34*x82 + x66*x82 + V{-0.0555555555555555};
auto x84 = V{0.666666666666667}*cell[2] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[8] + x44*x81 + x47*x81 + x83;
auto x85 = -V{0.037037037037037}*x11*x68 + x84;
auto x86 = V{0.111111111111111}*x14;
auto x87 = x11*(x79 - V{4.5}*x74*x74) + V{1};
auto x88 = -V{0.166666666666667}*cell[1] + V{0.833333333333333}*cell[3] - x34*x64 + x70;
auto x89 = -x63*x66 + x88;
auto x90 = x13 + x59;
auto x91 = -x90;
auto x92 = x11*(x19 + x79 - V{4.5}*x91*x91) + V{1};
auto x93 = V{0.037037037037037}*x11;
auto x94 = -V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[8] - V{0.037037037037037}*x11*x68 + x83;
auto x95 = x9*(-x44*x93 + x48 + x94);
auto x96 = -x12*(-x34*x53 - x40 - x51*x68 - x53*x66 + x57) + V{1};
auto x97 = x96*x96;
auto x98 = -V{1.5}*x97;
auto x99 = x11*(x46 + x98) + V{-1};
auto x100 = V{0.0231481481481481}*x11;
auto x101 = V{0.0277777777777778}*x11;
auto x102 = x20 + V{-4};
auto x103 = x102 + x26;
auto x104 = -x75;
auto x105 = x23 + x78 + V{-4};
auto x106 = x105 + x61;
auto x107 = x9*(x45 - x47*x93 + x94);
auto x108 = x11*(x19 + x41 + x98) + V{-1};
auto x109 = V{2}*x14 + V{-2};
auto x110 = V{0.0555555555555556}*x14;
cell[0] = x49*x9 - x62;
cell[1] = -x72*x9 - x73*x80;
cell[2] = -x85*x9 - x86*x87;
cell[3] = -x73*x92 - x89*x9;
cell[4] = x86*x99 - x95;
cell[5] = -(x73*(x11*(x106 + x19 - V{4.5}*x104*x104) + V{1}) + x9*(x100*x34 - x101*(x103 + x19 + x33) + x71));
cell[6] = -x86*(x11*(x105 - V{3}*x60) + V{1}) - x9*(V{0.0740740740740741}*x11*x68 - x56*(x102 + x67) + x84);
cell[7] = -(x73*(x11*(x106 + x30 - V{4.5}*x90*x90) + V{1}) + x9*(x100*x66 - x101*(x103 + x30 + x65) + x88));
cell[8] = -x107 + x108*x86;
cell.template getFieldPointer<descriptors::VELOCITY>()[0] = -V{1}*x12*(-x107 + V{0.111111111111111}*x108*x14 - x109*x72 - x109*x85 - x109*x89 - x110*x80 - x110*x92 - V{0.222222222222222}*x14*x87 + V{0.111111111111111}*x14*x99 + x49*x9 - x62 - x95 + V{1}) + V{1};
cell.template getFieldPointer<descriptors::VELOCITY>()[1] = x13;
return { x11, x22 + V{1}*x97 };
}
};

}

}
