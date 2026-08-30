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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x12 = x11 + V{-1};
auto x13 = x9 + V{-1};
auto x14 = V{1} / (x13);
auto x15 = cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1};
auto x16 = x14*x15;
auto x17 = V{1.38777878078145e-17}*x16;
auto x18 = V{3}*x10;
auto x19 = x10*x10;
auto x20 = x9*x9;
auto x21 = V{1.5}*x20;
auto x22 = -x21;
auto x23 = V{3}*x19 + x22 + V{1};
auto x24 = x18 + x23;
auto x25 = -x18;
auto x26 = x23 + x25;
auto x27 = V{0.0740740740740741}*x14;
auto x28 = V{1.5}*x19;
auto x29 = x28 + V{-1};
auto x30 = x21 + x29;
auto x31 = V{3}*x9;
auto x32 = V{3}*x20;
auto x33 = x29 + x31 - x32;
auto x34 = V{0.0555555555555555}*x16;
auto x35 = x10 + x9;
auto x36 = V{4.5}*(x35*x35);
auto x37 = x30 + x31;
auto x38 = x18 - x36 + x37;
auto x39 = x10 - x9;
auto x40 = V{4.5}*(x39*x39);
auto x41 = -x40;
auto x42 = x25 + x37 + x41;
auto x43 = V{0.111111111111111}*x14;
auto x44 = x15*x43;
auto x45 = V{1.33226762955019e-15}*cell[1] + V{6.66133814775094e-16}*cell[2] + V{8.88178419700125e-16}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x46 = V{0.444444444444444}*x14*x15*x30 + V{0.222222222222222}*x16*x33 - x24*x44 - x26*x44 + x34*x38 + x34*x42 + x45;
auto x47 = x24*x46;
auto x48 = x26*x46;
auto x49 = V{0.0277777777777778}*x14;
auto x50 = x42*x46;
auto x51 = V{0.00462962962962963}*x14;
auto x52 = x38*x46;
auto x53 = V{0.0833333333333334}*cell[4];
auto x54 = V{0.0833333333333334}*cell[8];
auto x55 = V{0.166666666666667}*cell[2];
auto x56 = V{1.73472347597681e-18}*x16;
auto x57 = V{0.00925925925925926}*x14;
auto x58 = -V{4.62592926927149e-18}*x14*x15*x33 - V{0.0185185185185185}*x14*x33*x46 + x24*x56 + x26*x56 + x47*x57 + x48*x57 + x53 + x54 + x55 + V{0.0555555555555556};
auto x59 = x12*(-V{0.833333333333333}*cell[1] + V{0.166666666666667}*cell[3] + V{1.87928376564154e-18}*x14*x15*x38 + V{1.01192202765314e-18}*x14*x15*x42 + V{0.0231481481481481}*x14*x42*x46 - x51*x52 - x58) + V{-0.0277777777777778};
auto x60 = V{0.0185185185185185}*x14;
auto x61 = V{0.00925925925925926}*x14;
auto x62 = V{1.44560289664734e-18}*x16;
auto x63 = x33*x46;
auto x64 = V{3.85185988877447e-34}*x16;
auto x65 = -V{0.333333333333334}*cell[1] - V{0.333333333333333}*cell[3] + V{-0.0555555555555556};
auto x66 = -V{0.666666666666667}*cell[2] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[8] + x65;
auto x67 = V{0.833333333333333}*cell[3];
auto x68 = -1/x13;
auto x69 = x15*x68;
auto x70 = V{3.46944695195361e-18}*x69;
auto x71 = V{0.037037037037037}*x68;
auto x72 = V{0.111111111111111}*x69;
auto x73 = -x31;
auto x74 = V{1} - x28;
auto x75 = x32 + x74;
auto x76 = x73 + x75;
auto x77 = V{0.0555555555555555}*x69;
auto x78 = x18 + x73;
auto x79 = x22 + x74;
auto x80 = x40 + x78 + x79;
auto x81 = -x38;
auto x82 = x24*x72 + x26*x72 - V{0.444444444444444}*x30*x69 + x45 + V{0.222222222222222}*x69*x76 + x77*x80 + x77*x81;
auto x83 = x24*x82;
auto x84 = x26*x82;
auto x85 = x76*x82;
auto x86 = V{0.00925925925925926}*x68;
auto x87 = x80*x82;
auto x88 = x81*x82;
auto x89 = V{1.44560289664734e-18}*x69;
auto x90 = x80*x89 + x81*x89 + x86*x87 + x86*x88;
auto x91 = -x12*(V{0.333333333333333}*cell[2] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[8] + x24*x70 + x26*x70 + x65 - V{0.037037037037037}*x68*x85 - V{2.31296463463574e-18}*x69*x76 + x71*x83 + x71*x84 + x90) + V{0.111111111111111};
auto x92 = V{3.85185988877447e-34}*x69;
auto x93 = V{0.0185185185185185}*x68;
auto x94 = V{1.73472347597681e-18}*x69;
auto x95 = V{0.00925925925925926}*x68;
auto x0 = -x12*(-V{2.66666666666667}*cell[1] - V{1.33333333333333}*cell[2] - V{2.66666666666667}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[8] + V{3.70074341541719e-17}*x14*x15*x33 + V{1.15648231731787e-17}*x14*x15*x38 + V{1.15648231731787e-17}*x14*x15*x42 + V{0.148148148148148}*x14*x33*x46 + V{0.0740740740740741}*x14*x38*x46 + V{0.0740740740740741}*x14*x42*x46 - x17*x24 - x17*x26 - x27*x47 - x27*x48 + V{-0.444444444444444}) + V{0.444444444444444}*x14*x30*x46 + V{-0.444444444444444};
auto x1 = x49*x50 + x59;
auto x2 = x12*(V{0.0740740740740741}*x14*x63 + V{1.15648231731787e-17}*x16*x33 - x24*x64 - x26*x64 + x38*x62 + x42*x62 + x47*x60 + x48*x60 + x50*x61 + x52*x61 + x66) + x43*x63 + V{-0.111111111111111};
auto x3 = x12*(V{0.166666666666666}*cell[1] + V{1.01192202765314e-18}*x14*x15*x38 + V{1.87928376564154e-18}*x14*x15*x42 + V{0.0231481481481481}*x14*x38*x46 - x50*x51 - x58 - x67) + x49*x52 + V{-0.0277777777777778};
auto x4 = -x43*x48 - x91;
auto x5 = x46*x49*(x30 + x41 + x78) + x59;
auto x6 = x12*(x24*x92 + x26*x92 + x66 + V{0.0740740740740741}*x68*x85 + V{1.15648231731787e-17}*x69*x76 - x83*x93 - x84*x93 + x90) - x43*x46*(x31 + x75) + V{-0.111111111111111};
auto x7 = x12*(V{0.166666666666666}*cell[1] + x24*x94 + x26*x94 - x53 - x54 - x55 - x67 + V{0.0185185185185185}*x68*x85 - V{0.00462962962962963}*x68*x87 + V{0.0231481481481481}*x68*x88 + V{4.62592926927149e-18}*x69*x76 + V{1.87928376564154e-18}*x69*x80 + V{1.01192202765314e-18}*x69*x81 + x83*x95 + x84*x95 + V{-0.0555555555555556}) - x46*x49*(x18 + x31 + x36 + x79) + V{-0.0277777777777778};
auto x8 = -x43*x47 - x91;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -V{1}*x14*x46, x19 + x20 };
}
};

}

}
