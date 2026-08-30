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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x12 = x11 + V{-1};
auto x13 = x10 + V{-1};
auto x14 = V{1} / (x13);
auto x15 = cell[0] + cell[2] + V{2}*cell[3] + V{2}*cell[4] + V{2}*cell[5] + cell[6] + V{1};
auto x16 = x14*x15;
auto x17 = V{3}*x9;
auto x18 = x9*x9;
auto x19 = V{3}*x18;
auto x20 = x10*x10;
auto x21 = V{1.5}*x20;
auto x22 = -x21;
auto x23 = x19 + x22 + V{1};
auto x24 = x17 + x23;
auto x25 = V{0.0555555555555555}*x16;
auto x26 = x10 - x9;
auto x27 = -x26;
auto x28 = V{4.5}*(x27*x27);
auto x29 = -x17;
auto x30 = V{3}*x10;
auto x31 = V{1.5}*x18;
auto x32 = x21 + V{-1};
auto x33 = x31 + x32;
auto x34 = x30 + x33;
auto x35 = x29 + x34;
auto x36 = -x28 + x35;
auto x37 = V{0.111111111111111}*x14;
auto x38 = x15*x37;
auto x39 = x17 - x19 + x32;
auto x40 = x10 + x9;
auto x41 = V{4.5}*(x40*x40);
auto x42 = x17 + x34 - x41;
auto x43 = -x30;
auto x44 = V{1} - x31;
auto x45 = V{3}*x20 + x44;
auto x46 = x43 + x45;
auto x47 = V{1.66533453693773e-16}*cell[2] + V{1.33226762955019e-15}*cell[3] + V{6.66133814775094e-16}*cell[4] + V{8.88178419700125e-16}*cell[5] + V{1.66533453693773e-16}*cell[6] + V{1.11022302462516e-16};
auto x48 = V{0.444444444444444}*x14*x15*x33 - V{0.222222222222222}*x16*x46 - x24*x38 + x25*x42 + x38*x39 + x47;
auto x49 = x25*x36 + x48;
auto x50 = x24*x49;
auto x51 = x46*x49;
auto x52 = V{0.0277777777777778}*x14;
auto x53 = -V{4.5}*x26*x26;
auto x54 = x17 + x43;
auto x55 = x35 + x53;
auto x56 = x25*x55 + x48;
auto x57 = V{0.00462962962962963}*x14;
auto x58 = x42*x49;
auto x59 = V{0.0833333333333334}*cell[2];
auto x60 = V{0.0833333333333334}*cell[6];
auto x61 = V{0.166666666666667}*cell[4];
auto x62 = -V{1.73472347597681e-18}*x14*x15*x39 - V{0.00925925925925926}*x14*x39*x49 + V{0.00925925925925926}*x14*x50 + V{0.0185185185185185}*x14*x51 + V{1.73472347597681e-18}*x16*x24 + V{4.62592926927149e-18}*x16*x46 + x59 + x60 + x61 + V{0.0555555555555556};
auto x63 = x12*(V{0.166666666666666}*cell[3] - V{0.833333333333333}*cell[5] + V{1.01192202765314e-18}*x14*x15*x36 + V{1.87928376564154e-18}*x14*x15*x42 + V{0.0231481481481481}*x14*x36*x49 - x57*x58 - x62) + V{-0.0277777777777778};
auto x64 = V{0.00925925925925926}*x14;
auto x65 = x36*x49;
auto x66 = V{1.44560289664734e-18}*x16;
auto x67 = V{3.46944695195361e-18}*x16;
auto x68 = V{0.037037037037037}*x14;
auto x69 = -V{0.333333333333334}*cell[3] - V{0.333333333333333}*cell[5] + V{-0.0555555555555556};
auto x70 = -V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[6] + x69;
auto x71 = V{0.833333333333333}*cell[3];
auto x72 = -1/x13;
auto x73 = x15*x72;
auto x74 = V{3.85185988877447e-34}*x73;
auto x75 = x23 + x29;
auto x76 = V{0.111111111111111}*x73;
auto x77 = V{0.0555555555555555}*x73;
auto x78 = x22 + x44;
auto x79 = x28 + x54 + x78;
auto x80 = -x42;
auto x81 = x24*x76 - V{0.444444444444444}*x33*x73 + V{0.222222222222222}*x46*x73 + x47 + x75*x76 + x77*x79 + x77*x80;
auto x82 = x46*x81;
auto x83 = V{0.0185185185185185}*x72;
auto x84 = x24*x81;
auto x85 = x75*x81;
auto x86 = V{0.00925925925925926}*x72;
auto x87 = x79*x81;
auto x88 = x80*x81;
auto x89 = V{1.44560289664734e-18}*x73;
auto x90 = x79*x89 + x80*x89 + x86*x87 + x86*x88;
auto x91 = -x12*(V{0.166666666666667}*cell[2] - V{0.666666666666667}*cell[4] + V{0.166666666666667}*cell[6] + x24*x74 + V{1.15648231731787e-17}*x46*x73 + x69 + V{0.0740740740740741}*x72*x82 + x74*x75 - x83*x84 - x83*x85 + x90) + V{0.111111111111111};
auto x92 = V{3.46944695195361e-18}*x73;
auto x93 = V{0.037037037037037}*x72;
auto x94 = V{1.73472347597681e-18}*x73;
auto x95 = V{0.00925925925925926}*x72;
auto x0 = -x12*(-V{0.666666666666667}*cell[2] - V{2.66666666666667}*cell[3] - V{1.33333333333333}*cell[4] - V{2.66666666666667}*cell[5] - V{0.666666666666667}*cell[6] + V{1.15648231731787e-17}*x14*x15*x36 + V{1.38777878078145e-17}*x14*x15*x39 + V{1.15648231731787e-17}*x14*x15*x42 + V{0.0740740740740741}*x14*x36*x49 + V{0.0740740740740741}*x14*x39*x49 + V{0.0740740740740741}*x14*x42*x49 - V{0.0740740740740741}*x14*x50 - V{0.148148148148148}*x14*x51 - V{1.38777878078145e-17}*x16*x24 - V{3.70074341541719e-17}*x16*x46 + V{-0.444444444444444}) + V{0.444444444444444}*x14*x33*x49 + V{-0.444444444444444};
auto x1 = x52*x56*(x33 + x53 + x54) + x63;
auto x2 = x12*(V{0.037037037037037}*x14*x51 + V{2.31296463463574e-18}*x16*x46 - x24*x67 + x36*x66 + x39*x49*x68 + x39*x67 + x42*x66 - x50*x68 + x58*x64 + x64*x65 + x70) + x37*x39*x56 + V{-0.111111111111111};
auto x3 = x12*(V{0.166666666666667}*cell[5] + V{1.87928376564154e-18}*x14*x15*x36 + V{1.01192202765314e-18}*x14*x15*x42 + V{0.0231481481481481}*x14*x42*x49 - x57*x65 - x62 - x71) + x42*x52*x56 + V{-0.0277777777777778};
auto x4 = -x37*x51 - x91;
auto x5 = x52*x55*x56 + x63;
auto x6 = x12*(x24*x92 - V{2.31296463463574e-18}*x46*x73 + x70 - V{0.037037037037037}*x72*x82 + x75*x92 + x84*x93 + x85*x93 + x90) - x37*x50 + V{-0.111111111111111};
auto x7 = x12*(V{0.166666666666667}*cell[5] + x24*x94 + V{4.62592926927149e-18}*x46*x73 - x59 - x60 - x61 - x71 + V{0.0185185185185185}*x72*x82 - V{0.00462962962962963}*x72*x87 + V{0.0231481481481481}*x72*x88 + V{1.87928376564154e-18}*x73*x79 + V{1.01192202765314e-18}*x73*x80 + x75*x94 + x84*x95 + x85*x95 + V{-0.0555555555555556}) - x49*x52*(x17 + x30 + x41 + x78) + V{-0.0277777777777778};
auto x8 = -x37*x49*(x30 + x45) - x91;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -V{1}*x14*x56, x18 + x20 };
}
};

}

}
