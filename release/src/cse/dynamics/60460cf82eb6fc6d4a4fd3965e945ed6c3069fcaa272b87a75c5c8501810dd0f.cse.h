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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x12 = x11 + V{-1};
auto x13 = V{1} / (x9 + V{1});
auto x14 = cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1};
auto x15 = x13*x14;
auto x16 = V{1.38777878078145e-17}*x15;
auto x17 = V{3}*x10;
auto x18 = x10*x10;
auto x19 = x9*x9;
auto x20 = V{1.5}*x19;
auto x21 = V{1} - x20;
auto x22 = V{3}*x18 + x21;
auto x23 = x17 + x22;
auto x24 = -x17;
auto x25 = x22 + x24;
auto x26 = V{0.0740740740740741}*x13;
auto x27 = V{0.0555555555555555}*x15;
auto x28 = x10 - x9;
auto x29 = -x28;
auto x30 = V{4.5}*(x29*x29);
auto x31 = V{3}*x9;
auto x32 = V{1.5}*x18;
auto x33 = x32 + V{-1};
auto x34 = x20 + x33;
auto x35 = x17 + x34;
auto x36 = -x31 + x35;
auto x37 = -x30 + x36;
auto x38 = V{0.444444444444444}*x13;
auto x39 = x14*x38;
auto x40 = V{0.111111111111111}*x13;
auto x41 = x14*x40;
auto x42 = V{3}*x19;
auto x43 = x31 - x32;
auto x44 = x42 + x43 + V{1};
auto x45 = x10 + x9;
auto x46 = V{4.5}*(x45*x45);
auto x47 = x21 + x43;
auto x48 = x17 + x46 + x47;
auto x49 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{0.222222222222222}*x15*x44 + x23*x41 + x25*x41 + x27*x48 + V{1.11022302462516e-16};
auto x50 = -x34*x39 + x49;
auto x51 = -x27*x37 + x50;
auto x52 = x23*x51;
auto x53 = x25*x51;
auto x54 = x48*x51;
auto x55 = x44*x51;
auto x56 = V{1.15648231731787e-17}*x15;
auto x57 = x37*x51;
auto x58 = V{0.0277777777777778}*x13;
auto x59 = -V{4.5}*x28*x28;
auto x60 = V{1.01192202765314e-18}*x15;
auto x61 = x24 + x30 + x47;
auto x62 = V{1.87928376564154e-18}*x15;
auto x63 = V{0.0231481481481481}*x13;
auto x64 = x27*x61 - x34*x39 + x49;
auto x65 = x61*x64;
auto x66 = V{0.00462962962962963}*x13;
auto x67 = x48*x64;
auto x68 = V{0.00925925925925926}*x13;
auto x69 = x23*x64;
auto x70 = x25*x64;
auto x71 = V{0.0185185185185185}*x13;
auto x72 = x44*x64;
auto x73 = V{1.73472347597681e-18}*x15;
auto x74 = -V{0.0833333333333334}*cell[4] - V{0.166666666666667}*cell[6] - V{0.0833333333333334}*cell[8] + V{4.62592926927149e-18}*x15*x44 + x23*x73 + x25*x73 + V{-0.0555555555555556};
auto x75 = x68*x69 + x68*x70 + x71*x72 + x74;
auto x76 = -x12*(-V{0.833333333333333}*cell[5] + V{0.166666666666667}*cell[7] + x48*x62 + x60*x61 + x63*x65 - x66*x67 + x75) + V{0.0277777777777778};
auto x77 = V{0.00925925925925926}*x13;
auto x78 = V{1.44560289664734e-18}*x15;
auto x79 = V{0.0740740740740741}*x13;
auto x80 = V{0.0185185185185185}*x13;
auto x81 = V{0.333333333333333}*cell[7];
auto x82 = V{0.333333333333334}*cell[5];
auto x83 = V{3.85185988877447e-34}*x15;
auto x84 = V{0.166666666666667}*cell[4] - V{0.666666666666667}*cell[6] + V{0.166666666666667}*cell[8] + V{1.15648231731787e-17}*x15*x44 + x23*x83 + x25*x83 + x48*x78 - x81 - x82 + V{-0.0555555555555556};
auto x85 = V{0.166666666666666}*cell[5] - V{0.833333333333333}*cell[7] + x48*x60;
auto x86 = -x27*(x36 + x59) + x50;
auto x87 = x57*x77;
auto x88 = x37*x78;
auto x89 = x12*(-V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[8] + V{3.46944695195361e-18}*x13*x14*x23 + V{3.46944695195361e-18}*x13*x14*x25 + V{1.44560289664734e-18}*x13*x14*x48 + V{0.037037037037037}*x13*x23*x51 + V{0.037037037037037}*x13*x25*x51 + V{0.00925925925925926}*x13*x48*x51 - V{0.037037037037037}*x13*x55 - V{2.31296463463574e-18}*x15*x44 - x81 - x82 - x87 - x88 + V{-0.0555555555555556}) + V{-0.111111111111111};
auto x0 = -x12*(-V{0.666666666666667}*cell[4] - V{2.66666666666667}*cell[5] - V{1.33333333333333}*cell[6] - V{2.66666666666667}*cell[7] - V{0.666666666666667}*cell[8] + V{0.148148148148148}*x13*x55 + V{3.70074341541719e-17}*x15*x44 + x16*x23 + x16*x25 + x26*x52 + x26*x53 + x26*x54 - x26*x57 - x37*x56 + x48*x56 + V{-0.444444444444444}) - x34*x38*x51 + V{-0.444444444444444};
auto x1 = -x51*x58*(x24 + x31 + x34 + x59) - x76;
auto x2 = x12*(x61*x78 + x65*x77 + x67*x77 - x69*x80 - x70*x80 + x72*x79 + x84) - x40*x51*(x31 + x33 - x42) + V{-0.111111111111111};
auto x3 = x12*(x61*x62 + x63*x67 - x65*x66 + x75 + x85) - x51*x58*(x31 + x35 - x46) + V{-0.0277777777777778};
auto x4 = x25*x40*x86 + x89;
auto x5 = -x57*x58 - x76;
auto x6 = x12*(-x52*x80 - x53*x80 + x54*x77 + x55*x79 + x84 - x87 - x88) + x40*x44*x86 + V{-0.111111111111111};
auto x7 = x12*(-x37*x62 + x52*x68 + x53*x68 + x54*x63 + x55*x71 + x57*x66 + x74 + x85) + x48*x58*x86 + V{-0.0277777777777778};
auto x8 = x23*x40*x86 + x89;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x13*x86, x18 + x19 };
}
};

}

}
