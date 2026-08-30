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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x12 = x11 + V{-1};
auto x13 = V{1} / (x10 + V{1});
auto x14 = cell[0] + V{2}*cell[1] + cell[2] + cell[6] + V{2}*cell[7] + V{2}*cell[8] + V{1};
auto x15 = x13*x14;
auto x16 = V{3}*x9;
auto x17 = x9*x9;
auto x18 = V{3}*x17;
auto x19 = x10*x10;
auto x20 = V{1.5}*x19;
auto x21 = x20 + V{-1};
auto x22 = x16 - x18 + x21;
auto x23 = V{3}*x10;
auto x24 = -x23;
auto x25 = x10 - x9;
auto x26 = V{4.5}*(x25*x25);
auto x27 = V{1.5}*x17;
auto x28 = x21 + x27;
auto x29 = x16 + x28;
auto x30 = x24 - x26 + x29;
auto x31 = V{0.0740740740740741}*x13;
auto x32 = V{0.111111111111111}*x13;
auto x33 = x14*x32;
auto x34 = V{0.444444444444444}*x13;
auto x35 = x14*x34;
auto x36 = V{0.0555555555555555}*x15;
auto x37 = V{1} - x20;
auto x38 = x18 + x37;
auto x39 = x16 + x38;
auto x40 = -x27;
auto x41 = x23 + x40;
auto x42 = V{3}*x19 + V{1};
auto x43 = x41 + x42;
auto x44 = x10 + x9;
auto x45 = V{4.5}*(x44*x44);
auto x46 = x37 + x41;
auto x47 = x16 + x45 + x46;
auto x48 = V{8.88178419700125e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.66533453693773e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + V{0.222222222222222}*x15*x43 + x33*x39 + x36*x47 + V{1.11022302462516e-16};
auto x49 = -x22*x33 - x28*x35 - x30*x36 + x48;
auto x50 = x22*x49;
auto x51 = x30*x49;
auto x52 = V{0.0277777777777778}*x13;
auto x53 = V{1.73472347597681e-18}*x15;
auto x54 = V{0.0231481481481481}*x13;
auto x55 = -x16;
auto x56 = x26 + x46 + x55;
auto x57 = x38 + x55;
auto x58 = -x28*x35 + x33*x57 + x36*x56 + x48;
auto x59 = x56*x58;
auto x60 = V{8.67361737988404e-19}*x15;
auto x61 = V{0.00462962962962963}*x13;
auto x62 = x47*x58;
auto x63 = V{0.00925925925925926}*x13;
auto x64 = x39*x58;
auto x65 = x57*x58;
auto x66 = V{0.0185185185185185}*x13;
auto x67 = x43*x58;
auto x68 = -V{0.0833333333333334}*cell[2] - V{0.0833333333333334}*cell[6] - V{0.166666666666667}*cell[8] + V{4.62592926927149e-18}*x15*x43 + x39*x53 + V{-0.0555555555555556};
auto x69 = x53*x57 + x63*x64 + x63*x65 + x66*x67 + x68;
auto x70 = -x12*(-V{0.833333333333333}*cell[1] + V{0.166666666666667}*cell[7] + x47*x53 + x54*x59 + x56*x60 - x61*x62 + x69) + V{0.0277777777777778};
auto x71 = V{0.333333333333333}*cell[7];
auto x72 = -x71;
auto x73 = V{0.333333333333334}*cell[1];
auto x74 = V{0.333333333333333}*cell[2];
auto x75 = V{0.333333333333333}*cell[6];
auto x76 = V{1.73472347597681e-18}*x15;
auto x77 = V{0.00925925925925926}*x13;
auto x78 = V{3.46944695195361e-18}*x15;
auto x79 = V{0.037037037037037}*x13;
auto x80 = V{0.037037037037037}*x13;
auto x81 = V{2.31296463463574e-18}*x15*x43;
auto x82 = V{0.166666666666667}*cell[1] - V{0.833333333333333}*cell[7] + x47*x60;
auto x83 = V{3.85185988877447e-34}*x15;
auto x84 = V{0.0185185185185185}*x13;
auto x85 = V{8.67361737988404e-19}*x15;
auto x86 = x47*x49;
auto x87 = x43*x49;
auto x88 = x39*x49;
auto x89 = x51*x77;
auto x90 = x12*(-V{0.333333333333333}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] - V{0.666666666666667}*cell[8] + V{0.0740740740740741}*x13*x87 + V{1.15648231731787e-17}*x15*x43 - x22*x83 - x30*x85 + x39*x83 + x47*x85 + x50*x84 + x72 + x77*x86 - x84*x88 - x89 + V{-0.0555555555555556}) + V{-0.111111111111111};
auto x91 = -x25;
auto x0 = -x12*(-V{2.66666666666667}*cell[1] - V{0.666666666666667}*cell[2] - V{0.666666666666667}*cell[6] - V{2.66666666666667}*cell[7] - V{1.33333333333333}*cell[8] + V{1.38777878078145e-17}*x13*x14*x39 + V{3.70074341541719e-17}*x13*x14*x43 + V{1.04083408558608e-17}*x13*x14*x47 + V{0.0740740740740741}*x13*x39*x49 + V{0.148148148148148}*x13*x43*x49 + V{0.0740740740740741}*x13*x47*x49 - V{1.38777878078145e-17}*x15*x22 - V{1.04083408558608e-17}*x15*x30 - x31*x50 - x31*x51 + V{-0.444444444444445}) - x28*x34*x49 + V{-0.444444444444444};
auto x1 = -x51*x52 - x70;
auto x2 = x12*(V{0.333333333333333}*cell[8] + x39*x78 + x47*x76 + x56*x76 + x57*x78 + x59*x77 + x62*x77 + x64*x79 + x65*x79 - x67*x80 + x72 - x73 - x74 - x75 - x81 + V{-0.0555555555555556}) - x32*x50 + V{-0.111111111111111};
auto x3 = x12*(x53*x56 + x54*x62 - x59*x61 + x69 + x82) - x49*x52*(x23 + x29 - x45) + V{-0.0277777777777778};
auto x4 = x32*x49*(x24 + x40 + x42) + x90;
auto x5 = -(x49*x52*(x23 + x28 + x55 - V{4.5}*x91*x91) + x70);
auto x6 = x12*(V{0.333333333333333}*cell[8] + V{3.46944695195361e-18}*x13*x14*x39 + V{1.73472347597681e-18}*x13*x14*x47 + V{0.037037037037037}*x13*x39*x49 + V{0.00925925925925926}*x13*x47*x49 - x22*x78 - x30*x76 - x50*x79 - x71 - x73 - x74 - x75 - x80*x87 - x81 - x89 + V{-0.0555555555555556}) + x32*x88 + V{-0.111111111111111};
auto x7 = x12*(-x22*x53 - x30*x53 - x50*x63 + x51*x61 + x54*x86 + x63*x88 + x66*x87 + x68 + x82) + x52*x86 + V{-0.0277777777777778};
auto x8 = x32*x87 + x90;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x13*x49, x17 + x19 };
}
};

}

}
