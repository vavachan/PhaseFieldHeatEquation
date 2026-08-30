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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<0, -1>, momenta::NoStress, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x12 = parameters.template get<descriptors::OMEGA>();
auto x11 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x10 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x9 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x13 = V{3}*x9;
auto x14 = V{3}*x10;
auto x15 = x14 + V{1};
auto x16 = x13 + x15;
auto x17 = V{0.0123456790123457}*x11;
auto x18 = -x13 + x15;
auto x19 = x13 + V{1};
auto x20 = -x14 + x19;
auto x21 = V{0.0493827160493827}*x11;
auto x22 = x14 + V{-1};
auto x23 = x13 + x22;
auto x24 = x13 + V{-1};
auto x25 = x12 + V{-1};
auto x26 = V{0.000771604938271605}*x11;
auto x27 = V{0.00308641975308642}*x11;
auto x28 = x11*x19;
auto x29 = V{3.08395284618099e-18}*cell[1] - V{7.70988211545248e-19}*cell[2] - V{3.08395284618099e-18}*cell[3] + x15*x27 + x16*x26 + x17 + x18*x26 + x20*x26 + V{0.00308641975308642}*x28;
auto x30 = -x22*x27 - x23*x26 - x24*x27 + x29;
auto x31 = x20*x30;
auto x32 = V{1.2335811384724e-17}*cell[1] - V{3.08395284618099e-18}*cell[2] - V{1.2335811384724e-17}*cell[3] + x15*x17 + x16*x27 + x17*x19 + x18*x27 + x20*x27 + x21;
auto x33 = -x17*x22 - x17*x24 - x23*x27 + x32;
auto x34 = x15*x33;
auto x35 = V{0.0833333333333333}*x34;
auto x36 = x22*x33;
auto x37 = V{0.0833333333333333}*x36;
auto x38 = V{0.0833333333333333}*cell[8];
auto x39 = V{0.0833333333333333}*cell[4];
auto x40 = x38 - x39;
auto x41 = V{0.166666666666667}*cell[2];
auto x42 = V{0.00925925925925926}*x11;
auto x43 = x19*x33;
auto x44 = x24*x33;
auto x45 = -V{0.00925925925925926}*x11*x19 + x24*x42 + x41 + V{0.0833333333333333}*x43 + V{0.0833333333333333}*x44 + V{0.0277777777777778};
auto x46 = x25*(-V{0.333333333333333}*cell[1] + V{1.85037170770859e-17}*cell[3] + V{0.00462962962962963}*x11*x18 + V{0.00462962962962963}*x11*x20 + V{0.166666666666667}*x18*x30 - V{0.166666666666667}*x31 + x35 + x37 - x40 - x45);
auto x47 = x18*x30;
auto x48 = V{0.666666666666667}*cell[2];
auto x49 = V{0.666666666666667}*cell[3];
auto x50 = -x24;
auto x51 = V{0.037037037037037}*x11;
auto x52 = -x23;
auto x53 = -x22;
auto x54 = x17*x50 + x17*x53 + x27*x52 + x32;
auto x55 = V{0.333333333333333}*x54;
auto x56 = x26*x52 + x27*x50 + x27*x53 + x29;
auto x57 = V{0.333333333333333}*x56;
auto x58 = x18*x57;
auto x59 = x20*x57;
auto x60 = x16*x42 - x16*x57 + x42*x52 + x52*x57;
auto x61 = V{0.666666666666667}*cell[1];
auto x62 = x18*x42;
auto x63 = x20*x42;
auto x64 = -x61 + x62 + x63;
auto x65 = x23*x30;
auto x66 = V{0.333333333333333}*cell[3];
auto x67 = V{0.00462962962962963}*x11;
auto x68 = V{0.0833333333333333}*x54;
auto x69 = V{0.166666666666667}*x56;
auto x70 = V{0.333333333333333}*cell[4];
auto x71 = V{0.666666666666667}*cell[3];
auto x72 = x61 - x62 - x63;
auto x73 = V{0.333333333333333}*x31;
auto x74 = V{0.333333333333333}*x47;
auto x75 = x16*x30;
auto x76 = -V{0.00925925925925926}*x11*x16 + x23*x42 + V{0.333333333333333}*x65 + V{0.333333333333333}*x75;
auto x77 = V{0.0277777777777778}*x11;
auto x78 = V{0.111111111111111}*x11;
auto x0 = V{4.93432455388958e-17}*cell[1] - V{1.2335811384724e-17}*cell[2] - V{4.93432455388958e-17}*cell[3] + V{0.197530864197531}*x11 + x15*x21 + x16*x17 + x17*x18 + x17*x20 - x17*x23 + x19*x21 - x21*x22 - x21*x24 + V{-0.444444444444444};
auto x1 = x46 + x47 + V{-0.0277777777777778};
auto x2 = x25*(-x19*x55 + V{0.037037037037037}*x28 - x48 - x49 + x50*x51 + x50*x55 + x58 - x59 + x60 + x64 + V{-0.111111111111111}) - x44 + V{-0.111111111111111};
auto x3 = x25*(-x15*x68 + x16*x67 - x16*x69 - x19*x68 + V{0.00925925925925926}*x28 + x40 - x41 + x42*x50 + x50*x68 + x52*x67 + x52*x69 + x53*x68 - x66 + V{-0.0277777777777778}) - x65 + V{-0.0277777777777778};
auto x4 = x25*(V{0.333333333333333}*cell[8] - x15*x55 + x53*x55 - x58 + x59 + x60 - x70 - x71 + x72) - x36 + V{-0.111111111111111};
auto x5 = x20*x30 - x46 + V{-0.0277777777777778};
auto x6 = x19*x33 - x25*(V{0.037037037037037}*x11*x19 - x24*x51 - V{0.333333333333333}*x43 - V{0.333333333333333}*x44 - x48 - x49 - x72 - x73 + x74 - x76 + V{-0.111111111111111}) + V{-0.111111111111111};
auto x7 = x16*x30 - x25*(V{0.00462962962962963}*x11*x16 - x23*x67 - x35 - x37 + x38 - x39 - x45 - V{0.166666666666667}*x65 - x66 - V{0.166666666666667}*x75) + V{-0.0277777777777778};
auto x8 = x15*x33 - x25*(V{0.333333333333333}*cell[8] - V{0.333333333333333}*x34 - V{0.333333333333333}*x36 - x64 - x70 - x71 + x73 - x74 - x76) + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1.11022302462516e-16}*cell[1] - V{2.77555756156289e-17}*cell[2] - V{1.11022302462516e-16}*cell[3] + V{0.444444444444444}*x11 + x15*x78 + x16*x77 + x18*x77 + x20*x77 - x22*x78 - x23*x77 - x24*x78 + V{0.111111111111111}*x28, x10*x10 + x9*x9 };
}
};

}

}
