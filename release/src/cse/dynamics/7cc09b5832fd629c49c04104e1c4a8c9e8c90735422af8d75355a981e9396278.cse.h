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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x7 = cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(0);
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x6 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x5 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x8 = cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(1);
auto x10 = V{1} / (x6 + V{1});
auto x11 = x10*x8;
auto x12 = V{3}*x5;
auto x13 = x12 + V{1};
auto x14 = cell[0] + cell[1] + cell[3] + V{2}*cell[4] - x8 + V{1};
auto x15 = -x14;
auto x16 = x10*x15;
auto x17 = V{0.0555555555555556}*x16;
auto x18 = V{3}*x6;
auto x19 = x18 + V{1};
auto x20 = x12 + V{-1};
auto x21 = x18 + V{-1};
auto x22 = V{0.0277777777777778}*x16;
auto x23 = V{0.0277777777777778}*x13;
auto x24 = x10*(V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[3] + V{0.111111111111111}*cell[4] + V{0.0555555555555556}) - V{0.0555555555555556}*x11 + V{1.85037170770859e-17};
auto x25 = -x16*x23 - x19*x22 + x20*x22 + x21*x22 + x24;
auto x26 = x20*x25;
auto x27 = x9 + V{-1};
auto x28 = -x20;
auto x29 = x10*x14;
auto x30 = V{0.0277777777777778}*x29;
auto x31 = -x21;
auto x32 = V{0.5}*x19*x30 + V{0.5}*x23*x29 + V{0.5}*x24 + V{0.5}*x28*x30 + V{0.5}*x30*x31;
auto x33 = V{0.0833333333333333}*x29;
auto x34 = x21*x25;
auto x35 = V{0.0833333333333333}*x16;
auto x36 = V{0.166666666666667}*x16;
auto x0 = V{0.0555555555555556}*x10*x15*x20 + V{0.0555555555555556}*x10*x15*x21 + x10*(V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[3] + V{0.222222222222222}*cell[4] + V{0.111111111111111}) - V{0.111111111111111}*x11 - x13*x17 - x17*x19 + V{-0.333333333333333};
auto x1 = -x26 + x27*(-x13*x32 + x13*x33 + x28*x32 - x28*x33 + V{0.5}*x7) + V{-0.166666666666667};
auto x2 = x27*(-x19*x32 + x19*x33 + x31*x32 - x31*x33 + V{0.5}*x8) - x34 + V{-0.166666666666667};
auto x3 = x13*x25 - x27*(-V{0.5}*x13*x25 - x13*x35 - x20*x35 - V{0.5}*x26 + V{0.5}*x7) + V{-0.166666666666667};
auto x4 = x19*x25 - x27*(-V{0.5}*x19*x25 - x19*x35 - x21*x35 - V{0.5}*x34 + V{0.5}*x8) + V{-0.166666666666667};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
return { x10*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.333333333333333}) - V{0.333333333333333}*x11 - x13*x36 - x19*x36 + x20*x36 + x21*x36 + V{1.11022302462516e-16}, x5*x5 + x6*x6 };
}
};

}

}
