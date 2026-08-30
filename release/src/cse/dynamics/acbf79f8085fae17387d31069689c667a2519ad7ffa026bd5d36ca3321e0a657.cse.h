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
struct CSE<AdvectionDiffusionCornerDynamics2D<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = parameters.template get<descriptors::OMEGA>();
auto x6 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x5 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x7 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x8 = V{3}*x5;
auto x9 = x8 + V{1};
auto x11 = V{0.0555555555555556}*x7;
auto x12 = V{3}*x6;
auto x13 = x12 + V{1};
auto x14 = x8 + V{-1};
auto x15 = x12 + V{-1};
auto x16 = V{0.0277777777777778}*x7;
auto x17 = V{0.166666666666667}*cell[0] + x13*x16 + x16*x9 + V{0.0555555555555556};
auto x18 = -x14*x16 - x15*x16 + x17;
auto x19 = x14*x18;
auto x20 = x10 + V{-1};
auto x21 = -x14;
auto x22 = -V{0.5}*x15*x16 + V{0.5}*x16*x21 + V{0.5}*x17;
auto x23 = V{0.0833333333333333}*x7;
auto x24 = V{1}*cell[3] - x23*x9 + V{0.166666666666667};
auto x25 = x13*x18;
auto x26 = x15*x18;
auto x27 = x20*(V{1}*cell[2] - x13*x23 + x15*x23 + V{0.5}*x25 + V{0.5}*x26 + V{0.166666666666667});
auto x28 = V{0.166666666666667}*x7;
auto x0 = V{0.333333333333333}*cell[0] + x11*x13 - x11*x14 - x11*x15 + x11*x9 + V{-0.222222222222222};
auto x1 = -x19 + x20*(x21*x22 - x21*x23 - x22*x9 + x24) + V{-0.166666666666667};
auto x2 = -x26 - x27 + V{-0.166666666666667};
auto x3 = x18*x9 - x20*(x14*x23 - V{0.5}*x18*x9 - V{0.5}*x19 + x24) + V{-0.166666666666667};
auto x4 = x25 + x27 + V{-0.166666666666667};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
return { cell[0] + x13*x28 - x14*x28 - x15*x28 + x28*x9 + V{0.333333333333333}, x5*x5 + x6*x6 };
}
};

}

}
