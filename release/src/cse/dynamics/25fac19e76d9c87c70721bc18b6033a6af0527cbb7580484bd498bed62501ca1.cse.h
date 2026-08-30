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
struct CSE<AdvectionDiffusionBoundariesDynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 2, 1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x14 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x8 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x7 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x11 = x14 + V{-1};
auto x12 = V{0.5}*cell[1];
auto x13 = V{0.5}*cell[4];
auto x15 = V{4}*x7;
auto x16 = x15 + V{1};
auto x17 = V{0.0625}*x10;
auto x18 = x15 + V{-1};
auto x19 = x11*(x12 - x13 + x16*x17 + x17*x18);
auto x20 = V{0.125}*x10;
auto x21 = V{0.5}*cell[2];
auto x22 = V{0.5}*cell[5];
auto x23 = V{4}*x8;
auto x24 = x23 + V{1};
auto x25 = x23 + V{-1};
auto x26 = x11*(x17*x24 + x17*x25 + x21 - x22);
auto x27 = V{4}*x9;
auto x28 = x27 + V{-1};
auto x29 = x27 + V{1};
auto x30 = V{0.5}*cell[0] + V{1}*cell[6] - V{0.5}*x10 + x12 + x13 - x17*x29 + x21 + x22 + V{0.5};
cell[0] = V{0.25}*x10 + V{-0.25};
cell[1] = -x18*x20 - x19 + V{-0.125};
cell[2] = -x20*x25 - x26 + V{-0.125};
cell[3] = x11*(-x17*x28 + x30) - x20*x28 + V{-0.125};
cell[4] = x16*x20 + x19 + V{-0.125};
cell[5] = x20*x24 + x26 + V{-0.125};
cell[6] = V{0.125}*x10*x29 - x11*(-x17*x28 + x30) + V{-0.125};
return { x10, x7*x7 + x8*x8 + x9*x9 };
}
};

}

}
