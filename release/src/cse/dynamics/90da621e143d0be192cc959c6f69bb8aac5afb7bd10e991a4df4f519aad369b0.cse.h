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
struct CSE<AdvectionDiffusionBoundariesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, 1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x20 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x19 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x23 = V{3}*x20;
auto x24 = x23 + V{-1};
auto x25 = x23 + V{1};
auto x26 = V{3}*x19;
auto x27 = V{0.0277777777777778}*x22;
auto x28 = V{3}*x21;
auto x29 = V{1} - x23;
auto x0 = cell[0];
auto x1 = cell[1];
auto x2 = -cell[11] - V{0.0555555555555556}*x22*x24 + V{0.0555555555555556}*x22*x25 + V{-0.111111111111111};
auto x3 = cell[3];
auto x4 = -cell[13] + V{0.0277777777777778}*x22*(x25 + x26) - x27*(x24 + x26) + V{-0.0555555555555556};
auto x5 = cell[5];
auto x6 = cell[6];
auto x7 = cell[7];
auto x8 = -cell[17] + V{0.0277777777777778}*x22*(x25 + x28) - x27*(x24 + x28) + V{-0.0555555555555556};
auto x9 = -cell[18] + V{0.0277777777777778}*x22*(x25 - x28) + V{0.0277777777777778}*x22*(x28 + x29) + V{-0.0555555555555556};
auto x10 = cell[10];
auto x11 = cell[11];
auto x12 = cell[12];
auto x13 = cell[13];
auto x14 = -cell[5] + V{0.0277777777777778}*x22*(x25 - x26) + V{0.0277777777777778}*x22*(x26 + x29) + V{-0.0555555555555556};
auto x15 = cell[15];
auto x16 = cell[16];
auto x17 = cell[17];
auto x18 = cell[18];
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { V{-1}, V{-1} };
}
};

}

}
