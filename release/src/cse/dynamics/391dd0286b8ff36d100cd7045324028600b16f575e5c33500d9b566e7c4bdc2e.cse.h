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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::CahnHilliardZerothOrder, collision::BGK, forcing::WellBalancedCahnHilliard>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x13 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0);
auto x10 = cell.template getFieldComponent<descriptors::SOURCE>(0);
auto x11 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x12 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = x13 + V{-1};
auto x15 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x16 = x9 + V{-1};
auto x17 = V{0.0277777777777778}*x10 + V{0.0277777777777778}*x13*x16;
auto x18 = V{0.0555555555555556}*x10 - V{0.111111111111111}*x13*x16;
auto x0 = -cell[0]*x14 - V{0.888888888888889}*x10 + x13*(x15 - V{0.555555555555556}*x9 + V{0.555555555555556});
auto x1 = -cell[1]*x14 + x17;
auto x2 = -cell[2]*x14 - x18;
auto x3 = -cell[3]*x14 + x17;
auto x4 = -cell[4]*x14 - x18;
auto x5 = -cell[5]*x14 + x17;
auto x6 = -cell[6]*x14 - x18;
auto x7 = -cell[7]*x14 + x17;
auto x8 = -cell[8]*x14 - x18;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x15 + V{1}, x11*x11 + x12*x12 };
}
};

}

}
