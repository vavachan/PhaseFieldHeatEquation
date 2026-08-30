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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::CahnHilliardZerothOrder, collision::BGK, forcing::WellBalancedCahnHilliard>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x20 = cell.template getFieldComponent<descriptors::SOURCE>(0);
auto x24 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x19 = cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0);
auto x22 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x25 = x24 + V{-1};
auto x26 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x27 = V{0.0277777777777778}*x20;
auto x28 = x19 + V{-1};
auto x29 = -V{0.0555555555555556}*x24*x28 + x27;
auto x30 = V{0.0277777777777778}*x24*x28 + x27;
auto x0 = -cell[0]*x25 - V{0.666666666666667}*x20 + x24*(-V{0.666666666666667}*x19 + x26 + V{0.666666666666667});
auto x1 = -cell[1]*x25 - x29;
auto x2 = -cell[2]*x25 - x29;
auto x3 = -cell[3]*x25 - x29;
auto x4 = -cell[4]*x25 + x30;
auto x5 = -cell[5]*x25 + x30;
auto x6 = -cell[6]*x25 + x30;
auto x7 = -cell[7]*x25 + x30;
auto x8 = -cell[8]*x25 + x30;
auto x9 = -cell[9]*x25 + x30;
auto x10 = -cell[10]*x25 - x29;
auto x11 = -cell[11]*x25 - x29;
auto x12 = -cell[12]*x25 - x29;
auto x13 = -cell[13]*x25 + x30;
auto x14 = -cell[14]*x25 + x30;
auto x15 = -cell[15]*x25 + x30;
auto x16 = -cell[16]*x25 + x30;
auto x17 = -cell[17]*x25 + x30;
auto x18 = -cell[18]*x25 + x30;
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
return { x26 + V{1}, x21*x21 + x22*x22 + x23*x23 };
}
};

}

}
