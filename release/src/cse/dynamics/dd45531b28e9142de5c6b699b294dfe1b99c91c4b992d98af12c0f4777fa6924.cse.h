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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x12 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0);
auto x13 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x14 = x12 + V{-1};
auto x15 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x16 = x15 + V{1};
auto x17 = x11*x11;
auto x18 = V{1.5}*x17;
auto x19 = x10*x10;
auto x20 = V{1.5}*x19;
auto x21 = x20 + V{-1};
auto x22 = x18 + x21;
auto x23 = x15 + V{1};
auto x24 = -V{3}*x13*x9 + x23;
auto x25 = V{0.0277777777777778}*x12;
auto x26 = V{3}*x11;
auto x27 = x10 - x11;
auto x28 = -x27;
auto x29 = V{3}*x10;
auto x30 = x22 + x29;
auto x31 = x24*x25;
auto x32 = V{0.111111111111111}*x12;
auto x33 = x24*x32;
auto x34 = -x29;
auto x35 = V{1} - x18;
auto x36 = V{3}*x19 + x35;
auto x37 = x10 + x11;
auto x38 = V{4.5}*(x37*x37);
auto x39 = V{3}*x17;
auto x40 = -x20 + x26;
auto x0 = -cell[0]*x14 + V{0.555555555555556}*x12*x24 - V{0.444444444444444}*x12*(x16*x22 + V{1});
auto x1 = -(cell[1]*x14 + x25*(x16*(-x26 + x30 - V{4.5}*x28*x28) + V{1}) + x31);
auto x2 = -cell[2]*x14 + V{0.111111111111111}*x12*(x16*(x34 + x36) + V{-1}) - x33;
auto x3 = -cell[3]*x14 - x25*(x16*(x26 + x30 - x38) + V{1}) - x31;
auto x4 = -cell[4]*x14 - x32*(x16*(x21 + x26 - x39) + V{1}) - x33;
auto x5 = -(cell[5]*x14 + x25*(x16*(x22 + x26 + x34 - V{4.5}*x27*x27) + V{1}) + x31);
auto x6 = -cell[6]*x14 + V{0.111111111111111}*x12*(x16*(x29 + x36) + V{-1}) - x33;
auto x7 = -cell[7]*x14 + V{0.0277777777777778}*x12*(x16*(x29 + x35 + x38 + x40) + V{-1}) - x31;
auto x8 = -cell[8]*x14 + V{0.111111111111111}*x12*(x16*(x39 + x40 + V{1}) + V{-1}) - x33;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x23, x17 + x19 };
}
};

}

}
