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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x12 = parameters.template get<statistics::AVERAGE_RHO>();
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x13 = x11 + V{-1};
auto x14 = x9*x9;
auto x15 = V{1.5}*x14;
auto x16 = x10*x10;
auto x17 = V{1.5}*x16;
auto x18 = x17 + V{-1};
auto x19 = x15 + x18;
auto x20 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x21 = -(x12 + V{-1})/(x20 + V{1}) + V{1};
auto x22 = x20 + V{1};
auto x23 = x21*x22;
auto x24 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x25 = V{3}*x10;
auto x26 = -x25;
auto x27 = x10 - x9;
auto x28 = V{3}*x9;
auto x29 = x19 + x28;
auto x30 = x26 + x29 - V{4.5}*x27*x27;
auto x31 = V{0.0277777777777778}*x23;
auto x32 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x33 = V{3}*x14;
auto x34 = x18 + x28 - x33;
auto x35 = x10 + x9;
auto x36 = V{4.5}*(x35*x35);
auto x37 = x25 + x29 - x36;
auto x38 = V{1} - x15;
auto x39 = V{3}*x16 + x38;
auto x40 = x26 + x39;
auto x41 = -x27;
auto x42 = x19 + x25 - x28 - V{4.5}*x41*x41;
auto x43 = -x17 + x28;
auto x44 = x33 + x43 + V{1};
auto x45 = x25 + x36 + x38 + x43;
auto x46 = x25 + x39;
auto x0 = -x13*(cell[0] + x19*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444}) - V{0.444444444444444}*x19*x23 + V{-0.444444444444444};
auto x1 = -x13*(cell[1] + x24*x30 + V{0.0277777777777778}) - x30*x31 + V{-0.0277777777777778};
auto x2 = -x13*(cell[2] + x32*x34 + V{0.111111111111111}) - V{0.111111111111111}*x23*x34 + V{-0.111111111111111};
auto x3 = -x13*(cell[3] + x24*x37 + V{0.0277777777777778}) - x31*x37 + V{-0.0277777777777778};
auto x4 = -x13*(cell[4] - x32*x40 + V{0.111111111111111}) + V{0.111111111111111}*x21*x22*x40 + V{-0.111111111111111};
auto x5 = -x13*(cell[5] + x24*x42 + V{0.0277777777777778}) - x31*x42 + V{-0.0277777777777778};
auto x6 = -x13*(cell[6] - x32*x44 + V{0.111111111111111}) + V{0.111111111111111}*x21*x22*x44 + V{-0.111111111111111};
auto x7 = -x13*(cell[7] - x24*x45 + V{0.0277777777777778}) + V{0.0277777777777778}*x21*x22*x45 + V{-0.0277777777777778};
auto x8 = -x13*(cell[8] - x32*x46 + V{0.111111111111111}) + V{0.111111111111111}*x21*x22*x46 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x12 + x20 + V{2}, x14 + x16 };
}
};

}

}
