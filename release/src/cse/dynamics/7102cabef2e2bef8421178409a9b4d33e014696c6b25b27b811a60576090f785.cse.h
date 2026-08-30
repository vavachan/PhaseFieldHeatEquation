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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x10 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x12 = parameters.template get<statistics::AVERAGE_RHO>();
auto x13 = x11 + V{-1};
auto x14 = x10 + V{1};
auto x15 = V{1} / (x14);
auto x16 = V{2}*cell[1] + V{2}*cell[7] + V{2}*cell[8] + V{1};
auto x17 = cell[0] + cell[2] + cell[6] + x16;
auto x18 = x15*x17;
auto x19 = V{0.444444444444444}*x18;
auto x20 = x9*x9;
auto x21 = V{1.5}*x20;
auto x22 = x10*x10;
auto x23 = V{1.5}*x22;
auto x24 = x23 + V{-1};
auto x25 = x21 + x24;
auto x26 = -x14*(x12 + V{-1})/x17 + V{1};
auto x27 = V{0.0277777777777778}*x18;
auto x28 = V{3}*x10;
auto x29 = -x28;
auto x30 = x10 - x9;
auto x31 = V{3}*x9;
auto x32 = x25 + x31;
auto x33 = x29 + x32 - V{4.5}*x30*x30;
auto x34 = x26*x27;
auto x35 = V{0.111111111111111}*x18;
auto x36 = V{3}*x20;
auto x37 = x24 + x31 - x36;
auto x38 = x10 + x9;
auto x39 = V{4.5}*(x38*x38);
auto x40 = x28 + x32 - x39;
auto x41 = V{1} - x21;
auto x42 = V{3}*x22 + x41;
auto x43 = x29 + x42;
auto x44 = -x30;
auto x45 = x25 + x28 - x31 - V{4.5}*x44*x44;
auto x46 = -x23 + x31;
auto x47 = x36 + x46 + V{1};
auto x48 = x28 + x39 + x41 + x46;
auto x49 = x28 + x42;
auto x0 = -x13*(cell[0] + x19*x25 + V{0.444444444444444}) - x19*x25*x26 + V{-0.444444444444444};
auto x1 = -x13*(cell[1] + x27*x33 + V{0.0277777777777778}) - x33*x34 + V{-0.0277777777777778};
auto x2 = -x13*(cell[2] + x35*x37 + V{0.111111111111111}) - x26*x35*x37 + V{-0.111111111111111};
auto x3 = -x13*(cell[3] + x27*x40 + V{0.0277777777777778}) - x34*x40 + V{-0.0277777777777778};
auto x4 = -x13*(cell[4] - x35*x43 + V{0.111111111111111}) + V{0.111111111111111}*x15*x17*x26*x43 + V{-0.111111111111111};
auto x5 = -x13*(cell[5] + x27*x45 + V{0.0277777777777778}) - x34*x45 + V{-0.0277777777777778};
auto x6 = -x13*(cell[6] - x35*x47 + V{0.111111111111111}) + V{0.111111111111111}*x15*x17*x26*x47 + V{-0.111111111111111};
auto x7 = -x13*(cell[7] - x27*x48 + V{0.0277777777777778}) + V{0.0277777777777778}*x15*x17*x26*x48 + V{-0.0277777777777778};
auto x8 = -x13*(cell[8] - x35*x49 + V{0.111111111111111}) + V{0.111111111111111}*x15*x17*x26*x49 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x12 + x15*(V{1}*cell[0] + V{1}*cell[2] + V{1}*cell[6] + x16) + V{1}, x20 + x22 };
}
};

}

}
