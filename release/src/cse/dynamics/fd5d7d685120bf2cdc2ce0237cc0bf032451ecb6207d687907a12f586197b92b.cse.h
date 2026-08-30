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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, -1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x12 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x10 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, -1>::VELOCITY>(0);
auto x11 = x12 + V{-1};
auto x13 = x10*x10;
auto x14 = V{1.5}*x13;
auto x15 = V{1} / (x9);
auto x16 = x15*(cell[0] + cell[2] + V{2}*cell[3] + V{2}*cell[4] + V{2}*cell[5] + cell[6] + V{1});
auto x17 = V{1} - x16;
auto x18 = x17*x17;
auto x19 = V{1.5}*x18;
auto x20 = x14 + x19;
auto x21 = V{0.0277777777777778}*x12;
auto x22 = V{3}*x10;
auto x23 = x10 + x16 + V{-1};
auto x24 = -x23;
auto x25 = x15*(V{3}*cell[0] + V{3}*cell[2] + V{6}*cell[3] + V{6}*cell[4] + V{6}*cell[5] + V{3}*cell[6] + V{3});
auto x26 = x25 + V{-4};
auto x27 = x20 + x26;
auto x28 = V{0.111111111111111}*x12;
auto x29 = -x22;
auto x30 = V{3}*x13 - x19 + V{1};
auto x31 = x10 + x17;
auto x32 = -x31;
auto x33 = x20 - x25 + V{2};
auto x34 = -x17;
auto x0 = -cell[0]*x11 - V{0.444444444444444}*x12*(x9*(x20 + V{-1}) + V{1});
auto x1 = -(cell[1]*x11 + x21*(x9*(x22 + x27 - V{4.5}*x24*x24) + V{1}));
auto x2 = -cell[2]*x11 + x28*(x9*(x29 + x30) + V{-1});
auto x3 = -(cell[3]*x11 + x21*(x9*(x22 + x33 - V{4.5}*x32*x32) + V{1}));
auto x4 = -(cell[4]*x11 + x28*(x9*(x33 - V{4.5}*x34*x34) + V{1}));
auto x5 = -(cell[5]*x11 + x21*(x9*(x29 + x33 - V{4.5}*x23*x23) + V{1}));
auto x6 = -cell[6]*x11 + x28*(x9*(x22 + x30) + V{-1});
auto x7 = -(cell[7]*x11 + x21*(x9*(x27 + x29 - V{4.5}*x31*x31) + V{1}));
auto x8 = -cell[8]*x11 - x28*(x9*(x14 - V{3}*x18 + x26) + V{1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x9, x13 + V{1}*x18 };
}
};

}

}
