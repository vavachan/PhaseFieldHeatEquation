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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x11 = x10 + V{-1};
auto x12 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x13 = x12 + V{1};
auto x14 = x12 + V{1};
auto x15 = V{1} / ((x14)*(x14));
auto x16 = V{1.5}*x15;
auto x17 = x9*x9;
auto x18 = cell[1] - cell[5];
auto x19 = -cell[4] + cell[8];
auto x20 = -cell[3] + cell[7] + x18 + x19;
auto x21 = x20*x20;
auto x22 = x17*x21;
auto x23 = x16*x22;
auto x24 = cell[2] - cell[6];
auto x25 = cell[3] - cell[7] + x18 + x24;
auto x26 = -x25;
auto x27 = x16*x17*(x26*x26) + V{-1};
auto x28 = x23 + x27;
auto x29 = V{0.0277777777777778}*x10;
auto x30 = V{3}*x9/x14;
auto x31 = x20*x30;
auto x32 = V{4.5}*x15;
auto x33 = V{2}*cell[1] - V{2}*cell[5] + x19 + x24;
auto x34 = x26*x30;
auto x35 = x28 + x34;
auto x36 = V{0.111111111111111}*x10;
auto x37 = x25*x30;
auto x38 = V{3}*x15;
auto x39 = x25*x25;
auto x40 = x17*x39;
auto x41 = V{1} - x23;
auto x42 = x38*x40 + x41;
auto x43 = -V{2}*cell[3] - cell[4] + V{2}*cell[7] + cell[8] - x24;
auto x44 = x17*x32*(x43*x43);
auto x45 = x22*x38;
auto x46 = -x33;
auto x47 = -x37;
auto x48 = -x16*x40 + x31;
auto x0 = -cell[0]*x11 - V{0.444444444444444}*x10*(x13*x28 + V{1});
auto x1 = -(cell[1]*x11 + x29*(x13*(-x17*x32*x33*x33 - x31 + x35) + V{1}));
auto x2 = -cell[2]*x11 + x36*(x13*(x37 + x42) + V{-1});
auto x3 = -cell[3]*x11 - x29*(x13*(x31 + x35 - x44) + V{1});
auto x4 = -cell[4]*x11 - x36*(x13*(x27 + x31 - x45) + V{1});
auto x5 = -(cell[5]*x11 + x29*(x13*(-x17*x32*x46*x46 + x28 + x31 - x34) + V{1}));
auto x6 = -cell[6]*x11 + x36*(x13*(x42 + x47) + V{-1});
auto x7 = -cell[7]*x11 + x29*(x13*(x41 + x44 + x47 + x48) + V{-1});
auto x8 = -cell[8]*x11 + x36*(x13*(x45 + x48 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x14, V{1}*x15*x17*(x21 + x39) };
}
};

}

}
