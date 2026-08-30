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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x13 = parameters.template get<descriptors::OMEGA>();
auto x12 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = x13 + V{-1};
auto x15 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x16 = x15 + V{1};
auto x17 = V{1} / (x13);
auto x18 = x17*x9;
auto x19 = x11 + x18;
auto x20 = x19*x19;
auto x21 = V{1.5}*x20;
auto x22 = x10*x17;
auto x23 = x12 + x22;
auto x24 = x23*x23;
auto x25 = V{1.5}*x24;
auto x26 = x25 + V{-1};
auto x27 = x21 + x26;
auto x28 = V{0.0277777777777778}*x13;
auto x29 = -x11 - x18 + x23;
auto x30 = V{1} - x25;
auto x31 = V{3}*x12 + V{3}*x22;
auto x32 = -x21 + x31;
auto x33 = V{3}*x11;
auto x34 = V{3}*x18;
auto x35 = -x33 - x34;
auto x36 = V{0.111111111111111}*x13;
auto x37 = V{3}*x20;
auto x38 = x33 + x34;
auto x39 = x19 + x23;
auto x40 = V{4.5}*(x39*x39);
auto x41 = x27 + x31;
auto x42 = V{3}*x24;
auto x43 = -x29;
auto x44 = x30 + x38;
auto x45 = x11 + V{0.5}*x9;
auto x46 = V{0.5}*x10 + x12;
auto x0 = -cell[0]*x14 - V{0.444444444444444}*x13*(x16*x27 + V{1});
auto x1 = -cell[1]*x14 + x28*(x16*(x30 + x32 + x35 + V{4.5}*(x29*x29)) + V{-1});
auto x2 = -cell[2]*x14 - x36*(x16*(x26 - x37 + x38) + V{1});
auto x3 = -cell[3]*x14 - x28*(x16*(x38 - x40 + x41) + V{1});
auto x4 = -cell[4]*x14 - x36*(x16*(x21 + x31 - x42 + V{-1}) + V{1});
auto x5 = -(cell[5]*x14 + x28*(x16*(x35 + x41 - V{4.5}*x43*x43) + V{1}));
auto x6 = -cell[6]*x14 + x36*(x16*(x37 + x44) + V{-1});
auto x7 = -cell[7]*x14 + x28*(x16*(x32 + x40 + x44) + V{-1});
auto x8 = -cell[8]*x14 + x36*(x16*(x32 + x42 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x15 + V{1}, x45*x45 + x46*x46 };
}
};

}

}
