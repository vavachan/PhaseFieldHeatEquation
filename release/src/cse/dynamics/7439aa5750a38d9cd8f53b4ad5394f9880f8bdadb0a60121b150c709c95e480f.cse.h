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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x13 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x12 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = x13 + V{-1};
auto x15 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x16 = x15 + V{1};
auto x17 = x12*x12;
auto x18 = V{1.5}*x17;
auto x19 = x11*x11;
auto x20 = V{1.5}*x19;
auto x21 = x20 + V{-1};
auto x22 = x18 + x21;
auto x23 = V{0.5}*x13 + V{-1};
auto x24 = x10*x12;
auto x25 = x11*x9;
auto x26 = x15 + V{1};
auto x27 = V{0.0277777777777778}*x13;
auto x28 = V{3}*x12;
auto x29 = x11 - x12;
auto x30 = -x29;
auto x31 = V{3}*x11;
auto x32 = x22 + x31;
auto x33 = V{9}*x11;
auto x34 = V{6}*x12;
auto x35 = x34 + V{3};
auto x36 = V{9}*x12;
auto x37 = V{6}*x11;
auto x38 = x23*x26;
auto x39 = V{0.0277777777777778}*x38;
auto x40 = V{0.111111111111111}*x13;
auto x41 = -x31;
auto x42 = V{1} - x18;
auto x43 = V{3}*x19 + x42;
auto x44 = V{0.333333333333333}*x24;
auto x45 = x37 + V{-3};
auto x46 = V{0.111111111111111}*x9;
auto x47 = x11 + x12;
auto x48 = V{4.5}*(x47*x47);
auto x49 = x34 + V{-3};
auto x50 = V{3}*x17;
auto x51 = V{0.111111111111111}*x10;
auto x52 = -V{0.333333333333333}*x25;
auto x53 = x37 + V{3};
auto x54 = -x20 + x28;
auto x0 = -cell[0]*x14 - V{0.444444444444444}*x13*(x16*x22 + V{1}) + V{1.33333333333333}*x23*x26*(x24 + x25);
auto x1 = -(cell[1]*x14 + x27*(x16*(-x28 + x32 - V{4.5}*x30*x30) + V{1}) + x39*(x10*(-x33 + x35) - x9*(x36 - x37 + V{3})));
auto x2 = -cell[2]*x14 + x38*(x44 - x45*x46) + x40*(x16*(x41 + x43) + V{-1});
auto x3 = -cell[3]*x14 - x27*(x16*(x28 + x32 - x48) + V{1}) - x39*(x10*(x33 + x49) + x9*(x36 + x45));
auto x4 = -cell[4]*x14 - x38*(x49*x51 + x52) - x40*(x16*(x21 + x28 - x50) + V{1});
auto x5 = -(cell[5]*x14 + x27*(x16*(x22 + x28 + x41 - V{4.5}*x29*x29) + V{1}) + x39*(-x10*(x33 - x34 + V{3}) + x9*(-x36 + x53)));
auto x6 = -cell[6]*x14 + x38*(x44 - x46*x53) + x40*(x16*(x31 + x43) + V{-1});
auto x7 = -cell[7]*x14 + V{0.0277777777777778}*x13*(x16*(x31 + x42 + x48 + x54) + V{-1}) - x39*(x10*(x33 + x35) + x9*(x36 + x53));
auto x8 = -cell[8]*x14 + V{0.111111111111111}*x13*(x16*(x50 + x54 + V{1}) + V{-1}) - x38*(x35*x51 + x52);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x26, x17 + x19 };
}
};

}

}
