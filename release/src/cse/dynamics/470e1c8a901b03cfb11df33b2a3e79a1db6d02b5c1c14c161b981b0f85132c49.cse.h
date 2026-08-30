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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q7<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x7 = cell.template getFieldComponent<descriptors::SOURCE>(0);
auto x9 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x8 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x12 = x11 + V{-1};
auto x13 = V{0.125}*x7;
auto x14 = V{0.5}*x11 + V{-1};
auto x15 = V{4}*x8;
auto x16 = V{0.125}*cell[0] + V{0.125}*cell[1] + V{0.125}*cell[2] + V{0.125}*cell[3] + V{0.125}*cell[4] + V{0.125}*cell[5] + V{0.125}*cell[6] + V{0.0625}*x7 + V{0.125};
auto x17 = x13*x14;
auto x18 = V{4}*x9;
auto x19 = V{4}*x10;
auto x0 = -cell[0]*x12 + x11*(V{0.25}*cell[0] + V{0.25}*cell[1] + V{0.25}*cell[2] + V{0.25}*cell[3] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.25}*cell[6] + x13) - V{0.25}*x14*x7;
auto x1 = -cell[1]*x12 - x11*(x16*(x15 + V{-1}) + V{0.125}) - x17;
auto x2 = -cell[2]*x12 - x11*(x16*(x18 + V{-1}) + V{0.125}) - x17;
auto x3 = -cell[3]*x12 - x11*(x16*(x19 + V{-1}) + V{0.125}) - x17;
auto x4 = -cell[4]*x12 + x11*(x16*(x15 + V{1}) + V{-0.125}) - x17;
auto x5 = -cell[5]*x12 + x11*(x16*(x18 + V{1}) + V{-0.125}) - x17;
auto x6 = -cell[6]*x12 + x11*(x16*(x19 + V{1}) + V{-0.125}) - x17;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
return { cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{0.5}*x7 + V{1}, x10*x10 + x8*x8 + x9*x9 };
}
};

}

}
