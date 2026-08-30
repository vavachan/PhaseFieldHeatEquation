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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x22 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x19 = cell.template getFieldComponent<descriptors::SOURCE>(0);
auto x23 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x24 = x23 + V{-1};
auto x25 = V{0.5}*x23 + V{-1};
auto x26 = x19*x25;
auto x27 = V{3}*x20;
auto x28 = x27 + V{-1};
auto x29 = V{0.0277777777777778}*x19;
auto x30 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + x29 + V{0.0555555555555556};
auto x31 = V{0.0555555555555556}*x26;
auto x32 = V{3}*x21;
auto x33 = x32 + V{-1};
auto x34 = V{3}*x22;
auto x35 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0138888888888889}*x19 + V{0.0277777777777778};
auto x36 = x25*x29;
auto x37 = -x27;
auto x38 = x32 + V{1};
auto x39 = x34 + V{1};
auto x40 = -x32;
auto x41 = x27 + V{1};
auto x42 = -x34;
auto x0 = -cell[0]*x24 + x23*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.166666666666667}*x19) - V{0.333333333333333}*x26;
auto x1 = -cell[1]*x24 - x23*(x28*x30 + V{0.0555555555555556}) - x31;
auto x2 = -cell[2]*x24 - x23*(x30*x33 + V{0.0555555555555556}) - x31;
auto x3 = -cell[3]*x24 - x23*(x30*(x34 + V{-1}) + V{0.0555555555555556}) - x31;
auto x4 = -cell[4]*x24 - x23*(x35*(x28 + x32) + V{0.0277777777777778}) - x36;
auto x5 = -cell[5]*x24 + x23*(x35*(x37 + x38) + V{-0.0277777777777778}) - x36;
auto x6 = -cell[6]*x24 - x23*(x35*(x28 + x34) + V{0.0277777777777778}) - x36;
auto x7 = -cell[7]*x24 + x23*(x35*(x37 + x39) + V{-0.0277777777777778}) - x36;
auto x8 = -cell[8]*x24 - x23*(x35*(x33 + x34) + V{0.0277777777777778}) - x36;
auto x9 = -cell[9]*x24 + x23*(x35*(x39 + x40) + V{-0.0277777777777778}) - x36;
auto x10 = -cell[10]*x24 + x23*(x30*x41 + V{-0.0555555555555556}) - x31;
auto x11 = -cell[11]*x24 + x23*(x30*x38 + V{-0.0555555555555556}) - x31;
auto x12 = -cell[12]*x24 + x23*(x30*x39 + V{-0.0555555555555556}) - x31;
auto x13 = -cell[13]*x24 + x23*(x35*(x32 + x41) + V{-0.0277777777777778}) - x36;
auto x14 = -cell[14]*x24 + x23*(x35*(x40 + x41) + V{-0.0277777777777778}) - x36;
auto x15 = -cell[15]*x24 + x23*(x35*(x34 + x41) + V{-0.0277777777777778}) - x36;
auto x16 = -cell[16]*x24 + x23*(x35*(x41 + x42) + V{-0.0277777777777778}) - x36;
auto x17 = -cell[17]*x24 + x23*(x35*(x34 + x38) + V{-0.0277777777777778}) - x36;
auto x18 = -cell[18]*x24 + x23*(x35*(x38 + x42) + V{-0.0277777777777778}) - x36;
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
return { cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{0.5}*x19 + V{1}, x20*x20 + x21*x21 + x22*x22 };
}
};

}

}
