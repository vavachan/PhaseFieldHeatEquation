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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x23 = x22 + V{-1};
auto x24 = (cell[0] + cell[10] + cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[1] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + V{1})/(x21 + V{1});
auto x25 = x19*x19;
auto x26 = V{1.5}*x25;
auto x27 = x20*x20;
auto x28 = V{1.5}*x27;
auto x29 = x21*x21;
auto x30 = V{1.5}*x29;
auto x31 = x28 + x30 + V{-1};
auto x32 = x26 + x31;
auto x33 = V{0.0555555555555556}*x22;
auto x34 = V{3}*x19;
auto x35 = V{3}*x20;
auto x36 = x26 + V{-1};
auto x37 = V{3}*x21;
auto x38 = V{0.0277777777777778}*x22;
auto x39 = x19 + x20;
auto x40 = x39*x39;
auto x41 = x32 + x34;
auto x42 = -x35;
auto x43 = x19 - x20;
auto x44 = -x43;
auto x45 = x19 + x21;
auto x46 = x45*x45;
auto x47 = -x37;
auto x48 = -x21;
auto x49 = x19 + x48;
auto x50 = -x49;
auto x51 = x20 + x21;
auto x52 = x51*x51;
auto x53 = x32 + x35;
auto x54 = x20 + x48;
auto x55 = -x54;
auto x56 = -x28;
auto x57 = V{1} - x30;
auto x58 = x56 + x57;
auto x59 = x34 + x58;
auto x60 = -x26;
auto x61 = x35 + x60;
auto x62 = x37 + x60;
auto x63 = -x34;
auto x64 = x32 + x37;
auto x0 = -cell[0]*x23 - V{0.333333333333333}*x22*(x24*x32 + V{1});
auto x1 = -cell[1]*x23 - x33*(-x24*(V{3}*x25 - x31 - x34) + V{1});
auto x2 = -cell[2]*x23 - x33*(-x24*(V{3}*x27 - x30 - x35 - x36) + V{1});
auto x3 = -cell[3]*x23 - x33*(-x24*(-x28 + V{3}*x29 - x36 - x37) + V{1});
auto x4 = -cell[4]*x23 - x38*(-x24*(-x35 + V{4.5}*x40 - x41) + V{1});
auto x5 = -(cell[5]*x23 + x38*(x24*(x41 + x42 - V{4.5}*x44*x44) + V{1}));
auto x6 = -cell[6]*x23 - x38*(-x24*(-x37 - x41 + V{4.5}*x46) + V{1});
auto x7 = -(cell[7]*x23 + x38*(x24*(x41 + x47 - V{4.5}*x50*x50) + V{1}));
auto x8 = -cell[8]*x23 - x38*(-x24*(-x37 + V{4.5}*x52 - x53) + V{1});
auto x9 = -(cell[9]*x23 + x38*(x24*(x47 + x53 - V{4.5}*x55*x55) + V{1}));
auto x10 = -cell[10]*x23 - x33*(-x24*(V{3}*x25 + x59) + V{1});
auto x11 = -cell[11]*x23 - x33*(-x24*(V{3}*x27 + x57 + x61) + V{1});
auto x12 = -cell[12]*x23 - x33*(-x24*(V{3}*x29 + x56 + x62 + V{1}) + V{1});
auto x13 = -cell[13]*x23 - x38*(-x24*(V{4.5}*x40 + x59 + x61) + V{1});
auto x14 = -(cell[14]*x23 + x38*(x24*(x53 + x63 - V{4.5}*x43*x43) + V{1}));
auto x15 = -cell[15]*x23 - x38*(-x24*(V{4.5}*x46 + x59 + x62) + V{1});
auto x16 = -(cell[16]*x23 + x38*(x24*(x63 + x64 - V{4.5}*x49*x49) + V{1}));
auto x17 = -cell[17]*x23 - x38*(-x24*(x37 + V{4.5}*x52 + x58 + x61) + V{1});
auto x18 = -(cell[18]*x23 + x38*(x24*(x42 + x64 - V{4.5}*x54*x54) + V{1}));
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
return { V{1}*x24, x25 + x27 + x29 };
}
};

}

}
