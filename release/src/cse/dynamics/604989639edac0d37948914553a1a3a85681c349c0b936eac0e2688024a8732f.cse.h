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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = parameters.template get<descriptors::OMEGA>();
auto x21 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x19 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x22 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(2);
auto x20 = x23 + V{-1};
auto x24 = x21*x21;
auto x25 = V{1.5}*x24;
auto x26 = x22*x22;
auto x27 = V{1.5}*x26;
auto x28 = x25 + x27;
auto x29 = V{1} / (x19);
auto x30 = x29*(cell[0] + V{2}*cell[10] + cell[11] + cell[12] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + V{1});
auto x31 = V{1} - x30;
auto x32 = -x31;
auto x33 = x32*x32;
auto x34 = V{1.5}*x33;
auto x35 = x34 + V{-1};
auto x36 = x28 + x35;
auto x37 = V{0.0555555555555556}*x23;
auto x38 = x31*x31;
auto x39 = x29*(V{3}*cell[0] + V{6}*cell[10] + V{3}*cell[11] + V{3}*cell[12] + V{6}*cell[13] + V{6}*cell[14] + V{6}*cell[15] + V{6}*cell[16] + V{3}*cell[17] + V{3}*cell[18] + V{3}*cell[2] + V{3}*cell[3] + V{3}*cell[8] + V{3}*cell[9] + V{3});
auto x40 = x28 + x34 + x39 + V{-4};
auto x41 = V{3}*x21;
auto x42 = V{3}*x24;
auto x43 = V{3}*x22;
auto x44 = V{3}*x26;
auto x45 = V{0.0277777777777778}*x23;
auto x46 = x30 + V{-1};
auto x47 = x21 + x46;
auto x48 = -x47;
auto x49 = -x41;
auto x50 = x21 + x31;
auto x51 = x22 + x46;
auto x52 = -x51;
auto x53 = -x43;
auto x54 = x22 + x31;
auto x55 = x21 + x22;
auto x56 = V{4.5}*(x55*x55);
auto x57 = x36 + x41;
auto x58 = x21 - x22;
auto x59 = -x58;
auto x60 = x28 - x39 + V{2};
auto x61 = V{1} - V{1.5}*x38;
auto x62 = -x27 + x41 + x61;
auto x63 = -x25 + x43;
auto x64 = x34 + x60;
auto x65 = -x50;
auto x66 = -x54;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x23*(x19*x36 + V{1});
auto x1 = -cell[1]*x20 - x37*(x19*(-V{4.5}*x38 + x40) + V{1});
auto x2 = -cell[2]*x20 - x37*(x19*(x27 + x35 + x41 - x42) + V{1});
auto x3 = -cell[3]*x20 - x37*(x19*(x25 + x35 + x43 - x44) + V{1});
auto x4 = -(cell[4]*x20 + x45*(x19*(x40 + x41 - V{4.5}*x48*x48) + V{1}));
auto x5 = -(cell[5]*x20 + x45*(x19*(x40 + x49 - V{4.5}*x50*x50) + V{1}));
auto x6 = -(cell[6]*x20 + x45*(x19*(x40 + x43 - V{4.5}*x52*x52) + V{1}));
auto x7 = -(cell[7]*x20 + x45*(x19*(x40 + x53 - V{4.5}*x54*x54) + V{1}));
auto x8 = -cell[8]*x20 - x45*(x19*(x43 - x56 + x57) + V{1});
auto x9 = -(cell[9]*x20 + x45*(x19*(x53 + x57 - V{4.5}*x59*x59) + V{1}));
auto x10 = -cell[10]*x20 - x37*(x19*(-V{3}*x33 + x60) + V{1});
auto x11 = -cell[11]*x20 + x37*(x19*(x42 + x62) + V{-1});
auto x12 = -cell[12]*x20 + x37*(x19*(x44 + x61 + x63) + V{-1});
auto x13 = -(cell[13]*x20 + x45*(x19*(x49 + x64 - V{4.5}*x47*x47) + V{1}));
auto x14 = -(cell[14]*x20 + x45*(x19*(x41 + x64 - V{4.5}*x65*x65) + V{1}));
auto x15 = -(cell[15]*x20 + x45*(x19*(x53 + x64 - V{4.5}*x51*x51) + V{1}));
auto x16 = -(cell[16]*x20 + x45*(x19*(x43 + x64 - V{4.5}*x66*x66) + V{1}));
auto x17 = -cell[17]*x20 + x45*(x19*(x56 + x62 + x63) + V{-1});
auto x18 = -(cell[18]*x20 + x45*(x19*(x36 + x43 + x49 - V{4.5}*x58*x58) + V{1}));
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
return { x19, x24 + x26 + V{1}*x38 };
}
};

}

}
