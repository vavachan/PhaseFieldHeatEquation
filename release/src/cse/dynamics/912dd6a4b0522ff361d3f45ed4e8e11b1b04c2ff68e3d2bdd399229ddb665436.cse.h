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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x12 = x11 + V{-1};
auto x13 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x14 = x13 + V{1};
auto x15 = x13 + V{1};
auto x16 = V{1} / (x15);
auto x17 = V{1}*cell[3];
auto x18 = V{1}*cell[1] - V{1}*cell[5];
auto x19 = V{0.5}*x10 + x16*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x17 + x18);
auto x20 = x19*x19;
auto x21 = V{1.5}*x20;
auto x22 = V{0.5}*x9;
auto x23 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x17 + x18;
auto x24 = -x23;
auto x25 = x16*x24 + x22;
auto x26 = x21 + V{-1} + V{1.5}*(x25*x25);
auto x27 = V{0.5}*x11 + V{-1};
auto x28 = V{3}*cell[3];
auto x29 = V{3}*cell[1] - V{3}*cell[5];
auto x30 = V{1.5}*x10 + x16*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x28 + x29);
auto x31 = x10*x30;
auto x32 = V{1.5}*x9;
auto x33 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x28 + x29;
auto x34 = -x16*x33;
auto x35 = x32 + x34;
auto x36 = x35*x9;
auto x37 = -x32;
auto x38 = x16*x33;
auto x39 = x16*x23;
auto x40 = x19 - x22;
auto x41 = V{4.5}*cell[3];
auto x42 = V{4.5}*cell[7];
auto x43 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x44 = V{4.5}*cell[2] - V{4.5}*cell[6] + x41 - x42 + x43;
auto x45 = x16*x44;
auto x46 = V{2.25}*x9;
auto x47 = V{2.25}*x10 + x16*(-V{4.5}*cell[4] + V{4.5}*cell[8] - x41 + x42 + x43);
auto x48 = -x46 + x47;
auto x49 = x22 - x39;
auto x50 = x49*x49;
auto x51 = -x21 - V{1.5}*x50 + V{1};
auto x52 = x30 + x51;
auto x53 = V{4.5}*x9;
auto x54 = V{9}*cell[3];
auto x55 = V{9}*cell[7];
auto x56 = V{9}*cell[1] - V{9}*cell[5];
auto x57 = V{9}*cell[2] - V{9}*cell[6] + x54 - x55 + x56;
auto x58 = x16*x57;
auto x59 = V{3}*x10;
auto x60 = V{6}*cell[3];
auto x61 = V{6}*cell[1] - V{6}*cell[5];
auto x62 = x16*(-V{6}*cell[4] + V{6}*cell[7] + V{6}*cell[8] - x60 + x61);
auto x63 = x59 + x62;
auto x64 = x63 + V{3};
auto x65 = V{6}*cell[2] - V{6}*cell[6] - V{6}*cell[7] + x60 + x61;
auto x66 = V{3}*x9;
auto x67 = -x66;
auto x68 = x67 + V{3};
auto x69 = V{4.5}*x10 + x16*(-V{9}*cell[4] + V{9}*cell[8] - x54 + x55 + x56);
auto x70 = x15*x27;
auto x71 = V{0.0277777777777778}*x70;
auto x72 = V{0.111111111111111}*x11;
auto x73 = -x44;
auto x74 = x16*x73 + x46;
auto x75 = x26 + x35;
auto x76 = x16*x65;
auto x77 = V{0.111111111111111}*x70;
auto x78 = V{0.0277777777777778}*x11;
auto x79 = x63 + V{-3};
auto x80 = x53 - x58;
auto x81 = x66 - x76;
auto x82 = x69 + V{-3};
auto x83 = x19*x47;
auto x84 = x26 + x30;
auto x85 = -x36;
auto x86 = -x45 + x46;
auto x87 = x32 - x38 + x51;
auto x88 = x81 + V{3};
auto x0 = -cell[0]*x12 - V{0.444444444444444}*x11*(x14*x26 + V{1}) + V{0.444444444444444}*x15*x27*(x31 + x36);
auto x1 = -cell[1]*x12 + V{0.0277777777777778}*x11*(x14*(x37 + x38 + x52 + (x39 + x40)*(x45 + x48)) + V{-1}) - x71*(x10*(-x53 + x58 + x64) - x9*(x16*x65 + x68 + x69));
auto x2 = -cell[2]*x12 - x72*(x14*(-x25*x74 + x75) + V{1}) - x77*(-x31 + x9*(-x68 - x76));
auto x3 = -cell[3]*x12 - x71*(x10*(x79 + x80) + x9*(x81 + x82)) - x78*(x14*(x30 + x75 - (x19 + x25)*(x47 + x74)) + V{1});
auto x4 = -cell[4]*x12 - x72*(x14*(-x83 + x84) + V{1}) - x77*(x10*x79 + x85);
auto x5 = -cell[5]*x12 - x71*(-x10*(-x16*x57 + x53 - x59 - x62 + V{3}) + x9*(-x67 - x76 - x82)) - x78*(x14*(-x34 + x37 + x84 - (x16*x24 - x40)*(x16*x73 - x48)) + V{1});
auto x6 = -cell[6]*x12 + V{0.111111111111111}*x11*(x14*(x49*x86 + x87) + V{-1}) - x77*(-x31 + x88*x9);
auto x7 = -cell[7]*x12 + V{0.0277777777777778}*x11*(x14*(x30 + x87 + (x19 + x49)*(x47 + x86)) + V{-1}) - x71*(x10*(x64 + x80) + x9*(x69 + x88));
auto x8 = -cell[8]*x12 + V{0.111111111111111}*x11*(x14*(x52 + x83) + V{-1}) - x77*(x10*x64 + x85);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x15, x20 + x50 };
}
};

}

}
