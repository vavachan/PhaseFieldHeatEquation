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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x9 = x11 + V{-1};
auto x10 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x10 + V{1};
auto x13 = x10 + V{1};
auto x14 = V{1} / ((x13)*(x13));
auto x15 = V{1.5}*x14;
auto x16 = cell[1] - cell[5];
auto x17 = -cell[4] + cell[8];
auto x18 = -cell[3] + cell[7] + x16 + x17;
auto x19 = x18*x18;
auto x20 = x15*x19;
auto x21 = cell[2] - cell[6];
auto x22 = cell[3] - cell[7] + x16 + x21;
auto x23 = -x22;
auto x24 = x15*(x23*x23) + V{-1};
auto x25 = x20 + x24;
auto x26 = cell[0]*x9 + V{0.444444444444444}*x11*(x12*x25 + V{1});
auto x27 = cell[1]*x9;
auto x28 = V{0.0277777777777778}*x11;
auto x29 = V{4.5}*x14;
auto x30 = V{2}*cell[1] - V{2}*cell[5] + x17 + x21;
auto x31 = -x30;
auto x32 = V{1} / (x13);
auto x33 = V{3}*cell[3];
auto x34 = V{3}*cell[7];
auto x35 = V{3}*cell[1] - V{3}*cell[5];
auto x36 = V{3}*cell[2] - V{3}*cell[6] + x33 - x34 + x35;
auto x37 = x32*x36;
auto x38 = V{1} - x20;
auto x39 = x32*(-V{3}*cell[4] + V{3}*cell[8] - x33 + x34 + x35);
auto x40 = x22*x22;
auto x41 = x15*x40;
auto x42 = x39 - x41;
auto x43 = x38 + x42;
auto x44 = x37 + x43;
auto x45 = x28*(x12*(x29*(x31*x31) + x44) + V{-1});
auto x46 = cell[2]*x9;
auto x47 = V{0.111111111111111}*x11;
auto x48 = V{3}*x14;
auto x49 = x38 + x40*x48;
auto x50 = x12*(x37 + x49) + V{-1};
auto x51 = cell[3]*x9;
auto x52 = V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8] + x21;
auto x53 = -x52;
auto x54 = x29*(x53*x53);
auto x55 = x28*(x12*(x25 - x32*x36 + x39 - x54) + V{1});
auto x56 = x51 + x55;
auto x57 = cell[4]*x9;
auto x58 = x19*x48;
auto x59 = x47*(x12*(x24 + x39 - x58) + V{1});
auto x60 = x57 + x59;
auto x61 = cell[5]*x9;
auto x62 = x30*x30;
auto x63 = -V{4.5}*x14*x62 + x20 + x37 + x39 + x41 + V{-1};
auto x64 = x28*(-x12*x63 + V{-1});
auto x65 = cell[6]*x9;
auto x66 = -x37;
auto x67 = x47*(x12*(x49 + x66) + V{-1});
auto x68 = cell[7]*x9;
auto x69 = x43 + x66;
auto x70 = x28*(x12*(x54 + x69) + V{-1});
auto x71 = cell[8]*x9;
auto x72 = x47*(x12*(x42 + x58 + V{1}) + V{-1});
auto x73 = -V{0.111111111111111}*x11*x50;
auto x74 = -V{1} / (-V{0.0277777777777778}*x11*(x12*(x29*x62 + x44) + V{-1}) - V{0.0277777777777778}*x11*(x12*(x29*(x52*x52) + x69) + V{-1}) + x26 + x27 + x28*(x12*x63 + V{1}) + x46 + x56 + x60 + x61 + x65 - x67 + x68 + x71 - x72 + x73 + V{-1});
auto x75 = V{1}*x27;
auto x76 = V{1}*x61;
auto x77 = -V{1}*cell[7]*x9 + V{1}*x51 + x55 + x70;
cell[0] = -x26;
cell[1] = -x27 + x45;
cell[2] = -x46 + x47*x50;
cell[3] = -x56;
cell[4] = -x60;
cell[5] = -x61 + x64;
cell[6] = -x65 + x67;
cell[7] = -x68 + x70;
cell[8] = -x71 + x72;
cell.template getFieldPointer<descriptors::VELOCITY>()[0] = -x74*(V{1}*cell[6]*x9 + x45 - V{1}*x46 - x64 - x67 - x73 - x75 + x76 - x77);
cell.template getFieldPointer<descriptors::VELOCITY>()[1] = -x74*(V{1}*cell[8]*x9 - x45 - V{1}*x57 - x59 + x64 - x72 + x75 - x76 - x77);
return { x13, V{1}*x14*(x19 + x40) };
}
};

}

}
