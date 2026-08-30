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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x19 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x22 = cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x25 = parameters.template get<descriptors::OMEGA>();
auto x24 = cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x26 = x25 + V{-1};
auto x27 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x28 = x27 + V{1};
auto x29 = V{1} / (x25);
auto x30 = x19*x29;
auto x31 = x22 + x30;
auto x32 = x31*x31;
auto x33 = V{1.5}*x32;
auto x34 = x20*x29;
auto x35 = x23 + x34;
auto x36 = x35*x35;
auto x37 = V{1.5}*x36;
auto x38 = x21*x29;
auto x39 = x24 + x38;
auto x40 = x39*x39;
auto x41 = V{1.5}*x40;
auto x42 = x37 + x41 + V{-1};
auto x43 = x33 + x42;
auto x44 = V{0.0555555555555556}*x25;
auto x45 = V{3}*x32;
auto x46 = V{3}*x22;
auto x47 = V{3}*x30;
auto x48 = x46 + x47;
auto x49 = V{3}*x36;
auto x50 = x33 + V{-1};
auto x51 = V{3}*x23;
auto x52 = V{3}*x34;
auto x53 = x51 + x52;
auto x54 = V{3}*x40;
auto x55 = V{3}*x24;
auto x56 = V{3}*x38;
auto x57 = x55 + x56;
auto x58 = V{0.0277777777777778}*x25;
auto x59 = x31 + x35;
auto x60 = V{4.5}*(x59*x59);
auto x61 = x43 + x48;
auto x62 = -x23 + x31 - x34;
auto x63 = -x62;
auto x64 = -x51 - x52;
auto x65 = x31 + x39;
auto x66 = V{4.5}*(x65*x65);
auto x67 = -x21*x29 - x24;
auto x68 = x31 + x67;
auto x69 = -x68;
auto x70 = -x55 - x56;
auto x71 = x35 + x39;
auto x72 = V{4.5}*(x71*x71);
auto x73 = x43 + x53;
auto x74 = x35 + x67;
auto x75 = -x74;
auto x76 = -x37;
auto x77 = V{1} - x41;
auto x78 = x76 + x77;
auto x79 = x48 + x78;
auto x80 = -x33;
auto x81 = x53 + x80;
auto x82 = x57 + x80;
auto x83 = -x46 - x47;
auto x84 = x43 + x57;
auto x85 = V{0.5}*x19 + x22;
auto x86 = V{0.5}*x20 + x23;
auto x87 = V{0.5}*x21 + x24;
auto x0 = -cell[0]*x26 - V{0.333333333333333}*x25*(x28*x43 + V{1});
auto x1 = -cell[1]*x26 - x44*(x28*(x42 - x45 + x48) + V{1});
auto x2 = -cell[2]*x26 - x44*(x28*(x41 - x49 + x50 + x53) + V{1});
auto x3 = -cell[3]*x26 - x44*(x28*(x37 + x50 - x54 + x57) + V{1});
auto x4 = -cell[4]*x26 - x58*(x28*(x53 - x60 + x61) + V{1});
auto x5 = -(cell[5]*x26 + x58*(x28*(x61 + x64 - V{4.5}*x63*x63) + V{1}));
auto x6 = -cell[6]*x26 - x58*(x28*(x57 + x61 - x66) + V{1});
auto x7 = -(cell[7]*x26 + x58*(x28*(x61 + x70 - V{4.5}*x69*x69) + V{1}));
auto x8 = -cell[8]*x26 - x58*(x28*(x57 - x72 + x73) + V{1});
auto x9 = -(cell[9]*x26 + x58*(x28*(x70 + x73 - V{4.5}*x75*x75) + V{1}));
auto x10 = -cell[10]*x26 + x44*(x28*(x45 + x79) + V{-1});
auto x11 = -cell[11]*x26 + x44*(x28*(x49 + x77 + x81) + V{-1});
auto x12 = -cell[12]*x26 + x44*(x28*(x54 + x76 + x82 + V{1}) + V{-1});
auto x13 = -cell[13]*x26 + x58*(x28*(x60 + x79 + x81) + V{-1});
auto x14 = -(cell[14]*x26 + x58*(x28*(x73 + x83 - V{4.5}*x62*x62) + V{1}));
auto x15 = -cell[15]*x26 + x58*(x28*(x66 + x79 + x82) + V{-1});
auto x16 = -(cell[16]*x26 + x58*(x28*(x83 + x84 - V{4.5}*x68*x68) + V{1}));
auto x17 = -cell[17]*x26 + x58*(x28*(x57 + x72 + x78 + x81) + V{-1});
auto x18 = -(cell[18]*x26 + x58*(x28*(x64 + x84 - V{4.5}*x74*x74) + V{1}));
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
return { x27 + V{1}, x85*x85 + x86*x86 + x87*x87 };
}
};

}

}
