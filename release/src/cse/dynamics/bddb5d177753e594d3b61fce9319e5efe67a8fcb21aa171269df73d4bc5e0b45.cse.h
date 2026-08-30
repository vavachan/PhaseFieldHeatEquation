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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x30 = parameters.template get<descriptors::OMEGA>();
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x31 = x30 + V{-1};
auto x32 = (cell[0] + V{2}*cell[11] + V{2}*cell[13] + cell[14] + cell[15] + V{2}*cell[16] + cell[17] + cell[18] + V{2}*cell[19] + cell[1] + V{2}*cell[21] + V{2}*cell[23] + V{2}*cell[25] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + V{1})/(x29 + V{1});
auto x33 = x27*x27;
auto x34 = V{1.5}*x33;
auto x35 = x28*x28;
auto x36 = V{1.5}*x35;
auto x37 = x29*x29;
auto x38 = V{1.5}*x37;
auto x39 = x36 + x38 + V{-1};
auto x40 = x34 + x39;
auto x41 = V{0.0740740740740741}*x30;
auto x42 = V{3}*x27;
auto x43 = V{3}*x28;
auto x44 = x34 + V{-1};
auto x45 = V{3}*x29;
auto x46 = V{0.0185185185185185}*x30;
auto x47 = x27 + x28;
auto x48 = x47*x47;
auto x49 = x40 + x42;
auto x50 = x43 + x49;
auto x51 = x27 - x28;
auto x52 = -x51;
auto x53 = -x43;
auto x54 = x49 + x53;
auto x55 = x27 + x29;
auto x56 = x55*x55;
auto x57 = -x45;
auto x58 = -x29;
auto x59 = x27 + x58;
auto x60 = -x59;
auto x61 = x28 + x29;
auto x62 = x61*x61;
auto x63 = x40 + x43;
auto x64 = x45 + x63;
auto x65 = x28 + x58;
auto x66 = -x65;
auto x67 = V{0.00462962962962963}*x30;
auto x68 = x29 + x47;
auto x69 = x68*x68;
auto x70 = x47 + x58;
auto x71 = -x70;
auto x72 = x29 + x51;
auto x73 = -x72;
auto x74 = -x42;
auto x75 = -x27 + x61;
auto x76 = -x36;
auto x77 = V{1} - x38;
auto x78 = x76 + x77;
auto x79 = -x34;
auto x80 = x43 + x79;
auto x81 = x45 + x78 + x80;
auto x82 = x42 + x78;
auto x83 = x45 + x79;
auto x84 = x80 + x82;
auto x85 = x82 + x83;
auto x86 = x40 + x45;
auto x87 = -x75;
auto x0 = -cell[0]*x31 - V{0.296296296296296}*x30*(x32*x40 + V{1});
auto x1 = -cell[1]*x31 - x41*(-x32*(V{3}*x33 - x39 - x42) + V{1});
auto x2 = -cell[2]*x31 - x41*(-x32*(V{3}*x35 - x38 - x43 - x44) + V{1});
auto x3 = -cell[3]*x31 - x41*(-x32*(-x36 + V{3}*x37 - x44 - x45) + V{1});
auto x4 = -cell[4]*x31 - x46*(-x32*(V{4.5}*x48 - x50) + V{1});
auto x5 = -(cell[5]*x31 + x46*(x32*(x54 - V{4.5}*x52*x52) + V{1}));
auto x6 = -cell[6]*x31 - x46*(-x32*(-x45 - x49 + V{4.5}*x56) + V{1});
auto x7 = -(cell[7]*x31 + x46*(x32*(x49 + x57 - V{4.5}*x60*x60) + V{1}));
auto x8 = -cell[8]*x31 - x46*(-x32*(V{4.5}*x62 - x64) + V{1});
auto x9 = -(cell[9]*x31 + x46*(x32*(x57 + x63 - V{4.5}*x66*x66) + V{1}));
auto x10 = -cell[10]*x31 - x67*(-x32*(-x45 - x50 + V{4.5}*x69) + V{1});
auto x11 = -(cell[11]*x31 + x67*(x32*(x50 + x57 - V{4.5}*x71*x71) + V{1}));
auto x12 = -(cell[12]*x31 + x67*(x32*(x45 + x54 - V{4.5}*x73*x73) + V{1}));
auto x13 = -(cell[13]*x31 + x67*(-x32*(x74 + x81 + V{4.5}*(x75*x75)) + V{1}));
auto x14 = -cell[14]*x31 - x41*(-x32*(V{3}*x33 + x82) + V{1});
auto x15 = -cell[15]*x31 - x41*(-x32*(V{3}*x35 + x77 + x80) + V{1});
auto x16 = -cell[16]*x31 - x41*(-x32*(V{3}*x37 + x76 + x83 + V{1}) + V{1});
auto x17 = -cell[17]*x31 - x46*(-x32*(V{4.5}*x48 + x84) + V{1});
auto x18 = -(cell[18]*x31 + x46*(x32*(x63 + x74 - V{4.5}*x51*x51) + V{1}));
auto x19 = -cell[19]*x31 - x46*(-x32*(V{4.5}*x56 + x85) + V{1});
auto x20 = -(cell[20]*x31 + x46*(x32*(x74 + x86 - V{4.5}*x59*x59) + V{1}));
auto x21 = -cell[21]*x31 - x46*(-x32*(V{4.5}*x62 + x81) + V{1});
auto x22 = -(cell[22]*x31 + x46*(x32*(x53 + x86 - V{4.5}*x65*x65) + V{1}));
auto x23 = -cell[23]*x31 - x67*(-x32*(x45 + V{4.5}*x69 + x84) + V{1});
auto x24 = -(cell[24]*x31 + x67*(-x32*(x57 + x84 + V{4.5}*(x70*x70)) + V{1}));
auto x25 = -(cell[25]*x31 + x67*(-x32*(x53 + x85 + V{4.5}*(x72*x72)) + V{1}));
auto x26 = -(cell[26]*x31 + x67*(x32*(x64 + x74 - V{4.5}*x87*x87) + V{1}));
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { V{1}*x32, x33 + x35 + x37 };
}
};

}

}
