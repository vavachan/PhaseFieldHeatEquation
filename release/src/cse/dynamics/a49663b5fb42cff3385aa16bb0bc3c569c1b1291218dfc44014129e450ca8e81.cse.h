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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = parameters.template get<descriptors::OMEGA>();
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x31 = x30 + V{-1};
auto x32 = x29 + V{-1};
auto x33 = cell[0] + V{2}*cell[10] + V{2}*cell[12] + cell[14] + cell[15] + cell[17] + cell[18] + cell[1] + V{2}*cell[20] + V{2}*cell[22] + V{2}*cell[24] + V{2}*cell[26] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + V{1};
auto x34 = -x33/x32;
auto x35 = x27*x27;
auto x36 = V{1.5}*x35;
auto x37 = x28*x28;
auto x38 = V{1.5}*x37;
auto x39 = x29*x29;
auto x40 = V{1.5}*x39;
auto x41 = x38 + x40 + V{-1};
auto x42 = x36 + x41;
auto x43 = V{0.0740740740740741}*x30;
auto x44 = V{3}*x27;
auto x45 = V{3}*x28;
auto x46 = x36 + V{-1};
auto x47 = V{3}*x29;
auto x48 = V{0.0185185185185185}*x30;
auto x49 = x27 + x28;
auto x50 = x49*x49;
auto x51 = x42 + x44;
auto x52 = x45 + x51;
auto x53 = x27 - x28;
auto x54 = -x53;
auto x55 = -x45;
auto x56 = x51 + x55;
auto x57 = x27 + x29;
auto x58 = x57*x57;
auto x59 = -x47;
auto x60 = -x29;
auto x61 = x27 + x60;
auto x62 = -x61;
auto x63 = x28 + x29;
auto x64 = x63*x63;
auto x65 = x42 + x45;
auto x66 = x47 + x65;
auto x67 = x28 + x60;
auto x68 = -x67;
auto x69 = V{0.00462962962962963}*x30;
auto x70 = x29 + x49;
auto x71 = x70*x70;
auto x72 = x49 + x60;
auto x73 = -x72;
auto x74 = x29 + x53;
auto x75 = -x74;
auto x76 = -x44;
auto x77 = -x27 + x63;
auto x78 = -x38;
auto x79 = V{1} - x40;
auto x80 = x78 + x79;
auto x81 = -x36;
auto x82 = x45 + x81;
auto x83 = x47 + x80 + x82;
auto x84 = x44 + x80;
auto x85 = x47 + x81;
auto x86 = x82 + x84;
auto x87 = x84 + x85;
auto x88 = x42 + x47;
auto x89 = -x77;
auto x0 = -cell[0]*x31 - V{0.296296296296296}*x30*(x34*x42 + V{1});
auto x1 = -cell[1]*x31 - x43*(-x34*(V{3}*x35 - x41 - x44) + V{1});
auto x2 = -cell[2]*x31 - x43*(-x34*(V{3}*x37 - x40 - x45 - x46) + V{1});
auto x3 = -cell[3]*x31 - x43*(-x34*(-x38 + V{3}*x39 - x46 - x47) + V{1});
auto x4 = -cell[4]*x31 - x48*(-x34*(V{4.5}*x50 - x52) + V{1});
auto x5 = -(cell[5]*x31 + x48*(x34*(x56 - V{4.5}*x54*x54) + V{1}));
auto x6 = -cell[6]*x31 - x48*(-x34*(-x47 - x51 + V{4.5}*x58) + V{1});
auto x7 = -(cell[7]*x31 + x48*(x34*(x51 + x59 - V{4.5}*x62*x62) + V{1}));
auto x8 = -cell[8]*x31 - x48*(-x34*(V{4.5}*x64 - x66) + V{1});
auto x9 = -(cell[9]*x31 + x48*(x34*(x59 + x65 - V{4.5}*x68*x68) + V{1}));
auto x10 = -cell[10]*x31 - x69*(-x34*(-x47 - x52 + V{4.5}*x71) + V{1});
auto x11 = -(cell[11]*x31 + x69*(x34*(x52 + x59 - V{4.5}*x73*x73) + V{1}));
auto x12 = -(cell[12]*x31 + x69*(x34*(x47 + x56 - V{4.5}*x75*x75) + V{1}));
auto x13 = -(cell[13]*x31 + x69*(-x34*(x76 + x83 + V{4.5}*(x77*x77)) + V{1}));
auto x14 = -cell[14]*x31 - x43*(-x34*(V{3}*x35 + x84) + V{1});
auto x15 = -cell[15]*x31 - x43*(-x34*(V{3}*x37 + x79 + x82) + V{1});
auto x16 = -cell[16]*x31 - x43*(-x34*(V{3}*x39 + x78 + x85 + V{1}) + V{1});
auto x17 = -cell[17]*x31 - x48*(-x34*(V{4.5}*x50 + x86) + V{1});
auto x18 = -(cell[18]*x31 + x48*(x34*(x65 + x76 - V{4.5}*x53*x53) + V{1}));
auto x19 = -cell[19]*x31 - x48*(-x34*(V{4.5}*x58 + x87) + V{1});
auto x20 = -(cell[20]*x31 + x48*(x34*(x76 + x88 - V{4.5}*x61*x61) + V{1}));
auto x21 = -cell[21]*x31 - x48*(-x34*(V{4.5}*x64 + x83) + V{1});
auto x22 = -(cell[22]*x31 + x48*(x34*(x55 + x88 - V{4.5}*x67*x67) + V{1}));
auto x23 = -cell[23]*x31 - x69*(-x34*(x47 + V{4.5}*x71 + x86) + V{1});
auto x24 = -(cell[24]*x31 + x69*(-x34*(x59 + x86 + V{4.5}*(x72*x72)) + V{1}));
auto x25 = -(cell[25]*x31 + x69*(-x34*(x55 + x87 + V{4.5}*(x74*x74)) + V{1}));
auto x26 = -(cell[26]*x31 + x69*(x34*(x66 + x76 - V{4.5}*x89*x89) + V{1}));
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
return { -V{1}*x33/x32, x35 + x37 + x39 };
}
};

}

}
