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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = parameters.template get<descriptors::OMEGA>();
auto x24 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x22 = cell.template getFieldComponent<descriptors::FORCE>(2);
auto x19 = cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0);
auto x25 = x23 + V{-1};
auto x26 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x27 = x26 + V{1};
auto x28 = x20*x20;
auto x29 = V{1.5}*x28;
auto x30 = x21*x21;
auto x31 = V{1.5}*x30;
auto x32 = x22*x22;
auto x33 = V{1.5}*x32;
auto x34 = x31 + x33 + V{-1};
auto x35 = x29 + x34;
auto x36 = x26 + V{1};
auto x37 = -V{3}*x19*x24 + x36;
auto x38 = V{0.0555555555555556}*x23;
auto x39 = V{3}*x20;
auto x40 = V{3}*x28;
auto x41 = x37*x38;
auto x42 = V{3}*x21;
auto x43 = V{3}*x30;
auto x44 = x29 + V{-1};
auto x45 = V{3}*x22;
auto x46 = V{3}*x32;
auto x47 = V{0.0277777777777778}*x23;
auto x48 = x20 + x21;
auto x49 = V{4.5}*(x48*x48);
auto x50 = x35 + x39;
auto x51 = x37*x47;
auto x52 = -x42;
auto x53 = x20 - x21;
auto x54 = -x53;
auto x55 = x20 + x22;
auto x56 = V{4.5}*(x55*x55);
auto x57 = -x45;
auto x58 = -x22;
auto x59 = x20 + x58;
auto x60 = -x59;
auto x61 = x21 + x22;
auto x62 = V{4.5}*(x61*x61);
auto x63 = x35 + x42;
auto x64 = x21 + x58;
auto x65 = -x64;
auto x66 = -x31;
auto x67 = V{1} - x33;
auto x68 = x66 + x67;
auto x69 = x39 + x68;
auto x70 = -x29;
auto x71 = x42 + x70;
auto x72 = x45 + x70;
auto x73 = -x39;
auto x74 = x35 + x45;
auto x0 = -cell[0]*x25 + V{0.666666666666667}*x23*x37 - V{0.333333333333333}*x23*(x27*x35 + V{1});
auto x1 = -cell[1]*x25 - x38*(x27*(x34 + x39 - x40) + V{1}) - x41;
auto x2 = -cell[2]*x25 - x38*(x27*(x33 + x42 - x43 + x44) + V{1}) - x41;
auto x3 = -cell[3]*x25 - x38*(x27*(x31 + x44 + x45 - x46) + V{1}) - x41;
auto x4 = -cell[4]*x25 - x47*(x27*(x42 - x49 + x50) + V{1}) - x51;
auto x5 = -(cell[5]*x25 + x47*(x27*(x50 + x52 - V{4.5}*x54*x54) + V{1}) + x51);
auto x6 = -cell[6]*x25 - x47*(x27*(x45 + x50 - x56) + V{1}) - x51;
auto x7 = -(cell[7]*x25 + x47*(x27*(x50 + x57 - V{4.5}*x60*x60) + V{1}) + x51);
auto x8 = -cell[8]*x25 - x47*(x27*(x45 - x62 + x63) + V{1}) - x51;
auto x9 = -(cell[9]*x25 + x47*(x27*(x57 + x63 - V{4.5}*x65*x65) + V{1}) + x51);
auto x10 = -cell[10]*x25 + V{0.0555555555555556}*x23*(x27*(x40 + x69) + V{-1}) - x41;
auto x11 = -cell[11]*x25 + V{0.0555555555555556}*x23*(x27*(x43 + x67 + x71) + V{-1}) - x41;
auto x12 = -cell[12]*x25 + V{0.0555555555555556}*x23*(x27*(x46 + x66 + x72 + V{1}) + V{-1}) - x41;
auto x13 = -cell[13]*x25 + V{0.0277777777777778}*x23*(x27*(x49 + x69 + x71) + V{-1}) - x51;
auto x14 = -(cell[14]*x25 + x47*(x27*(x63 + x73 - V{4.5}*x53*x53) + V{1}) + x51);
auto x15 = -cell[15]*x25 + V{0.0277777777777778}*x23*(x27*(x56 + x69 + x72) + V{-1}) - x51;
auto x16 = -(cell[16]*x25 + x47*(x27*(x73 + x74 - V{4.5}*x59*x59) + V{1}) + x51);
auto x17 = -cell[17]*x25 + V{0.0277777777777778}*x23*(x27*(x45 + x62 + x68 + x71) + V{-1}) - x51;
auto x18 = -(cell[18]*x25 + x47*(x27*(x52 + x74 - V{4.5}*x64*x64) + V{1}) + x51);
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
return { x36, x28 + x30 + x32 };
}
};

}

}
