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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x14 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x13 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x16 = parameters.template get<multiphase::RHO_LIQUID>();
auto x15 = parameters.template get<multiphase::RHO_VAPOR>();
auto x11 = -x16;
auto x12 = x11 + x15;
auto x17 = V{1} / (x12);
auto x18 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x19 = x18 + V{1};
auto x20 = x13*(x11 + x19);
auto x21 = x17*x20;
auto x22 = x14*(-x15 + x19);
auto x23 = -x22/x12;
auto x24 = x21 + x23 + V{-1};
auto x25 = x20 - x22;
auto x26 = x17*x25;
auto x27 = x18 + V{1};
auto x28 = V{1} / (x19);
auto x29 = V{1}*cell[3];
auto x30 = V{1}*cell[1] - V{1}*cell[5];
auto x31 = V{0.5}*x10 + x28*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x29 + x30);
auto x32 = x31*x31;
auto x33 = V{1.5}*x32;
auto x34 = V{0.5}*x9;
auto x35 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x29 + x30;
auto x36 = -x35;
auto x37 = x28*x36 + x34;
auto x38 = x33 + V{-1} + V{1.5}*(x37*x37);
auto x39 = V{3}*cell[3];
auto x40 = V{3}*cell[1] - V{3}*cell[5];
auto x41 = V{1.5}*x10 + x28*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x39 + x40);
auto x42 = x10*x41;
auto x43 = V{1.5}*x9;
auto x44 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x39 + x40;
auto x45 = -x28*x44;
auto x46 = x43 + x45;
auto x47 = x46*x9;
auto x48 = V{0.5}*x21 + V{0.5}*x23 + V{-1};
auto x49 = V{4.5}*x9;
auto x50 = V{9}*cell[3];
auto x51 = V{9}*cell[7];
auto x52 = V{9}*cell[1] - V{9}*cell[5];
auto x53 = V{9}*cell[2] - V{9}*cell[6] + x50 - x51 + x52;
auto x54 = x28*x53;
auto x55 = V{3}*x10;
auto x56 = V{6}*cell[3];
auto x57 = V{6}*cell[1] - V{6}*cell[5];
auto x58 = x28*(-V{6}*cell[4] + V{6}*cell[7] + V{6}*cell[8] - x56 + x57);
auto x59 = x55 + x58;
auto x60 = x59 + V{3};
auto x61 = V{6}*cell[2] - V{6}*cell[6] - V{6}*cell[7] + x56 + x57;
auto x62 = V{3}*x9;
auto x63 = -x62;
auto x64 = x63 + V{3};
auto x65 = V{4.5}*x10 + x28*(-V{9}*cell[4] + V{9}*cell[8] - x50 + x51 + x52);
auto x66 = x19*x48;
auto x67 = V{0.0277777777777778}*x66;
auto x68 = -x43;
auto x69 = x28*x44;
auto x70 = x28*x35;
auto x71 = x31 - x34;
auto x72 = V{4.5}*cell[3];
auto x73 = V{4.5}*cell[7];
auto x74 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x75 = V{4.5}*cell[2] - V{4.5}*cell[6] + x72 - x73 + x74;
auto x76 = x28*x75;
auto x77 = V{2.25}*x9;
auto x78 = V{2.25}*x10 + x28*(-V{4.5}*cell[4] + V{4.5}*cell[8] - x72 + x73 + x74);
auto x79 = -x77 + x78;
auto x80 = x34 - x70;
auto x81 = x80*x80;
auto x82 = -x33 - V{1.5}*x81 + V{1};
auto x83 = x41 + x82;
auto x84 = V{0.111111111111111}*x26;
auto x85 = -x75;
auto x86 = x28*x85 + x77;
auto x87 = x38 + x46;
auto x88 = x28*x61;
auto x89 = V{0.111111111111111}*x66;
auto x90 = V{0.0277777777777778}*x26;
auto x91 = x59 + V{-3};
auto x92 = x49 - x54;
auto x93 = x62 - x88;
auto x94 = x65 + V{-3};
auto x95 = x31*x78;
auto x96 = x38 + x41;
auto x97 = -x47;
auto x98 = x93 + V{3};
auto x99 = -x76 + x77;
auto x100 = x43 - x69 + x82;
auto x0 = -cell[0]*x24 + V{0.444444444444444}*x19*x48*(x42 + x47) - V{0.444444444444444}*x26*(x27*x38 + V{1});
auto x1 = -cell[1]*x24 + V{0.0277777777777778}*x17*x25*(x27*(x68 + x69 + x83 + (x70 + x71)*(x76 + x79)) + V{-1}) - x67*(x10*(-x49 + x54 + x60) - x9*(x28*x61 + x64 + x65));
auto x2 = -cell[2]*x24 - x84*(x27*(-x37*x86 + x87) + V{1}) - x89*(-x42 + x9*(-x64 - x88));
auto x3 = -cell[3]*x24 - x67*(x10*(x91 + x92) + x9*(x93 + x94)) - x90*(x27*(x41 + x87 - (x31 + x37)*(x78 + x86)) + V{1});
auto x4 = -cell[4]*x24 - x84*(x27*(-x95 + x96) + V{1}) - x89*(x10*x91 + x97);
auto x5 = -cell[5]*x24 - x67*(-x10*(-x28*x53 + x49 - x55 - x58 + V{3}) + x9*(-x63 - x88 - x94)) - x90*(x27*(-x45 + x68 + x96 - (x28*x36 - x71)*(x28*x85 - x79)) + V{1});
auto x6 = -cell[6]*x24 + V{0.111111111111111}*x17*x25*(x27*(x100 + x80*x99) + V{-1}) - x89*(-x42 + x9*x98);
auto x7 = -cell[7]*x24 + V{0.0277777777777778}*x17*x25*(x27*(x100 + x41 + (x31 + x80)*(x78 + x99)) + V{-1}) - x67*(x10*(x60 + x92) + x9*(x65 + x98));
auto x8 = -cell[8]*x24 + V{0.111111111111111}*x17*x25*(x27*(x83 + x95) + V{-1}) - x89*(x10*x60 + x97);
cell.template getFieldPointer<descriptors::OMEGA>()[0] = x26;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x19, x32 + x81 };
}
};

}

}
