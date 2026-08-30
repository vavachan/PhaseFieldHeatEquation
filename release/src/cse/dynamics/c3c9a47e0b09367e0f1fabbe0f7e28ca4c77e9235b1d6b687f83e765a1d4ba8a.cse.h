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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x12 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x13 = cell[2] + cell[3];
auto x14 = cell[7] + cell[8];
auto x15 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + x13 + x14;
auto x16 = x15 + V{1};
auto x17 = V{1} / (x16);
auto x18 = cell[1] - cell[5];
auto x19 = -cell[6] - cell[7] + x13 + x18;
auto x20 = x17*(x15 + V{1});
auto x21 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x22 = V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8] - x17*x19*x19 - x19*x20*x9 + x21;
auto x23 = -cell[3] - cell[4] + x14 + x18;
auto x24 = -V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8] + x10*x20*x23 - x17*x23*x23 + x21;
auto x25 = V{1}*cell[3];
auto x26 = V{1}*cell[1];
auto x27 = x20*(-x10*x19 + x23*x9);
auto x28 = x17*x19*x23;
auto x29 = V{1}*cell[7];
auto x30 = -V{1}*cell[5];
auto x31 = x29 + x30;
auto x32 = V{1} / (V{2.52268963608289}*util::sqrt(x17*(x12*x12)*util::sqrt((-x25 + x26 - V{0.5}*x27 - V{1}*x28 - x31)*(V{2}*cell[1] - V{2}*cell[3] + V{2}*cell[5] - V{2}*cell[7] - V{1}*x27 - V{2}*x28) + V{1}*(x22*x22) + V{1}*(x24*x24)) + V{0.0392836979096202}/((x11)*(x11))) + V{0.5}/x11);
auto x33 = V{1}*cell[8];
auto x34 = V{1}*cell[4];
auto x35 = V{0.5}*x10 + x17*(-x25 + x26 + x31 + x33 - x34);
auto x36 = x35*x35;
auto x37 = V{1.5}*x36;
auto x38 = V{0.5}*x9;
auto x39 = V{1}*cell[2] - V{1}*cell[6] + x25 + x26 - x29 + x30;
auto x40 = -x39;
auto x41 = x17*x40 + x38;
auto x42 = x37 + V{-1} + V{1.5}*(x41*x41);
auto x43 = V{1} - x32;
auto x44 = V{3}*cell[3];
auto x45 = V{3}*cell[1] - V{3}*cell[5];
auto x46 = V{1.5}*x10 + x17*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x44 + x45);
auto x47 = x10*x46;
auto x48 = V{1.5}*x9;
auto x49 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x44 + x45;
auto x50 = -x17*x49;
auto x51 = x48 + x50;
auto x52 = V{1} - V{0.5}*x32;
auto x53 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x54 = -x48;
auto x55 = x17*x49;
auto x56 = x17*x39;
auto x57 = x35 - x38;
auto x58 = V{4.5}*cell[3];
auto x59 = V{4.5}*cell[7];
auto x60 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x61 = V{4.5}*cell[2] - V{4.5}*cell[6] + x58 - x59 + x60;
auto x62 = x17*x61;
auto x63 = V{2.25}*x9;
auto x64 = V{2.25}*x10 + x17*(-V{4.5}*cell[4] + V{4.5}*cell[8] - x58 + x59 + x60);
auto x65 = -x63 + x64;
auto x66 = x38 - x56;
auto x67 = x66*x66;
auto x68 = -x37 - V{1.5}*x67 + V{1};
auto x69 = x46 + x68;
auto x70 = V{4.5}*x9;
auto x71 = V{9}*cell[3];
auto x72 = V{9}*cell[7];
auto x73 = V{9}*cell[1] - V{9}*cell[5];
auto x74 = V{9}*cell[2] - V{9}*cell[6] + x71 - x72 + x73;
auto x75 = x17*x74;
auto x76 = V{3}*x10;
auto x77 = V{6}*cell[7];
auto x78 = V{6}*cell[3];
auto x79 = V{6}*cell[1] - V{6}*cell[5];
auto x80 = x17*(-V{6}*cell[4] + V{6}*cell[8] + x77 - x78 + x79);
auto x81 = x76 + x80;
auto x82 = x81 + V{3};
auto x83 = V{3}*x9;
auto x84 = x17*(V{6}*cell[2] - V{6}*cell[6] - x77 + x78 + x79);
auto x85 = -x83 + x84;
auto x86 = x85 + V{3};
auto x87 = V{4.5}*x10 + x17*(-V{9}*cell[4] + V{9}*cell[8] - x71 + x72 + x73);
auto x88 = V{0.0277777777777778}*x16;
auto x89 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x90 = -x61;
auto x91 = x17*x90 + x63;
auto x92 = x42 + x51;
auto x93 = V{0.111111111111111}*x16;
auto x94 = x81 + V{-3};
auto x95 = x70 - x75;
auto x96 = x83 - x84;
auto x97 = x87 + V{-3};
auto x98 = x35*x64;
auto x99 = x42 + x46;
auto x100 = x48 - x55;
auto x101 = -x100*x9;
auto x102 = -x62 + x63;
auto x103 = x100 + x68;
auto x104 = x96 + V{3};
auto x0 = V{1}*cell[0]*x43 - V{0.444444444444444}*x16*x52*(x47 + x51*x9) - x32*(x42*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444});
auto x1 = x26*x43 + x32*(x53*(x54 + x55 + x69 + (x56 + x57)*(x62 + x65)) + V{-0.0277777777777778}) + x52*x88*(x10*(-x70 + x75 + x82) - x9*(x86 + x87));
auto x2 = V{1}*cell[2]*x43 - x32*(x89*(-x41*x91 + x92) + V{0.111111111111111}) - x52*x93*(x47 + x86*x9);
auto x3 = x25*x43 - x32*(x53*(x46 + x92 - (x35 + x41)*(x64 + x91)) + V{0.0277777777777778}) + x52*x88*(x10*(x94 + x95) + x9*(x96 + x97));
auto x4 = -x32*(x89*(-x98 + x99) + V{0.111111111111111}) + x34*x43 + x52*x93*(x10*x94 + x101);
auto x5 = V{1}*cell[5]*x43 - x32*(x53*(-x50 + x54 + x99 - (x17*x40 - x57)*(x17*x90 - x65)) + V{0.0277777777777778}) - x52*x88*(x10*(-x17*x74 + x70 - x76 - x80 + V{3}) - x9*(-x85 - x97));
auto x6 = V{1}*cell[6]*x43 + x32*(x89*(x102*x66 + x103) + V{-0.111111111111111}) - x52*x93*(-x104*x9 + x47);
auto x7 = x29*x43 + x32*(x53*(x103 + x46 + (x102 + x64)*(x35 + x66)) + V{-0.0277777777777778}) + x52*x88*(x10*(x82 + x95) + x9*(x104 + x87));
auto x8 = x32*(x89*(x69 + x98) + V{-0.111111111111111}) + x33*x43 + x52*x93*(x10*x82 + x101);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x16, x36 + x67 };
}
};

}

}
