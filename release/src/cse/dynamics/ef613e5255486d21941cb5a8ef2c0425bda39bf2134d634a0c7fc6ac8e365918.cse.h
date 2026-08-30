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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, guoZhao::GuoZhaoSecondOrder, collision::BGK, guoZhao::GuoZhaoForcing<momenta::GuoZhaoForced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x14 = cell.template getFieldComponent<descriptors::K>(0);
auto x9 = cell.template getFieldComponent<descriptors::BODY_FORCE>(0);
auto x15 = cell.template getFieldComponent<descriptors::NU>(0);
auto x16 = parameters.template get<descriptors::OMEGA>();
auto x11 = cell.template getFieldComponent<descriptors::EPSILON>(0);
auto x10 = cell.template getFieldComponent<descriptors::BODY_FORCE>(1);
auto x12 = x16 + V{-1};
auto x13 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x17 = x13 + V{1};
auto x18 = V{1} / (x11);
auto x19 = V{0.5}*x11;
auto x20 = x15/x14;
auto x21 = x19*x20 + V{1};
auto x22 = util::sqrt(x21*x21);
auto x23 = V{1} / ((x21 + x22)*(x21 + x22));
auto x24 = x13 + V{1};
auto x25 = V{1} / (x24);
auto x26 = V{1}*cell[7];
auto x27 = V{1}*cell[3];
auto x28 = V{1}*cell[1] - V{1}*cell[5];
auto x29 = -V{1}*cell[4] + V{1}*cell[8] + x26 - x27 + x28;
auto x30 = x25*x29;
auto x31 = x10*x19 + x30;
auto x32 = x31*x31;
auto x33 = V{6}*x32;
auto x34 = V{1}*cell[2] - V{1}*cell[6] - x26 + x27 + x28;
auto x35 = x25*x34;
auto x36 = x19*x9;
auto x37 = -x35 + x36;
auto x38 = x37*x37;
auto x39 = -x18*x23*(x33 + V{6}*x38) + V{1};
auto x40 = x24*(V{0.5}*x16 + V{-1});
auto x41 = V{1} / (V{0.25}*x11*x20 + V{0.5}*x22 + V{0.5});
auto x42 = V{1}*x20;
auto x43 = x10 - x31*x41*x42;
auto x44 = -x37*x41*x42 + x9;
auto x45 = V{1.5}*x11;
auto x46 = x45*x9;
auto x47 = -V{3}*x35 + x46;
auto x48 = x41*x47;
auto x49 = V{18}*x18;
auto x50 = x31 - x36;
auto x51 = x35 + x50;
auto x52 = x10*x45 + V{3}*x30;
auto x53 = x41*x52;
auto x54 = x39 + x53;
auto x55 = V{0.0277777777777778}*x40;
auto x56 = x41*x52;
auto x57 = -x56;
auto x58 = V{3}*x11;
auto x59 = V{9}*x35;
auto x60 = V{4.5}*x11;
auto x61 = x60*x9;
auto x62 = x10*x60 + V{9}*x30;
auto x63 = -x61 + x62;
auto x64 = x41*(x59 + x63) + x58;
auto x65 = -x25*x34;
auto x66 = x46 + V{3}*x65;
auto x67 = x41*x66;
auto x68 = V{0.111111111111111}*x16;
auto x69 = x36 + x65;
auto x70 = x69*x69;
auto x71 = x18*x23*(x33 + V{6}*x70) + V{-1};
auto x72 = x41*x66 + x71;
auto x73 = -x59 + x61;
auto x74 = x58 + x67;
auto x75 = -x41*x43*(V{0.166666666666667}*x10*x11 + V{0.333333333333333}*x25*x29);
auto x76 = x31 + x69;
auto x77 = V{9}*x65;
auto x78 = -x41*(x61 + x62 + x77);
auto x79 = x56 + x58;
auto x80 = x23*x32*x49;
auto x81 = -x41*x44*(V{0.166666666666667}*x11*x9 - V{0.333333333333333}*x25*x34);
auto x82 = x50 - x65;
auto x83 = x39 + x48;
auto x84 = x41*(x63 - x77);
auto x85 = x41*x47;
auto x86 = x58 - x85;
auto x87 = x31 + x37;
auto x88 = x41*(x62 + x73);
auto x89 = x57 + x58;
auto x0 = -cell[0]*x12 + V{0.444444444444444}*x16*(x17*x39 + V{-1}) + V{1.33333333333333}*x40*x41*(x31*x43 + x37*x44);
auto x1 = -(cell[1]*x12 - V{0.0277777777777778}*x16*(x17*(x23*x49*(x51*x51) - x48 + x54) + V{-1}) + x55*(x43*(x57 + x64) - x44*(x64 + x67)));
auto x2 = -cell[2]*x12 - x40*(V{0.111111111111111}*x44*(x41*x73 - x74) + x75) - x68*(x17*(-x23*x49*x70 + x72) + V{1});
auto x3 = -(cell[3]*x12 + V{0.0277777777777778}*x16*(x17*(-x23*x49*x76*x76 + x53 + x72) + V{1}) + x55*(x43*(-x78 - x79) + x44*(-x74 - x78)));
auto x4 = -cell[4]*x12 - x40*(V{0.111111111111111}*x43*(x41*x62 - x79) + x81) - x68*(x17*(x53 + x71 - x80) + V{1});
auto x5 = -(cell[5]*x12 - V{0.0277777777777778}*x16*(x17*(x23*x49*(x82*x82) - x53 + x83) + V{-1}) + x55*(x43*(-x79 + x84) + x44*(x58 - x84 - x85)));
auto x6 = -cell[6]*x12 + V{0.111111111111111}*x16*(x17*(x23*x38*x49 + x83) + V{-1}) - x40*(V{0.111111111111111}*x44*(x41*x73 + x86) + x75);
auto x7 = -(cell[7]*x12 - V{0.0277777777777778}*x16*(x17*(x23*x49*(x87*x87) + x53 + x83) + V{-1}) + x55*(x43*(x88 + x89) + x44*(x86 + x88)));
auto x8 = -cell[8]*x12 + V{0.111111111111111}*x16*(x17*(x54 + x80) + V{-1}) - x40*(V{0.111111111111111}*x43*(x41*x62 + x89) + x81);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x24, V{4}*x23*(x32 + x38) };
}
};

}

}
