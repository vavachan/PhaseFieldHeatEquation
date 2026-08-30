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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FreeEnergyInletOutletDensity, momenta::FreeEnergyInletOutletMomentum<1, -1>, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately>, equilibria::FreeEnergy, collision::FreeEnergyInletOutlet<1, -1>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(1);
auto x9 = cell.template getFieldComponent<descriptors::FORCE>(0);
auto x13 = parameters.template get<descriptors::OMEGA>();
auto x11 = x13 + V{-1};
auto x12 = cell[0]*x11;
auto x14 = x10*x10;
auto x15 = x9*(V{1.5}*x14 + V{-1}) + V{1};
auto x16 = V{0.444444444444444}*x13*x15;
auto x17 = V{0.111111111111111}*x13;
auto x18 = x15*x17;
auto x19 = cell[2]*x11;
auto x20 = cell[3]*x11;
auto x21 = cell[4]*x11;
auto x22 = cell[5]*x11;
auto x23 = cell[6]*x11;
auto x24 = V{0.0277777777777778}*x13;
auto x25 = x9*(-V{3}*x10 + V{3}*x14 + V{1}) + V{-1};
auto x26 = -x24*x25;
auto x27 = x13*x9;
auto x28 = V{0.166666666666667}*x12 + x18 + V{0.166666666666667}*x19 + V{0.166666666666667}*x20 + V{0.166666666666667}*x21 + V{0.166666666666667}*x22 + V{0.166666666666667}*x23 + x26 - V{0.0277777777777778}*x27 + V{0.166666666666667}*x9 + V{-0.166666666666667};
auto x29 = x17*x9;
auto x30 = x18 + x29;
auto x31 = x24*x9 + x26;
auto x32 = -x17*x25;
cell[0] = -x12 + V{0.555555555555556}*x13*x9 - x16;
cell[1] = x28;
cell[2] = -x19 - x30;
cell[3] = -x20 - x31;
cell[4] = -x21 - x29 - x32;
cell[5] = -x22 - x31;
cell[6] = -x23 - x30;
cell[7] = x28;
cell[8] = V{0.666666666666667}*x12 + x16 + V{0.666666666666667}*x19 + V{0.666666666666667}*x20 + V{0.666666666666667}*x21 + V{0.666666666666667}*x22 + V{0.666666666666667}*x23 - V{0.111111111111111}*x27 + x32 + V{0.666666666666667}*x9 + V{-0.666666666666667};
return { x9, x14 };
}
};

}

}
