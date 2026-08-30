/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_OPERATOR_SCOPE_H
#define CORE_OPERATOR_SCOPE_H

namespace olb {

/// Block-wide operator application scopes
/**
 * Declares how the actual OPERATOR::apply template wants to be called.
 **/
enum struct OperatorScope {
  /// Per-cell application, i.e. OPERATOR::apply is passed a CELL concept implementation
  PerCell,
  /// Per-block application, i.e. OPERATOR::apply is passed a ConcreteBlockLattice
  PerBlock,
  /// Per-cell application with parameters, i.e. OPERATOR::apply is passed a CELL concept implementation and parameters
  PerCellWithParameters,
};

/// Returns human-readable name of scope
std::string getName(OperatorScope scope) {
  switch (scope) {
  case OperatorScope::PerCell:
    return "PerCell";
  case OperatorScope::PerBlock:
    return "PerBlock";
  case OperatorScope::PerCellWithParameters:
    return "PerCellWithParameters";
  default:
    throw std::invalid_argument("OperatorScope value is not defined");
  }
}

}

#endif
