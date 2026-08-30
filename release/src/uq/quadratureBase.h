/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
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

// quadrature_base.h
#ifndef QUADRATURE_BASE_H
#define QUADRATURE_BASE_H

#include "filesIO.h"

namespace olb {

namespace uq {

namespace Quadrature {

enum class QuadratureMethod { HouseholderQR, WilkinsonShiftQR };

template <typename T>
class QuadratureBase {
public:
  virtual ~QuadratureBase()                        = default;
  virtual const std::vector<T>& getPoints() const  = 0;
  virtual const std::vector<T>& getWeights() const = 0;
};

} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // QUADRATURE_BASE_H
