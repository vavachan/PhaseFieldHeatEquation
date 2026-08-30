/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert
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

#ifndef ANALYTICAL_PHASE_FIELD_F_2D_H
#define ANALYTICAL_PHASE_FIELD_F_2D_H

#include<vector>
#include<cmath>
#include<string>
#include"math.h"

#include "functors/genericF.h"
#include "analyticalF.h"
#include "core/superLattice2D.h"

/** \file
  This file contains functors needed to describe multiphase phenomena such
  as interfacial profiles and Young-Laplace pressure fields for initialization
  as well as comparison to analytical solutions.
*/


namespace olb {

// hyperbolic tangent between values 0 and factor
template <typename T>
class CircularInterface2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> _center;
  T _interfaceWidth;
  T _radius;
  T _factor;
  bool _bubble;

public:
  CircularInterface2D(std::vector<T> center, T radius, T interfaceWidth, T _factor, bool bubble);
  bool operator()(T output[], const T input[]) override;
};

// smooth Young-Laplace pressure over interface
template <typename T>
class LaplacePressure2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> _center;
  T _interfaceWidth;
  T _radius;
  T _surfaceTension;

public:
  LaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension);
  bool operator()(T output[], const T x[]) override;
};

// shifted Young-Laplace pressure over interface
template <typename T>
class ShiftedLaplacePressure2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> _center;
  T _interfaceWidth;
  T _radius;
  T _surfaceTension;
  T _pressureShift;

public:
  ShiftedLaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension, T pressureShift);
  bool operator()(T output[], const T x[]) override;
};

} // end namespace olb
#endif
