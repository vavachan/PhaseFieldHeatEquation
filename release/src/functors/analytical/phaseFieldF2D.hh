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

#ifndef ANALYTICAL_PHASE_FIELD_F_2D_HH
#define ANALYTICAL_PHASE_FIELD_F_2D_HH

#include<vector>
#include<cmath>
#include<string>

#include "phaseFieldF2D.h"
#include "functors/genericF.h"
#include "analyticalF.h"
#include "functors/lattice/superBaseF2D.h"
#include "geometry/superGeometry.h"

#include "core/superLattice2D.h"
#include "utilities/vectorHelpers.h"  // for normalize


namespace olb {

template <typename T>
CircularInterface2D<T>::CircularInterface2D(std::vector<T> center, T radius, T interfaceWidth, T factor, bool bubble) : AnalyticalF2D<T,T>(1)
{
  this->getName() = "CircularInterface2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _factor = factor;
  _bubble = bubble;
}

template <typename T>
bool CircularInterface2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  if (_bubble) output[0] = T(_factor*(1.+tanh(2*d/_interfaceWidth))/2.);
  else output[0] = T(_factor*(1.-tanh(2*d/_interfaceWidth))/2.);
  return true;
}


template <typename T>
LaplacePressure2D<T>::LaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension) : AnalyticalF2D<T,T>(1)
{
  this->getName() = "LaplacePressure2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _surfaceTension = surfaceTension;
}

template <typename T>
bool LaplacePressure2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  T tanh_d = tanh(2*d/_interfaceWidth);
  output[0] = _surfaceTension/(_radius*4)*tanh_d*(tanh_d*tanh_d-3.);
  return true;
}


template <typename T>
ShiftedLaplacePressure2D<T>::ShiftedLaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension, T pressureShift): AnalyticalF2D<T,T>(1)
{
  this->getName() = "ShiftedLaplacePressure2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _surfaceTension = surfaceTension;
  _pressureShift = pressureShift;
}

template <typename T>
bool ShiftedLaplacePressure2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  T tanh_d = tanh(2*d/_interfaceWidth);
  output[0] = _surfaceTension/(_radius*4)*tanh_d*(tanh_d*tanh_d-3.) + _pressureShift;
  return true;
}


} // end namespace olb
#endif
