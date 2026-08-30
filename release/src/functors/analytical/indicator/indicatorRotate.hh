/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Michael Grinschewski
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

#ifndef INDICATOR_ROTATE_HH
#define INDICATOR_ROTATE_HH

#include <vector>
#include "utilities/omath.h"
#include <cassert>
#include <sstream>

#include "indicatorF3D.h"
#include "indicatorF2D.h"
#include "utilities/vectorHelpers.h"

namespace olb{

// Constructor for 2D Indicator
template <typename S, unsigned D>
IndicatorRotate<S,D>::IndicatorRotate(Vector<S,2> rotationPoint, S rotationAngle, IndicatorF<S,D>& indicator)
:  _rotationPoint(rotationPoint), _rotationAngle(rotationAngle), _indicator(indicator), _min(indicator.getMin()), _max(indicator.getMax())
{}

// Constructor for 3D Indicator
template <typename S, unsigned D>
IndicatorRotate<S,D>::IndicatorRotate(Vector<S,3> rotationPoint, Vector<S,3> rotationAxis, S rotationAngle, IndicatorF<S,D>& indicator)
: _rotationPoint(rotationPoint), _rotationAxis(rotationAxis), _rotationAngle(rotationAngle), _indicator(indicator), _min(indicator.getMin()), _max(indicator.getMax())
{}

// Rotate method for 2D
template <typename S, unsigned D>
bool IndicatorRotate<S,D>::rotate(S output[], const S input[], Vector<S,2> rotationPoint, S rotationAngle )
{
  // Unpack input coordinates
  S x = input[0], y = input[1];
  S x_r = _rotationPoint[0], y_r = _rotationPoint[1];
  // Translate point to origin
  Vector<S, 2> p = {x - x_r, y - y_r};
  // Precompute trigonometric functions
  S cos_theta = util::cos(_rotationAngle);
  S sin_theta = util::sin(_rotationAngle);
  // Compute rotation around the origin
  Vector<S, 2> p_rot = {p[0] * cos_theta - p[1] * sin_theta, p[1] * cos_theta + p[0] * sin_theta};
  // Translate point back to its proper place
  output[0] = p_rot[0] + x_r;
  output[1] = p_rot[1] + y_r;
  return true;
}

// 3D rotate function
template <typename S, unsigned D>
bool IndicatorRotate<S,D>::rotate(S output[], const S input[], Vector<S,3> rotationPoint, Vector<S,3> rotationAxis, S rotationAngle )
{
  // Unpack input coordinates
  S x = input[0], y = input[1], z = input[2];
  S x_r = _rotationPoint[0], y_r = _rotationPoint[1], z_r = _rotationPoint[2];
  S a = _rotationAxis[0], b = _rotationAxis[1], c = _rotationAxis[2];

  // Translate point to origin
  Vector<S, 3> p = {x - x_r, y - y_r, z - z_r};

  // Normalize the rotation axis
  S norm = util::sqrt(a * a + b * b + c * c);
  Vector<S, 3> u = {a / norm, b / norm, c / norm};

  // Precompute trigonometric functions
  S cos_theta = util::cos(_rotationAngle);
  S sin_theta = util::sin(_rotationAngle);

  // Compute cross product u × p
  Vector<S, 3> cross_product = crossProduct(u,p);

  // Compute dot product u ⋅ p
  S dot_product = u[0] * p[0] + u[1] * p[1] + u[2] * p[2];

  // Rodrigues' rotation formula
  Vector<S, 3> p_rot = {
    p[0] * cos_theta + cross_product[0] * sin_theta + u[0] * dot_product * (1 - cos_theta),
    p[1] * cos_theta + cross_product[1] * sin_theta + u[1] * dot_product * (1 - cos_theta),
    p[2] * cos_theta + cross_product[2] * sin_theta + u[2] * dot_product * (1 - cos_theta)
  };

  output[0] = p_rot[0] + x_r;
  output[1] = p_rot[1] + y_r;
  output[2] = p_rot[2] + z_r;

  return true;
}

// getMin function split by dimensions
template <typename S, unsigned D>
Vector<S,D>& IndicatorRotate<S,D>::getMin()
{
  if constexpr (D == 2){
    S minRot[2] {};
    const S inMin[2] = {_indicator.getMin()[0], _indicator.getMin()[1]};
    IndicatorRotate<S,D>::rotate(minRot, inMin, _rotationPoint, _rotationAngle);
    S maxRot[2] {};
    const S inMax[2] = {_indicator.getMax()[0], _indicator.getMax()[1]};
    IndicatorRotate<S,D>::rotate(maxRot, inMax, _rotationPoint, _rotationAngle);
    _min[0] = util::min(minRot[0],maxRot[0]);
    _min[1] = util::min(minRot[1],maxRot[1]);
    return _min;
  }
  else if constexpr (D == 3){
    S minRot[3] {};
    const S inMin[3] = {_indicator.getMin()[0], _indicator.getMin()[1], _indicator.getMin()[2]};
    IndicatorRotate<S,D>::rotate(minRot, inMin, _rotationPoint, _rotationAxis, _rotationAngle);
    S maxRot[3] {};
    const S inMax[3] = {_indicator.getMax()[0], _indicator.getMax()[1], _indicator.getMax()[2]};
    IndicatorRotate<S,D>::rotate(maxRot, inMax, _rotationPoint, _rotationAxis, _rotationAngle);
    _min[0] = util::min(minRot[0],maxRot[0]);
    _min[1] = util::min(minRot[1],maxRot[1]);
    _min[2] = util::min(minRot[2],maxRot[2]);
    return _min;
  }
}

template <typename S, unsigned D>
Vector<S,D>& IndicatorRotate<S,D>::getMax()
{
  if constexpr (D == 2){
    S minRot[2] {};
    const S inMin[2] = {_indicator.getMin()[0], _indicator.getMin()[1]};
    IndicatorRotate<S,D>::rotate(minRot, inMin, _rotationPoint, _rotationAngle);
    S maxRot[2] {};
    const S inMax[2] = {_indicator.getMax()[0], _indicator.getMax()[1]};
    IndicatorRotate<S,D>::rotate(maxRot, inMax, _rotationPoint, _rotationAngle);
    _max[0] = util::max(minRot[0],maxRot[0]);
    _max[1] = util::max(minRot[1],maxRot[1]);
    return _max;
  }
  else if constexpr (D == 3){
    S minRot[3] {};
    const S inMin[3] = {_indicator.getMin()[0], _indicator.getMin()[1], _indicator.getMin()[2]};
    IndicatorRotate<S,D>::rotate(minRot, inMin, _rotationPoint, _rotationAxis, _rotationAngle);
    S maxRot[3] {};
    const S inMax[3] = {_indicator.getMax()[0], _indicator.getMax()[1], _indicator.getMax()[2]};
    IndicatorRotate<S,D>::rotate(maxRot, inMax, _rotationPoint, _rotationAxis, _rotationAngle);
    _max[0] = util::max(minRot[0],maxRot[0]);
    _max[1] = util::max(minRot[1],maxRot[1]);
    _max[2] = util::max(minRot[2],maxRot[2]);
    return _max;
  }
}

// operator function
template <typename S, unsigned D>
bool IndicatorRotate<S,D>::operator() (bool output[], const S input[] )
{
  if constexpr (D == 2){
    S inputTranslated[2] {};
    IndicatorRotate<S,D>::rotate(inputTranslated, input, _rotationPoint, _rotationAngle);
    _indicator.operator()(output, inputTranslated);
    return output[0];
  }
  else if constexpr (D == 3){
    S inputTranslated[3] {};
    IndicatorRotate<S,D>::rotate(inputTranslated, input, _rotationPoint, _rotationAxis, _rotationAngle);
    _indicator.operator()(output, inputTranslated);
    return output[0];
  }
}

}

#endif
