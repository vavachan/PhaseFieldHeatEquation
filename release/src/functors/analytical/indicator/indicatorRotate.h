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

#ifndef INDICATOR_ROTATE_H
#define INDICATOR_ROTATE_H

#include <vector>

#include "indicatorBaseF3D.h"
#include "indicatorBaseF2D.h"
#include "io/xmlReader.h"
#include "utilities/functorPtr.h"

/* This file contains Content for Rotation of Indicators */

namespace olb {

template <typename S, unsigned D>
class IndicatorRotate : public IndicatorF<S,D> {
private:
  Vector<S,D> _rotationPoint;
  Vector<S,D> _rotationAxis;
  S _rotationAngle;
  IndicatorF<S,D>& _indicator;
  Vector<S,D> _min;
  Vector<S,D> _max;
public:
  // Constructs 2D Indicator that is rotated around the rotation point by the specified angle
  IndicatorRotate(Vector<S,2> rotationPoint, S rotationAngle, IndicatorF<S,D>& indicator);
  // Constructs 3D Indicator that is rotated around the rotation axis through the point by the specified angle
  IndicatorRotate(Vector<S,3> rotationPoint, Vector<S,3> rotationAxis, S rotationAngle, IndicatorF<S,D>& indicator);
  bool rotate(S output[], const S input[], Vector<S,2> rotationPoint, S rotationAngle );
  bool rotate(S output[], const S input[], Vector<S,3> rotationPoint, Vector<S,3> rotationAxis, S rotationAngle );
  bool operator() (bool output[], const S input[]) override;
  Vector<S,D>& getMin() override;
  Vector<S,D>& getMax() override;
};

}

#endif
