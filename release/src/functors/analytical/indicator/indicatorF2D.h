/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Cyril Masquelier, Jan Marquardt, Mathias J. Krause,
 *  Davide Dapelo
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

#ifndef INDICATOR_F_2D_H
#define INDICATOR_F_2D_H

#include<vector>
#include "indicatorBaseF2D.h"
#include "io/xmlReader.h"

#include "core/blockData.h"
#include "core/unitConverter.h"
#include "indicatorBaseF3D.h"

/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Cuboid
 - Circle
 - Triangle

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {

/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
/// theta rotates cuboid around its center, theta in radian measure
template <typename S>
class IndicatorF2DfromIndicatorF3D : public IndicatorF2D<S> {
private:
  IndicatorF3D<S>& _indicator3D;
public:
  IndicatorF2DfromIndicatorF3D(IndicatorF3D<S>& indicator3D);

  bool operator() (bool output[], const S input[]) override;
};

/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
/// theta rotates cuboid around its center, theta in radian measure
template <typename S>
class IndicatorCuboid2D : public IndicatorF2D<S> {
private:
  const Vector<S,2> _center;
  const S _xLength;
  const S _yLength;
  const S _theta;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid2D(Vector<S,2> extend, Vector<S,2> origin, S theta=0);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid2D(S xlength, S ylength, Vector<S,2> center= {S(), S()}, S theta=0);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]) override;
  Vector<S,2> const& getCenter() const;
  S const getxLength() const;
  S const getyLength() const;
  S signedDistance(const Vector<S,2>& input) override;
};

/// indicator function for a 2D NACA airfoil
template <typename S>
class IndicatorAirfoil2D : public IndicatorF2D<S> {
private:
  const Vector<S,2> _center;      // Center point of the airfoil
  const S _chordLength;           // Length from leading to trailing edge
  const S _camber;   // Thickness as percentage of chord length
  const S _camberPos;                // Max Camber
  const S _thicknessPercentage;             // Max Camber Position

  // NACA coefficients for thickness distribution formula
  static constexpr S _coeffs[5] = {0.2969, -0.1260, -0.3516, 0.2843, -0.1036};
  static constexpr S _exponents[5] = {0.5, 1.0, 2.0, 3.0, 4.0};

  // Function to compute the local thickness at a given point along the camber line
  S computeThickness(S x) const;

  S computeCamber(S x) const;

  // Function to transform a point considering rotation
  void transformPoint(S& x, S& y) const;

public:
  IndicatorAirfoil2D(Vector<S,2> center, S chordLength, S camber, S camberPos, S thicknessPercentage);
  Vector<S,2> const& getCenter() const;
  S const getChordLength() const;
  S const getThicknessPercentage() const;
  bool operator() (bool output[], const S input[]) override;
  S signedDistance(const Vector<S,2>& input) override;
};

/// indicator function for a 2D circle
template <typename S>
class IndicatorCircle2D : public IndicatorF2D<S> {
private:
  const Vector<S,2> _center;
  const S _radius;
  const S _radius2;
public:
  IndicatorCircle2D(Vector<S,2> center, S radius);
  Vector<S,2> const& getCenter() const;
  S const getRadius() const;
  bool operator() (bool output[], const S input[]) override;
  bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction,  int iC=-1) override;
  bool normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC=-1) override;
  S signedDistance(const Vector<S,2>& input) override;
  using IndicatorF2D<S>::distance;
};


/// indicator function for a 2D triangle
template <typename S>
class IndicatorTriangle2D : public IndicatorF2D<S> {
private:
  const Vector<S,2> _a;
  const Vector<S,2> _b;
  const Vector<S,2> _c;
public:
  IndicatorTriangle2D(Vector<S,2> a, Vector<S,2> b, Vector<S,2> c);
  Vector<S,2> const& getVertexA() const;
  Vector<S,2> const& getVertexB() const;
  Vector<S,2> const& getVertexC() const;
  bool operator() (bool output[], const S input[]) override;
  S signedDistance(const Vector<S,2>& input) override;
};

/// indicator function for a 2D equilateral triangle
template <typename S>
class IndicatorEquiTriangle2D : public IndicatorF2D<S> {
private:
  const Vector<S,2> _center;
  const S _radius;
  const Vector<S,2> _a;
  const Vector<S,2> _b;
  const Vector<S,2> _c;

public:
  IndicatorEquiTriangle2D(Vector<S,2> center, S radius);
  Vector<S,2> const& getVertexA() const;
  Vector<S,2> const& getVertexB() const;
  Vector<S,2> const& getVertexC() const;
  Vector<S,2> const& getCenter() const;
  S const getRadius() const;
  bool operator() (bool output[], const S input[]) override;
  S signedDistance(const Vector<S,2>& input) override;
};

/// indicator from VTIreader
template <typename S>
class IndicatorBlockData2D : public IndicatorF2D<S> {
private:
  BlockData<3,S,S>& _blockData;
  S const _deltaR;
  bool const _invert;

public:
  IndicatorBlockData2D(BlockData<3,S,S>& blockData,
    Vector<S,3> extend, Vector<S,3> origin, S deltaR, bool invert);
  S signedDistance(const Vector<S,2>& input) override;
};

///indicator from VTIreader
template <typename S>
class IndicatorBlockData2Dvti : public IndicatorF2D<S> {
private:
  BlockData<2,S,S>& _blockData;
  S const _deltaR;
  bool const _invert;

public:
  IndicatorBlockData2Dvti(BlockData<2,S,S>& blockData,
    Vector<S,2> extend, Vector<S,2> origin, S deltaR, bool invert);
  S signedDistance(const Vector<S,2>& input) override;
};

/// Indicator function creating an layer around an input indicator (for positive \p layerSize) or
/// reducing the input indicator by a layer (for negative \p layerSize).
/// \param[in] indicatorF - some indicator (e.g. IndicatorCircle2D)
/// \param[in] layerSize - size of the layer (can be negative) [physical units]
template <typename S>
class IndicatorLayer2D : public IndicatorF2D<S> {
private:
  IndicatorF2D<S>& _indicatorF;
  S _layerSize;
  bool _isPositive;
public:
  IndicatorLayer2D(IndicatorF2D<S>& indicatorF, S layerSize);
  bool operator() (bool output[], const S input[]) override;
};


template <typename T>
class IndicatorSDF2D : public IndicatorF2D<T> {
private:
  std::function<T(Vector<T, 2>)> _f;

public:
  IndicatorSDF2D(std::function<T(Vector<T, 2>)> f);

  bool operator()(bool output[], const T input[]) override;
};

/////////creatorFunctions//////////////////////
template <typename S>
IndicatorCuboid2D<S>* createIndicatorCuboid2D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCircle2D<S>* createIndicatorCircle2D(XMLreader const& params, bool verbose=false);


}

#endif
