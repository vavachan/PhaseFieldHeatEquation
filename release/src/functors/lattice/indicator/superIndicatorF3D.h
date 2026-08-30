/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2018 Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_INDICATOR_F_3D_H
#define SUPER_INDICATOR_F_3D_H

#include <list>

#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "superIndicatorBaseF3D.h"
#include "geometry/superGeometry.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {


/// SuperIndicatorF3D from IndicatorF3D
template <typename T>
class SuperIndicatorFfromIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<IndicatorF3D<T>> _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromIndicatorF3D(FunctorPtr<IndicatorF3D<T>>&& indicatorF,
                                  SuperGeometry<T,3>&           geometry);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// SuperIndicatorF3D from SmoothIndicatorF3D
/**
 * Returns true iff SmoothIndicatorF3D output is not near zero.
 **/
template <typename T, bool HLBM>
class SuperIndicatorFfromSmoothIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>> _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromSmoothIndicatorF3D(FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>>&& indicatorF,
                                        SuperGeometry<T,3>&                     geometry);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator functor from material numbers
/**
 * Material number comparison is implemented in BlockIndicatorMaterial3D.
 **/
template <typename T>
class SuperIndicatorMaterial3D : public SuperIndicatorF3D<T> {
public:
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials Vector of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry<T,3>& geometry, std::vector<int> materials);
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials List of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry<T,3>& geometry, std::list<int> materials);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};

/// threshold for external field
template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperIndicatorFieldThreshold3D : public SuperIndicatorF3D<T> {
private:
SuperLattice<T, DESCRIPTOR>& _superLattice;
T _thresholdValue;
std::string _condition;

public:
  SuperIndicatorFieldThreshold3D(SuperGeometry<T,3>& geometry, SuperLattice<T, DESCRIPTOR>& superLattice, std::vector<int> materials, T thresholdValue, std::string condition);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};

/// Indicator extended by a layer
template <typename T>
class SuperIndicatorLayer3D : public SuperIndicatorF3D<T> {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;

public:
  SuperIndicatorLayer3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};

/// Indicator identity functor
/**
 * Proxies a given indicator functor for simplified memory management.
 * i.e. mixing non-owning indicator references and owning indicator pointers
 *      in functor compositions.
 **/
template <typename T>
class SuperIndicatorIdentity3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperIndicatorIdentity3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};

/// Indicator intersection functor
template <typename T>
class SuperIndicatorMultiplication3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _f;
  FunctorPtr<SuperIndicatorF3D<T>> _g;
public:
  SuperIndicatorMultiplication3D(FunctorPtr<SuperIndicatorF3D<T>>&& f,
                                 FunctorPtr<SuperIndicatorF3D<T>>&& g);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};

/// Indicator substraction functor
template <typename T>
class SuperIndicatorSubstraction3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _f;
  FunctorPtr<SuperIndicatorF3D<T>> _g;
public:
  SuperIndicatorSubstraction3D(FunctorPtr<SuperIndicatorF3D<T>>&& f,
                               FunctorPtr<SuperIndicatorF3D<T>>&& g);

  using SuperIndicatorF3D<T>::operator();

  bool operator() (bool output[], const int input[]) override;
};

/// Indicator identifying neighbors of boundary cells
/**
 * Checks if neighbor is part of a boundary via passed SuperIndicatorMaterial.
 **/
template <typename T>
class SuperIndicatorBoundaryNeighbor3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
  int _overlap;
public:
  SuperIndicatorBoundaryNeighbor3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF, int overlap);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};
} // namespace olb

#endif
