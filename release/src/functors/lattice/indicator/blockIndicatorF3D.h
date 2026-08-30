/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_INDICATOR_F_3D_H
#define BLOCK_INDICATOR_F_3D_H

#include <list>

#include "blockIndicatorBaseF3D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {

/// BlockIndicatorF3D from IndicatorF3D
template <typename T>
class BlockIndicatorFfromIndicatorF3D : public BlockIndicatorF3D<T> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param indicatorF    Indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromIndicatorF3D(IndicatorF3D<T>&    indicatorF,
                                  BlockGeometry<T,3>& blockGeometry);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator helper for people that don't like writing boilerplate
template <typename T>
class BlockIndicatorFfromCallableF3D : public BlockIndicatorF3D<T> {
protected:
  std::function<bool(LatticeR<3>)> _callableF;

public:
  template <typename F>
  BlockIndicatorFfromCallableF3D(BlockGeometry<T,3>& blockGeometry, F&& f)
  : BlockIndicatorF3D<T>(blockGeometry)
  , _callableF(f)
  { }

  using BlockIndicatorF3D<T>::operator();

  bool operator() (bool output[], const int input[]) override {
    output[0] = _callableF(input);
    return true;
  }

  Vector<int,3> getMin() override {
    return std::numeric_limits<int>::min();
  }
  Vector<int,3> getMax() override {
    return std::numeric_limits<int>::max();
  }
};

/// BlockIndicatorF3D from SmoothIndicatorF3D
template <typename T, bool HLBM>
class BlockIndicatorFfromSmoothIndicatorF3D : public BlockIndicatorF3D<T> {
protected:
  SmoothIndicatorF3D<T,T,HLBM>& _indicatorF;
public:
  /**
   * \param indicatorF    Smooth indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromSmoothIndicatorF3D(SmoothIndicatorF3D<T,T,HLBM>& indicatorF,
                                        BlockGeometry<T,3>&           blockGeometry);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};


/// Block indicator functor from material numbers
template <typename T>
class BlockIndicatorMaterial3D : public BlockIndicatorF3D<T> {
protected:
  const std::vector<int> _materials;
public:
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param materials     Material number vector
   **/
  BlockIndicatorMaterial3D(BlockGeometry<T,3>& blockGeometry,
                           std::vector<int>    materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param materials     Material number list
   **/
  BlockIndicatorMaterial3D(BlockGeometry<T,3>& blockGeometry,
                           std::list<int>      materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param material      Material number
   **/
  BlockIndicatorMaterial3D(BlockGeometry<T,3>& blockGeometry,
                           int                 material);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns true iff indicated domain subset is empty
  bool isEmpty() override;
  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator threshold for external field
template <typename T, typename DESCRIPTOR, typename FIELD>
class BlockIndicatorFieldThreshold3D : public BlockIndicatorF3D<T> {
protected:
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  const std::vector<int> _materials;
  T _thresholdValue;
  std::string _condition;
public:
  BlockIndicatorFieldThreshold3D(BlockGeometry<T,3>& blockGeometry, BlockLattice<T,DESCRIPTOR>& blockLattice, std::vector<int> materials, T thresholdValue, std::string condition);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
  /// Returns true iff indicated domain subset is empty
  bool isEmpty() override;
  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator extended by a layer
template <typename T>
class BlockIndicatorLayer3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _indicatorF;

public:
  /**
   * \param indicatorF Block indicator to be extended by a layer
   **/
  BlockIndicatorLayer3D(BlockIndicatorF3D<T>& indicatorF);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator identity
template <typename T>
class BlockIndicatorIdentity3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Block indicator to be proxied
   **/
  BlockIndicatorIdentity3D(BlockIndicatorF3D<T>& indicatorF);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator intersection
template <typename T>
class BlockIndicatorMultiplication3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _f;
  BlockIndicatorF3D<T>& _g;
public:
  /**
   * \param f Block indicator to be multiplied
   * \param g Block indicator to be multiplied
   **/
  BlockIndicatorMultiplication3D(BlockIndicatorF3D<T>& f, BlockIndicatorF3D<T>& g);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator substraction
template <typename T>
class BlockIndicatorSubstraction3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _f;
  BlockIndicatorF3D<T>& _g;
public:
  BlockIndicatorSubstraction3D(BlockIndicatorF3D<T>& f, BlockIndicatorF3D<T>& g);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

/// Block indicator identifying neighbors of boundary cells
template <typename T>
class BlockIndicatorBoundaryNeighbor3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _indicatorF;
  int _overlap;
public:
  BlockIndicatorBoundaryNeighbor3D(BlockIndicatorF3D<T>& indicatorF, int overlap);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

} // namespace olb

#endif
