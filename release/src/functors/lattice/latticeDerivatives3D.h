/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2015 Patrick Nathen, Mathias J. Krause
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

#ifndef LATTICE_DERIVATIVES_3D_H
#define LATTICE_DERIVATIVES_3D_H

#include <list>

#include "blockBaseF3D.h"
#include "superBaseF3D.h"
#include "core/unitConverter.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"


namespace olb {

/// functor to get pointwise finite difference Dissipation on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T>
class BlockFiniteDifference3D : public BlockF3D<T> {
private:
  BlockGeometry<T,3>& _blockGeometry;
  BlockF3D<T>& _blockFunctor;
  std::list<int>& _matNumber;
  int _targetDim;
  int _n[3];
public:
  BlockFiniteDifference3D(BlockGeometry<T,3>& blockGeometry,
    BlockF3D<T>& blockFunctor,
    std::list<int>& matNumber);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T>
class SuperFiniteDifference3D : public SuperF3D<T> {
private:
  SuperGeometry<T,3>& _sGeometry;
  SuperF3D<T>& _sFunctor;
  std::list<int>& _matNumber;
public:
  SuperFiniteDifference3D(SuperGeometry<T,3>& sGeometry,
    SuperF3D<T>& sFunctor,
    std::list<int>& matNumber);
};

template <typename T, typename DESCRIPTOR>
class BlockPhysFiniteDifference3D : public BlockF3D<T> {
private:
  BlockF3D<T>& _blockFinDiff;
  int _targetDim;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockPhysFiniteDifference3D(BlockF3D<T>& blockFunctor,
    const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperPhysFiniteDifference3D : public SuperF3D<T> {
private:
  SuperFiniteDifference3D<T> _sFinDiff;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperPhysFiniteDifference3D(SuperGeometry<T,3>& sGeometry,
    SuperF3D<T>& sFunctor,
    std::list<int>& matNumber,
    const UnitConverter<T,DESCRIPTOR>& converter);
};


/// functor to get pointwise finite difference Laplacian operator
// for each component of given functor
// uses second order scheme (first order at the boundaries)
template <typename T>
class BlockLaplacian3D : public BlockF3D<T> {
private:
  BlockGeometry<T,3>& _blockGeometry;
  BlockF3D<T>& _blockFunctor;
  int _n[3];
  bool _forthOrder;
public:
  BlockLaplacian3D(BlockGeometry<T,3>& blockGeometry,
    BlockF3D<T>& blockFunctor,
    bool forthOrder);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise finite difference Laplacian operator
template <typename T>
class SuperLaplacian3D : public SuperF3D<T> {
private:
  SuperGeometry<T,3>& _sGeometry;
  SuperF3D<T>& _sFunctor;
public:
  SuperLaplacian3D(SuperGeometry<T,3>& sGeometry,
    SuperF3D<T>& sFunctor,
    bool forthOrder=true);
};

/// functor to get pointwise finite difference Laplacian operator
// scaled with h^-2
template <typename T, typename DESCRIPTOR>
class BlockPhysLaplacian3D : public BlockF3D<T> {
private:
  BlockF3D<T>&                       _blockLaplacian;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  const T                            _factor;
public:
  BlockPhysLaplacian3D(BlockF3D<T>& blockFunctor,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool forthOrder);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise finite difference Laplacian operator
// scaled with h^-2
template <typename T, typename DESCRIPTOR>
class SuperPhysLaplacian3D : public SuperF3D<T> {
private:
  SuperLaplacian3D<T> _laplacian;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperPhysLaplacian3D(SuperGeometry<T,3>& sGeometry,
    SuperF3D<T>& sFunctor,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool forthOrder=true);
};

} // end namespace olb

#endif
