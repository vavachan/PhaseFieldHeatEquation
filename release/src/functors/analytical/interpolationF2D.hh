/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2018 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef INTERPOLATION_F_2D_HH
#define INTERPOLATION_F_2D_HH

#include "interpolationF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbm.h"

namespace olb {

/// Bilinear interpolation for rectangular lattice with dimensions delta[i];
/// This functor performs bilinear interpolation for 2D data.
template <typename T, typename W>
SpecialAnalyticalFfromBlockF2D<T,W>::SpecialAnalyticalFfromBlockF2D(
  BlockF2D<W>& f, Cuboid2D<T>& cuboid,
  Vector<T,2> delta, T scale)
  : AnalyticalF2D<T,W>(f.getTargetDim()), _f(f), _cuboid(cuboid), _delta(delta), _scale(scale)
{
  this->getName() = "fromBlockF";
}

template <typename T, typename W>
bool SpecialAnalyticalFfromBlockF2D<T,W>::operator()(W output[], const T physC[])
{
  Vector<T,2> origin = _cuboid.getOrigin();

  // Scale physC in all 2 dimensions
  Vector<T,2> physCv;
  for (int i=0; i<2; i++) {
    physCv[i] = origin[i] + (physC[i] - origin[i]) * ( _cuboid.getDeltaR() / _delta[i] );
  }

  int latticeR[2];
  for (int i=0; i<2; i++) {
    latticeR[i] = util::max((int)util::floor( (physCv[i] - origin[i]) / _cuboid.getDeltaR() ), 0);
  }
  Vector<T,2> physRiC;
  Vector<W,2> d, e;
  W output_tmp[3]; // Assuming maximum of 3 components
  Vector<T,2> latticeRv;

  for (int i=0; i<2; i++) {
    latticeRv[i] = (T) latticeR[i];
  }
  physRiC = origin + latticeRv * _cuboid.getDeltaR();
  T dr = 1. / _cuboid.getDeltaR();

  // Compute weights
  d = (physCv - physRiC) * dr;
  e = 1. - d;

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
    output_tmp[iD] = W();
  }

  // Perform bilinear interpolation
  // Corner (0,0)
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0] * e[1];
  }

  // Corner (0,1)
  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0] * d[1];
  }

  // Corner (1,0)
  if (_cuboid.getNy() != 1) {
    latticeR[1]--;
  }
  if (_cuboid.getNx() != 1) {
    latticeR[0]++;
  }
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0] * e[1];
  }

  // Corner (1,1)
  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0] * d[1];
  }

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] *= _scale;
  }

  return true;
}

template <typename T, typename W>
AnalyticalFfromBlockF2D<T,W>::AnalyticalFfromBlockF2D(
  BlockF2D<W>& f, Cuboid2D<T>& cuboid)
  : AnalyticalF2D<T,W>(f.getTargetDim()),
    _f(f), _cuboid(cuboid)
{
  this->getName() = "fromBlockF";
}

template <typename T, typename W>
bool AnalyticalFfromBlockF2D<T,W>::operator()(W output[], const T physC[])
{
  int latticeC[2];
  Vector<T,2> physCv(physC);
  LatticeR<2> latticeR = _cuboid.getFloorLatticeR(physCv);

  auto& block = _f.getBlockStructure();
  auto padding = std::min(1, block.getPadding());

  if (LatticeR<2>(latticeR) >= -padding && LatticeR<2>(latticeR) < block.getExtent()+padding-1) {
    const int& locX = latticeR[0];
    const int& locY = latticeR[1];

    Vector<T,2> physRiC = _cuboid.getPhysR({locX, locY});

    // compute weights
    Vector<W,2> d = (physCv - physRiC) * (1. / _cuboid.getDeltaR());
    Vector<W,2> e = 1. - d;

    T output_tmp[_f.getTargetDim()];
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output_tmp[iD] = T();
    }

    latticeC[0] = locX;
    latticeC[1] = locY;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * e[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX;
    latticeC[1] = locY + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * d[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * e[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * d[1];
      output_tmp[iD] = T();
    }

    return true;
  }
  else {
    return false;
  }
}

template <typename T, typename W>
AnalyticalFfromSuperF2D<T,W>::AnalyticalFfromSuperF2D(SuperF2D<T>& f,
    bool communicateToAll, bool communicateOverlap)
  : AnalyticalF2D<T,W>(f.getTargetDim()),
    _communicateToAll(communicateToAll),
    _communicateOverlap(communicateOverlap),
    _f(f),
    _cuboidDecomposition(_f.getSuperStructure().getCuboidDecomposition())
{
  this->getName() = "fromSuperF("+ f.getName()+")";

  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new AnalyticalFfromBlockF2D<T>(_f.getBlockF(iC),
                                     _cuboidDecomposition.get(load.glob(iC)))
    );
  }
}

template <typename T, typename W>
bool AnalyticalFfromSuperF2D<T,W>::operator() (T output[], const T physC[])
{
  const auto targetDim = _f.getTargetDim();
  for (int iD = 0; iD < targetDim; ++iD) {
    output[iD] = W();
  }
  Vector<T,2> physCV(physC);
  auto latticeR = _cuboidDecomposition.getLatticeR(physCV);
  if (!latticeR) {
    return false;
  }

  if (_communicateOverlap) {
    _f.getSuperStructure().communicate();
  }

  int dataFound = 0;

  const LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    if (_blockF[iC]->operator()(output, physC)) {
      ++dataFound;
    }
  }

  if (_communicateToAll) {
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(dataFound, MPI_SUM);
    for (int iD = 0; iD < targetDim; ++iD) {
      singleton::mpi().reduceAndBcast(output[iD], MPI_SUM);
    }
#endif
    for (int iD = 0; iD < targetDim; ++iD) {
      output[iD]/=dataFound;
    }
  }
  else {
    if (dataFound!=0) {
      for (int iD = 0; iD < targetDim; ++iD) {
        output[iD]/=dataFound;
      }
    }
  }

  if (dataFound>0) {
    return true;
  }
  return false;
}

template <typename T, typename W>
int AnalyticalFfromSuperF2D<T,W>::getBlockFSize() const
{
  OLB_ASSERT(_blockF.size() < UINT32_MAX,
             "it is safe to cast std::size_t to int");
  return _blockF.size();
}

template <typename T, typename W>
AnalyticalFfromBlockF2D<T,W>& AnalyticalFfromSuperF2D<T,W>::getBlockF(int iCloc)
{
  OLB_ASSERT(size_t(iCloc) < _blockF.size() && iCloc >= 0,
             "block functor index within bounds");
  return *(_blockF[iCloc]);
}


} // end namespace olb

#endif
