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

#ifndef SUPER_INDICATOR_F_3D_HH
#define SUPER_INDICATOR_F_3D_HH

#include <numeric>

#include "superIndicatorF3D.h"
#include "blockIndicatorF3D.h"
#include "core/util.h"

namespace olb {

template <typename T>
SuperIndicatorFfromIndicatorF3D<T>::SuperIndicatorFfromIndicatorF3D(
  FunctorPtr<IndicatorF3D<T>>&& indicatorF, SuperGeometry<T,3>& geometry)
  : SuperIndicatorF3D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromIndicatorF3D<T>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
  }
}

template <typename T>
bool SuperIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  auto physR = this->_superStructure.getCuboidDecomposition().getPhysR(input);
  return _indicatorF(output, physR.data());
}


template <typename T, bool HLBM>
SuperIndicatorFfromSmoothIndicatorF3D<T, HLBM>::SuperIndicatorFfromSmoothIndicatorF3D(
  FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>>&& indicatorF,
  SuperGeometry<T,3>&                     geometry)
  : SuperIndicatorF3D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
  }
}

template <typename T, bool HLBM>
bool SuperIndicatorFfromSmoothIndicatorF3D<T, HLBM>::operator() (bool output[], const int input[])
{
  T inside[1];
  auto physR = this->_superStructure.getCuboidDecomposition().getPhysR(input);
  _indicatorF(inside, physR.data());
  return !util::nearZero(inside[0]);
}

template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry<T,3>& geometry, std::vector<int> materials)
  : SuperIndicatorF3D<T>(geometry)
{
  geometry.updateStatistics(false);
  const std::string matString = std::accumulate(
                                  materials.begin()+1,
                                  materials.end(),
                                  std::to_string(materials[0]),
  [](const std::string& a, int b) {
    return a + '_' + std::to_string(b);
  });
  this->getName() = "SuperIndicator_on_Material_" + matString;

  for (int iC = 0; iC < this->_superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMaterial3D<T>(this->_superGeometry.getBlockGeometry(iC),
                                      materials)
    );
  }
}

template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry<T,3>& geometry, std::list<int> materials)
  : SuperIndicatorMaterial3D(geometry,
                             std::vector<int>(materials.begin(), materials.end()))
{ }

template <typename T>
bool SuperIndicatorMaterial3D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    // query material number of appropriate block indicator
    return this->getBlockF(load.loc(input[0]))(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
SuperIndicatorFieldThreshold3D<T,DESCRIPTOR,FIELD>::SuperIndicatorFieldThreshold3D(
  SuperGeometry<T,3>& geometry, SuperLattice<T, DESCRIPTOR>& superLattice, std::vector<int> materials, T thresholdValue, std::string condition)
  : SuperIndicatorF3D<T>(geometry),
    _superLattice(superLattice),
    _thresholdValue(thresholdValue),
    _condition(condition)
{
  geometry.updateStatistics(false);
  for (int iC=0; iC < this->_superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFieldThreshold3D<T,DESCRIPTOR,FIELD>(this->_superGeometry.getBlockGeometry(iC), this->_superLattice.getBlock(iC), materials, _thresholdValue, _condition)
    );
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
bool SuperIndicatorFieldThreshold3D<T,DESCRIPTOR,FIELD>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    // query material number of appropriate block indicator
    return this->getBlockF(load.loc(input[0]))(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T>
SuperIndicatorLayer3D<T>::SuperIndicatorLayer3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIndicatorF3D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorLayer3D<T>(_indicatorF->getBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorLayer3D<T>::operator()(bool output[], const int input[])
{
  _indicatorF(output, input);
  for (int iPop=1; iPop < descriptors::D3Q27<>::q; ++iPop) {
    bool tmpOutput{};
    Vector<int,4> tmpInput(input);
    auto c_i = descriptors::c<descriptors::D3Q27<>>(iPop);
    tmpInput += Vector<int,4>{0, c_i[0], c_i[1], c_i[2]};
    _indicatorF(&tmpOutput, tmpInput.data());
    output[0] |= tmpOutput;
  }
  return true;
}

template <typename T>
SuperIndicatorIdentity3D<T>::SuperIndicatorIdentity3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIndicatorF3D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorIdentity3D<T>(_indicatorF->getBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorIdentity3D<T>::operator()(bool output[], const int input[])
{
  return _indicatorF(output, input);
}

template <typename T>
SuperIndicatorMultiplication3D<T>::SuperIndicatorMultiplication3D(
  FunctorPtr<SuperIndicatorF3D<T>>&& f,
  FunctorPtr<SuperIndicatorF3D<T>>&& g)
  : SuperIndicatorF3D<T>(f->getSuperGeometry()),
    _f(std::move(f)), _g(std::move(g))
{
  this->getName() = _f->getName() + " * " + _g->getName();

  for (int iC = 0; iC < _f->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMultiplication3D<T>(
        _f->getBlockIndicatorF(iC),
        _g->getBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorMultiplication3D<T>::operator()(
  bool output[], const int input[])
{
  _f(output, input);
  if (output[0]) {
    _g(output, input);
  }
  return output[0];
}

template <typename T>
SuperIndicatorSubstraction3D<T>::SuperIndicatorSubstraction3D(
  FunctorPtr<SuperIndicatorF3D<T>>&& f,
  FunctorPtr<SuperIndicatorF3D<T>>&& g)
  : SuperIndicatorF3D<T>(f->getSuperGeometry()),
    _f(std::move(f)), _g(std::move(g))
{
  this->getName() = _f->getName() + " - " + _g->getName();

  for (int iC = 0; iC < _f->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorSubstraction3D<T>(
        _f->getBlockIndicatorF(iC),
        _g->getBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorSubstraction3D<T>::operator()(
  bool output[], const int input[])
{
  _f(output, input);
  if (output[0]) {
    bool tmp { };
    _g(&tmp, input);
    output[0] &= !tmp;
  }
  return output[0];
}

template <typename T>
SuperIndicatorBoundaryNeighbor3D<T>::SuperIndicatorBoundaryNeighbor3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF, int overlap)
  : SuperIndicatorF3D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF)),
    _overlap(overlap)
{
  this->getName() = "SuperIndicatorBoundaryNeighbor_for_BoundaryMaterial_" + _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorBoundaryNeighbor3D<T>(_indicatorF->getBlockIndicatorF(iC), _overlap));
  }
}

template <typename T>
bool SuperIndicatorBoundaryNeighbor3D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (load.isLocal(input[0])) {
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

} // namespace olb

#endif
