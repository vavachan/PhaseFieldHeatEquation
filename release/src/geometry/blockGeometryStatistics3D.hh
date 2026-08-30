/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny
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

/** \file
 * Representation of a statistic for a 3D geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_STATISTICS_3D_HH
#define BLOCK_GEOMETRY_STATISTICS_3D_HH

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include "utilities/omath.h"

#include "geometry/blockGeometry.h"
#include "geometry/blockGeometryStatistics3D.h"

namespace olb {

template<typename T>
BlockGeometryStatistics3D<T>::BlockGeometryStatistics3D(BlockGeometry<T,3>* blockGeometry)
  : _blockGeometry(blockGeometry),
    clout(std::cout,"BlockGeometryStatistics3D")
{
  _statisticsUpdateNeeded = true;
}

template<typename T>
bool& BlockGeometryStatistics3D<T>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T>
bool const & BlockGeometryStatistics3D<T>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}


template<typename T>
void BlockGeometryStatistics3D<T>::update(bool verbose)
{
  const_this = const_cast<const BlockGeometryStatistics3D<T>*>(this);

  if (getStatisticsStatus() ) {
    _material2n.clear();
    _blockGeometry->forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
      takeStatistics(iX,iY,iZ);
    });

    _nMaterials=int();
    std::map<int, std::size_t>::iterator iter;
    for (iter = _material2n.begin(); iter != _material2n.end(); iter++) {
      _nMaterials++;
    }

    if (verbose) {
      clout << "updated" << std::endl;
    }
    getStatisticsStatus()=false;
  }
}


template<typename T>
int BlockGeometryStatistics3D<T>::getNmaterials()
{
  update();
  return const_this->getNmaterials();
}

template<typename T>
int BlockGeometryStatistics3D<T>::getNmaterials() const
{
  return _nMaterials;
}

template<typename T>
std::size_t BlockGeometryStatistics3D<T>::getNvoxel(int material)
{
  update();
  return const_this->getNvoxel(material);
}

template<typename T>
std::size_t BlockGeometryStatistics3D<T>::getNvoxel(int material) const
{
  try {
    return _material2n.at(material);
  }
  catch (std::out_of_range& ex) {
    return 0;
  }
}

template<typename T>
std::map<int, std::size_t> BlockGeometryStatistics3D<T>::getMaterial2n()
{
  update();
  return const_this->getMaterial2n();
}

template<typename T>
std::map<int, std::size_t> BlockGeometryStatistics3D<T>::getMaterial2n() const
{
  return _material2n;
}

template<typename T>
std::size_t BlockGeometryStatistics3D<T>::getNvoxel()
{
  update();
  return const_this->getNvoxel();
}

template<typename T>
std::size_t BlockGeometryStatistics3D<T>::getNvoxel() const
{
  std::size_t total = 0;
  for (const auto& material : _material2n) {
    total += material.second;
  }
  return total;
}
template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMinLatticeR(int material)
{
  update();
  return const_this->getMinLatticeR(material);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMinLatticeR(int material) const
{
  try {
    return _material2min.at(material);
  }
  catch (std::out_of_range& ex) {
    std::vector<int> null;
    return null;
  }
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMaxLatticeR(int material)
{
  update();
  return const_this->getMaxLatticeR(material);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getMaxLatticeR(int material) const
{
  try {
    return _material2max.at(material);
  }
  catch (std::out_of_range& ex) {
    std::vector<int> null;
    return null;
  }
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getMinPhysR(int material) const
{
  std::vector<T> tmp(3,T());
  Vector<T,3> physR;
  // _blockGeometry->getPhysR(&(tmp[0]), &(getMinLatticeR(material)[0]));
  _blockGeometry->getPhysR(physR, &(getMinLatticeR(material)[0]));
  tmp[0] = physR[0];
  tmp[1] = physR[1];
  tmp[2] = physR[2];
  return tmp;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getMaxPhysR(int material) const
{
  std::vector<T> tmp(3,T());
  Vector<T,3> physR;
  // _blockGeometry->getPhysR(&(tmp[0]), &(getMaxLatticeR(material)[0]));
  _blockGeometry->getPhysR(physR, &(getMaxLatticeR(material)[0]));
  tmp[0] = physR[0];
  tmp[1] = physR[1];
  tmp[2] = physR[2];
  return tmp;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getLatticeExtend(int material)
{
  update();
  return const_this->getLatticeExtend(material);
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getLatticeExtend(int material) const
{
  try {
    std::vector<T> extend;
    for (int iDim = 0; iDim < 3; iDim++) {
      extend.push_back(_material2max.at(material)[iDim] - _material2min.at(material)[iDim]);
    }
    return extend;
  }
  catch (std::out_of_range& ex) {
    std::vector<T, std::allocator<T>> null;
    return null;
  }
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysExtend(int material)
{
  update();
  return const_this->getPhysExtend(material);
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysExtend(int material) const
{
  std::vector<T> extend;
  for (int iDim = 0; iDim < 3; iDim++) {
    extend.push_back(getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim]);
  }
  return extend;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysRadius(int material)
{
  update();
  return const_this->getPhysRadius(material);
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getPhysRadius(int material) const
{
  std::vector<T> radius;
  for (int iDim=0; iDim<3; iDim++) {
    radius.push_back((getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.);
  }
  return radius;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getCenterPhysR(int material)
{
  update();
  return const_this->getCenterPhysR(material);
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::getCenterPhysR(int material) const
{
  std::vector<T> center;
  for (int iDim=0; iDim<3; iDim++) {
    center.push_back(getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim]);
  }
  return center;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getType(int iX, int iY, int iZ)
{
  update();
  return const_this->getType(iX, iY, iZ);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getType(const int* input) const
{
  return const_this->getType(input[0], input[1], input[2]);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getType(int iX, int iY, int iZ,
                                                       BlockIndicatorF3D<T>& fluidI,
                                                       BlockIndicatorF3D<T>& outsideI) const
{
  const auto [normalType, normal] = computeBoundaryTypeAndNormal(fluidI, outsideI, {iX,iY,iZ});
  return {static_cast<int>(normalType), normal[0], normal[1], normal[2]};
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::getType(int iX, int iY, int iZ) const
{
  BlockIndicatorMaterial3D<T> fluidI(*_blockGeometry, 1);
  BlockIndicatorMaterial3D<T> outsideI(*_blockGeometry, 0);
  return getType(iX,iY,iZ,fluidI,outsideI);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeNormal(int iX, int iY, int iZ)
{
  update();
  return const_this->computeNormal(iX, iY, iZ);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeNormal(int iX, int iY, int iZ) const
{
  std::vector<int> normal (3,int(0));

  if (iX != 0) {
    if (_blockGeometry->getMaterial({iX - 1, iY, iZ}) == 1) {
      normal[0] = -1;
    }
  }
  if (iX != _nX - 1) {
    if (_blockGeometry->getMaterial({iX + 1, iY, iZ}) == 1) {
      normal[0] = 1;
    }
  }
  if (iY != 0) {
    if (_blockGeometry->getMaterial({iX, iY - 1, iZ}) == 1) {
      normal[1] = -1;
    }
  }
  if (iY != _nY - 1) {
    if (_blockGeometry->getMaterial({iX, iY + 1, iZ}) == 1) {
      normal[1] = 1;
    }
  }
  if (iZ != 0) {
      if (_blockGeometry->getMaterial({iX, iY, iZ - 1}) == 1) {
      normal[2] = -1;
    }
  }
  if (iZ != _nZ - 1) {
      if (_blockGeometry->getMaterial({iX, iY, iZ + 1}) == 1) {
      normal[2] = 1;
    }
  }
  return normal;
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::computeNormal(int material)
{

  update();
  return const_this->computeNormal(material);
}

template<typename T>
std::vector<T> BlockGeometryStatistics3D<T>::computeNormal(int material) const
{
  std::vector<T> normal (3,int(0));
  std::vector<int> minC = getMinLatticeR(material);
  std::vector<int> maxC = getMaxLatticeR(material);
  for (int iX = minC[0]; iX<=maxC[0]; iX++) {
    for (int iY = minC[1]; iY<=maxC[1]; iY++) {
      for (int iZ = minC[2]; iZ<=maxC[2]; iZ++) {
        if (_blockGeometry->getMaterial({iX,iY,iZ}) == material) {
          auto n = computeNormal(iX,iY,iZ);
          normal[0]+=n[0];
          normal[1]+=n[1];
          normal[2]+=n[2];
        }
      }
    }
  }
  T norm = util::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if (norm>0.) {
    normal[0]/=norm;
    normal[1]/=norm;
    normal[2]/=norm;
  }
  return normal;
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeDiscreteNormal(int material, T maxNorm)
{
  update();
  return const_this->computeDiscreteNormal(material, maxNorm);
}

template<typename T>
std::vector<int> BlockGeometryStatistics3D<T>::computeDiscreteNormal(int material, T maxNorm) const
{
  std::vector<T> normal = computeNormal(material);
  std::vector<int> discreteNormal(3,int(0));

  T smallestAngle = T(0);
  for (int iX = -1; iX<=1; iX++) {
    for (int iY = -1; iY<=1; iY++) {
      for (int iZ = -1; iZ<=1; iZ++) {
        T norm = util::sqrt(iX*iX+iY*iY+iZ*iZ);
        if (norm>0.&& norm<maxNorm) {
          T angle = (iX*normal[0] + iY*normal[1] + iZ*normal[2])/norm;
          if (angle>=smallestAngle) {
            smallestAngle=angle;
            discreteNormal[0] = iX;
            discreteNormal[1] = iY;
            discreteNormal[2] = iZ;
          }
        }
      }
    }
  }
  return discreteNormal;
}

template<typename T>
bool BlockGeometryStatistics3D<T>::check(int material, int iX, int iY,
    int iZ, unsigned offsetX, unsigned offsetY, unsigned offsetZ)
{
  update();
  return const_this->check(material, iX, iY, iZ, offsetX, offsetY, offsetZ);
}

template<typename T>
bool BlockGeometryStatistics3D<T>::check(int material, int iX, int iY,
    int iZ, unsigned offsetX, unsigned offsetY, unsigned offsetZ) const
{
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      for (int iOffsetZ = -offsetZ; iOffsetZ <= (int) offsetZ; ++iOffsetZ) {
        if (_blockGeometry->getMaterial({iX + iOffsetX, iY + iOffsetY,
          iZ + iOffsetZ}) != material) {
          found = false;
        }
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStatistics3D<T>::find(int material, unsigned offsetX,
                                        unsigned offsetY, unsigned offsetZ, int& foundX, int& foundY,
                                        int& foundZ)
{
  update();
  return const_this->find(material, offsetX, offsetY, offsetZ, foundX, foundY, foundZ);
}

template<typename T>
bool BlockGeometryStatistics3D<T>::find(int material, unsigned offsetX,
                                        unsigned offsetY, unsigned offsetZ, int& foundX, int& foundY,
                                        int& foundZ) const
{
  bool found = false;
  for (foundX = 0; foundX < _nX; foundX++) {
    for (foundY = 0; foundY < _nY; foundY++) {
      for (foundZ = 0; foundZ < _nZ; foundZ++) {
        found = check(material, foundX, foundY, foundZ, offsetX,
                      offsetY, offsetZ);
        if (found) {
          return found;
        }
      }
    }
  }
  return found;
}

template<typename T>
void BlockGeometryStatistics3D<T>::print()
{
  update();
  return const_this->print();
}

template<typename T>
void BlockGeometryStatistics3D<T>::print() const
{
  try {
    for (const auto& material : _material2n) {
      clout << "materialNumber=" << material.first
            << "; count=" << material.second
            << "; minLatticeR=(" << _material2min.at(material.first)[0] <<","
            << _material2min.at(material.first)[1] <<","<< _material2min.at(material.first)[2] <<")"
            << "; maxLatticeR=(" << _material2max.at(material.first)[0] <<","
            << _material2max.at(material.first)[1] <<","<< _material2max.at(material.first)[2] <<")"
            << std::endl;
    }
  }
  catch (std::out_of_range& ex) {
  }
}

template<typename T>
void BlockGeometryStatistics3D<T>::takeStatistics(int iX, int iY, int iZ)
{
  int type = _blockGeometry->getMaterial({iX, iY, iZ});
  if (_material2n.count(type) == 0) {
    _material2n[type] = 1;
    std::vector<int> minCo;
    std::vector<int> maxCo;
    minCo.push_back(iX);
    minCo.push_back(iY);
    minCo.push_back(iZ);
    _material2min[type] = minCo;
    maxCo.push_back(iX);
    maxCo.push_back(iY);
    maxCo.push_back(iZ);
    _material2max[type] = maxCo;

  }
  else {
    _material2n[type]++;
    if (iX < _material2min[type][0]) {
      _material2min[type][0] = iX;
    }
    if (iY < _material2min[type][1]) {
      _material2min[type][1] = iY;
    }
    if (iZ < _material2min[type][2]) {
      _material2min[type][2] = iZ;
    }
    if (iX > _material2max[type][0]) {
      _material2max[type][0] = iX;
    }
    if (iY > _material2max[type][1]) {
      _material2max[type][1] = iY;
    }
    if (iZ > _material2max[type][2]) {
      _material2max[type][2] = iZ;
    }
  }
}

} // namespace olb

#endif
