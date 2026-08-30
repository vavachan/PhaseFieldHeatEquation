/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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
 * Input in STL format -- header file. - nope
 */

#ifndef BLOCK_LATTICE_STL_READER_HH
#define BLOCK_LATTICE_STL_READER_HH


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "octree.hh"
#include "blockLatticeSTLreader.h"


// All OpenLB code is contained in this namespace.
namespace olb {



/*
 * BlockLatticeSTLreader functions
 */
template<typename T>
BlockLatticeSTLreader<T>::BlockLatticeSTLreader(CuboidDecomposition3D<T>& cbg3d, LoadBalancer<T>& hlb, const std::string fName, T voxelSize, T stlSize,
                        int method, bool verbose, T overlap, T max)
  : _cuboidDecomposition (cbg3d),
    _loadBalancer(hlb),
    _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName(fName),
    _mesh(fName, stlSize),
    _verbose(verbose),
    clout(std::cout, "BlockLatticeSTLreader")
{
  this->getName() = "BlockLatticeSTLreader";

  if (_verbose) {
    clout << "Voxelizing ..." << std::endl;
  }

  Vector<T,3> extension = _mesh.getMax() - _mesh.getMin();
  if ( util::nearZero(max) ) {
    max = util::max(extension[0], util::max(extension[1], extension[2])) + _voxelSize;
  }
  int j = 0;
  for (; _voxelSize * util::pow(2, j) < max; j++)
    ;
  Vector<T,3> center;
  T radius = _voxelSize * util::pow(2, j - 1);

  /// Find center of tree and move by _voxelSize/4.
  for (unsigned i = 0; i < 3; i++) {
    center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize / 4.;
  }

  /// Create tree
  _tree = new Octree<T>(center, radius, &_mesh, j, _overlap);

  /// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel inside the STL.
  for (int i = 0; i < 3; i++) {
    this->_myMin[i] = center[i] + _voxelSize / 2.;
    this->_myMax[i] = center[i] - _voxelSize / 2.;
  }
  for (int i = 0; i < 3; i++) {
    while (this->_myMin[i] > _mesh.getMin()[i]) {
      this->_myMin[i] -= _voxelSize;
    }
    while (this->_myMax[i] < _mesh.getMax()[i]) {
      this->_myMax[i] += _voxelSize;
    }
    this->_myMax[i] -= _voxelSize;
    this->_myMin[i] += _voxelSize;
  }

  /// Indicate nodes of the tree. (Inside/Outside)
  switch (method) {
  case 1:
    indicate1();
    break;
  case 3:
    indicate3();
    break;
  case 4:
    indicate2_Xray();
    break;
  case 5:
    indicate2_Yray();
    break;
  default:
    indicate2();
    break;
  }

  if (_verbose) {
    print();
  }

  clout <<"creating of the lists for signed distance function"<< std::endl;
  _trianglesInCuboidList.resize(_loadBalancer.size());
  //_cuboidDecomposition.size();
  for( int iC=0; iC < _loadBalancer.size();iC++)
  {
    //_trainglesInCuboidList.emplace_back();
    int globC = _loadBalancer.glob(iC);
    auto& cuboid = _cuboidDecomposition.get(globC);
  for (STLtriangle<T>& triangle : _mesh.getTriangles())
  {
    Vector<T,3> center = triangle.getCenter() ;
/* <<<<<<< HEAD
    if (cuboid.intersects(center[0], center[1], center[2], 3))
======= */
    //TODO!!! For the signed distance function to work correctly, the overlapp (6) in checkInters must be higher than the real overlapp to capture all neighbouring triangles in the overlap correctly
    if (cuboid.isInside(center, 6))

    {
      _trianglesInCuboidList[iC].push_back(triangle);

    }
  }
  }
 // clout << "create neighbouring triangels removed" << std::endl;
 createNeighbouringTriangleInCuboidVector();
 clout <<'\n' << '\n' << "WARNING! The sensitivity in signed distance function could be insufficient for your case!!!" << '\n' << std::endl;

  if (_verbose) {

    clout << "Voxelizing ... OK" << std::endl;
  }
}

template<typename T>
void BlockLatticeSTLreader<T>::createNeighbouringTriangleInCuboidVector()
{
  T sensitivity = 1.e-30;
  _neighbouringTriangleInCuboidList.resize(_loadBalancer.size());
  clout << "createNeighbouringTriangleInCuboidVector called" << std::endl;
  for(int iC=0; iC < _loadBalancer.size();iC++)
  {
  _neighbouringTriangleInCuboidList[iC].resize(_trianglesInCuboidList[iC].size());


  //searches for the neighbouring triangles in the cuboid
    for (int ii=0;ii<_neighbouringTriangleInCuboidList[iC].size(); ii++)
    {
      std::vector<STLtriangle<T>> vec;

      for(int jj=ii+1;jj<_trianglesInCuboidList[iC].size();jj++)
      {

      auto& i = (_trianglesInCuboidList[iC])[ii];
      auto& j = (_trianglesInCuboidList[iC])[jj];
      STLtriangle<T> pi = i;
          //searches if there are common edges or vortex in the triangle pair
      if(norm(i.point[0]-j.point[0]) < sensitivity || norm(i.point[1]-j.point[1]) < sensitivity || norm(i.point[2]-j.point[2]) < sensitivity ||
      norm(i.point[0]-j.point[1]) < sensitivity || norm(i.point[1]-j.point[2]) < sensitivity || norm(i.point[2]-j.point[0]) < sensitivity ||
      norm(i.point[0]-j.point[2]) < sensitivity || norm(i.point[1]-j.point[0]) < sensitivity || norm(i.point[2]-j.point[1]) < sensitivity )
         {
              STLtriangle<T> pj= j;
              (_neighbouringTriangleInCuboidList[iC][ii]).push_back(pj);
              (_neighbouringTriangleInCuboidList[iC][jj]).push_back(pi);
         }
      }
    }
  }
  clout << "createneighbouringTriangleInCuboidVector finished" << std::endl;
}

/*
 * BlockLatticeSTLreader functions
 */
template<typename T>
BlockLatticeSTLreader<T>::BlockLatticeSTLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize,
                        int method, bool verbose, T overlap, T max)
  : _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName("meshPoints.stl"),
    _mesh(meshPoints, stlSize),
    _verbose(verbose),
    clout(std::cout, "BlockLatticeSTLreader")
{
  this->getName() = "BlockLatticeSTLreader";

  if (_verbose) {
    clout << "Voxelizing ..." << std::endl;
  }

  Vector<T,3> extension = _mesh.getMax() - _mesh.getMin();
  if ( util::nearZero(max) ) {
    max = util::max(extension[0], util::max(extension[1], extension[2])) + _voxelSize;
  }
  int j = 0;
  for (; _voxelSize * util::pow(2, j) < max; j++)
    ;
  Vector<T,3> center;
  T radius = _voxelSize * util::pow(2, j - 1);

  /// Find center of tree and move by _voxelSize/4.
  for (unsigned i = 0; i < 3; i++) {
    center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize / 4.;
  }

  /// Create tree

  _tree = new Octree<T>(center, radius, &_mesh, j, _overlap);

  /// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel inside the STL.
  for (int i = 0; i < 3; i++) {
    this->_myMin[i] = center[i] + _voxelSize / 2.;
    this->_myMax[i] = center[i] - _voxelSize / 2.;
  }
  for (int i = 0; i < 3; i++) {
    while (this->_myMin[i] > _mesh.getMin()[i]) {
      this->_myMin[i] -= _voxelSize;
    }
    while (this->_myMax[i] < _mesh.getMax()[i]) {
      this->_myMax[i] += _voxelSize;
    }
    this->_myMax[i] -= _voxelSize;
    this->_myMin[i] += _voxelSize;
  }
  //automaticly choose the method with minimum extension in its direction

  /*if(extension[0] == std::min_element(extension.begin(), extension.end())){
    method = 4;
  }
  else if(extension[1] == std::min_element(extension.begin(), extension.end())){
    method = 5;
  }
  else if(extension[2] == std::min_element(extension.begin(), extension.end())){
    method = 0;
  }
  */


  // Indicate nodes of the tree. (Inside/Outside)
  switch (method) {
  case 1:
    indicate1();
    break;
  case 3:
    indicate3();
    break;
  case 4:
    indicate2_Xray();
    break;
  case 5:
    indicate2_Yray();
    break;
  default:
    indicate2();
    break;
  }

  if (_verbose) {
    print();
  }
  if (_verbose) {
    clout << "Voxelizing ... OK" << std::endl;
  }
}

template<typename T>
BlockLatticeSTLreader<T>::~BlockLatticeSTLreader()
{
  delete _tree;
}

/*
 *  Old indicate function (slower, more stable)
 *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
 *  with STL for each ray. Odd number of intersection means inside (Majority vote).
 */

template<typename T>
void BlockLatticeSTLreader<T>::indicate1()
{
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();
  Vector<T,3> dir, pt, s;

  int intersections = 0;
  int inside = 0;
  Octree<T>* node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  for (; it != leafs.end(); ++it) {
    inside = 0;

    pt = (*it)->getCenter();
    intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);
    (*it)->setInside(inside > 1);
  }
}

/*
 *  New indicate function (faster, less stable)
 *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes on the fly.
 */
template<typename T>
void BlockLatticeSTLreader<T>::indicate2()
{
  T rad = _tree->getRadius();
  Vector<T,3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T,3> pt = rayPt;
  Vector<T,3> rayDir;
  rayDir[0] = 0.;
  rayDir[1] = 0.;
  rayDir[2] = 1.;
  //Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T>* node = nullptr;
  unsigned short rayInside = 0;
  Vector<T,3> nodeInters;
  while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[2] = node->getCenter()[2] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[2] = node->getCenter()[2] - node->getRadius();
      rayInside = 0;
      while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[2] = rayPt[2];
      pt[1] += _voxelSize;
    }
    pt[1] = rayPt[1];
    pt[0] += _voxelSize;
  }
}



/*
 *  New indicate function (faster, less stable)
 *  Define ray in X-direction for each Voxel in YZ-layer. Indicate all nodes on the fly.
 */

template<typename T>
void BlockLatticeSTLreader<T>::indicate2_Xray()
{
  T rad = _tree->getRadius();
  Vector<T,3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T,3> pt = rayPt;
  Vector<T,3> rayDir;
  rayDir[0] = 1.;
  rayDir[1] = 0.;
  rayDir[2] = 0.;
  //Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T>* node = nullptr;
  unsigned short rayInside = 0;
  Vector<T,3> nodeInters;
  while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[0] = node->getCenter()[0] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[0] = node->getCenter()[0] - node->getRadius();
      rayInside = 0;
      while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[0] = rayPt[0];
      pt[1] += _voxelSize;
    }
    pt[1] = rayPt[1];
    pt[2] += _voxelSize;
  }
}

/*
 *  New indicate function (faster, less stable)
 *  Define ray in Y-direction for each Voxel in XZ-layer. Indicate all nodes on the fly.
 */

template<typename T>
void BlockLatticeSTLreader<T>::indicate2_Yray()
{
  T rad = _tree->getRadius();
  Vector<T,3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T,3> pt = rayPt;
  Vector<T,3> rayDir;
  rayDir[0] = 0.;
  rayDir[1] = 1.;
  rayDir[2] = 0.;
  //Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T>* node = nullptr;
  unsigned short rayInside = 0;
  Vector<T,3> nodeInters;
  while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[1] = node->getCenter()[1] - node->getRadius();
    rayInside = 0;
    while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[1] = node->getCenter()[1] - node->getRadius();
      rayInside = 0;
      while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[1] = rayPt[1];
      pt[0] += _voxelSize;
    }
    pt[0] = rayPt[0];
    pt[2] += _voxelSize;
  }
}

/*
 *  Double ray approach: two times (X-, Y-, Z-direction) for each leaf.
 *  Could be use to deal with double layer triangles and face intersections.
 */
template<typename T>
void BlockLatticeSTLreader<T>::indicate3()
{
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();

  Vector<T,3> dir, pt, s;
  Octree<T>* node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  int intersections;
  int sum_intersections;

  for (; it != leafs.end(); ++it) {
    pt = (*it)->getCenter();
    intersections = 0;
    sum_intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// X- dir
    intersections = 0;
    s = pt;  // + step;
    dir[0] = -1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] > _mesh.getMin()[0] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = -1;
    dir[2] = 0;
    while (s[1] > _mesh.getMin()[1] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = -1;
    while (s[2] > _mesh.getMin()[2] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }
    (*it)->setInside(sum_intersections > 5);
  }
}

template<typename T>
bool BlockLatticeSTLreader<T>::operator() (bool output[], const T input[])
{
  output[0] = false;
  T coords = _tree->getRadius();
  Vector<T,3> c(_tree->getCenter());
  if (c[0] - coords < input[0] && input[0] < c[0] + coords && c[1] - coords < input[1]
      && input[1] < c[1] + coords && c[2] - coords < input[2] && input[2] < c[2] + coords) {
    std::vector<T> tmp(input, input + 3);
    output[0] = _tree->find(tmp)->getInside();
  }
  return true;
}

template<typename T>
bool BlockLatticeSTLreader<T>::distance(T& distance, const Vector<T,3>& origin,
                            const Vector<T,3>& direction, int iC)
{
  Octree<T>* node = nullptr;
  Vector<T,3> dir(direction);
  dir = normalize(dir);
  Vector<T,3> extends = _mesh.getMax() - _mesh.getMin();
  Vector<T,3> pt(origin);
  Vector<T,3> q;
  Vector<T,3> s;
  Vector<T,3> center = _mesh.getMin() + 1 / 2. * extends;
  T step = _voxelSize / 1000., a = 0;

  for (int i = 0; i < 3; i++) {
    extends[i] /= 2.;
  }

  if (!(_mesh.getMin()[0] < origin[0] && origin[0] < _mesh.getMax()[0]
        && _mesh.getMin()[1] < origin[1] && origin[1] < _mesh.getMax()[1]
        && _mesh.getMin()[2] < origin[2] && origin[2] < _mesh.getMax()[2])) {
    T t = T(), d = T();
    bool foundQ = false;

    if (dir[0] > 0) {
      d = _mesh.getMin()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];

      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }
    else if (dir[0] < 0) {
      d = _mesh.getMax()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[1] > 0 && !foundQ) {
      d = _mesh.getMin()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }
    else if (dir[1] < 0 && !foundQ) {
      d = _mesh.getMax()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[2] > 0 && !foundQ) {
      d = _mesh.getMin()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    }
    else if (dir[2] < 0 && !foundQ) {
      d = _mesh.getMax()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    }

    if (!foundQ) {
      return false;
    }
  }

  while ((util::fabs(pt[0] - center[0]) < extends[0])
         && (util::fabs(pt[1] - center[1]) < extends[1])
         && (util::fabs(pt[2] - center[2]) < extends[2])) {
    node = _tree->find(pt);
    if (node->closestIntersection(Vector<T,3>(origin), dir, q, a)) {
      Vector<T,3> vek(q - Vector<T,3>(origin));
      distance = norm(vek);
      return true;
    }
    else {
      Octree<T>* tmpNode = _tree->find(pt);
      tmpNode->intersectRayNode(pt, dir, s);
      for (int i = 0; i < 3; i++) {
        pt[i] = s[i] + step * dir[i];
      }
    }
  }

  if (_verbose) {
    clout << "Returning false" << std::endl;
  }
  return false;
}

template<typename T>
Vector<T,3> BlockLatticeSTLreader<T>::findNormalOnSurface(const PhysR<T,3>& pt)
{
  // Check if the position is on the corner of a triangle
  unsigned countTriangles = 0;
  Vector<T,3> normal(T(0));

  for (const STLtriangle<T>& triangle : _mesh.getTriangles()) {
    if (triangle.isPointInside(pt)) {
      ++countTriangles;
      normal+=triangle.getNormal();
    }
  }
  if (countTriangles > 0) {
    return normal / countTriangles;
  }

  // if nothing was found return (0,0,0) to indicate that nothing was found
  return normal;
}

template<typename T>
Vector<T,3> BlockLatticeSTLreader<T>::evalSurfaceNormal(const Vector<T,3>& origin)
{
  Vector<T,3> normal(0.);
  PhysR<T,3> closestPoint;
  T distance = std::numeric_limits<T>::max();
  STLtriangle<T> t;
  for (const STLtriangle<T>& triangle : _mesh.getTriangles()) {
    PhysR<T,3> const pointOnTriangle = triangle.closestPtPointTriangle(origin);
    //TODO! There could be a problem with concave surfaces on STL when there is a bigger distance there could be no point on the triangle and function does not work!
    if (triangle.isPointInside(pointOnTriangle))
    {
    PhysR<T,3> const currDistance = pointOnTriangle - origin;
    T currDistanceNorm = norm(currDistance);
    if (util::nearZero(currDistanceNorm)) {
      std::cout << "util near zero" << std::endl;
      return findNormalOnSurface(origin);
    }
    else if (distance > currDistanceNorm) {
      normal = currDistance;
      distance = currDistanceNorm;
      closestPoint = pointOnTriangle;
      t = triangle;
    }
    }
  }

  if (!util::nearZero(norm(normal))) {
    normal = normalize(normal);
  }
  //FP not sure about, how it should work, without that it simply points to the STL!
 /* else {
    return normal;
  }*/
/*
  // Possible change of the sign so that the normal fits to the SDF logic
  if(distance < _voxelSize) {
    bool isInsideInnerPoint;
    this->operator()(&isInsideInnerPoint, (closestPoint-_voxelSize*normal).data());
    bool isInsideOuterPoint;
    this->operator()(&isInsideOuterPoint, (closestPoint+_voxelSize*normal).data());
    normal *= (isInsideInnerPoint && !isInsideOuterPoint ? 1 : -1);
  }
  else {
    bool isInside;
    this->operator()(&isInside, origin.data());
    normal *= (isInside ? 1 : -1);
  }*/

  if(util::dotProduct3D(normal, t.getNormal())> 0)
    {

      distance =-1.*util::fabs(distance);//inside
    }
    else
    {
      distance =util::fabs(distance);//inside
    }
  return normal;
}
/*
template<typename T>
T BlockLatticeSTLreader<T>::signedDistance(int iC, const Vector<T,3>& input )
{
  T distance = std::numeric_limits<T>::max();
  auto& triangles = _trianglesInCuboidList[iC];
  STLtriangle<T> _triangle;
  Vector<T,3> normal;
  for (const STLtriangle<T>& triangle : triangles)
  {
    PhysR<T,3> pointOnTriangle = triangle.closestPtPointTriangle(input);
    T distance2 = util::min(distance, norm(pointOnTriangle - input));
    //TODO! There could be a problem with concave surfaces on STL when there is a bigger distance there could be no point on the triangle and function does not work!
    if(distance2 < distance && triangle.isPointInside(pointOnTriangle))
    {
      distance =distance2;
      _triangle = triangle;
      normal = normalize(pointOnTriangle - input);
    }
  }
  if(util::dotProduct3D(normal, _triangle.getNormal())> 0)
  {

      distance =-1.*util::fabs(distance);
  }
  else
  {
      distance =util::fabs(distance);
  }
  return distance;
}

template<typename T>
T BlockLatticeSTLreader<T>::signedDistance(const Vector<T,3>& input)
{
  int iC = -1;
  int globC = _cuboidDecomposition.get_iC(input, 0);
  if (_loadBalancer.isLocal(globC)) {
    iC = _loadBalancer.loc(globC);
  } else {
    globC = _cuboidDecomposition.get_iC(input, 3);
    if (_loadBalancer.isLocal(globC)) {
      iC = _loadBalancer.loc(globC);
    } else {
      throw std::runtime_error("Distance origin must be located within local cuboid (overlap)");
    }
  }
  return signedDistance(iC, input);
}
*/


///new signed distance function, sign based on the dot product computation
template<typename T>
T BlockLatticeSTLreader<T>::signedDistance(int iC, const Vector<T,3>& input )
{

  using namespace util;
  T distance = std::numeric_limits<T>::max();
  auto& triangles = _trianglesInCuboidList[iC];
  STLtriangle<T> _triangle;
  int triangle_index=-1;
  Vector<T,3> normal;

  for (int i=0;i<triangles.size();i++)
  {
    auto& triangle =  triangles[i];
    PhysR<T,3> pointOnTriangle = triangle.closestPtPointTriangle(input);
    T distance2 = util::min(distance, norm(pointOnTriangle - input));
    if(distance2 < distance  && triangle.isPointInside(pointOnTriangle) )
    {
      distance =distance2;
      _triangle = triangle;
      triangle_index = i;
      normal = normalize(pointOnTriangle - input);
    }
  }
  Vector<T,3> distancesToEdges = {}, VortexPoint = {}, EdgePoint1 = {}, EdgePoint2 = {};
  _triangle.getPointToEdgeDistances ((_triangle.closestPtPointTriangle(input)).data(), distancesToEdges);
  Vector<T,3> pseudoNormal = normal;//used for edge and vortex point type
  //different threatment for edge points
  if(_triangle.isEdgePoint(distancesToEdges, EdgePoint1, EdgePoint2))
  {
    for(auto& i: _neighbouringTriangleInCuboidList[iC][triangle_index])
    {
      if( (EdgePoint1 == i.point[0] || EdgePoint1 == i.point[1] || EdgePoint1 == i.point[2]) && (EdgePoint2 == i.point[0] || EdgePoint2 == i.point[1] || EdgePoint2 == i.point[2]))
      {
        pseudoNormal = normalize(M_PI*_triangle.getNormal() + M_PI*i.getNormal());
        break;
      }
     /* if ( i == (_neighbouringTriangleInCuboidList[iC][triangle_index][_neighbouringTriangleInCuboidList[iC][triangle_index].size()-1]))
        std::cout << "edge not found" << std::endl <<std::endl << std::endl;*/
    }
    //set the sign
  if(util::dotProduct3D(normal, pseudoNormal)> 0)
  {

      distance =-1.*util::fabs(distance);
  }
  else
  {
      distance =util::fabs(distance);
  }
  }
  //different threatment for vortex point
  else if(_triangle.isVortexPoint(distancesToEdges, VortexPoint))
  {
    int triangle_counter = 0;
    Vector<T,3> vortex1, vortex2, edge1, edge2;
    T dotP;
    if(VortexPoint == _triangle.point[0])
    {
      vortex1 = _triangle.point[1];
      vortex2 = _triangle.point[2];
    }
    else
    {
    if(VortexPoint == _triangle.point[1])
    {
      vortex1 = _triangle.point[0];
      vortex2 = _triangle.point[2];
    }
    else
    {
    if(VortexPoint == _triangle.point[2])
    {
      vortex1 = _triangle.point[0];
      vortex2 = _triangle.point[1];
    }
    }
    }
    edge1 = normalize(vortex1 - VortexPoint);
    edge2 = normalize(vortex2 - VortexPoint);
    dotP = dotProduct3D(edge1, edge2);

       //some problems with 0 and pi angle in acos
    if(dotP < -0.999999)
          pseudoNormal = M_PI*_triangle.getNormal();
        else if(dotP > 0.999999)
        {pseudoNormal = 0;}
        else
        {
    pseudoNormal = acos(dotP)*_triangle.getNormal();
    triangle_counter++;
        }

    for(auto& i: _neighbouringTriangleInCuboidList[iC][triangle_index])
    {
      if( VortexPoint ==  i.point[0] || VortexPoint ==  i.point[1] || VortexPoint ==  i.point[2])
      {
        Vector<T,3> vortex1, vortex2, edge1, edge2;
        T dotP;
        if(VortexPoint == i.point[0])
        {
          vortex1 = i.point[1];
          vortex2 = i.point[2];
        }
        else
        {
        if(VortexPoint == i.point[1])
        {
          vortex1 = i.point[0];
          vortex2 = i.point[2];
        }
        else
        {
        if(VortexPoint == i.point[2])
        {
          vortex1 = i.point[0];
          vortex2 = i.point[1];
        }
        else
          std::cout << "something went wrong" << std::endl;
        }
        }
        edge1 = normalize(vortex1 - VortexPoint);
        edge2 = normalize(vortex2 - VortexPoint);
        dotP = dotProduct3D(edge1, edge2);
        //std::cout << "inside if" << std::endl;
        //some problems with 0 and pi angle in acos
        if(dotP < -0.999999)
        { pseudoNormal += M_PI*i.getNormal();}
        else if(dotP > 0.999999)
        {;}
        else
        {
        pseudoNormal += acos(dotP)*i.getNormal();
         triangle_counter++;
        }
      }

    }
    pseudoNormal = normalize(pseudoNormal);
  /*   std::cout << "normal " << normal  << "pseudonormal " << pseudoNormal << std::endl << std::endl;
    std::cout << "dotProduct3D is " << util::dotProduct3D(normal, pseudoNormal) << std::endl;
    std::cout << "number of triangle " << triangle_counter << std::endl; */

  if(util::dotProduct3D(normal, pseudoNormal)> 0)
  {

      distance =-1.*util::fabs(distance);
  }
  else
  {
      distance =util::fabs(distance);
  }

  }
  else
  {
  if(util::dotProduct3D(normal, _triangle.getNormal())> 0)
  {

      distance =-1.*util::fabs(distance);
  }
  else
  {
      distance =util::fabs(distance);
  }
  }
  return distance;
}

template<typename T>
T BlockLatticeSTLreader<T>::signedDistance(const Vector<T,3>& input)
{
  int iC = -1;
  auto globC = _cuboidDecomposition.getC(input, 0);
  if (_loadBalancer.isLocal(*globC)) {
    iC = _loadBalancer.loc(*globC);
  } else {
    globC = _cuboidDecomposition.getC(input, 3);
    if (_loadBalancer.isLocal(*globC)) {
      iC = _loadBalancer.loc(*globC);
    } else {
      throw std::runtime_error("Distance origin must be located within local cuboid (overlap)");
    }
  }
  return signedDistance(iC, input);
}


template <typename T>
Vector<T,3> BlockLatticeSTLreader<T>::surfaceNormal(const Vector<T,3>& pos, const T meshSize)
{
  return evalSurfaceNormal(pos);
}

template<typename T>
void BlockLatticeSTLreader<T>::print()
{
  _mesh.print();
  _tree->print();
  clout << "voxelSize=" << _voxelSize << "; stlSize=" << _stlSize << std::endl;
  clout << "minPhysR(VoxelMesh)=(" << this->_myMin[0] << "," << this->_myMin[1]
        << "," << this->_myMin[2] << ")";
  clout << "; maxPhysR(VoxelMesh)=(" << this->_myMax[0] << ","
        << this->_myMax[1] << "," << this->_myMax[2] << ")" << std::endl;
}

template<typename T>
void BlockLatticeSTLreader<T>::writeOctree()
{
  _tree->write(_fName);
}

template<typename T>
void BlockLatticeSTLreader<T>::writeSTL(std::string stlName)
{
  if (stlName == "") {
    _mesh.write(_fName);
  }
  else {
    _mesh.write(stlName);
  }
}

template<typename T>
void BlockLatticeSTLreader<T>::setNormalsOutside()
{
  unsigned int noTris = _mesh.triangleSize();
  Vector<T,3> center;
  //Octree<T>* node = nullptr;
  for (unsigned int i = 0; i < noTris; i++) {
    center = _mesh.getTri(i).getCenter();
    if (_tree->find(
          center + _mesh.getTri(i).normal * util::sqrt(3.) * _voxelSize)->getInside()) {
      //      cout << "Wrong direction" << std::endl;
      Vector<T,3> pt(_mesh.getTri(i).point[0]);
      _mesh.getTri(i).point[0] = _mesh.getTri(i).point[2];
      _mesh.getTri(i).point[2] = pt;
      _mesh.getTri(i).init();
      //      _mesh.getTri(i).getNormal()[0] *= -1.;
      //      _mesh.getTri(i).getNormal()[1] *= -1.;
      //      _mesh.getTri(i).getNormal()[2] *= -1.;
    }
  }
}

template<typename T>
void BlockLatticeSTLreader<T>::setBoundaryInsideNodes()
{
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  for (auto it = leafs.begin(); it != leafs.end(); ++it) {
    if ((*it)->getBoundaryNode()) {
      (*it)->setInside(true);
    }
  }
}

}  // namespace olb

#endif
