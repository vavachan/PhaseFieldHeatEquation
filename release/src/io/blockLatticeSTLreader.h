/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause
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
 * Input in STL format -- header file.
 */

#ifndef BLOCK_LATTICE_STL_READER_H
#define BLOCK_LATTICE_STL_READER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <limits>

#include "communication/loadBalancer.h"
#include "geometry/cuboidDecomposition.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"
#include "octree.h"
#include "core/vector.h"


// All OpenLB code is contained in this namespace.
namespace olb {




template<typename T>
class BlockLatticeSTLreader : public IndicatorF3D<T> {
private:
  /*
   *  Old indicate function (slower, more stable)
   *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
   *  with STL for each ray. Odd number of intersection means inside (Majority vote).
   */
  void indicate1();
  /*
   *  New indicate function (faster, less stable)
   *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes on the fly.
   */
  void indicate2();
  /*
  *  New indicate function (faster, less stable)
  *  Define ray in X-direction for each Voxel in YZ-layer. Indicate all nodes on the fly.
  */
  void indicate2_Xray();
  /*
  *  New indicate function (faster, less stable)
  *  Define ray in Y-direction for each Voxel in XZ-layer. Indicate all nodes on the fly.
  */
  void indicate2_Yray();
  /*
   *  Double ray approach: two times (X-, Y-, Z-direction) for each leaf.
   *  Could be use to deal with double layer triangles and face intersections.
   */

  void indicate3();

  /// Finds normal for points on the surface (do not use for points that aren't on the surface!)
  Vector<T,3> findNormalOnSurface(const PhysR<T,3>& pt);

  /// Finds surface normal
  Vector<T,3> evalSurfaceNormal(const Vector<T,3>& origin);

  /// Size of the smallest voxel
  T _voxelSize;
  /// Factor to get Si unit (m), i.e. "0.001" means mm
  T _stlSize;
  /// Overlap increases Octree radius by _overlap
  T _overlap;
  /// Pointer to tree
  Octree<T>* _tree;
  /// The filename
  const std::string _fName;
  /// The mesh
  STLmesh<T> _mesh;
  /// Variable for output
  bool _verbose;
  /// The OstreamManager
  mutable OstreamManager clout;

  CuboidDecomposition3D<T>& _cuboidDecomposition;

  LoadBalancer<T>& _loadBalancer;

  std::vector<std::vector<STLtriangle<T>>> _trianglesInCuboidList;

  std::vector<std::vector<std::vector<STLtriangle<T>>>> _neighbouringTriangleInCuboidList;

public:
  /**
   * Constructs a new BlockLatticeSTLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */

  BlockLatticeSTLreader(CuboidDecomposition3D<T>& cbg3d, LoadBalancer<T>& hlb, const std::string fName, T voxelSize, T stlSize=1, int method=2,
            bool verbose = false, T overlap=0., T max=0.);
  /**
   * Constructs a new BlockLatticeSTLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */
  BlockLatticeSTLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize=1, int method=2,
            bool verbose = false, T overlap=0., T max=0.);

  ~BlockLatticeSTLreader() override;
  /// Returns whether node is inside or not.
  bool operator() (bool output[], const T input[]) override;

  /// Computes distance to closest triangle intersection
  bool distance(T& distance,const Vector<T,3>& origin, const Vector<T,3>& direction, int iC=-1) override;

  /// Computes signed distance to closest triangle in direction of the surface normal in local cuboid locC
  T signedDistance(int locC, const Vector<T,3>& input);
  /// Computes signed distance to closest triangle in direction of the surface normal
  T signedDistance(const Vector<T,3>& input) override;

  /*  /// Computes signed distance to closest triangle in direction of the surface normal in local cuboid locC
  T signedDistance(int locC, const Vector<T,3>& input, STLtriangle<T>& _triangle);
  /// Computes signed distance to closest triangle in direction of the surface normal
  T signedDistance(const Vector<T,3>& input, STLtriangle<T>& _triangle);// override;*/

  /// Finds and returns normal of the closest surface (triangle)
  Vector<T,3> surfaceNormal(const Vector<T,3>& pos, const T meshSize=0) override;

  /// Prints console output
  void print();

  /// Writes STL mesh in Si units
  void writeSTL(std::string stlName="");

  /// Writes Octree
  void writeOctree();

  /// Rearranges normals of triangles to point outside of geometry
  void setNormalsOutside();

  /// Every octree leaf intersected by the STL will be part of the inside nodes.
  /// Artificially enlarges all details that would otherwise be cut off by the voxelSize.
  void setBoundaryInsideNodes();

  /// Returns tree
  inline Octree<T>* getTree() const
  {
    return _tree;
  };


  /// Returns mesh
  inline STLmesh<T>& getMesh()
  {
    return _mesh;
  };
  ///Needed for the signed Distance function: Creates a list of neighbouring triangles for each triangle
  void createNeighbouringTriangleInCuboidVector();
};
}  // namespace olb

#endif
