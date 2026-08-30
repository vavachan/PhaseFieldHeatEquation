/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2024 Thomas Henn, Mathias J. Krause, Christoph Gaul
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

#ifndef STL_READER_H
#define STL_READER_H

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
class Octree;

template<typename T>
struct STLtriangle {
  /** Test intersection between ray and triangle
   * \param pt Raypoint
   * \param dir Direction
   * \param q Point of intersection (if intersection occurs)
   * \param alpha Explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
   *            q = pt + alpha * dir
   * \param rad It's complicated. Imagine you have a sphere with radius rad moving a long the ray. Then q becomes the first point of the sphere to touch the triangle.
   */
  bool testRayIntersect(const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha, const T& rad = T(), bool print = false);
  Vector<T,3> closestPtPointTriangle(const Vector<T,3>& pt) const;

  /// A triangle contains 3 Points
  std::vector<Vector<T,3> > point;

  /// normal of triangle
  Vector<T,3> normal;

  /// variables explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
  /// precomputed for speedup
  Vector<T,3> uBeta, uGamma;
  T d, kBeta, kGamma;

public:
  /// Constructor constructs
  STLtriangle():point(3, Vector<T,3>()), normal(T()), uBeta(T()), uGamma(T()), d(T()), kBeta(T()), kGamma(T()) {};
  /// CopyConstructor copies
  STLtriangle(STLtriangle<T> const& tri):point(tri.point), normal(tri.normal), uBeta(tri.uBeta), uGamma(tri.uGamma), d(tri.d), kBeta(tri.kBeta), kGamma(tri.kGamma) {};
  /// Operator= equals
  STLtriangle<T>& operator=(STLtriangle<T> const& tri)
  {
    point = tri.point;
    normal = tri.normal;
    uBeta = tri.uBeta;
    uGamma = tri.uGamma;
    d = tri.d;
    kBeta = tri.kBeta;
    kGamma = tri.kGamma;
    return *this;
  };

  ~STLtriangle() {};

  /// Initializes triangle and precomputes member variables.
  void init();
  /// Return write access to normal
  inline Vector<T,3>& getNormal()
  {
    return normal;
  }
  /// Return read access to normal
  inline const Vector<T,3>& getNormal() const
  {
    return normal;
  }
  /// Returns center
  Vector<T,3> getCenter();
  /// Returns Pt0-Pt1
  std::vector<T> getE0();
  /// Returns Pt0-Pt2
  std::vector<T> getE1();
  /// Check whether a point is inside a triangle
  bool isPointInside(const PhysR<T,3>& pt) const;
  ///Returns true if the point is on a edge (smaller than sensitivity) and gives the perpendicular distances to the edges in output
  ///https://www.math.kit.edu/ianm2/lehre/am22016s/media/distance-harvard.pdf
  bool getPointToEdgeDistances (const Vector<T,3>& input, Vector<T,3> & output, T sensitivity = 1.e-15);
  ///Returns true if is near edge (smaller than sensitivity) and not near vortex and saves in P1 and P2 the vortex points of the edge
  bool isEdgePoint (const Vector<T,3> & input,  Vector<T,3>& P1, Vector<T,3> & P2, T sensitivity = 1.e-15);
  ///Returns true if is near vortex (smaller than sensitivity) and saves in P the vortex points
  bool isVortexPoint (const Vector<T,3> & input,Vector<T,3>& P, T sensitivity = 1.e-15);

  /// Returns the closest point to physR on the triangle
  Vector<T,3> closestPointTo(Vector<T,3> physR);

};

template<typename T>
class STLmesh {
  /// Computes distance squared betwenn p1 and p2
  T distPoints(Vector<T,3>& p1, Vector<T,3>& p2);
  /// Filename
  const std::string _fName;
  /// Vector of Triangles
  std::vector<STLtriangle<T> > _triangles;
  /// Min and Max points of axis aligned bounding box coordinate in SI units
  Vector<T,3> _min, _max;
  /// largest squared length of edge of all triangles
  T _maxDist2;
  /// OstreamManager
  mutable OstreamManager clout;

public:
  /**
   * Constructs a new STLmesh from a file
   * \param Filename - Filename
   * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   */
  STLmesh(std::string, T stlSize = 1.);

  /**
   * Constructs a new STLmesh from a file
   * \param Filename - Filename
   * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   */
  STLmesh(const std::vector<std::vector<T>> meshPoints, T stlSize = 1.);

  /// Returns reference to a triangle
  inline STLtriangle<T>& getTri(unsigned int i)
  {
    return _triangles[i];
  }
  /// Returns reference to all triangles
  inline std::vector<STLtriangle<T> >& getTriangles()
  {
    return _triangles;
  }
  /// Returns number of triangles
  inline unsigned int triangleSize() const
  {
    return _triangles.size();
  }
  /// Returns _min
  inline Vector<T,3>& getMin()
  {
    return _min;
  };
  /// Returns _max
  inline Vector<T,3>& getMax()
  {
    return _max;
  };
  /// Returns maxDist squared
  inline float maxDist2() const
  {
    return _maxDist2;
  }
  /// Prints console output
  void print(bool full = false);
  /// Writes STL mesh in Si units
  void write(std::string fName);
  /// Compute intersection between Ray and set of triangles; returns true if intersection is found
  bool testRayIntersect(const std::set<unsigned int>& tris, const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha);
};

/**
 * @brief Enum class that specifies the mode to use for computing the sign of the signed distance.
 *
 * This enum class specifies the mode to use for computing the sign of the signed distance.
 * It can take one of the following values:
 * - SignMode::EXACT: The sign is computed exactly using the winding number close to the surface (slow).
 * - SignMode::CACHED: The sign is computed from cache (fast, but less accurate in vicinity of the surface).
 */
enum class SignMode { EXACT, CACHED };


enum class RayMode : int {
  FastRayZ  = 1,  /// Indicate function with ray in Z-direction(faster, less stable). Default option
  Robust    = 2,  /// Old indicate function (slower, more stable)
  FastRayX  = 3,  /// Indicate function with ray in X-direction(faster, less stable).
  FastRayY  = 4,  /// Indicate function with ray in Y-direction(faster, less stable).
  DoubleRay = 5   /// Indicate function with double ray
};


template<typename T>
class STLreader : public IndicatorF3D<T> {
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

  /// Iterates over triangles close to passed point using the octree.
  /// If the point is outside of the mesh, it iterates over all triangles as fallback.
  template <typename F>
  void iterateOverCloseTriangles(const PhysR<T,3>& pt, F func, Octree<T>* leafNode = nullptr);

  /// Evaluates the normal for points on the surface (do not use for points that aren't on the surface!)
  /// Due to rounding errors it's possible that points that were found via STLtriangle.closestPtPointTriangle return false on STLtriangle.isPointInside.
  /// Therefore, we provide a fallback normal that should be set to the normal of the triangle the cloest point was located on.
  Vector<T,3> evalNormalOnSurface(const PhysR<T,3>& pt, const Vector<T,3>& fallbackNormal);

  /// Finds surface normal
  template <SignMode SIGNMODE>
  Vector<T,3> evalSurfaceNormal(const Vector<T,3>& origin);

  /**
   * @brief Computes the sign for the signed distance of a point to the surface of the geometry.
   *
   * @tparam SIGNMODE The mode to use for computing the sign.
   * @param pt The point for which the sign should be computed.
   * @param distance The distance of the point to the surface of the geometry. This parameter is marked as maybe_unused as it's only used in EXACT mode.
   * @return short The sign of the signed distance. Returns -1 if the point is inside the geometry, 1 otherwise.
   *
   * This function computes the sign for the signed distance of a point to the surface of the geometry.
   * The mode used for computing the sign is determined by the template parameter SIGNMODE.
   * If SIGNMODE is SignMode::EXACT, the sign is computed exactly using the distance parameter.
   * If SIGNMODE is SignMode::CACHE, the sign is computed from cache (less accurate).
   */
  template <SignMode SIGNMODE>
  short evalSignForSignedDistance(const Vector<T,3>& pt, [[maybe_unused]] const T distance, Vector<T,3>  vecdist = {} , STLtriangle<T> stlT = STLtriangle<T>() );

template<SignMode SIGNMODE>
Vector<T,3> evalSurfaceNormalForPseudoNormal(const Vector<T,3>& origin, Vector<T,3> & outputPointOnSurface);


  short evalSignForSignedDistanceFromNormal(const Vector<T,3>& normal, const Vector<T,3>& distance);
  /// Evaluate sign for signed distance
  /// Computes the sign following 10.1109/TVCG.2005.49 (distance corresponds to p - c)
  /// p: the position that is evaluated
  /// c: the closest point on the surface
  /// Note: This method is less robust.
  short evalSignForSignedDistanceFromPseudonormal(const Vector<T,3>& pseudonormal, const Vector<T,3>& distance);
  /// Evaluate sign for signed distance
  /// Computest the sign following 10.1145/2461912.2461916
  /// Note: This method is more robust.
  short evalSignForSignedDistanceFromWindingNumber(const Vector<T,3>& pt);
  /// Evaluate sign for signed distance using cached information
  /// Less accurate
  short evalSignForSignedDistanceFromCache(const Vector<T,3>& pt);

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

public:
  /**
   * Constructs a new STLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */

  STLreader(const std::string fName, T voxelSize, T stlSize=1, RayMode method = RayMode::FastRayZ,
            bool verbose = false, T overlap=0., T max=0.);
  /**
   * Constructs a new STLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */
  STLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize=1, RayMode method= RayMode::FastRayZ,
            bool verbose = false, T overlap=0., T max=0.);

  ~STLreader() override;
  /// Returns whether node is inside or not.
  bool operator() (bool output[], const T input[]) override;

  /// Returns whether node is inside the top-level octree or not.
  bool isInsideRootTree(const T input[]);

  /// Returns the closest point in the bounding box
  /// If `input` is already inside, then it returns `input`.
  Vector<T,3> closestPointInBoundingBox(const Vector<T,3>& input);

  /// Computes distance to closest triangle intersection
  bool distance(T& distance,const Vector<T,3>& origin, const Vector<T,3>& direction, int iC=-1) override;

  T distanceToClosestSurfacePoint(Vector<T,3> physR) {
    T closest = std::numeric_limits<T>::max();
    if (auto tree = _tree->find(physR)) {
      auto& triangles = tree->getTriangles();
      for (unsigned int iTriangle : triangles) {
        auto& triangle = _mesh.getTri(iTriangle);
        auto closeR = triangle.closestPointTo(physR);
        auto dist = norm(closeR - physR);
        if (dist < closest) {
          closest = dist;
        }
      }
    }
    return closest;
  }

  Vector<T,3> normalAtClosestSurfacePoint(Vector<T,3> physR) {
    T closest = std::numeric_limits<T>::max();
    int iT = -1;
    if (auto tree = _tree->find(physR)) {
      auto& triangles = tree->getTriangles();
      for (unsigned int iTriangle : triangles) {
        auto& triangle = _mesh.getTri(iTriangle);
        auto closeR = triangle.closestPointTo(physR);
        auto dist = norm(closeR - physR);
        if (dist < closest) {
          closest = dist;
          iT = iTriangle;
        }
      }
    }
    if (iT > 0) {
      auto& triangle = _mesh.getTri(iT);
      return triangle.getNormal();
    } else {
      return 0;
    }
  }

  /// Computes signed distance to closest triangle in direction of the surface normal
  /// Using the cached information (faster, but less accurate)
  T signedDistance(const Vector<T,3>& input) override;
  /// Computes exact signed distance to closest triangle in direction of the surface normal
  /// Much slower, but more accurate
  T signedDistanceExact(const Vector<T,3>& input) override;
  /// Computes signed distance to closest triangle in direction of the surface normal
  template <SignMode SIGNMODE>
  T signedDistance(const Vector<T,3>& input);

  /// Finds and returns normal of the closest surface (triangle)
  /// Using the cached information (faster, but less accurate)
  Vector<T,3> surfaceNormal(const Vector<T,3>& pos, const T meshSize=0) override;
  /// Finds and returns normal of the closest surface (triangle)
  /// Much slower, but more accurate
  Vector<T,3> surfaceNormalExact(const Vector<T,3>& pos, const T meshSize=0) override;
  /// Finds and returns normal of the closest surface (triangle)
  template <SignMode SIGNMODE>
  Vector<T,3> surfaceNormal(const Vector<T,3>& pos, const T meshSize=0);

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

};

}  // namespace olb

#endif
