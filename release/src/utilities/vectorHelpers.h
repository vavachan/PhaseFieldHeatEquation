/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Lukas Baron, Mathias J. Krause
 *                2023 Julius Jessberger
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

#ifndef VECTOR_HELPERS_H
#define VECTOR_HELPERS_H

#include <assert.h>
#include <vector>
#include <string>
#include <sstream>
#include <limits>

#include "io/ostreamManager.h"
#include "utilities/omath.h"
#include "core/vector.h"
#include "utilities/aDiff.h"

namespace olb {

template<typename T, unsigned Size> class Vector;

namespace util {

template <class T, unsigned DIM> class ADf;

template <class T, unsigned DIM> inline ADf<T,DIM> sqrt (const ADf<T,DIM>& a);

template<typename S>
using StdVector = std::vector<S,std::allocator<S>>;

template <class T>
inline void copyN(T c[], const T a[], const unsigned dim) any_platform
{
  for (unsigned i=0; i<dim; i++) {
    c[i] = a[i];
  }
}

template <class S, class T>
inline void copyN(S c[], const T a[], const unsigned dim) any_platform
{
  for (unsigned i=0; i<dim; i++) {
    c[i] = a[i];
  }
}

template <class T>
inline void copy3(T c[], const T a[])
{
  for (unsigned i=0; i<3; i++) {
    c[i] = a[i];
  }
}


template <typename T>
std::vector<T> fromVector3(const Vector<T,3>& vec)
{
  std::vector<T> v;
  v.push_back(vec[0]);
  v.push_back(vec[1]);
  v.push_back(vec[2]);
  return v;
}
template <typename T>
std::vector<T> fromVector2(const Vector<T,2>& vec)
{
  std::vector<T> v;
  v.push_back(vec[0]);
  v.push_back(vec[1]);
  return v;
}

/// l2 norm to the power of 2 of a vector of arbitrary length
template <typename T>
T norm2(const T* a, unsigned size)
{
  T v{};
  for (unsigned iD=0; iD<size; ++iD) {
    v += a[iD]*a[iD];
  }
  return v;
}

/// l2 norm to the power of 2 of a vector of arbitrary length
template <typename T>
T norm2(const std::vector<T>& a)
{
  return norm2(a.data(), a.size());
}

/// l2 norm of a vector of arbitrary length
template <typename T>
T norm(const T* a, unsigned size)
{
  return util::sqrt(norm2(a, size));
}

/// l2 norm of a vector of arbitrary length
template <typename T>
T norm(const std::vector<T>& a)
{
  return norm(a.data(), a.size());
}

/// dot product, only valid in 3d
template <typename T>
T dotProduct3D(const Vector<T,3>& a, const Vector<T,3>& b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/// dot product, only valid in 2d
template <typename T>
T dotProduct2D(const Vector<T,2>& a, const Vector<T,2>& b)
{
  return a[0]*b[0] + a[1]*b[1];
}

/// dot product
template <typename T,unsigned D>
T dotProduct(const Vector<T,D>& a, const Vector<T,D>& b)
{
  if constexpr (D==2){
    return dotProduct2D(a, b);
  } else {
    return dotProduct3D(a, b);
  }
}

template <typename T, unsigned D>
Vector<T,D> normalize(const Vector<T,D>& a)
{
  return a / norm(a);
}

/// returns a normalized vector, works for arbitrary lengths
template <typename T>
std::vector<T> normalize(const std::vector<T>& a)
{
  std::vector<T> out(a);
  T scale = norm(a);
  assert(scale>0);
  for (unsigned int iDim=0; iDim<a.size(); iDim++) {
    out[iDim] /= scale;
  }
  return out;
}

/// applies floor to each component of a vector
template <typename T, unsigned Size>
Vector<T,Size> floor(const Vector<T,Size>& a) any_platform
{
  Vector<T,Size> out;
  for (unsigned int iDim=0; iDim < Size; ++iDim) {
    out[iDim] = util::floor(a[iDim]);
  }
  return out;
}

/// applies ceil to each component of a vector
template <typename T, unsigned Size>
Vector<T,Size> ceil(const Vector<T,Size>& a) any_platform
{
  Vector<T,Size> out;
  for (unsigned int iDim=0; iDim < Size; ++iDim) {
    out[iDim] = util::ceil(a[iDim]);
  }
  return out;
}

/// applies fmod to each component of a vector
template <typename T, typename S, unsigned Size>
Vector<T,Size> fmod(const Vector<T,Size>& a, S b)
{
  Vector<T,Size> out;
  for (unsigned int iDim=0; iDim < Size; ++iDim) {
    out[iDim] = util::fmod(a[iDim], b);
  }
  return out;
}

/// computes the average of all elements
template <typename T, unsigned Size>
T average(const Vector<T,Size>& a)
{
  T sum = a[0];
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    sum += a[iDim];
  }
  return sum/Size;
}

/// finds maximum element of all elements
template <typename T, unsigned Size>
T max_element(const Vector<T,Size>& a)
{
  T max = a[0];
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    max = std::max(max,a[iDim]);
  }
  return max;
}

/// finds minimum element of all elements
template <typename T, unsigned Size>
T min_element(const Vector<T,Size>& a)
{
  T min = a[0];
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    min = std::min(min,a[iDim]);
  }
  return min;
}

/// finds position of maximum element of all elements
template <typename T, unsigned Size>
unsigned maxElementPos(const Vector<T,Size>& a)
{
  unsigned maxPos = 0;
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    if (a[iDim]>a[maxPos]){
      maxPos = iDim;
    }
  }
  return maxPos;
}

/// finds position of minimum element of all elements
template <typename T, unsigned Size>
unsigned minElementPos(const Vector<T,Size>& a)
{
  unsigned minPos = 0;
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    if (a[iDim]<a[minPos]){
      minPos = iDim;
    }
  }
  return minPos;
}

/// finds maximum element of all absolute elements
template <typename T, unsigned Size>
T maxElementAbs(const Vector<T,Size>& a)
{
  T maxAbs = a[0];
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    if (abs(a[iDim])>abs(maxAbs)){
      maxAbs = a[iDim];
    }
  }
  return maxAbs;
}

/// finds position of maximum element of all absolute elements
template <typename T, unsigned Size>
unsigned maxElementAbsPos(const Vector<T,Size>& a)
{
  unsigned maxAbsPos = 0;
  for (unsigned int iDim=1; iDim < Size; ++iDim) {
    if (abs(a[iDim])>abs(a[maxAbsPos])){
      maxAbsPos = iDim;
    }
  }
  return maxAbsPos;
}

/// Calculates angles between two 2D vectors
template <typename T, bool ensureAngularBounds=true>
T angleBetweenVectors(const Vector<T,2>& a, const Vector<T,2>& b)
{
  if constexpr(ensureAngularBounds){
    return std::fmod(util::atan2(b[1]*a[0]-b[0]*a[1], a[0]*b[0]+a[1]*b[1]), M_PI);
  } else {
    return util::atan2(b[1]*a[0]-b[0]*a[1], a[0]*b[0]+a[1]*b[1]);
  }
}

/// Calculates angles between two 3D vectors
template <typename T, bool ensureAngularBounds=true>
Vector<T,3> angleBetweenVectors(const Vector<T,3>& a, const Vector<T,3>& b)
{
  Vector<T,3> angles;
  angles[0] = angleBetweenVectors<T,ensureAngularBounds>(Vector<T,2>(a[1], a[2]), Vector<T,2>(b[1], b[2]));
  angles[1] = angleBetweenVectors<T,ensureAngularBounds>(Vector<T,2>(a[0], a[2]), Vector<T,2>(b[0], b[2]));
  angles[2] = angleBetweenVectors<T,ensureAngularBounds>(Vector<T,2>(a[0], a[1]), Vector<T,2>(b[0], b[1]));
  return angles;
}

template <typename T>
constexpr T determinant(const Vector<T,2>& v0, const Vector<T,2>& v1)
{
  return crossProduct2D(v0, v1);
}

template <typename T>
constexpr T determinant(const Vector<T,3>& v0, const Vector<T,3>& v1,
  const Vector<T,3>& v2)
{
  return v0[0]*v1[1]*v2[2] + v1[0]*v2[1]*v0[2] + v2[0]*v0[1]*v1[2]
    - v0[2]*v1[1]*v2[0] - v1[2]*v2[1]*v0[0] - v2[2]*v0[1]*v1[0];
}


/// Solve (a0 a1) x = rhs
// Returns +-inf if a0, a1 are linearly dependent and rhs is not (no solution)
// Returns nan if a0, a1, rhs are linearly dependent (infinitely many solutions)
template <typename T>
Vector<T,2> solveLinearSystem(const Vector<T,2>& a0,
  const Vector<T,2>& a1, const Vector<T,2>& rhs)
{
  const T det = crossProduct2D(a0, a1);
  return Vector<T,2> (crossProduct2D(rhs, a1) / det, crossProduct2D(a0, rhs) / det);
}

template <typename T>
Vector<T,3> solveLinearSystem_help(const Vector<T,3>& a0,
  const Vector<T,3>& a1, const Vector<T,3>& a2, const Vector<T,3>& rhs,
  T det)
{
  Vector<T,3> res;
  res[0] = determinant(rhs, a1, a2) / det;
  res[1] = determinant(a0, rhs, a2) / det;
  res[2] = determinant(a0, a1, rhs) / det;
  return res;
}

/// Solve (a0 a1 a2) x = rhs
// Returns +-inf if a0, a1, a2 are linearly dependent and rhs is not (no solution)
// Returns nan if a0, a1, a2, rhs are linearly dependent (infinitely many solutions)
// Cf. https://danceswithcode.net/engineeringnotes/linear_equations/linear_equations.html
// for formula
template <typename T>
Vector<T,3> solveLinearSystem(const Vector<T,3>& a0,
  const Vector<T,3>& a1, const Vector<T,3>& a2, const Vector<T,3>& rhs)
{
  const T det = determinant(a0, a1, a2);
  return solveLinearSystem_help(a0, a1, a2, rhs, det);
}


/*
/// algorithm by Möller–Trumbore (TODO add ref), implemented by Lucas Cruz and Mathias J. Krause
/// returns true if there is an intersection of a triangle given by (point0, point1, point1) and a ray given by its origin and direction and computes the distance
template <typename T>
bool triangleIntersectionWithNormalDirection(const std::vector<T>& point0,
    const std::vector<T>& point1, const std::vector<T>& point2,
    const std::vector<T>& origin, const std::vector<T>& normalDirection,
    T& distance)
{
  T EPSILON = std::numeric_limits<T>::epsilon();
  std::vector<T> e1, e2;
  std::vector<T> P, Q, TT;
  T det, inv_det;
  T t, u, v;
  e1 = point1 - point0;
  e2 = point2 - point0;
  P = crossProduct3D(normalDirection, e2);
  det = dotProduct3D(P, e1);
  if (det > -EPSILON && det < EPSILON) {
    return false;
  }
  inv_det = T(1) / det;
  TT = origin - point0;
  u = dotProduct3D(TT, P)*inv_det;
  if (u < T() || u > T(1)) {
    return false;
  }
  Q = crossProduct3D(TT, e1);
  v = dotProduct3D(normalDirection, Q) * inv_det;
  if (v < T() || u + v  > T(1)) {
    return false;
  }
  t = dotProduct3D(e2, Q)*inv_det;
  if (t > EPSILON) {
    distance = t;
    return true;
  }
  return false;
}

template <typename T>
bool triangleIntersection(const std::vector<T>& point0, const std::vector<T>& point1, const std::vector<T>& point2, const std::vector<T>& origin, const std::vector<T>& direction, T& distance)
{
  std::vector<T> normalDirection(normalize(direction) );
  return triangleIntersectionWithNormalDirection(point0, point1, point2, origin, normalDirection, distance );
}
*/
template <typename T>
std::vector<T> assign(T a, T b)
{
  std::vector<T> v1;
  v1.push_back(a);
  v1.push_back(b);
  return v1;
}

template <typename T>
std::vector<T> assign(T a, T b, T c)
{
  std::vector<T> v1;
  v1.push_back(a);
  v1.push_back(b);
  v1.push_back(c);
  return v1;
}


template<typename U>
void print(U data, const std::string& name="", OstreamManager clout = OstreamManager(std::cout,"print"),
           const char delimiter=',')
{
  static_assert(!std::is_integral<U>::value && !std::is_floating_point<U>::value, "passed integral or floating_point value to function print()");
  if (name != "") {
    clout << name << " = ";
  }
  for ( auto& element : data ) {
    clout << std::fixed << element << delimiter << ' ';
  }
  clout << std::endl;
}

/// Check, if object is contained in iteratable container c.
// U must provide equality 'operator==' and this has to fit the elements of c.
template<typename C, typename U>
bool isContained(const C& c, U object) {
  return (std::find(c.begin(), c.end(), object) != c.end());
}

/// Creates a container of type C.
// See below for template specializations.
template<typename C>
struct ContainerCreator { };

template<typename T>
struct ContainerCreator<std::vector<T>> {
  using C = std::vector<T>;

  static constexpr C create(std::size_t size) {
    return C(size);
  }
};

template<typename T, std::size_t SIZE>
struct ContainerCreator<std::array<T,SIZE>> {
  using C = std::array<T,SIZE>;

  static constexpr C create(std::size_t size) {
    OLB_PRECONDITION(SIZE == size);
    return C{};
  }
};

template<typename T, unsigned SIZE>
struct ContainerCreator<Vector<T,SIZE>> {
  using C = Vector<T,SIZE>;

  static constexpr C create(std::size_t size) {
    OLB_PRECONDITION(SIZE == size);
    return C{};
  }
};

/// @brief Compute serial index of symmetric tensor
/// @tparam DIM Spatial dimension
/// @param i First spatial index
/// @param j Second spatial index
/// @return Serialized index, as it is used in computeStress
template<unsigned DIM>
unsigned serialSymmetricTensorIndex(unsigned i, unsigned j) any_platform
{
  OLB_PRECONDITION(DIM > 1);
  OLB_PRECONDITION(DIM < 4);
  if constexpr (DIM==2) {
    return i+j;
  }
  else if constexpr (DIM==3) {
    const unsigned res = i+j;
    if ((i>0) && (j>0)) {
      return res+1;
    } else {
      return res;
    }
  }
  else {
    // other dimensions are not implemented
    exit(1);
  }
}

} // namespace util

} // namespace olb

#endif
