/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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
 * Some arithmetic helper functions.
 */

#ifndef GEOMETRIC_OPERATIONS_H
#define GEOMETRIC_OPERATIONS_H

#include "core/vector.h"
#include "dimensionConverter.h"
#include "matrix.h"

namespace olb {

namespace util {

// angle conversions
template <typename T, unsigned D>
decltype(Vector<decltype(util::sqrt(T())), D>())
degreeToRadian(const Vector<T, D>& angle)
{
  constexpr BaseType<decltype(util::sqrt(T()))> conversionFactor = M_PI / 180.;
  return conversionFactor * angle;
}

template <typename T>
decltype(util::sqrt(T())) degreeToRadian(T angle)
{
  return degreeToRadian(Vector<T, 1>(angle))[0];
}

template <typename T, unsigned D>
decltype(Vector<decltype(util::sqrt(T())), D>())
radianToDegree(const Vector<T, D>& angle)
{
  constexpr BaseType<decltype(util::sqrt(T()))> conversionFactor = 180. / M_PI;
  return conversionFactor * angle;
}

template <typename T>
decltype(util::sqrt(T())) radianToDegree(T angle)
{
  return radianToDegree(Vector<T, 1>(angle))[0];
}

template <typename T, unsigned D>
Vector<T, utilities::dimensions::convert<D>::matrix> calculateRotationMatrix(
    const Vector<T, utilities::dimensions::convert<D>::rotation>& angle) any_platform
{
  Vector<T, utilities::dimensions::convert<D>::matrix> rotationMatrix;

  if constexpr (D == 2) {
    T const cos = util::cos(angle[0]);
    T const sin = util::sin(angle[0]);

    // row 1
    rotationMatrix[0] = cos;
    rotationMatrix[1] = -sin;
    // row 2
    rotationMatrix[2] = sin;
    rotationMatrix[3] = cos;
  }
  else {
    T const cos[3] = {util::cos(angle[0]), util::cos(angle[1]),
                      util::cos(angle[2])};
    T const sin[3] = {util::sin(angle[0]), util::sin(angle[1]),
                      util::sin(angle[2])};

    // |x0| / 0 1 2 \   |x1|
    // |y0| | 3 4 5 | = |x1|
    // |z0| \ 6 7 8 /   |x1|

    // row 1
    rotationMatrix[0] = cos[1] * cos[2];
    rotationMatrix[1] = sin[0] * sin[1] * cos[2] - cos[0] * sin[2];
    rotationMatrix[2] = cos[0] * sin[1] * cos[2] + sin[0] * sin[2];
    // row 2
    rotationMatrix[3] = cos[1] * sin[2];
    rotationMatrix[4] = sin[0] * sin[1] * sin[2] + cos[0] * cos[2];
    rotationMatrix[5] = cos[0] * sin[1] * sin[2] - sin[0] * cos[2];
    // row 3
    rotationMatrix[6] = -sin[1];
    rotationMatrix[7] = sin[0] * cos[1];
    rotationMatrix[8] = cos[0] * cos[1];
  }

  return rotationMatrix;
}

template <typename T>
Vector<T,4> inverseRotation(const Vector<T,4>& m) any_platform {
  return Vector{m[0], m[2],
                m[1], m[3]};
}

template <typename T>
Vector<T,9> inverseRotation(const Vector<T,9>& m) any_platform {
  return Vector{m[0], m[3], m[6],
                m[1], m[4], m[7],
                m[2], m[5], m[8]};
}

/// Compute v=M b
template <typename T, unsigned D>
Vector<T,D> matrixVectorProduct(
  const Vector<T,utilities::dimensions::convert<D>::matrix>& M,
  const Vector<T,D>& b) any_platform
{
  return Vector<T,D>([&](unsigned iD) -> T {
    T element{};
    for (unsigned jD=0; jD < D; ++jD) {
      element += M[iD*D + jD] * b[jD];
    }
    return element;
  });
}

/// Compute C=AB for square matrices
template <typename T, unsigned D>
Vector<T,utilities::dimensions::convert<D>::matrix> matrixMatrixProduct(
  const Vector<T,utilities::dimensions::convert<D>::matrix>& A,
  const Vector<T,utilities::dimensions::convert<D>::matrix>& B) any_platform
{
  Vector<T,utilities::dimensions::convert<D>::matrix> C{};
  for (unsigned i=0; i < D; ++i) {
    for (unsigned j=0; j < D; ++j) {
      for (unsigned k=0; k < D; ++k) {
        C[i*D + j] += A[i*D + k] * B[k*D + j];
      }
    }
  }
  return C;
}

template <typename T, unsigned D>
Vector<T, utilities::dimensions::convert<D>::matrix> invertRotationMatrix(
    const Vector<T, utilities::dimensions::convert<D>::matrix>& rotationMatrix)
{
  if constexpr (D == 2) {
    return Vector<T, 4>(rotationMatrix[0], rotationMatrix[2], rotationMatrix[1],
                        rotationMatrix[3]);
  }
  else {

    //Individual Entries ( source: https://doi.org/10.1016/j.compfluid.2018.02.027 )
    // |x0| / 0 3 6 \   |x1|
    // |y0| | 1 4 7 | = |x1|
    // |z0| \ 2 5 8 /   |x1|

    return Vector<T, 9>(rotationMatrix[0], rotationMatrix[3], rotationMatrix[6],
                        rotationMatrix[1], rotationMatrix[4], rotationMatrix[7],
                        rotationMatrix[2], rotationMatrix[5],
                        rotationMatrix[8]);
  }
  __builtin_unreachable();
}

template <typename T, unsigned D>
Vector<T, utilities::dimensions::convert<D>::matrix>
calculateInverseRotationMatrix(
    const Vector<T, utilities::dimensions::convert<D>::rotation>& angle)
{
  return invertRotationMatrix<T, D>(calculateRotationMatrix<T, D>(angle));
}

/// Rotates the input around the rotationCenter with a given rotationMatrix
template <typename T, unsigned D,
          bool OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN = false>
Vector<T, D> executeRotation(
    const Vector<T, D>&                                         input,
    const Vector<T, utilities::dimensions::convert<D>::matrix>& rotationMatrix,
    const Vector<T, D>& rotationCenter = Vector<T, D>(0.))
{
  const Vector<T, D> dist = input - rotationCenter;
  Vector<T, D>       rotated;

  if constexpr (D == 2) {
    rotated =
        Vector<T, 2>(dist[0] * rotationMatrix[0] + dist[1] * rotationMatrix[1],
                     dist[0] * rotationMatrix[2] + dist[1] * rotationMatrix[3]);
  }
  else {
    rotated =
        Vector<T, 3>(rotationMatrix[0] * dist[0] + rotationMatrix[1] * dist[1] +
                         rotationMatrix[2] * dist[2],
                     rotationMatrix[3] * dist[0] + rotationMatrix[4] * dist[1] +
                         rotationMatrix[5] * dist[2],
                     rotationMatrix[6] * dist[0] + rotationMatrix[7] * dist[1] +
                         rotationMatrix[8] * dist[2]);
  }

  if constexpr (!OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN) {
    return rotationCenter + rotated;
  }
  else {
    return rotated;
  }

  __builtin_unreachable();
}

/// Rotates the input around the rotationCenter with a given rotationMatrix in the opposite direction
template <typename T, unsigned D,
          bool OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN = false>
Vector<T, D> invertRotation(
    const Vector<T, D>&                                         input,
    const Vector<T, utilities::dimensions::convert<D>::matrix>& rotationMatrix,
    const Vector<T, D>& rotationCenter = Vector<T, D>(0.))
{
  const Vector<T, utilities::dimensions::convert<D>::matrix> invRotationMatrix =
      invertRotationMatrix<T, D>(rotationMatrix);
  return executeRotation(input, invRotationMatrix, rotationCenter);
}

/// Calculate local velocity
template <typename T, unsigned D>
constexpr Vector<T, D> calculateLocalVelocity(
    const Vector<T, D>& rotationCenter, const Vector<T, D>& velocity,
    const Vector<T, utilities::dimensions::convert<D>::rotation>&
                        angularVelocity,
    const Vector<T, D>& position)
{
  if constexpr (D == 2) {
    // two dimensions: u = U + w x r = (Ux, Uy, 0) + (0,0,w) x (X,Y,0) = (Ux, Uy, 0) + (-w*Y, w*X, 0)
    return Vector<T, 2>(
        velocity[0] - angularVelocity[0] * (position[1] - rotationCenter[1]),
        velocity[1] + angularVelocity[0] * (position[0] - rotationCenter[0]));
  }
  else {
    // three dimensions: u = U + w x r = (Ux, Uy, Uz) + (wx,wy,wz) x (X,Y,Z) = (Ux, Uy, Uz) + (wy*Z-wz*Y, wz*X-wx*Z, wx*Y-wy*X)
    return velocity +
           crossProduct3D(angularVelocity, position - rotationCenter);
  }
  __builtin_unreachable();
}

/// Rotate moment of inertia (mofi)
template <typename T>
constexpr Matrix<T, 3, 3> rotateMofi(const Vector<T, 3>& mofi,
                                     const Vector<T, 9>& rotationMatrix)
{
  // I' = R(angle) * I * R(angle)^T
  // TODO: The inertia tensor is symmetric - make use of it
  T data[3][3] = {
      {mofi[0],    T {},    T {}},
      {   T {}, mofi[1],    T {}},
      {   T {},    T {}, mofi[2]}
  };
  Matrix<T, 3, 3>       inertiaTensor(data);
  const Matrix<T, 3, 3> _rotationMatrix(rotationMatrix);
  inertiaTensor =
      _rotationMatrix * inertiaTensor * (_rotationMatrix.transpose());

  return inertiaTensor;
}

} // namespace util

} // namespace olb

#endif
