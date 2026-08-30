/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2024 Thomas Henn, Mathias J. Krause, Jonathan Jeppener-Haltenhoff, Christoph Gaul
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

#ifndef STL_READER_HH
#define STL_READER_HH


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "octree.hh"
#include "stlReader.h"

#ifdef FEATURE_VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
void STLtriangle<T>::init()
{
  Vector<T,3> A = point[0];
  Vector<T,3> B = point[1];
  Vector<T,3> C = point[2];
  Vector<T,3> b, c;
  T bb = 0., bc = 0., cc = 0.;

  for (int i = 0; i < 3; i++) {
    b[i] = B[i] - A[i];
    c[i] = C[i] - A[i];
    bb += b[i] * b[i];
    bc += b[i] * c[i];
    cc += c[i] * c[i];
  }

  normal[0] = b[1] * c[2] - b[2] * c[1];
  normal[1] = b[2] * c[0] - b[0] * c[2];
  normal[2] = b[0] * c[1] - b[1] * c[0];

  T norm = util::sqrt(
             util::pow(normal[0], 2) + util::pow(normal[1], 2) + util::pow(normal[2], 2));
  normal[0] /= norm;
  normal[1] /= norm;
  normal[2] /= norm;

  T D = 1.0 / (cc * bb - bc * bc);
  T bbD = bb * D;
  T bcD = bc * D;
  T ccD = cc * D;

  kBeta = 0.;
  kGamma = 0.;
  d = 0.;

  for (int i = 0; i < 3; i++) {
    uBeta[i] = b[i] * ccD - c[i] * bcD;
    uGamma[i] = c[i] * bbD - b[i] * bcD;
    kBeta -= A[i] * uBeta[i];
    kGamma -= A[i] * uGamma[i];
    d += A[i] * normal[i];
  }
}

template<typename T>
Vector<T,3> STLtriangle<T>::getCenter()
{
  Vector<T,3> center( T(0) );

  center[0] = (point[0][0] + point[1][0]
               + point[2][0]) / 3.;
  center[1] = (point[0][1] + point[1][1]
               + point[2][1]) / 3.;
  center[2] = (point[0][2] + point[1][2]
               + point[2][2]) / 3.;

  return center;
}

template<typename T>
std::vector<T> STLtriangle<T>::getE0()
{
  Vector<T,3> vec;
  vec[0] = point[0][0] - point[1][0];
  vec[1] = point[0][1] - point[1][1];
  vec[2] = point[0][2] - point[1][2];
  return vec;
}

template<typename T>
std::vector<T> STLtriangle<T>::getE1()
{
  Vector<T,3> vec;
  vec[0] = point[0][0] - point[2][0];
  vec[1] = point[0][1] - point[2][1];
  vec[2] = point[0][2] - point[2][2];
  return vec;
}

template<typename T>
bool STLtriangle<T>::isPointInside(const PhysR<T,3>& pt) const
{
  // tests with T=double and T=float show that the epsilon must be increased
  const T epsilon = std::numeric_limits<BaseType<T>>::epsilon()*T(10);

  const T beta = pt * uBeta + kBeta;
  const T gamma = pt * uGamma + kGamma;

  // check if approximately equal
  if ( util::nearZero(norm(pt - (point[0] + beta*(point[1]-point[0]) + gamma*(point[2]-point[0]))), epsilon) ) {
    const T alpha = T(1) - beta - gamma;
    return (beta >= T(0) || util::nearZero(beta, epsilon))
           && (gamma >= T(0) || util::nearZero(gamma, epsilon))
           && (alpha >= T(0) || util::nearZero(alpha, epsilon));
  }
  return false;
}

/* Schnitttest nach
 * http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf
 *
 * Creative Commons Namensnennung 3.0 Deutschland
 * http://creativecommons.org/licenses/by/3.0/de/
 *
 * P. Hofmann, 22. August 2010
 *
 */
template<typename T>
bool STLtriangle<T>::testRayIntersect(const Vector<T,3>& pt,
                                      const Vector<T,3>& dir,
                                      Vector<T,3>& q, T& alpha, const T& rad,
                                      bool print)
{
  T rn = 0.;
  Vector<T,3> testPt = pt + rad * normal;
  /* Vector<T,3> help; */

  for (int i = 0; i < 3; i++) {
    rn += dir[i] * normal[i];
  }
#ifdef OLB_DEBUG

  if (print) {
    std::cout << "Pt: " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
  }
  if (print)
    std::cout << "testPt: " << testPt[0] << " " << testPt[1] << " " << testPt[2]
              << std::endl;
  if (print)
    std::cout << "PosNeg: "
              << normal[0] * testPt[0] + normal[1] * testPt[1] + normal[2] * testPt[2]
              - d << std::endl;
  if (print)
    std::cout << "Normal: " << normal[0] << " " << normal[1] << " " << normal[2]
              << std::endl;
#endif

  // Schnitttest Flugrichtung -> Ebene
  if (util::fabs(rn) < std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 1" << std::endl;
    }
#endif
    return false;
  }
  alpha = d - testPt[0] * normal[0] - testPt[1] * normal[1] - testPt[2] * normal[2];
  //  alpha -= testPt[i] * normal[i];
  alpha /= rn;

  // Abstand Partikel Ebene
  if (alpha < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 2" << std::endl;
    }
#endif
    return false;
  }
  for (int i = 0; i < 3; i++) {
    q[i] = testPt[i] + alpha * dir[i];
  }
  T beta = kBeta;
  for (int i = 0; i < 3; i++) {
    beta += uBeta[i] * q[i];
  }
#ifdef OLB_DEBUG
  T dist = util::sqrt(
             util::pow(q[0] - testPt[0], 2) + util::pow(q[1] - testPt[1], 2)
             + util::pow(q[2] - testPt[2], 2));
#endif

  // Schnittpunkt q in der Ebene?
  if (beta < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG

    if (print) {
      std::cout << "FALSE 3 BETA " << beta << " DIST " << dist << std::endl;
    }
#endif
    return false;
  }
  T gamma = kGamma;
  for (int i = 0; i < 3; i++) {
    gamma += uGamma[i] * q[i];
  }
  if (gamma < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 4 GAMMA " << gamma << " DIST " << dist << std::endl;
    }
#endif
    return false;
  }
  if (1. - beta - gamma < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print)
      std::cout << "FALSE 5 VAL " << 1 - beta - gamma << " DIST " << dist
                << std::endl;
#endif
    return false;
  }
#ifdef OLB_DEBUG
  if (print) {
    std::cout << "TRUE" << " GAMMA " << gamma << " BETA " << beta << std::endl;
  }
#endif
  return true;
}

/**
 * computes closest Point in a triangle to another point.
 * source: Real-Time Collision Detection. Christer Ericson. ISBN-10: 1558607323
 */
template<typename T>
Vector<T,3> STLtriangle<T>::closestPtPointTriangle(
  const Vector<T,3>& pt) const
{

  const T nEps = -std::numeric_limits<T>::epsilon();
  const T Eps = std::numeric_limits<T>::epsilon();

  Vector<T,3> ab = point[1] - point[0];
  Vector<T,3> ac = point[2] - point[0];
  Vector<T,3> bc = point[2] - point[1];

  T snom = (pt - point[0])*ab;
  T sdenom = (pt - point[1])*(point[0] - point[1]);

  T tnom = (pt - point[0])*ac;
  T tdenom = (pt - point[2])*(point[0] - point[2]);

  if (snom < nEps && tnom < nEps) {
    return point[0];
  }

  T unom = (pt - point[1])*bc;
  T udenom = (pt - point[2])*(point[1] - point[2]);

  if (sdenom < nEps && unom < nEps) {
    return point[1];
  }
  if (tdenom < nEps && udenom < nEps) {
    return point[2];
  }

  T vc = normal*crossProduct3D(point[0] - pt, point[1] - pt);

  if (vc < nEps && snom > Eps && sdenom > Eps) {
    return point[0] + snom / (snom + sdenom) * ab;
  }

  T va = normal*crossProduct3D(point[1] - pt, point[2] - pt);

  if (va < nEps && unom > Eps && udenom > Eps) {
    return point[1] + unom / (unom + udenom) * bc;
  }

  T vb = normal*crossProduct3D(point[2] - pt, point[0] - pt);

  if (vb < nEps && tnom > Eps && tdenom > Eps) {
    return point[0] + tnom / (tnom + tdenom) * ac;
  }

  T u = va / (va + vb + vc);
  T v = vb / (va + vb + vc);
  T w = 1. - u - v;

  return u * point[0] + v * point[1] + w * point[2];
}

template<typename T>
bool STLtriangle<T>::getPointToEdgeDistances(const Vector<T,3>& input, Vector<T,3>& output, T sensitivity)
{
  auto P1 = point[0];
  auto P2 = point[1];
  output[0] = norm(crossProduct3D(P2-input,P1-P2))/norm(P1-P2);
  P1 = point[0];
  P2 = point[2];
  output[1] = norm(crossProduct3D(P2-input,P1-P2))/norm(P1-P2);
  P1 = point[1];
  P2 = point[2];
  output[2] = norm(crossProduct3D(P2-input,P1-P2))/norm(P1-P2);
  return true;
}

template<typename T>
bool STLtriangle<T>::isEdgePoint (const Vector<T,3> & input,  Vector<T,3>& P1, Vector<T,3> & P2,   T sensitivity )
{
  using namespace std;
  if(input[0] < sensitivity && input[1] >= sensitivity && input[2] >= sensitivity )
  {
    P1 = point[0];
    P2 = point[1];
    //cout << "edge point " << P1 << " " << P2 << std::endl;
    return true;
  }
  if(input[0] >= sensitivity && input[1] < sensitivity && input[2] >= sensitivity )
  {
    P1 = point[0];
    P2 = point[2];
    //cout << "edge point" << P1 << " " << P2 << std::endl;
    return true;
  }
  if(input[0] >= sensitivity && input[1] >= sensitivity && input[2] < sensitivity )
  {
    P1 = point[1];
    P2 = point[2];
   //cout << "edge point" << P1 << " " << P2 << std::endl;
    return true;
  }
  return false;
}

template<typename T>
bool STLtriangle<T>::isVortexPoint (const Vector<T,3> & input,Vector<T,3>& P, T sensitivity )
{
  olb::OstreamManager clout = OstreamManager(std::cout, "STLtriangle");
  using namespace std;
  if (input[0] < sensitivity && input[1] < sensitivity && input[2] < sensitivity )
  {
    cout << "Error isVortexPoint! Possible reduce sensitivity!" << std::endl;
    return false;
  }
  if(input[0] < sensitivity && input[1] < sensitivity)
  {
    P = point[0];
    //cout << "vortex point " << P << std::endl;
    return true;
  }
  if(input[0] < sensitivity && input[2] < sensitivity)
  {
    P = point[1];
    //cout << "vortex point " << P << std::endl;
    return true;
  }
  if(input[1] < sensitivity && input[2] < sensitivity)
  {
    P = point[2];
    //cout << "vortex point " << P << std::endl;
    return true;
  }
  return false;
}

template<typename T>
Vector<T,3> STLtriangle<T>::closestPointTo(Vector<T,3> physR) {
  const auto AB = point[1] - point[0];
  const auto AC = point[2] - point[0];
  const auto AP = physR     - point[0];

  const auto d1 = AB * AP;
  const auto d2 = AC * AP;
  if (d1 <= T{0} && d2 <= T{0}) {
    return point[0];
  }

  auto BP = physR - point[1];
  auto d3 = AB * BP;
  auto d4 = AC * BP;
  if (d3 >= T{0} && d4 <= d3) {
    return point[1];
  }

  auto vc = d1 * d4 - d3 * d2;
  if (vc <= T{0} && d1 >= T{0} && d3 <= T{0}) {
    auto v = d1 / (d1 - d3);
    return point[0] + v * AB;
  }

  auto CP = physR - point[2];
  auto d5 = AB * CP;
  auto d6 = AC * CP;
  if (d6 >= T{0} && d5 <= d6) {
    return point[2];
  }

  auto vb = d5 * d2 - d1 * d6;
  if (vb <= T{0} && d2 >= T{0} && d6 <= T{0}) {
    auto w = d2 / (d2 - d6);
    return point[0] + w * AC;
  }

  auto va = d3 * d6 - d5 * d4;
  if (va <= T{0} && (d4 - d3) >= T{0} && (d5 - d6) >= T{0}) {
    auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return point[1] + w * (point[2] - point[1]);
  }

  auto denom = 1.0 / (va + vb + vc);
  auto v = vb * denom;
  auto w = vc * denom;
  return point[0] + AB * v + AC * w;
}

template<typename T>
STLmesh<T>::STLmesh(std::string fName, T stlSize)
  : _fName(fName),
    _min(T()),
    _max(T()),
    _maxDist2(0),
    clout(std::cout, "STLmesh")
{
  std::ifstream f(fName.c_str(), std::ios::in);
  _triangles.reserve(10000);
  if (!f.good()) {
    throw std::runtime_error("STL File not valid.");
  }
  char buf[6];
  buf[5] = 0;
  f.read(buf, 5);
  const std::string asciiHeader = "solid";
  if (std::string(buf) == asciiHeader) {
    f.seekg(0, std::ios::beg);
    if (f.good()) {
      std::string s0, s1;
      int i = 0;
      while (!f.eof()) {
        f >> s0;
        if (s0 == "facet") {
          STLtriangle<T> tri;
          f >> s1 >> tri.normal[0] >> tri.normal[1] >> tri.normal[2];
          f >> s0 >> s1;
          f >> s0 >> tri.point[0][0] >> tri.point[0][1]
            >> tri.point[0][2];
          f >> s0 >> tri.point[1][0] >> tri.point[1][1]
            >> tri.point[1][2];
          f >> s0 >> tri.point[2][0] >> tri.point[2][1]
            >> tri.point[2][2];
          f >> s0;
          f >> s0;
          for (int k = 0; k < 3; k++) {
            tri.point[0][k] *= stlSize;
            tri.point[1][k] *= stlSize;
            tri.point[2][k] *= stlSize;
          }
          if (i == 0) {
            _min = T();
            _max = T();

            _min[0] = tri.point[0][0];
            _min[1] = tri.point[0][1];
            _min[2] = tri.point[0][2];

            _max[0] = tri.point[0][0];
            _max[1] = tri.point[0][1];
            _max[2] = tri.point[0][2];

            _min[0] = util::min(_min[0], tri.point[1][0]);
            _min[1] = util::min(_min[1], tri.point[1][1]);
            _min[2] = util::min(_min[2], tri.point[1][2]);

            _max[0] = util::max(_max[0], tri.point[1][0]);
            _max[1] = util::max(_max[1], tri.point[1][1]);
            _max[2] = util::max(_max[2], tri.point[1][2]);

            _min[0] = util::min(_min[0], tri.point[2][0]);
            _min[1] = util::min(_min[1], tri.point[2][1]);
            _min[2] = util::min(_min[2], tri.point[2][2]);

            _max[0] = util::max(_max[0], tri.point[2][0]);
            _max[1] = util::max(_max[1], tri.point[2][1]);
            _max[2] = util::max(_max[2], tri.point[2][2]);

          }
          else {
            _min[0] = util::min(_min[0], tri.point[0][0]);
            _min[1] = util::min(_min[1], tri.point[0][1]);
            _min[2] = util::min(_min[2], tri.point[0][2]);

            _max[0] = util::max(_max[0], tri.point[0][0]);
            _max[1] = util::max(_max[1], tri.point[0][1]);
            _max[2] = util::max(_max[2], tri.point[0][2]);

            _min[0] = util::min(_min[0], tri.point[1][0]);
            _min[1] = util::min(_min[1], tri.point[1][1]);
            _min[2] = util::min(_min[2], tri.point[1][2]);

            _max[0] = util::max(_max[0], tri.point[1][0]);
            _max[1] = util::max(_max[1], tri.point[1][1]);
            _max[2] = util::max(_max[2], tri.point[1][2]);

            _min[0] = util::min(_min[0], tri.point[2][0]);
            _min[1] = util::min(_min[1], tri.point[2][1]);
            _min[2] = util::min(_min[2], tri.point[2][2]);

            _max[0] = util::max(_max[0], tri.point[2][0]);
            _max[1] = util::max(_max[1], tri.point[2][1]);
            _max[2] = util::max(_max[2], tri.point[2][2]);
          }

          i++;
          tri.init();
          _triangles.push_back(tri);

          _maxDist2 = util::max(distPoints(tri.point[0], tri.point[1]),
                                _maxDist2);
          _maxDist2 = util::max(distPoints(tri.point[2], tri.point[1]),
                                _maxDist2);
          _maxDist2 = util::max(distPoints(tri.point[0], tri.point[2]),
                                _maxDist2);
        }
        else if (s0 == "endsolid") {
          break;
        }
      }
    }
  }
  else {
    f.close();
    f.open(fName.c_str(), std::ios::in | std::ios::binary);
    char comment[80];
    f.read(comment, 80);

    if (!f.good()) {
      throw std::runtime_error("STL File not valid.");
    }

    comment[79] = 0;
    int32_t nFacets;
    f.read(reinterpret_cast<char *>(&nFacets), sizeof(int32_t));

    if (!f.good()) {
      throw std::runtime_error("STL File not valid.");
    }

    float v[12];
    std::uint16_t uint16;
    for (int32_t i = 0; i < nFacets; ++i) {
      for (unsigned int j = 0; j < 12; ++j) {
        f.read(reinterpret_cast<char *>(&v[j]), sizeof(float));
      }
      f.read(reinterpret_cast<char *>(&uint16), sizeof(std::uint16_t));
      STLtriangle<T> tri;
      tri.normal[0] = v[0];
      tri.normal[1] = v[1];
      tri.normal[2] = v[2];
      tri.point[0][0] = v[3];
      tri.point[0][1] = v[4];
      tri.point[0][2] = v[5];
      tri.point[1][0] = v[6];
      tri.point[1][1] = v[7];
      tri.point[1][2] = v[8];
      tri.point[2][0] = v[9];
      tri.point[2][1] = v[10];
      tri.point[2][2] = v[11];

      for (int k = 0; k < 3; k++) {
        tri.point[0][k] *= stlSize;
        tri.point[1][k] *= stlSize;
        tri.point[2][k] *= stlSize;
      }
      if (i == 0) {
        _min[0] = tri.point[0][0];
        _min[1] = tri.point[0][1];
        _min[2] = tri.point[0][2];

        _max[0] = tri.point[0][0];
        _max[1] = tri.point[0][1];
        _max[2] = tri.point[0][2];

        _min[0] = util::min(_min[0], (T) tri.point[1][0]);
        _min[1] = util::min(_min[1], (T) tri.point[1][1]);
        _min[2] = util::min(_min[2], (T) tri.point[1][2]);

        _max[0] = util::max(_max[0], (T) tri.point[1][0]);
        _max[1] = util::max(_max[1], (T) tri.point[1][1]);
        _max[2] = util::max(_max[2], (T) tri.point[1][2]);

        _min[0] = util::min(_min[0], (T) tri.point[2][0]);
        _min[1] = util::min(_min[1], (T) tri.point[2][1]);
        _min[2] = util::min(_min[2], (T) tri.point[2][2]);

        _max[0] = util::max(_max[0], (T) tri.point[2][0]);
        _max[1] = util::max(_max[1], (T) tri.point[2][1]);
        _max[2] = util::max(_max[2], (T) tri.point[2][2]);

      }
      else {
        _min[0] = util::min(_min[0], (T) tri.point[0][0]);
        _min[1] = util::min(_min[1], (T) tri.point[0][1]);
        _min[2] = util::min(_min[2], (T) tri.point[0][2]);

        _max[0] = util::max(_max[0], (T) tri.point[0][0]);
        _max[1] = util::max(_max[1], (T) tri.point[0][1]);
        _max[2] = util::max(_max[2], (T) tri.point[0][2]);

        _min[0] = util::min(_min[0], (T) tri.point[1][0]);
        _min[1] = util::min(_min[1], (T) tri.point[1][1]);
        _min[2] = util::min(_min[2], (T) tri.point[1][2]);

        _max[0] = util::max(_max[0], (T) tri.point[1][0]);
        _max[1] = util::max(_max[1], (T) tri.point[1][1]);
        _max[2] = util::max(_max[2], (T) tri.point[1][2]);

        _min[0] = util::min(_min[0], (T) tri.point[2][0]);
        _min[1] = util::min(_min[1], (T) tri.point[2][1]);
        _min[2] = util::min(_min[2], (T) tri.point[2][2]);

        _max[0] = util::max(_max[0], (T) tri.point[2][0]);
        _max[1] = util::max(_max[1], (T) tri.point[2][1]);
        _max[2] = util::max(_max[2], (T) tri.point[2][2]);
      }
      tri.init();
      _triangles.push_back(tri);

      _maxDist2 = util::max(distPoints(tri.point[0], tri.point[1]), _maxDist2);
      _maxDist2 = util::max(distPoints(tri.point[2], tri.point[1]), _maxDist2);
      _maxDist2 = util::max(distPoints(tri.point[0], tri.point[2]), _maxDist2);
    }
  }
  f.close();
}

template<typename T>
STLmesh<T>::STLmesh(const std::vector<std::vector<T>> meshPoints, T stlSize)
  : _fName("meshPoints.stl"),
    _min(T()),
    _max(T()),
    _maxDist2(0),
    clout(std::cout, "STLmesh")
{
  _triangles.reserve(10000);
  for (size_t i = 0; i < meshPoints.size() / 3; i++) {
    STLtriangle<T> tri;
    tri.point[0][0] = meshPoints[i*3 + 0][0];
    tri.point[0][1] = meshPoints[i*3 + 0][1];
    tri.point[0][2] = meshPoints[i*3 + 0][2];

    tri.point[1][0] = meshPoints[i*3 + 1][0];
    tri.point[1][1] = meshPoints[i*3 + 1][1];
    tri.point[1][2] = meshPoints[i*3 + 1][2];

    tri.point[2][0] = meshPoints[i*3 + 2][0];
    tri.point[2][1] = meshPoints[i*3 + 2][1];
    tri.point[2][2] = meshPoints[i*3 + 2][2];
    for (int k = 0; k < 3; k++) {
      tri.point[0][k] *= stlSize;
      tri.point[1][k] *= stlSize;
      tri.point[2][k] *= stlSize;
    }
    if (i == 0) {
      _min*=T();
      _max*=T();

      _min[0] = tri.point[0][0];
      _min[1] = tri.point[0][1];
      _min[2] = tri.point[0][2];

      _max[0] = tri.point[0][0];
      _max[1] = tri.point[0][1];
      _max[2] = tri.point[0][2];

      _min[0] = util::min(_min[0], (T) tri.point[1][0]);
      _min[1] = util::min(_min[1], (T) tri.point[1][1]);
      _min[2] = util::min(_min[2], (T) tri.point[1][2]);

      _max[0] = util::max(_max[0], (T) tri.point[1][0]);
      _max[1] = util::max(_max[1], (T) tri.point[1][1]);
      _max[2] = util::max(_max[2], (T) tri.point[1][2]);

      _min[0] = util::min(_min[0], (T) tri.point[2][0]);
      _min[1] = util::min(_min[1], (T) tri.point[2][1]);
      _min[2] = util::min(_min[2], (T) tri.point[2][2]);

      _max[0] = util::max(_max[0], (T) tri.point[2][0]);
      _max[1] = util::max(_max[1], (T) tri.point[2][1]);
      _max[2] = util::max(_max[2], (T) tri.point[2][2]);

    }
    else {
      _min[0] = util::min(_min[0], (T) tri.point[0][0]);
      _min[1] = util::min(_min[1], (T) tri.point[0][1]);
      _min[2] = util::min(_min[2], (T) tri.point[0][2]);

      _max[0] = util::max(_max[0], (T) tri.point[0][0]);
      _max[1] = util::max(_max[1], (T) tri.point[0][1]);
      _max[2] = util::max(_max[2], (T) tri.point[0][2]);

      _min[0] = util::min(_min[0], (T) tri.point[1][0]);
      _min[1] = util::min(_min[1], (T) tri.point[1][1]);
      _min[2] = util::min(_min[2], (T) tri.point[1][2]);

      _max[0] = util::max(_max[0], (T) tri.point[1][0]);
      _max[1] = util::max(_max[1], (T) tri.point[1][1]);
      _max[2] = util::max(_max[2], (T) tri.point[1][2]);

      _min[0] = util::min(_min[0], (T) tri.point[2][0]);
      _min[1] = util::min(_min[1], (T) tri.point[2][1]);
      _min[2] = util::min(_min[2], (T) tri.point[2][2]);

      _max[0] = util::max(_max[0], (T) tri.point[2][0]);
      _max[1] = util::max(_max[1], (T) tri.point[2][1]);
      _max[2] = util::max(_max[2], (T) tri.point[2][2]);
    }

    tri.init();
    _triangles.push_back(tri);

    _maxDist2 = util::max(distPoints(tri.point[0], tri.point[1]),
                          _maxDist2);
    _maxDist2 = util::max(distPoints(tri.point[2], tri.point[1]),
                          _maxDist2);
    _maxDist2 = util::max(distPoints(tri.point[0], tri.point[2]),
                          _maxDist2);
  }
}

template<typename T>
T STLmesh<T>::distPoints(Vector<T,3>& p1, Vector<T,3>& p2)
{
  return util::pow(T(p1[0] - p2[0]), 2)
         + util::pow(T(p1[1] - p2[1]), 2)
         + util::pow(T(p1[2] - p2[2]), 2);
}

template<typename T>
void STLmesh<T>::print(bool full)
{
  if (full) {
    int i = 0;
    clout << "Triangles: " << std::endl;
    typename std::vector<STLtriangle<T> >::iterator it = _triangles.begin();

    for (; it != _triangles.end(); ++it) {
      clout << i++ << ": " << it->point[0][0] << " " << it->point[0][1]
            << " " << it->point[0][2] << " | " << it->point[1][0] << " "
            << it->point[1][1] << " " << it->point[1][2] << " | "
            << it->point[2][0] << " " << it->point[2][1] << " "
            << it->point[2][2] << std::endl;
    }
  }
  clout << "nTriangles=" << _triangles.size() << "; maxDist2=" << _maxDist2
        << std::endl;
  clout << "minPhysR(StlMesh)=(" << getMin()[0] << "," << getMin()[1] << ","
        << getMin()[2] << ")";
  clout << "; maxPhysR(StlMesh)=(" << getMax()[0] << "," << getMax()[1] << ","
        << getMax()[2] << ")" << std::endl;
}

template<typename T>
void STLmesh<T>::write(std::string fName)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    std::string fullName = singleton::directories().getVtkOutDir() + fName
                           + ".stl";
    std::ofstream f(fullName.c_str());
    f << "solid ascii " << fullName << "\n";

    for (unsigned int i = 0; i < _triangles.size(); i++) {
      f << "facet normal " << _triangles[i].normal[0] << " "
        << _triangles[i].normal[1] << " " << _triangles[i].normal[2] << "\n";
      f << "    outer loop\n";
      f << "        vertex " << _triangles[i].point[0][0] << " "
        << _triangles[i].point[0][1] << " " << _triangles[i].point[0][2]
        << "\n";
      f << "        vertex " << _triangles[i].point[1][0] << " "
        << _triangles[i].point[1][1] << " " << _triangles[i].point[1][2]
        << "\n";
      f << "        vertex " << _triangles[i].point[2][0] << " "
        << _triangles[i].point[2][1] << " " << _triangles[i].point[2][2]
        << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f << "endsolid\n";
    f.close();
  }
  /*if (_verbose)*/clout << "Write ... OK" << std::endl;
}

template<typename T>
bool STLmesh<T>::testRayIntersect(const std::set<unsigned int>& tris, const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha)
{
  std::set<unsigned int>::iterator it = tris.begin();
  for (; it != tris.end(); ++it) {
    if (_triangles[*it].testRayIntersect(pt, dir, q, alpha) && alpha < 1) {
      return true;
    }
  }
  return false;
}


/*
 * STLReader functions
 */
template<typename T>
STLreader<T>::STLreader(const std::string fName, T voxelSize, T stlSize,
                        RayMode method, bool verbose, T overlap, T max)
  : _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName(fName),
    _mesh(fName, stlSize),
    _verbose(verbose),
    clout(std::cout, "STLreader")
{
  this->getName() = "STLreader";

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
  case RayMode::Robust:
    indicate1();
    break;
  case RayMode::DoubleRay:
    indicate3();
    break;
  case RayMode::FastRayX:
    indicate2_Xray();
    break;
  case RayMode::FastRayY:
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

/*
 * STLReader functions
 */
template<typename T>
STLreader<T>::STLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize,
                        RayMode method, bool verbose, T overlap, T max)
  : _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName("meshPoints.stl"),
    _mesh(meshPoints, stlSize),
    _verbose(verbose),
    clout(std::cout, "STLreader")
{
  this->getName() = "STLreader";

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
  case RayMode::Robust:
    indicate1();
    break;
  case RayMode::DoubleRay:
    indicate3();
    break;
  case RayMode::FastRayX:
    indicate2_Xray();
    break;
  case RayMode::FastRayY:
    indicate2_Yray();
    break;
  default:
    indicate2();
    break;
  }

  setNormalsOutside();

  if (_verbose) {
    print();
  }
  if (_verbose) {
    clout << "Voxelizing ... OK" << std::endl;
  }
}

template<typename T>
STLreader<T>::~STLreader()
{
  delete _tree;
}

/*
 *  Old indicate function (slower, more stable)
 *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
 *  with STL for each ray. Odd number of intersection means inside (Majority vote).
 */

template<typename T>
void STLreader<T>::indicate1()
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
void STLreader<T>::indicate2()
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
void STLreader<T>::indicate2_Xray()
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
void STLreader<T>::indicate2_Yray()
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
void STLreader<T>::indicate3()
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
bool STLreader<T>::operator() (bool output[], const T input[])
{
  output[0] = false;
  if (isInsideRootTree(input)) {
    std::vector<T> tmp(input, input + 3);
    output[0] = _tree->find(tmp)->getInside();
  }
  return true;
}


template<typename T>
bool STLreader<T>::isInsideRootTree(const T input[])
{
  T coords = _tree->getRadius();
  Vector<T,3> c(_tree->getCenter());
  return c[0] - coords < input[0] && input[0] < c[0] + coords && c[1] - coords < input[1]
      && input[1] < c[1] + coords && c[2] - coords < input[2] && input[2] < c[2] + coords;
}


template<typename T>
Vector<T,3> STLreader<T>::closestPointInBoundingBox(const Vector<T,3>& input)
{
  T coords = _tree->getRadius();
  Vector<T,3> c(_tree->getCenter());
  Vector<T,3> closestPoint;
  for(int i = 0; i < 3; ++i) {
    closestPoint[i] = util::max(c[i] - coords, util::min(c[i] + coords, input[i]));
  }
  return closestPoint;
}


template<typename T>
bool STLreader<T>::distance(T& distance, const Vector<T,3>& origin,
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
      if (tmpNode) {
        tmpNode->intersectRayNode(pt, dir, s);
        for (int i = 0; i < 3; i++) {
          pt[i] = s[i] + step * dir[i];
        }
      }
    }
  }

  return false;
}

template<typename T>
template<typename F>
void STLreader<T>::iterateOverCloseTriangles(const PhysR<T,3>& pt, F func, Octree<T>* leafNode) {
  // Find the leaf node in the Octree that contains the input point
  leafNode = _tree->find(pt);

  if (leafNode && !leafNode->getTriangles().empty()) {
    const std::vector<unsigned int>& triangleIndices = leafNode->getTriangles();
    for (unsigned int idx : triangleIndices) {
      const STLtriangle<T>& triangle = _mesh.getTri(idx);
      func(triangle);
    }
  }
  else {
    for (const STLtriangle<T>& triangle : _mesh.getTriangles()) {
      func(triangle);
    }
  }
};

template<typename T>
Vector<T,3> STLreader<T>::evalNormalOnSurface(const PhysR<T,3>& pt, const Vector<T,3>& fallbackNormal)
{
  // Check if the position is on the corner of a triangle
  unsigned countTriangles = 0;
  Vector<T,3> normal(T(0));

  iterateOverCloseTriangles(pt, [&](const STLtriangle<T>& triangle){
      // TODO: Calculate angle-weighted psuedonormal (see 10.1109/TVCG.2005.49) in case the point lies on corners
      // Edges correspond to an unweighted average (as calculated below) anyway
      if (triangle.isPointInside(pt)) {
        ++countTriangles;
        normal+=triangle.getNormal();
      }
    });

  if (countTriangles > 0) {
    return normal / countTriangles;
  }

  // If the provided point isn't located on the surface, return the predefined fallback
  return fallbackNormal;
}

template<typename T>
template<SignMode SIGNMODE>
Vector<T,3> STLreader<T>::evalSurfaceNormal(const Vector<T,3>& origin)
{
  Vector<T,3> normal(0.);
  Vector<T,3> closestPointOnSurface(0.);
  const STLtriangle<T>* closestTriangle = nullptr;
  T distance = std::numeric_limits<T>::max();
  const PhysR<T,3> pt(closestPointInBoundingBox(origin));

  iterateOverCloseTriangles(pt, [&](const STLtriangle<T>& triangle) {
      PhysR<T,3> const pointOnTriangle = triangle.closestPtPointTriangle(pt);
      PhysR<T,3> const currDistance = pt - pointOnTriangle;
      T currDistanceNorm = norm_squared(currDistance);
      if (distance > currDistanceNorm) {
        normal = currDistance;
        distance = currDistanceNorm;
        closestPointOnSurface = pointOnTriangle;
        closestTriangle = &triangle;
      }
    });

  distance = util::sqrt(distance);
  if (!util::nearZero(distance)) {
    normal = normal/distance;
    normal *= evalSignForSignedDistance<SIGNMODE>(origin, distance);
  }
  else {
    normal = evalNormalOnSurface(closestPointOnSurface, closestTriangle->getNormal());
    // The below isn't necessary because all normals are set to point outside by default (see constructor)

    /* const T tmpStepSize = _voxelSize * T{0.01}; */
    /* const short signBefore = evalSignForSignedDistance(origin - tmpStepSize * normal, *closestTriangle); */
    /* const short signAfter  = evalSignForSignedDistance(origin + tmpStepSize * normal, *closestTriangle); */

    /* if(signBefore > 0 || signAfter < 0) { */
    /*   normal *= T{-1}; */
    /* } */
  }

  return normal;
}


template<typename T>
template<SignMode SIGNMODE>
Vector<T,3> STLreader<T>::evalSurfaceNormalForPseudoNormal(const Vector<T,3>& origin, Vector<T,3> & outputPointOnSurface)
{
  Vector<T,3> normal(0.);
  Vector<T,3> closestPointOnSurface(0.);
  const STLtriangle<T>* closestTriangle = nullptr;
  T distance = std::numeric_limits<T>::max();
  const PhysR<T,3> pt(closestPointInBoundingBox(origin));

  iterateOverCloseTriangles(pt, [&](const STLtriangle<T>& triangle) {
      PhysR<T,3> const pointOnTriangle = triangle.closestPtPointTriangle(pt);
      PhysR<T,3> const currDistance = pt - pointOnTriangle;
      T currDistanceNorm = norm_squared(currDistance);
      if (distance > currDistanceNorm) {
        normal = currDistance;
        distance = currDistanceNorm;
        closestPointOnSurface = pointOnTriangle;
        closestTriangle = &triangle;
      }
    });

/*   distance = util::sqrt(distance);
  if (!util::nearZero(distance)) {
    normal = normal/distance;
    normal *= evalSignForSignedDistance<SIGNMODE>(origin, distance);
  }
  else { */
    normal = evalNormalOnSurface(closestPointOnSurface, closestTriangle->getNormal());
    outputPointOnSurface = closestPointOnSurface;

    // The below isn't necessary because all normals are set to point outside by default (see constructor)

    /* const T tmpStepSize = _voxelSize * T{0.01}; */
    /* const short signBefore = evalSignForSignedDistance(origin - tmpStepSize * normal, *closestTriangle); */
    /* const short signAfter  = evalSignForSignedDistance(origin + tmpStepSize * normal, *closestTriangle); */

    /* if(signBefore > 0 || signAfter < 0) { */
    /*   normal *= T{-1}; */
    /* } */
  //}

  return normal;
}


template <typename T>
template <SignMode SIGNMODE>
short STLreader<T>::evalSignForSignedDistance(
    const Vector<T,3>& pt, [[maybe_unused]] const T distance, Vector<T,3>  vecdist, STLtriangle<T> stlT)
{
  short sign;

  if constexpr (SIGNMODE == SignMode::EXACT)
  {

    if(distance < _voxelSize) {

    sign = evalSignForSignedDistanceFromWindingNumber(pt);
    //sign = evalSignForSignedDistanceFromNormal(stlT.getNormal(), vecdist);
   /* Vector<T,3> closestPointOnSurface (0.);
    Vector<T,3> normal = evalSurfaceNormalForPseudoNormal<SignMode::EXACT>(pt,closestPointOnSurface);
      sign = evalSignForSignedDistanceFromNormal(normal, pt - closestPointOnSurface);*/
    }
    else
    {
    sign = evalSignForSignedDistanceFromCache(pt);
    }
  }
  else {
  sign = evalSignForSignedDistanceFromCache(pt);
  }
  return sign;
}

template<typename T>
short STLreader<T>::evalSignForSignedDistanceFromPseudonormal(const Vector<T,3>& pseudonormal, const Vector<T,3>& distance)
{
  const T projection = pseudonormal * distance;

  if(projection > 0) {
    return 1;
  }
  else if (projection < 0 ) {
    return -1;
  }
  return 0;
}


//aproximation not working for concave surfaces
template<typename T>
short STLreader<T>::evalSignForSignedDistanceFromNormal(const Vector<T,3>& normal, const Vector<T,3>& distance)
{
  const T projection = normal * distance;

  if(projection > 0) {
    return 1;
  }
  else if (projection < 0 ) {
    return -1;
  }
  return 0;
}




template<typename T>
short STLreader<T>::evalSignForSignedDistanceFromWindingNumber(const Vector<T,3>& pt)
{
  T windingNumber{0};

  for (const STLtriangle<T>& triangle : _mesh.getTriangles()) {
    const PhysR<T,3> a = triangle.point[0] - pt;
    const PhysR<T,3> b = triangle.point[1] - pt;
    const PhysR<T,3> c = triangle.point[2] - pt;

    const T aNorm = norm(a);
    const T bNorm = norm(b);
    const T cNorm = norm(c);

    const T numerator = a * crossProduct3D(b, c);
    const T denominator = aNorm * bNorm * cNorm + (a*b) * cNorm + (b*c) * aNorm + (c*a) * bNorm;

    windingNumber += util::atan2(numerator,denominator);
  }
  windingNumber *= 2;

  if(int(util::round(windingNumber))%2 == 0) {
    return 1;
  }
  return -1;
}


template <typename T>
T STLreader<T>::signedDistance(const Vector<T, 3>& input)
{
  return signedDistance<SignMode::CACHED>(input);
}

template <typename T>
T STLreader<T>::signedDistanceExact(const Vector<T, 3>& input)
{
  return signedDistance<SignMode::EXACT>(input);
}




template <typename T>
template <SignMode SIGNMODE>
T STLreader<T>::signedDistance(const Vector<T, 3>& input)
{
  T distanceNorm = std::numeric_limits<T>::max();
  //Vector<T,3>  distVec  (std::numeric_limits<T>::max(),std::numeric_limits<T>::max(),std::numeric_limits<T>::max());
  T extraDistance{};
  //STLtriangle<T> stlT;

  const auto evalDistance = [&](const Vector<T,3>& pt) {
    iterateOverCloseTriangles(pt, [&](const STLtriangle<T>& triangle) {
        const PhysR<T, 3> currPointOnTriangle =
            triangle.closestPtPointTriangle(pt);
        const Vector<T, 3> currDistance     = pt - currPointOnTriangle;
        T                  currDistanceNorm = norm_squared(currDistance);
        distanceNorm = util::min(distanceNorm, currDistanceNorm);
       /* if(distanceNorm -currDistanceNorm < std::numeric_limits<T>::min())
       {distVec = currDistanceNorm;
       stlT = triangle;}*/
      });
    return util::sqrt(distanceNorm);
  };

  if(!isInsideRootTree(input.data())) {
    const PhysR<T,3> ptInsideTree = closestPointInBoundingBox(input);
    const T distance = evalDistance(ptInsideTree);
    extraDistance = norm(ptInsideTree - input);
    return distance + extraDistance;
  }
  const T distance = evalDistance(input);

  return evalSignForSignedDistance<SIGNMODE>(input, distance) * distance;
}


template<typename T>
short STLreader<T>::evalSignForSignedDistanceFromCache(const Vector<T,3>& pt)
{
  bool isInside;
  this->operator()(&isInside, pt.data());
  return (isInside ? -1 : 1);
}



template <typename T>
Vector<T,3> STLreader<T>::surfaceNormal(const Vector<T,3>& pos, const T meshSize)
{
  return surfaceNormal<SignMode::CACHED>(pos, meshSize);
}

template <typename T>
Vector<T,3> STLreader<T>::surfaceNormalExact(const Vector<T,3>& pos, const T meshSize)
{
  return surfaceNormal<SignMode::EXACT>(pos, meshSize);
}

template <typename T>
template <SignMode SIGNMODE>
Vector<T,3> STLreader<T>::surfaceNormal(const Vector<T,3>& pos, const T meshSize)
{
  return evalSurfaceNormal<SIGNMODE>(pos);
}

template<typename T>
void STLreader<T>::print()
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
void STLreader<T>::writeOctree()
{
  _tree->write(_fName);
}

template<typename T>
void STLreader<T>::writeSTL(std::string stlName)
{
  if (stlName == "") {
    _mesh.write(_fName);
  }
  else {
    _mesh.write(stlName);
  }
}

template<typename T>
void STLreader<T>::setNormalsOutside()
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
void STLreader<T>::setBoundaryInsideNodes()
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
