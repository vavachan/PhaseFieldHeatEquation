/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, 2020 Adrian Kummerlaender, 2024 Yuji Shimojima
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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef UTIL_H
#define UTIL_H

#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include "core/baseType.h"
#include "core/vector.h"
#include "meta.h"
#include "utilities/omath.h"
#include "utilities/oalgorithm.h"
#include "descriptor/functions.h"

namespace olb {

namespace util {

// *INDENT-OFF*

template<typename T> T norm(const std::vector<T>& a);

template <typename T>
inline int sign(T val) any_platform
{
  return (0 < val) - (val < 0);
}

template <typename T>
inline bool aligned_to_x(const std::vector<T>& vec)
{
  return (vec[0]!=0 and vec[1]==0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_y(const std::vector<T>& vec)
{
  return (vec[0]==0 and vec[1]!=0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_z(const std::vector<T>& vec)
{
  return (vec[0]==0 and vec[1]==0 and vec[2]!=0);
}

template <typename T>
inline bool aligned_to_grid(const std::vector<T>& vec)
{
  return (aligned_to_x<T>(vec) or
          aligned_to_y<T>(vec) or
          aligned_to_z<T>(vec));
}


inline bool intersect (
  int x0, int x1, int y0, int y1,
  int x0_, int x1_, int y0_, int y1_,
  int& newX0, int& newX1, int& newY0, int& newY1 )
{
  newX0 = util::max(x0,x0_);
  newY0 = util::max(y0,y0_);

  newX1 = util::min(x1,x1_);
  newY1 = util::min(y1,y1_);

  return newX1>=newX0 && newY1>=newY0;
}

inline bool intersect (
  int x0, int x1, int y0, int y1, int z0, int z1,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int& newX0, int& newX1, int& newY0, int& newY1, int& newZ0, int& newZ1 )
{
  newX0 = util::max(x0,x0_);
  newY0 = util::max(y0,y0_);
  newZ0 = util::max(z0,z0_);

  newX1 = util::min(x1,x1_);
  newY1 = util::min(y1,y1_);
  newZ1 = util::min(z1,z1_);

  return newX1>=newX0 && newY1>=newY0 && newZ1>=newZ0;
}

inline bool contained(int x, int y,
                      int x0, int x1, int y0, int y1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1;
}

inline bool contained(int x, int y, int z,
                      int x0, int x1, int y0, int y1, int z0, int z1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1 &&
         z>=z0 && z<=z1;
}


template<typename T>
any_platform T sqr(T arg) any_platform
{
  return arg*arg;
}

template<typename T>
any_platform T cube(T arg)
{
  return arg*arg*arg;
}

/// Compute norm square of a d-dimensional vector
template<typename ARRAY_LIKE, unsigned D>
auto normSqr(const ARRAY_LIKE& u) any_platform
{
  auto uSqr = decltype(u[0]){};
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<unsigned D, typename ARRAY_LIKE>
auto norm(const ARRAY_LIKE& u) any_platform
{
  return sqrt(normSqr<ARRAY_LIKE,D>(u));
}

template<typename T, unsigned D>
T normSqr(const T* u) any_platform
{
  auto uSqr = T{};
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<typename T>
T normSqr(std::initializer_list<T> data)
{
  return std::inner_product(data.begin(), data.end(), data.begin(), T(0));
}


/// Compute norm square of a d-dimensional vector
template<typename T, unsigned D, typename IMPL>
T normSqr(const ScalarVector<T,D,IMPL>& u) any_platform
{
  T uSqr = T();
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}



template<typename T, int d>
T scalarProduct(const Vector<T,d>& u1, const Vector<T,d>& u2)  any_platform
{
  T prod = T();
  if (u1.size() == u2.size()) {
    for (unsigned iD=0; iD < u1.size(); ++iD) {
      prod += u1[iD]*u2[iD];
    }
  }
  return prod;
}

template<typename T>
T scalarProduct(const std::vector<T>& u1, const std::vector<T>& u2)
{
  T prod = T();
  if (u1.size() == u2.size()) {
    for (unsigned iD=0; iD<u1.size(); ++iD) {
      prod += u1[iD]*u2[iD];
    }
  }
  return prod;
}

/// Compute number of elements of a symmetric d-dimensional tensor
template <typename DESCRIPTORBASE>
struct TensorVal {
  static constexpr int n = (DESCRIPTORBASE::d*(DESCRIPTORBASE::d+1))/2; ///< result stored in n
};

/// Return array of population indices where c[iVel] == value
template <typename DESCRIPTOR, unsigned iVel, int value>
constexpr auto populationsContributingToVelocity() any_platform {
  constexpr auto velocity_matches_value = [](unsigned iPop) -> bool {
    return descriptors::c<DESCRIPTOR>(iPop,iVel) == value;
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_value));
};

/// Return array of population indices where c[iPop][iD] == NORMAL[iD]
template <typename DESCRIPTOR, int... NORMAL>
constexpr auto populationsContributingToDirection() any_platform {
  constexpr auto velocity_matches_direction = [](unsigned iPop) -> bool {
    return meta::indexed_pack_contains<NORMAL...>([iPop](unsigned iD, int x) -> bool {
      return x != 0 && x == descriptors::c<DESCRIPTOR>(iPop,iD);
    });
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_direction));
};

template <typename DESCRIPTORBASE>
int findVelocity(const int v[DESCRIPTORBASE::d]) any_platform
{
  for (int iPop=0; iPop<DESCRIPTORBASE::q; ++iPop) {
    bool fit = true;
    for (int iD=0; iD<DESCRIPTORBASE::d; ++iD) {
      if (descriptors::c<DESCRIPTORBASE>(iPop,iD) != v[iD]) {
        fit = false;
        break;
      }
    }
    if (fit) {
      return iPop;
    }
  }
  return DESCRIPTORBASE::q;
}

/// Compute opposites of wall-incoming population indices
template <typename DESCRIPTOR, int direction, int orientation>
constexpr auto subIndexOutgoing() any_platform {
  // Should be expressed in terms of populationsContributingToVelocity
  // but std::array / std::index_sequence distinction makes this ugly
  // To be revisitited
  constexpr auto velocity_matches_value = [](unsigned iPop) {
    return descriptors::c<DESCRIPTOR>(iPop,direction) == orientation;
  };
  constexpr auto opposite_population = [](unsigned iPop) {
    return descriptors::opposite<DESCRIPTOR>(iPop);
  };
  return meta::array_from_index_sequence(
    meta::map_index_sequence(opposite_population,
                             descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_value)));
}

template <typename DESCRIPTOR, int direction, int orientation>
constexpr auto subIndexOutgoingRemaining() any_platform {
  constexpr auto velocity_matches = [](unsigned iPop) {
    for (unsigned jPop : subIndexOutgoing<DESCRIPTOR,direction,orientation>()) {
      if (iPop == jPop) {
        return false;
      }
    }
    return true;
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches));
}

template <typename DESCRIPTOR, int plane, int normal1, int normal2>
constexpr auto subIndexOutgoing3DonEdges() any_platform {
  constexpr auto velocity_matches = [](unsigned iPop) {
           if constexpr (plane == 0) {
      return ( descriptors::c<DESCRIPTOR>(iPop,0) * 0
             + descriptors::c<DESCRIPTOR>(iPop,1) * (normal1 == 1 ? 1 : -1)
             + descriptors::c<DESCRIPTOR>(iPop,2) * (normal2 == 1 ? 1 : -1)) < 0;
    } else if constexpr (plane == 1) {
      return ( descriptors::c<DESCRIPTOR>(iPop,0) * (normal1 == 1 ? 1 : -1)
             + descriptors::c<DESCRIPTOR>(iPop,1) * 0
             + descriptors::c<DESCRIPTOR>(iPop,2) * (normal2 == 1 ? 1 : -1)) < 0;
    } else if constexpr (plane == 2) {
      return ( descriptors::c<DESCRIPTOR>(iPop,0) * (normal1 == 1 ? 1 : -1)
             + descriptors::c<DESCRIPTOR>(iPop,1) * (normal2 == 1 ? 1 : -1)
             + descriptors::c<DESCRIPTOR>(iPop,2) * 0) < 0;
    } else {
      return false;
    }
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches));
}

// For 2D Corners
template <typename DESCRIPTOR, int normalX, int normalY>
constexpr auto subIndexOutgoing2DonCorners() any_platform {
  constexpr auto velocity_matches = [](unsigned iPop) {
    return ( descriptors::c<DESCRIPTOR>(iPop,0) * normalX
           + descriptors::c<DESCRIPTOR>(iPop,1) * normalY) < 0;
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches));
}

// For 3D Corners
template <typename DESCRIPTOR, int normalX, int normalY, int normalZ>
constexpr auto subIndexOutgoing3DonCorners() any_platform {
  constexpr auto velocity_matches = [](unsigned iPop) {
    return ( descriptors::c<DESCRIPTOR>(iPop,0) * normalX
           + descriptors::c<DESCRIPTOR>(iPop,1) * normalY
           + descriptors::c<DESCRIPTOR>(iPop,2) * normalZ) < 0;
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches));
}

///finds all the remaining indexes of a lattice given some other indexes
template <typename DESCRIPTORBASE>
std::vector<int> remainingIndexes(const std::vector<int> &indices)
{
  std::vector<int> remaining;
  for (int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop) {
    bool found = false;
    for (unsigned jPop = 0; jPop < indices.size(); ++jPop) {
      if (indices[jPop] == iPop) {
        found = true;
      }
    }
    if (!found) {
      remaining.push_back(iPop);
    }
  }
  return remaining;
}


/// Util Function for Wall Model of Malaspinas
/// get link with smallest angle to a vector
template <typename T, typename DESCRIPTOR>
int get_nearest_link(const std::vector<T>& vec)
{
  T max=-1;
  int max_index = 0;
  for (int iQ=1; iQ<DESCRIPTOR::q; ++iQ) {
    std::vector<T> c_i(DESCRIPTOR::c[iQ], DESCRIPTOR::c[iQ]+3);
    T tmp = util::scalarProduct<T>(c_i, vec)/util::norm(c_i);
    if (tmp > max) {
      max = tmp;
      max_index = iQ;
    }
  }
  return max_index;
}


namespace tensorIndices2D {
enum { xx=0, xy=1, yy=2 };
}

namespace tensorIndices3D {
enum { xx=0, xy=1, xz=2, yy=3, yz=4, zz=5 };
}

/// compute lattice density from lattice pressure
template <typename T, typename DESCRIPTOR>
any_platform T densityFromPressure( T latticePressure )
{
  // rho = p / c_s^2 + 1
  return latticePressure * descriptors::invCs2<T,DESCRIPTOR>() + 1.0;
}

/// compute lattice pressure from lattice density
template <typename T, typename DESCRIPTOR>
any_platform T pressureFromDensity( T latticeDensity )
{
  // p = (rho - 1) * c_s^2
  return (latticeDensity - 1.0) / descriptors::invCs2<T,DESCRIPTOR>();
}

/// return true if a is close to zero
template <typename T>
inline bool nearZero(T a) any_platform
{
  if (a==T()) {
    return true;
  }
  T EPSILON = std::numeric_limits<T>::epsilon();
  if (a > -EPSILON && a < EPSILON) {
    return true;
  }
  else {
    return false;
  }
}

template <typename T>
inline bool nearZero(T a, T epsilon) any_platform
{
  if (a > -epsilon && a < epsilon) {
    return true;
  }
  else {
    return false;
  }
}

template<typename T, typename U=T, typename W=T>
inline bool approxEqual(T a, U b, W epsilon) any_platform
{
  if (a==b) {
    return true;
  }
  return nearZero<T>(a - b, epsilon);
}

template<typename T, typename U=T>
inline bool approxEqual(T a, U b) any_platform
{
  if (a==b) {
    return true;
  }
  if (nearZero(a) && nearZero(b)) {
    return true;
  }
  T EPSILON = std::numeric_limits<T>::epsilon()*4.*util::fabs(a);
  return approxEqual(a,b,EPSILON);
}

template <unsigned D, typename F>
void forEachOffset(F fn) any_platform {
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      if constexpr (D == 3) {
        for (int k = -1; k <= 1; ++k) { fn(i, j, k); }
      }
      else { fn(i, j); }
    }
  }
}

}  // namespace util

template <typename T>
std::enable_if_t<std::is_arithmetic<T>::type::value, T> abs(T x) any_platform
{
  return util::fabs(x);
}

// *INDENT-OFF*

}  // namespace olb

#endif
