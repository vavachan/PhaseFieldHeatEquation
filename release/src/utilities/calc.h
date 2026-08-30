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

#ifndef CALC_H
#define CALC_H

#include "omath.h"

namespace olb {

namespace util {

/// Variant of fmod (floating point modulo) that always returns positive values
template<typename T, typename S>
inline auto fmod_pos(T a, S b) {
  const auto res = util::fmod(a, b);
  if (res < 0) {
    return (res + b);
  } else {
    return res;
  }
}

/// Accurate summation of floating point numbers with the Kahan algorithm
// cf. documentation of class below
template<typename T>
void kahanSum(T output[2], T summand)
{
  const T y {summand - output[1]};
  const T t {output[0] + y};
  output[1] = (t - output[0]) - y;
  output[0] = t;
}
template<typename T>
void kahanSum(T& output0, T& output1, T summand)
{
  const T y {summand - output1};
  const T t {output0 + y};
  output1 = (t - output0) - y;
  output0 = t;
}

/** Accurate summation of floating point numbers with the Kahan algorithm
 * Reduces round-off effects which arise if the total sum is significantly
 * larger than the single summands. Works very well for adding many numbers of
 * the same sign. Does not really help if sign changes, e.g. if the summands
 * are randomly distributed around 0.
 * Cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
template<typename T>
class KahanSummator
{
private:
  // 0th component = sum
  // 1st component tracks the cancellation errors
  T result[2] = {0,0};

public:
  void add(T summand) {
    kahanSum<T>(result, summand);
  }

  T getSum() const {
    return result[0];
  }

  void initialize(T value = 0.) {
    result[0] = value;
    result[1] = T();
  }
};

/// @brief Compute cross sum of an integer
/// @tparam K integer type
/// @param k an integer
/// @return cross sum
// cf. https://www.c-plusplus.net/forum/topic/267232/quersumme
template<typename K>
K crossSum(K k) {
  if (k < 0) {
    k = -k;
  }
  K res (0);
  while (k > 0) {
    res += k % 10;
    k /= 10;
  }
  return res;
}

} // namespace util

} // namespace olb

#endif
