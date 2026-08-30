/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
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

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>

namespace olb {

namespace uq {

// Enumeration to specify the distribution type
enum class DistributionType {
  Uniform,
  Normal,
  // MultivariateNormal
};

// Struct to hold distribution information
template <typename T>
struct Distribution {
  DistributionType type;

  // Parameters for univariate distributions
  T param1 = 0.0; // For Uniform: lower bound, for Normal: mean
  T param2 = 1.0; // For Uniform: upper bound, for Normal: standard deviation

  // Constructor for univariate distributions
  Distribution(DistributionType type, T param1 = 0.0, T param2 = 1.0)
      : type(type)
      , param1(param1)
      , param2(param2)
  {}
};

// Factory functions for cleaner syntax
template <typename T>
inline Distribution<T> uniform(T min, T max)
{
  return Distribution(DistributionType::Uniform, min, max);
}

template <typename T>
inline Distribution<T> normal(T mean, T stddev)
{
  return Distribution(DistributionType::Normal, mean, stddev);
}

// Joint distribution as std::vector<Distribution>
template <typename T>
inline std::vector<Distribution<T>> joint(const std::vector<Distribution<T>>& dists)
{
  return dists;
}

// Affine transformation function
template <typename T>
T affine(T x, const Distribution<T>& dist)
{
  switch (dist.type) {
  case DistributionType::Uniform:
    return 0.5 * (dist.param2 - dist.param1) * x + 0.5 * (dist.param1 + dist.param2);
  case DistributionType::Normal:
    return dist.param1 + dist.param2 * x;
  default:
    throw std::runtime_error("Unsupported distribution type for affine transformation.");
  }
}

} // namespace uq

} // namespace olb

#endif // DISTRIBUTION_H
