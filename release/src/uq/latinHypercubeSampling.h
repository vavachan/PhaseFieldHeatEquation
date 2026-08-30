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

#ifndef LATIN_HYPERCUBE_SAMPLING_H
#define LATIN_HYPERCUBE_SAMPLING_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

namespace olb {

namespace uq {

template <typename T>
class LatinHypercubeSampling {
public:
  // Constructor: numSamples is the number of samples, numDimensions is the dimensionality of the space
  LatinHypercubeSampling(int numSamples, int numDimensions)
      : numSamples(numSamples)
      , numDimensions(numDimensions)
  {
    // Initialize the random generator
    rng.seed(std::chrono::system_clock::now().time_since_epoch().count());
  }

  // Generate Latin Hypercube samples
  std::vector<std::vector<T>> generateSamples()
  {
    // Create a matrix to store samples
    std::vector<std::vector<T>> samples(numSamples, std::vector<T>(numDimensions));

    // Step 1: Divide the range [0, 1] into numSamples intervals for each dimension
    for (int i = 0; i < numDimensions; ++i) {
      // Generate evenly spaced intervals
      std::vector<T> intervals(numSamples);
      for (int j = 0; j < numSamples; ++j) {
        intervals[j] = (j + randomUniform(0.0, 1.0)) / numSamples;
      }

      // Step 2: Randomly permute the intervals to ensure a random distribution in each dimension
      std::shuffle(intervals.begin(), intervals.end(), rng);

      // Step 3: Assign the permuted intervals to the samples
      for (int j = 0; j < numSamples; ++j) {
        samples[j][i] = intervals[j];
      }
    }

    return samples;
  }

  // Get the number of samples
  int getSamplesNumber() const { return numSamples; }

private:
  int          numSamples;
  int          numDimensions;
  std::mt19937 rng; // Random number generator

  // Generate a random number in the range [low, high)
  T randomUniform(T low, T high)
  {
    std::uniform_real_distribution<T> dist(low, high);
    return dist(rng);
  }
};

} // namespace uq

} // namespace olb

#endif // LATIN_HYPERCUBE_SAMPLING_H
