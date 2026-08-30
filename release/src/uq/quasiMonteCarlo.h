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

#ifndef QUASI_MONTE_CARLO_H
#define QUASI_MONTE_CARLO_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "distribution.h" // For Distribution and DistributionType

namespace olb {

namespace uq {

enum class GeneratorType { Sobol, Halton };

template <typename T>
class QuasiMonteCarlo {
public:
  // Constructor for multiple uniform distributions per dimension
  QuasiMonteCarlo(std::size_t numSamples, std::size_t randomNumberDimension,
                  const std::vector<Distribution<T>>& distributions,
                  const std::string& dir_file = "new-joe-kuo-6.21201", GeneratorType generator = GeneratorType::Sobol)
      : numSamples(numSamples)
      , randomNumberDimension(randomNumberDimension)
      , distributions(distributions)
      , sobol_dir_file(dir_file)
      , generatorType(generator)
  {
    if (distributions.size() != randomNumberDimension) {
      throw std::invalid_argument("Number of distributions must match the random number dimension.");
    }
    for (const auto& dist : distributions) {
      if (dist.type != DistributionType::Uniform) {
        throw std::invalid_argument("QuasiMonteCarlo now only supports Uniform distributions.");
      }
    }
  }

  // Constructor for a single uniform distribution applied to all dimensions
  QuasiMonteCarlo(std::size_t numSamples, std::size_t randomNumberDimension, Distribution<T> distribution,
                  const std::string& dir_file = "new-joe-kuo-6.21201", GeneratorType generator = GeneratorType::Sobol)
      : numSamples(numSamples)
      , randomNumberDimension(randomNumberDimension)
      , sobol_dir_file(dir_file)
      , generatorType(generator)
  {
    if (distribution.type != DistributionType::Uniform) {
      throw std::invalid_argument("QuasiMonteCarlo now only supports Uniform distributions.");
    }
    distributions = std::vector<Distribution<T>>(randomNumberDimension, distribution);
  }

  // Generate samples
  std::vector<std::vector<T>> generateSamples()
  {
    if (numSamples == 0 || randomNumberDimension == 0) {
      throw std::runtime_error("Number of samples and random number dimension must be greater than zero.");
    }

    std::vector<std::vector<T>> samples(numSamples, std::vector<T>(randomNumberDimension));

    if (generatorType == GeneratorType::Halton) {
      auto haltonPoints = halton_points(numSamples, randomNumberDimension);
      mapToUniformDistributions(samples, haltonPoints);
    }
    else if (generatorType == GeneratorType::Sobol) {
      auto sobolPoints =
          sobol_points(static_cast<unsigned>(numSamples), static_cast<unsigned>(randomNumberDimension), sobol_dir_file);
      mapToUniformDistributions(samples, sobolPoints);
    }
    else {
      throw std::runtime_error("Unsupported generator type.");
    }

    return samples;
  }

  // Get the number of samples
  std::size_t getSamplesNumber() const { return numSamples; }

private:
  std::size_t                  numSamples;
  std::size_t                  randomNumberDimension;
  std::vector<Distribution<T>> distributions;
  std::string                  sobol_dir_file;
  GeneratorType                generatorType;

  void mapToUniformDistributions(std::vector<std::vector<T>>& samples, const std::vector<std::vector<T>>& points)
  {
    samples.resize(numSamples, std::vector<T>(randomNumberDimension));
    for (std::size_t i = 0; i < numSamples; ++i) {
      for (std::size_t j = 0; j < randomNumberDimension; ++j) {
        const Distribution<T>& dist = distributions[j];
        samples[i][j]               = dist.param1 + (dist.param2 - dist.param1) * points[i][j];
      }
    }
  }

  // Helper function to generate Sobol sequence points
  std::vector<std::vector<T>> sobol_points(unsigned N, unsigned D, const std::string& dir_file)
  {
    std::ifstream infile(dir_file);
    if (!infile) {
      throw std::runtime_error("Input file containing direction numbers cannot be found!");
    }

    // Read header line
    std::string buffer;
    std::getline(infile, buffer);

    // L = max number of bits needed
    unsigned L = static_cast<unsigned>(std::ceil(std::log(static_cast<T>(N)) / std::log(2.0)));

    // C[i] = index from the right of the first zero bit of i
    std::vector<unsigned> C(N);
    C[0] = 1;
    for (unsigned i = 1; i < N; i++) {
      C[i]           = 1;
      unsigned value = i;
      while (value & 1) {
        value >>= 1;
        C[i]++;
      }
    }

    // POINTS[i][j] = the jth component of the ith point
    std::vector<std::vector<T>> POINTS(N, std::vector<T>(D));

    // Initialize direction numbers V
    std::vector<std::vector<unsigned>> V(D, std::vector<unsigned>(L + 1));

    // Read direction numbers for each dimension
    for (unsigned j = 0; j < D; ++j) {
      if (j == 0) {
        // First dimension
        for (unsigned i = 1; i <= L; ++i) {
          V[j][i] = 1U << (32 - i); // All m's = 1
        }
      }
      else {
        // Read in parameters from file
        unsigned d, s, a;
        infile >> d >> s >> a;
        std::vector<unsigned> m(s + 1);
        for (unsigned i = 1; i <= s; ++i) {
          infile >> m[i];
        }

        // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
        if (L <= s) {
          for (unsigned i = 1; i <= L; ++i) {
            V[j][i] = m[i] << (32 - i);
          }
        }
        else {
          for (unsigned i = 1; i <= s; ++i) {
            V[j][i] = m[i] << (32 - i);
          }
          for (unsigned i = s + 1; i <= L; ++i) {
            V[j][i] = V[j][i - s] ^ (V[j][i - s] >> s);
            for (unsigned k = 1; k <= s - 1; ++k) {
              V[j][i] ^= (((a >> (s - 1 - k)) & 1U) * V[j][i - k]);
            }
          }
        }
      }
    }

    // Evaluate X[0] to X[N-1], scaled by pow(2,32)
    std::vector<unsigned> X(D, 0);
    for (unsigned i = 0; i < N; ++i) {
      if (i > 0) {
        for (unsigned j = 0; j < D; ++j) {
          X[j] ^= V[j][C[i - 1]];
        }
      }
      for (unsigned j = 0; j < D; ++j) {
        POINTS[i][j] = static_cast<T>(X[j]) / std::pow(2.0, 32);
      }
    }

    return POINTS;
  }

  std::vector<std::vector<T>> halton_points(std::size_t N, std::size_t D)
  {
    std::vector<std::vector<T>> points(N, std::vector<T>(D));
    for (std::size_t d = 0; d < D; ++d) {
      unsigned base = getPrime(d);
      for (std::size_t i = 0; i < N; ++i) {
        points[i][d] = halton_single_point(i + 1, base);
      }
    }
    return points;
  }

  T halton_single_point(std::size_t index, unsigned base)
  {
    T result = 0.0;
    T f      = 1.0 / base;
    while (index > 0) {
      result += f * (index % base);
      index /= base;
      f /= base;
    }
    return result;
  }

  unsigned getPrime(std::size_t index)
  {
    static const unsigned primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    if (index < sizeof(primes) / sizeof(primes[0])) {
      return primes[index];
    }
    else {
      throw std::invalid_argument("Prime index out of range for Halton sequence.");
    }
  }
};

} // namespace uq

} // namespace olb

#endif // QUASI_MONTE_CARLO_H
