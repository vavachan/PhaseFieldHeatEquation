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

#ifndef UNCERTAINTY_QUANTIFICATION_H
#define UNCERTAINTY_QUANTIFICATION_H

#include <vector>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

#include "distribution.h" // Include the Distribution definitions

// Include the generalized polynomial chaos header
#include "generalizedPolynomialChaos.h"

// Include the sampling methods
#include "monteCarlo.h"
#include "quasiMonteCarlo.h"
#include "latinHypercubeSampling.h"

namespace olb {

namespace uq {

// Enumeration to specify the uncertainty quantification method
enum class UQMethod { GPC, MonteCarlo, QuasiMonteCarlo, LatinHypercubeSampling };

template <typename T>
class UncertaintyQuantification {
public:
  // Constructor
  UncertaintyQuantification(UQMethod uqMethod);

  // Delete copy constructor and copy assignment operator to make the class non-copyable
  UncertaintyQuantification(const UncertaintyQuantification&)            = delete;
  UncertaintyQuantification& operator=(const UncertaintyQuantification&) = delete;

  // Default move constructor and move assignment operator
  UncertaintyQuantification(UncertaintyQuantification&&)            = default;
  UncertaintyQuantification& operator=(UncertaintyQuantification&&) = default;

  // Destructor
  ~UncertaintyQuantification() = default;

  // Initialization functions
  void initializeGPC(std::size_t order, std::size_t nq, Distribution<T> distribution,
                     Quadrature::QuadratureMethod quadratureMethod = Quadrature::QuadratureMethod::WilkinsonShiftQR);
  void initializeGPC(std::size_t order, std::size_t nq, const std::vector<Distribution<T>>& distributions,
                     Quadrature::QuadratureMethod quadratureMethod = Quadrature::QuadratureMethod::WilkinsonShiftQR);

  void initializeMonteCarlo(std::size_t numSamples, const std::vector<Distribution<T>>& distributions,
                            unsigned int seed);
  void initializeMonteCarlo(std::size_t numSamples, Distribution<T> distribution, unsigned int seed);

  void initializeQuasiMonteCarlo(std::size_t numSamples, Distribution<T> distribution,
                                 const std::string& dir_file  = "new-joe-kuo-6.21201",
                                 GeneratorType      generator = GeneratorType::Sobol);
  void initializeQuasiMonteCarlo(std::size_t numSamples, const std::vector<Distribution<T>>& distributions,
                                 const std::string& dir_file  = "new-joe-kuo-6.21201",
                                 GeneratorType      generator = GeneratorType::Sobol);

  void initializeLatinHypercubeSampling(std::size_t numSamples, std::size_t randomNumberDimension);

  // Function to get sampling points
  std::vector<std::vector<T>> getSamplingPoints();
  std::size_t                 getSamplesNumber();

  // Statistical moments
  T mean(const std::vector<T>& samples);
  T std(const std::vector<T>& samples);

  // Other methods common to all UQ methods
  // ...

  // getters
  std::unique_ptr<GeneralizedPolynomialChaos<T>> getOps() { return std::move(ops); };

private:
  UQMethod                     uqMethod;
  std::size_t                  numSamples;            // Number of samples (for MC, QMC, LHS)
  std::size_t                  randomNumberDimension; // Dimensionality of the random input
  std::vector<std::vector<T>>  points;                // Sampling points
  std::vector<Distribution<T>> distributions;

  // GPC-specific members
  std::size_t                                    order; // Order of polynomials (for GPC)
  std::size_t                                    No;    // total order of polynomials system (for GPC)
  std::size_t                                    nq;    // Number of quadrature points per dimension (for GPC)
  std::unique_ptr<GeneralizedPolynomialChaos<T>> ops;
  std::vector<std::vector<T>>                    weights;           // Weights for quadrature (GPC)
  std::vector<T>                                 weightsMultiplied; // Combined weights (GPC)
  std::vector<std::vector<std::size_t>>          multiIndices;      // Multi-indices (GPC)
  std::vector<std::shared_ptr<Polynomials::PolynomialBasis<T>>> polynomialBases;
  Quadrature::QuadratureMethod                                  quadratureMethod; // Quadrature method (GPC)

  // Sampling method instances
  std::unique_ptr<MonteCarlo<T>>             monteCarlo;
  std::unique_ptr<QuasiMonteCarlo<T>>        quasiMonteCarlo;
  std::unique_ptr<LatinHypercubeSampling<T>> lhs;

  // Random number generator
  std::mt19937 rng; // Random number generator

  // Sample generation functions are now in respective classes

  // Example evaluation function can be defined here or elsewhere
};

} // namespace uq

} // namespace olb

#endif // UNCERTAINTY_QUANTIFICATION_H
