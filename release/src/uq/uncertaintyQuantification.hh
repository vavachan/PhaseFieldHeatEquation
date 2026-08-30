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

#ifndef UNCERTAINTY_QUANTIFICATION_HH
#define UNCERTAINTY_QUANTIFICATION_HH

#include "uncertaintyQuantification.h"
#include <numeric>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

namespace olb {

namespace uq {

// Constructor
template <typename T>
UncertaintyQuantification<T>::UncertaintyQuantification(UQMethod uqMethod)
    : uqMethod(uqMethod)
    , rng(std::random_device {}())
{
  // Initialize variables to default values
  order                 = 0;
  nq                    = 0;
  numSamples            = 0;
  randomNumberDimension = 0;
  quadratureMethod      = Quadrature::QuadratureMethod::HouseholderQR;
  ops                   = nullptr; // Initialize the unique_ptr to nullptr
  monteCarlo            = nullptr;
  quasiMonteCarlo       = nullptr;
  lhs                   = nullptr;
}

// Initialization function for GPC
template <typename T>
void UncertaintyQuantification<T>::initializeGPC(std::size_t order, std::size_t nq, Distribution<T> distribution,
                                                 Quadrature::QuadratureMethod quadratureMethod)
{
  this->order                 = order;
  this->nq                    = nq;
  this->quadratureMethod      = quadratureMethod;
  this->randomNumberDimension = 1;
  this->distributions         = {distribution};

  // Initialize the GeneralizedPolynomialChaos object
  ops = std::make_unique<GeneralizedPolynomialChaos<T>>(order, nq, distributions, quadratureMethod);
  // Get quadrature points and weights from ops
  ops->getPointsAndWeights(points, weights);

  No = ops->getPolynomialsOrder();

  // Compute the combined weights
  numSamples        = ops->getQuadraturePointsNumber();
  weightsMultiplied = ops->getWeightsMultiplied();
  // Get multi-indices from ops
  multiIndices = ops->getMultiIndices();
}

template <typename T>
void UncertaintyQuantification<T>::initializeGPC(std::size_t order, std::size_t nq,
                                                 const std::vector<Distribution<T>>& distributions,
                                                 Quadrature::QuadratureMethod        quadratureMethod)
{
  this->order                 = order;
  this->nq                    = nq;
  this->quadratureMethod      = quadratureMethod;
  this->randomNumberDimension = distributions.size();
  this->distributions         = distributions;

  // Initialize the GeneralizedPolynomialChaos object
  ops = std::make_unique<GeneralizedPolynomialChaos<T>>(order, nq, distributions, quadratureMethod);
  // Get quadrature points and weights from ops
  ops->getPointsAndWeights(points, weights);

  No = ops->getPolynomialsOrder();

  // Compute the combined weights
  numSamples        = ops->getQuadraturePointsNumber();
  weightsMultiplied = ops->getWeightsMultiplied();
  // Get multi-indices from ops
  multiIndices = ops->getMultiIndices();
}

// Initialization function for Monte Carlo with vector of distributions
template <typename T>
void UncertaintyQuantification<T>::initializeMonteCarlo(std::size_t                         numSamples,
                                                        const std::vector<Distribution<T>>& distributions,
                                                        unsigned int                        seed)
{
  this->numSamples            = numSamples;
  this->randomNumberDimension = distributions.size();
  this->distributions         = distributions;

  // Create MonteCarlo instance
  monteCarlo = std::make_unique<MonteCarlo<T>>(numSamples, distributions, seed);
}

// Overload for a single distribution applied to all dimensions
template <typename T>
void UncertaintyQuantification<T>::initializeMonteCarlo(std::size_t numSamples, Distribution<T> distribution,
                                                        unsigned int seed)
{
  this->numSamples            = numSamples;
  this->randomNumberDimension = 1;
  this->distributions         = {distribution};

  // Create MonteCarlo instance
  monteCarlo = std::make_unique<MonteCarlo<T>>(numSamples, distribution, seed);
}

// Initialization function for Quasi-Monte Carlo
template <typename T>
void UncertaintyQuantification<T>::initializeQuasiMonteCarlo(std::size_t numSamples, Distribution<T> distribution,
                                                             const std::string& dir_file, GeneratorType generator)
{
  this->numSamples            = numSamples;
  this->randomNumberDimension = 1;
  this->distributions         = {distribution};

  // Create QuasiMonteCarlo instance
  quasiMonteCarlo =
      std::make_unique<QuasiMonteCarlo<T>>(numSamples, randomNumberDimension, distribution, dir_file, generator);
}

template <typename T>
void UncertaintyQuantification<T>::initializeQuasiMonteCarlo(std::size_t                         numSamples,
                                                             const std::vector<Distribution<T>>& distributions,
                                                             const std::string& dir_file, GeneratorType generator)
{
  this->numSamples            = numSamples;
  this->randomNumberDimension = distributions.size();
  this->distributions         = distributions;

  // Create QuasiMonteCarlo instance
  quasiMonteCarlo =
      std::make_unique<QuasiMonteCarlo<T>>(numSamples, randomNumberDimension, distributions, dir_file, generator);
}

// Initialization function for Latin Hypercube Sampling
template <typename T>
void UncertaintyQuantification<T>::initializeLatinHypercubeSampling(std::size_t numSamples,
                                                                    std::size_t randomNumberDimension)
{
  this->numSamples            = numSamples;
  this->randomNumberDimension = randomNumberDimension;

  // Create LatinHypercubeSampling instance
  lhs = std::make_unique<LatinHypercubeSampling<T>>(static_cast<int>(numSamples),
                                                    static_cast<int>(randomNumberDimension));
}

// Function to get sampling points
template <typename T>
std::vector<std::vector<T>> UncertaintyQuantification<T>::getSamplingPoints()
{
  switch (uqMethod) {
  case UQMethod::GPC:
    if (!ops) {
      throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
    }
    return ops->getStochasticCollocationSample();
    break;
  case UQMethod::MonteCarlo:
    if (!monteCarlo) {
      throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
    }
    return monteCarlo->generateSamples();
    break;
  case UQMethod::QuasiMonteCarlo:
    if (!quasiMonteCarlo) {
      throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
    }
    return quasiMonteCarlo->generateSamples();
    break;
  case UQMethod::LatinHypercubeSampling:
    if (!lhs) {
      throw std::runtime_error(
          "Latin Hypercube Sampling has not been initialized. Call initializeLatinHypercubeSampling() first.");
    }
    return lhs->generateSamples();
    break;
  default:
    throw std::runtime_error("Invalid UQ method.");
  }
}

template <typename T>
std::size_t UncertaintyQuantification<T>::getSamplesNumber()
{
  switch (uqMethod) {
  case UQMethod::GPC:
    if (!ops) {
      throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
    }
    return ops->getQuadraturePointsNumber();
    break;
  case UQMethod::MonteCarlo:
    if (!monteCarlo) {
      throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
    }
    return monteCarlo->getSamplesNumber();
    break;
  case UQMethod::QuasiMonteCarlo:
    if (!quasiMonteCarlo) {
      throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
    }
    return quasiMonteCarlo->getSamplesNumber();
    break;
  case UQMethod::LatinHypercubeSampling:
    if (!lhs) {
      throw std::runtime_error(
          "Latin Hypercube Sampling has not been initialized. Call initializeLatinHypercubeSampling() first.");
    }
    return lhs->getSamplesNumber();
    break;
  default:
    throw std::runtime_error("Invalid UQ method.");
  }
}

// Example function evaluation (to be customized)
template <typename T>
T UncertaintyQuantification<T>::mean(const std::vector<T>& input)
{
  T mean = 0.0;
  switch (uqMethod) {
  case UQMethod::GPC: {
    if (!ops) {
      throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
    }
    std::vector<T> chaos(No, 0.0);
    ops->randomToChaos(input, chaos);
    mean = ops->mean(chaos);
    break;
  }
  case UQMethod::MonteCarlo: {
    if (!monteCarlo) {
      throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
    }
    T sum = std::accumulate(input.begin(), input.end(), 0.0);
    mean  = sum / input.size();
    break;
  }
  case UQMethod::QuasiMonteCarlo: {
    if (!quasiMonteCarlo) {
      throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
    }
    T sum = std::accumulate(input.begin(), input.end(), 0.0);
    mean  = sum / input.size();
    break;
  }
  default:
    throw std::runtime_error("Invalid UQ method.");
  }
  return mean;
}

template <typename T>
T UncertaintyQuantification<T>::std(const std::vector<T>& input)
{
  T std = 0.0;
  switch (uqMethod) {
  case UQMethod::GPC: {
    if (!ops) {
      throw std::runtime_error("GPC has not been initialized. Call initializeGPC() first.");
    }
    std::vector<T> chaos(No, 0.0);
    ops->randomToChaos(input, chaos);
    std = ops->std(chaos);
    break;
  }
  case UQMethod::MonteCarlo: {
    if (!monteCarlo) {
      throw std::runtime_error("Monte Carlo has not been initialized. Call initializeMonteCarlo() first.");
    }
    T avg   = mean(input);
    T sumSq = std::accumulate(input.begin(), input.end(), 0.0, [avg](T acc, T val) {
      return acc + (val - avg) * (val - avg);
    });
    std     = std::sqrt(sumSq / (input.size() - 1));
    break;
  }
  case UQMethod::QuasiMonteCarlo: {
    if (!quasiMonteCarlo) {
      throw std::runtime_error("Quasi-Monte Carlo has not been initialized. Call initializeQuasiMonteCarlo() first.");
    }
    T avg   = mean(input);
    T sumSq = std::accumulate(input.begin(), input.end(), 0.0, [avg](T acc, T val) {
      return acc + (val - avg) * (val - avg);
    });
    std     = std::sqrt(sumSq / (input.size() - 1));
    break;
  }
  default:
    throw std::runtime_error("Invalid UQ method.");
  }
  return std;
}

} // namespace uq

} // namespace olb

#endif // UNCERTAINTY_QUANTIFICATION_HH
