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

// generalized_polynomial_chaos.hh
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_HH
#define GENERALIZED_POLYNOMIAL_CHAOS_HH

#include "generalizedPolynomialChaos.h"

#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>

// Include the polynomial basis and quadrature headers
// #include "quadrature.h"

// Namespace aliases
template <typename T>
using LegendreBasis = olb::uq::Polynomials::LegendreBasis<T>;
template <typename T>
using HermiteBasis = olb::uq::Polynomials::HermiteBasis<T>;
// using Quadrature = Quadrature::Quadrature;

namespace olb {

namespace uq {

// Constructor
template <typename T>
GeneralizedPolynomialChaos<T>::GeneralizedPolynomialChaos(std::size_t _order, std::size_t _nq,
                                                          const std::vector<Distribution<T>>& _distributions,
                                                          Quadrature::QuadratureMethod        _quadratureMethod)
    : pointsWeightsMethod(0)
    , No(0)
    , nq(_nq)
    , totalNq(0)
    , order(_order)
    , randomNumberDimension(0)
    , quadratureMethod(_quadratureMethod)
    , distributions(_distributions) // Initialize distributions
{
  // Map distributions to polynomial bases
  polynomialBases.clear();
  for (const auto& dist : distributions) {
    polynomialBases.push_back(createPolynomialBasis(dist));
  }

  // Set randomNumberDimension to the size of distributions
  randomNumberDimension = distributions.size();

  // Calculate multi-indices
  calculateMultiIndices(randomNumberDimension, order, inds);
  No = inds.size();

  // Initialize polynomial coefficients, quadratures, and matrices
  initializeQuadratures();
  initializeMatrices();

  // Evaluate polynomials at quadrature points
  evaluatePhiRan();

  // Compute tensors
  computeTensors();
}

// Initialize quadratures
template <typename T>
void GeneralizedPolynomialChaos<T>::initializeQuadratures()
{
  points.resize(randomNumberDimension);
  weights.resize(randomNumberDimension);

  totalNq = std::pow(nq, randomNumberDimension);

  for (std::size_t i = 0; i < randomNumberDimension; ++i) {
    auto quadrature = polynomialBases[i]->getQuadrature(nq, quadratureMethod);
    points[i]       = quadrature->getPoints();
    weights[i]      = quadrature->getWeights();
  }
}

// Initialize matrices
template <typename T>
void GeneralizedPolynomialChaos<T>::initializeMatrices()
{
  phiRan.resize(totalNq * No, 0.0);
  phiRan_T.resize(totalNq * No, 0.0);
  t2Product.resize(No, 0.0);
  t2Product_inv.resize(No, 0.0);
  t3Product.resize(No * No * No, 0.0);

  // Generate pointsWeightsIndexList
  pointsWeightsIndexList.resize(totalNq, std::vector<std::size_t>(randomNumberDimension));
  for (std::size_t i = 0; i < totalNq; ++i) {
    pointsWeightsIndexList[i] = findIndex(i, randomNumberDimension, nq);
  }

  // Compute weightsMultiplied and pointsTensor
  weightsMultiplied.resize(totalNq, 1.0);
  pointsTensor.resize(totalNq, std::vector<T>(randomNumberDimension));
  for (std::size_t k = 0; k < totalNq; ++k) {
    pointsTensor[k].resize(randomNumberDimension);
    for (std::size_t dim = 0; dim < randomNumberDimension; ++dim) {
      std::size_t idx = pointsWeightsIndexList[k][dim];
      weightsMultiplied[k] *= weights[dim][idx];
      pointsTensor[k][dim] = points[dim][idx];
    }
  }
}

// Initialize polynomial coefficients
template <typename T>
void GeneralizedPolynomialChaos<T>::initializePolynomialCoefficients()
{
  // std::cout << "Initializing polynomial coefficients..." << std::endl;
  coefficients.resize(randomNumberDimension);
  for (std::size_t phi_i = 0; phi_i < randomNumberDimension; ++phi_i) {
    auto basis = std::static_pointer_cast<LegendreBasis>(polynomialBases[phi_i]);

    coefficients[phi_i].resize(No);
    for (std::size_t i = 0; i < No; ++i) {
      coefficients[phi_i][i] = basis->computeCoefficients(i);
    }
  }
}

// Evaluate n_order polynomial at point k
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate(std::size_t n_order, std::size_t k)
{
  T result = 1.0;
  for (std::size_t i = 0; i < randomNumberDimension; ++i) {
    result *= evaluate(inds[n_order][i], k, i);
  }
  return result;
}

// Evaluate n_order polynomial at point k and dimension phi_i
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate(std::size_t n_order, std::size_t k, std::size_t phi_i)
{
  T x = points[phi_i][pointsWeightsIndexList[k][phi_i]];
  return evaluate(n_order, x, phi_i);
}

// Evaluate n_order polynomial at given multi-index
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate(std::size_t n_order, const std::vector<std::size_t>& idx)
{
  T result = 1.0;
  for (std::size_t i = 0; i < randomNumberDimension; ++i) {
    result *= evaluate(inds[n_order][i], points[i][idx[i]], i);
  }
  return result;
}

// Evaluate polynomial basis at given order, point x, and dimension phi_i
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate(std::size_t n_order, T x, std::size_t phi_i)
{
  if (phi_i >= polynomialBases.size()) {
    throw std::out_of_range("Invalid dimension index phi_i.");
  }
  return polynomialBases[phi_i]->evaluatePolynomial(n_order, x);
}

// Evaluate the polynomial at kth point up to order_max
template <typename T>
T GeneralizedPolynomialChaos<T>::evaluate_polynomial(std::size_t order_max, std::size_t k)
{
  T sum = 0.0;
  for (std::size_t i = 0; i <= order_max; ++i) {
    sum += evaluate(i, k);
  }
  return sum;
}

// Evaluate phiRan matrix
template <typename T>
void GeneralizedPolynomialChaos<T>::evaluatePhiRan()
{
  // std::cout << "Evaluating phiRan matrix..." << std::endl;
  for (std::size_t k = 0; k < totalNq; ++k) {
    for (std::size_t i = 0; i < No; ++i) {
      phiRan[k * No + i] = evaluate(i, pointsWeightsIndexList[k]);
      // std::cout << phiRan[k * No + i] << " ";
      phiRan_T[i * totalNq + k] = phiRan[k * No + i];
    }
    // std::cout << std::endl;
  }
}

// Helper functions
template <typename T>
void GeneralizedPolynomialChaos<T>::calculateMultiIndices(std::size_t d, std::size_t n,
                                                          std::vector<std::vector<std::size_t>>& indices)
{

  std::vector<std::size_t> index(d, 0);

  std::function<void(std::size_t, std::size_t, std::size_t)> recursiveFunction = [&](std::size_t pos, std::size_t sum,
                                                                                     std::size_t maxOrder) {
    if (pos == d - 1) {
      index[pos] = maxOrder - sum;
      indices.push_back(index);
      return;
    }

    for (std::size_t i = 0; i <= maxOrder - sum; ++i) {
      index[pos] = i;
      recursiveFunction(pos + 1, sum + i, maxOrder);
    }
  };

  for (std::size_t order = 0; order <= n; ++order) {
    recursiveFunction(0, 0, order);
  }
}

template <typename T>
std::vector<std::size_t> GeneralizedPolynomialChaos<T>::findIndex(std::size_t idx, std::size_t dimension,
                                                                  std::size_t nq)
{
  if (dimension == 1) {
    return {idx};
  }

  // General case for dimension > 1
  std::vector<std::size_t> index(dimension);
  for (std::size_t i = dimension; i-- > 0;) { // Loop from dimension-1 to 0
    index[i] = idx % nq;
    idx /= nq;
  }
  return index;
}

// Compute tensors (t2Product and t3Product)
template <typename T>
void GeneralizedPolynomialChaos<T>::computeTensors()
{
  if (loadSaveT2T3ProductMatrix) {
    const std::string directoryT2Product = "./t2Product/";
    if (!directoryExists(directoryT2Product)) {
      createDirectory(directoryT2Product);
    }
    const std::string directoryT3Product = "./t3Product/";
    if (!directoryExists(directoryT3Product)) {
      createDirectory(directoryT3Product);
    }
    const std::string t2ProductFile = directoryT2Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" +
                                      std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";
    const std::string t3ProductFile = directoryT3Product + "dims_" + std::to_string(randomNumberDimension) + "_order_" +
                                      std::to_string(order) + "_nq_" + std::to_string(nq) + ".bin";

    bool filesExist = fileExists(t2ProductFile) && fileExists(t3ProductFile);
    if (filesExist) {
      readVector1D(t2ProductFile, t2Product);
      for (std::size_t i = 0; i < No; ++i) {
        t2Product_inv[i] = 1.0 / t2Product[i];
      }
      readVector1D(t3ProductFile, t3Product);
    }
    else {
      for (std::size_t i = 0; i < No; ++i) {
        for (std::size_t m = 0; m < totalNq; ++m) {
          t2Product[i] += phiRan[m * No + i] * phiRan[m * No + i] * weightsMultiplied[m];
        }
      }

      for (std::size_t i = 0; i < No; ++i) {
        t2Product_inv[i] = 1.0 / t2Product[i];
      }

      for (std::size_t i = 0; i < No; ++i) {
        for (std::size_t j = 0; j < No; ++j) {
          for (std::size_t k = 0; k < No; ++k) {
            T sum = 0.0;
            for (std::size_t m = 0; m < totalNq; ++m) {
              sum += phiRan[m * No + i] * phiRan[m * No + j] * phiRan[m * No + k] * weightsMultiplied[m];
            }
            t3Product[i * No * No + j * No + k] = sum;
          }
        }
      }
      saveVector1D(t2ProductFile, t2Product);
      saveVector1D(t3ProductFile, t3Product);
    }
  }
  else {
    for (std::size_t i = 0; i < No; ++i) {
      for (std::size_t m = 0; m < totalNq; ++m) {
        t2Product[i] += phiRan[m * No + i] * phiRan[m * No + i] * weightsMultiplied[m];
      }
    }

    for (std::size_t i = 0; i < No; ++i) {
      t2Product_inv[i] = 1.0 / t2Product[i];
    }

    for (std::size_t i = 0; i < No; ++i) {
      for (std::size_t j = 0; j < No; ++j) {
        for (std::size_t k = 0; k < No; ++k) {
          T sum = 0.0;
          for (std::size_t m = 0; m < totalNq; ++m) {
            sum += phiRan[m * No + i] * phiRan[m * No + j] * phiRan[m * No + k] * weightsMultiplied[m];
          }
          t3Product[i * No * No + j * No + k] = sum;
        }
      }
    }
  }
}

// Transformation functions
template <typename T>
void GeneralizedPolynomialChaos<T>::chaosToRandom(const std::vector<T>& chaosCoefficients,
                                                  std::vector<T>&       randomVariables)
{
  randomVariables.resize(totalNq, 0.0);

  for (std::size_t k = 0; k < totalNq; ++k) {
    auto startIt       = phiRan.begin() + k * No;
    randomVariables[k] = std::inner_product(chaosCoefficients.begin(), chaosCoefficients.end(), startIt, 0.0);
  }
}

template <typename T>
void GeneralizedPolynomialChaos<T>::randomToChaos(const std::vector<T>& randomVariables,
                                                  std::vector<T>&       chaosCoefficients)
{
  chaosCoefficients.resize(No, 0.0);
  std::vector<T> weightedRandomVariables(totalNq);

  // Compute weighted random variables
  for (std::size_t k = 0; k < totalNq; ++k) {
    weightedRandomVariables[k] = weightsMultiplied[k] * randomVariables[k];
  }

  // Compute chaos coefficients
  for (std::size_t i = 0; i < No; ++i) {
    auto startIt = phiRan_T.begin() + i * totalNq;
    chaosCoefficients[i] =
        std::inner_product(weightedRandomVariables.begin(), weightedRandomVariables.end(), startIt, 0.0);
    chaosCoefficients[i] *= t2Product_inv[i];
  }
}

// Chaos operations
template <typename T>
void GeneralizedPolynomialChaos<T>::chaosProduct(const std::vector<T>& chaos1, const std::vector<T>& chaos2,
                                                 std::vector<T>& product)
{
  product.resize(No, 0.0);
  std::vector<T> precomputedProductsFlat(No * No);

  for (std::size_t j = 0; j < No; ++j) {
    for (std::size_t k = 0; k < No; ++k) {
      precomputedProductsFlat[j * No + k] = chaos1[j] * chaos2[k];
    }
  }

  for (std::size_t i = 0; i < No; ++i) {
    T sum = 0.0;
    for (std::size_t j = 0; j < No; ++j) {
      for (std::size_t k = 0; k < No; ++k) {
        std::size_t flatIndex = i * No * No + j * No + k;
        sum += precomputedProductsFlat[j * No + k] * t3Product[flatIndex];
      }
    }
    product[i] = sum * t2Product_inv[i];
  }
}

template <typename T>
void GeneralizedPolynomialChaos<T>::chaosSum(const std::vector<T>& chaos1, const std::vector<T>& chaos2,
                                             std::vector<T>& sum)
{
  sum.resize(No);
  for (std::size_t i = 0; i < No; ++i) {
    sum[i] = chaos1[i] + chaos2[i];
  }
}

// Statistical moments
template <typename T>
T GeneralizedPolynomialChaos<T>::mean(const std::vector<T>& chaosCoefficients)
{
  return chaosCoefficients[0];
}

template <typename T>
T GeneralizedPolynomialChaos<T>::std(const std::vector<T>& chaosCoefficients)
{
  T variance = 0.0;
  for (std::size_t i = 1; i < No; ++i) {
    variance += t2Product[i] * chaosCoefficients[i] * chaosCoefficients[i];
  }
  return std::sqrt(variance);
}

template <typename T>
void GeneralizedPolynomialChaos<T>::convert2affinePCE(const Distribution<T>& distribution, std::vector<T>& chaos)
{
  switch (distribution.type) {
  case DistributionType::Uniform: {
    T a1     = 0.5 * (distribution.param1 + distribution.param2);
    T a2     = 0.5 * (distribution.param2 - distribution.param1);
    chaos[0] = a1;
    chaos[1] = a2;
  }
  case DistributionType::Normal: {
    chaos[0] = distribution.param1;
    chaos[1] = distribution.param2;
  }
  // Add cases for other distributions
  default:
    throw std::runtime_error("Unsupported distribution type for GPC.");
  }
}

// Getters
template <typename T>
std::size_t GeneralizedPolynomialChaos<T>::getPolynomialsOrder() const
{
  return No;
}

template <typename T>
std::size_t GeneralizedPolynomialChaos<T>::getQuadraturePointsNumber() const
{
  return totalNq;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getPointsAndWeights(std::vector<std::vector<T>>& points,
                                                        std::vector<std::vector<T>>& weights)
{
  points  = this->points;
  weights = this->weights;
}

template <typename T>
std::vector<std::vector<T>> GeneralizedPolynomialChaos<T>::getStochasticCollocationSample()
{
  std::vector<std::vector<T>> samples(totalNq, std::vector<T>(randomNumberDimension));

  for (std::size_t j = 0; j < randomNumberDimension; ++j) {
    for (std::size_t i = 0; i < totalNq; ++i) {
      samples[i][j] = affine(points[j][pointsWeightsIndexList[i][j]], distributions[j]);
    }
  }

  return samples;
}

template <typename T>
std::vector<T> GeneralizedPolynomialChaos<T>::getWeightsMultiplied() const
{
  return weightsMultiplied;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getTensors(std::vector<T>& t2Product, std::vector<T>& t2Product_inv,
                                               std::vector<T>& t3Product)
{
  t2Product     = this->t2Product;
  t2Product_inv = this->t2Product_inv;
  t3Product     = this->t3Product;
}

template <typename T>
std::shared_ptr<Polynomials::PolynomialBasis<T>>
GeneralizedPolynomialChaos<T>::getPolynomialBasis(std::size_t dimension) const
{
  if (dimension >= randomNumberDimension) {
    throw std::out_of_range("Dimension is out of bounds");
  }

  // Cast the void pointer back to the correct polynomial basis type
  return polynomialBases[dimension];
}

template <typename T>
std::vector<std::vector<std::size_t>> GeneralizedPolynomialChaos<T>::getMultiIndices() const
{
  return inds;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getPhiRan(std::vector<T>& phiRan)
{
  phiRan = this->phiRan;
}

template <typename T>
void GeneralizedPolynomialChaos<T>::getCoefficients(std::vector<std::vector<std::vector<T>>>& polynomialCoeffs)
{
  polynomialCoeffs = this->coefficients;
}

} // namespace uq

} // namespace olb

#endif // GENERALIZED_POLYNOMIAL_CHAOS_HH
