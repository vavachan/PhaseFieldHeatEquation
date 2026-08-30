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

// generalized_polynomial_chaos.h
#ifndef GENERALIZED_POLYNOMIAL_CHAOS_H
#define GENERALIZED_POLYNOMIAL_CHAOS_H

#include <vector>
#include <memory>
#include <cmath>
#include <string>

#include "filesIO.h"

#include "distribution.h"

// Include the polynomial basis and quadrature headers
#include "polynomial.h"

// #include "quadrature.h"

namespace olb {

namespace uq {

template <typename T>
class GeneralizedPolynomialChaos {
public:
  // Constructor
  GeneralizedPolynomialChaos(std::size_t order, std::size_t nq, const std::vector<Distribution<T>>& distributions,
                             Quadrature::QuadratureMethod quadratureMethod);

  // Evaluation functions
  T evaluate(std::size_t n_order, std::size_t k);
  T evaluate(std::size_t n_order, std::size_t k, std::size_t phi_i);
  T evaluate(std::size_t n_order, const std::vector<std::size_t>& idx);
  T evaluate(std::size_t n_order, T x, std::size_t phi_i);
  T evaluate_polynomial(std::size_t order_max, std::size_t k);

  // Compute phiRan matrix
  void evaluatePhiRan();

  // Compute tensors
  void computeTensors();

  // Transformation functions
  void chaosToRandom(const std::vector<T>& chaosCoefficients, std::vector<T>& randomVariables);
  void randomToChaos(const std::vector<T>& randomVariables, std::vector<T>& chaosCoefficients);

  // Chaos operations
  void chaosProduct(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& product);
  void chaosSum(const std::vector<T>& chaos1, const std::vector<T>& chaos2, std::vector<T>& sum);

  // Statistical moments
  T mean(const std::vector<T>& chaosCoefficients);
  T std(const std::vector<T>& chaosCoefficients);

  void convert2affinePCE(const Distribution<T>& distribution, std::vector<T>& chaos);

  // Getters
  std::size_t getPolynomialsOrder() const;
  std::size_t getQuadraturePointsNumber() const;
  void        getPointsAndWeights(std::vector<std::vector<T>>& points, std::vector<std::vector<T>>& weights);
  std::vector<std::vector<T>> getStochasticCollocationSample();
  void           getTensors(std::vector<T>& t2Product, std::vector<T>& t2Product_inv, std::vector<T>& t3Product);
  std::vector<T> getWeightsMultiplied() const;
  // Template function to get the polynomial basis at a specific dimension (i)
  std::shared_ptr<Polynomials::PolynomialBasis<T>> getPolynomialBasis(std::size_t i) const;
  std::vector<std::vector<std::size_t>>            getMultiIndices() const;

  void getPhiRan(std::vector<T>& phiRan);
  void getCoefficients(std::vector<std::vector<std::vector<T>>>& polynomialCoeffs);

private:
  std::size_t                              pointsWeightsMethod;
  std::size_t                              No;      // Number of polynomials
  std::size_t                              nq;      // Number of quadrature points per dimension
  std::size_t                              totalNq; // Total number of quadrature points
  std::size_t                              order;
  std::size_t                              randomNumberDimension;
  std::vector<std::vector<std::size_t>>    inds;              // Multi-indices
  std::vector<std::vector<T>>              points;            // Points for each dimension
  std::vector<std::vector<T>>              weights;           // Weights for each dimension
  std::vector<std::vector<T>>              pointsTensor;      // Tensor product of points
  std::vector<T>                           weightsMultiplied; // Combined weights
  std::vector<std::vector<std::size_t>>    pointsWeightsIndexList;
  std::vector<std::vector<std::vector<T>>> coefficients; // Coefficients of polynomials

  Quadrature::QuadratureMethod quadratureMethod;

  std::vector<T> phiRan;   // Evaluated polynomials at quadrature points
  std::vector<T> phiRan_T; // Transpose of phiRan
  std::vector<T> t2Product;
  std::vector<T> t2Product_inv;
  std::vector<T> t3Product;

  bool loadSaveT2T3ProductMatrix = false; // Flag to load/save T2 and T3 product matrices

  // Distributions for each dimension
  std::vector<Distribution<T>> distributions;

  // Polynomial bases for each dimension
  std::vector<std::shared_ptr<Polynomials::PolynomialBasis<T>>> polynomialBases;

  // Initialization functions
  void initializeQuadratures();
  void initializeMatrices();

  void initializePolynomialCoefficients();

  // Helper functions
  std::vector<std::size_t> findIndex(std::size_t idx, std::size_t dimension, std::size_t nq);
  void calculateMultiIndices(std::size_t d, std::size_t n, std::vector<std::vector<std::size_t>>& indices);

  std::shared_ptr<Polynomials::PolynomialBasis<T>> createPolynomialBasis(const Distribution<T>& dist)
  {
    switch (dist.type) {
    case DistributionType::Uniform:
      return std::make_shared<Polynomials::LegendreBasis<T>>();
    case DistributionType::Normal:
      return std::make_shared<Polynomials::HermiteBasis<T>>();
    // Add cases for other distributions
    default:
      throw std::runtime_error("Unsupported distribution type for GPC.");
    }
  }
};

} // namespace uq

} // namespace olb

#endif // GENERALIZED_POLYNOMIAL_CHAOS_H
