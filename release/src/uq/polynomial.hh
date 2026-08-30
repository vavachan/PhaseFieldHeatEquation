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

#ifndef POLYNOMIAL_HH
#define POLYNOMIAL_HH
#include "polynomial.h"
#include "quadrature.h"

namespace olb {

namespace uq {

namespace Polynomials {

// Helper function to compute Legendre coefficients recursively
template <typename T>
static std::vector<T> computeLegendreCoefficients(std::size_t n)
{
  if (n == 0) {
    return {1.0}; // P_0(x) = 1
  }
  else if (n == 1) {
    return {0.0, 1.0}; // P_1(x) = x
  }
  else {
    // Recurrence relation: P_n(x) = ((2n - 1)/n) * x * P_{n-1}(x) - ((n - 1)/n) * P_{n-2}(x)
    std::vector<T> Pn_minus1 = computeLegendreCoefficients<T>(n - 1);
    std::vector<T> Pn_minus2 = computeLegendreCoefficients<T>(n - 2);
    std::vector<T> Pn(n + 1, 0.0);

    T a = (2.0 * n - 1.0) / n;
    T b = (n - 1.0) / n;

    // Multiply P_{n-1} by x (shift coefficients)
    std::vector<T> xPn_minus1(n + 1, 0.0);
    for (std::size_t i = 0; i < Pn_minus1.size(); ++i) {
      xPn_minus1[i + 1] = Pn_minus1[i];
    }

    // Compute Pn = a * x * P_{n-1} - b * P_{n-2}
    for (std::size_t i = 0; i <= n; ++i) {
      T term1 = a * xPn_minus1[i];
      T term2 = (i < Pn_minus2.size()) ? b * Pn_minus2[i] : 0.0;
      Pn[i]   = term1 - term2;
    }

    return Pn;
  }
}

// Compute coefficients for the Legendre polynomial of order n
template <typename T>
std::vector<T> LegendreBasis<T>::computeCoefficients(std::size_t n) const
{
  return computeLegendreCoefficients<T>(n);
}

// Construct the Jacobi matrix of size n for Legendre polynomials
template <typename T>
std::vector<std::vector<T>> LegendreBasis<T>::constructJacobiMatrix(std::size_t n) const
{
  if (n <= 1) {
    throw std::invalid_argument("The order (n) must be greater than 1.");
  }

  std::vector<std::vector<T>> J(n, std::vector<T>(n, 0.0));

  // Fill the Jacobi matrix according to the recurrence relation for Legendre polynomials
  for (std::size_t i = 1; i < n; ++i) {
    T a         = i / std::sqrt(4.0 * i * i - 1.0);
    J[i][i - 1] = a;
    J[i - 1][i] = a;
  }

  return J;
}

template <typename T>
std::shared_ptr<Quadrature::QuadratureBase<T>>
LegendreBasis<T>::getQuadrature(std::size_t nq, Quadrature::QuadratureMethod method) const
{
  return std::make_shared<Quadrature::Quadrature<T, LegendreBasis<T>>>(nq, method);
}

// Helper function to compute probabilist's Hermite coefficients recursively
template <typename T>
static std::vector<T> computeHermiteCoefficients(std::size_t n)
{
  // Base cases
  if (n == 0) {
    return {1.0}; // He_0(x) = 1
  }
  else if (n == 1) {
    return {0.0, 1.0}; // He_1(x) = x
  }
  else {
    // Initialize the coefficients for He_0(x) and He_1(x)
    std::vector<T> Hn_minus_two = {1.0};      // He_0(x)
    std::vector<T> Hn_minus_one = {0.0, 1.0}; // He_1(x)
    std::vector<T> Hn;                        // To store coefficients of He_n(x)

    // Iteratively compute He_n(x) using the recurrence relation
    for (std::size_t k = 2; k <= n; ++k) {
      // Initialize Hn with zeros, size is k + 1
      Hn.assign(k + 1, 0.0);

      // Compute He_n = x * He_{n-1}(x) - (k - 1) * He_{n-2}(x)

      // Term 1: x * He_{n-1}(x)
      for (std::size_t i = 0; i < Hn_minus_one.size(); ++i) {
        Hn[i + 1] += Hn_minus_one[i];
      }

      // Term 2: - (k - 1) * He_{n-2}(x)
      for (std::size_t i = 0; i < Hn_minus_two.size(); ++i) {
        Hn[i] -= (k - 1) * Hn_minus_two[i];
      }

      // Update for next iteration
      Hn_minus_two = Hn_minus_one;
      Hn_minus_one = Hn;
    }
    return Hn_minus_one; // Coefficients of He_n(x)
  }
}

// Compute coefficients for the Hermite polynomial of order n
template <typename T>
std::vector<T> HermiteBasis<T>::computeCoefficients(std::size_t n) const
{
  return computeHermiteCoefficients<T>(n);
}

// Construct the Jacobi matrix of size n for Hermite polynomials
template <typename T>
std::vector<std::vector<T>> HermiteBasis<T>::constructJacobiMatrix(std::size_t n) const
{
  if (n <= 1) {
    throw std::invalid_argument("The order (n) must be greater than 1.");
  }

  std::vector<std::vector<T>> J(n, std::vector<T>(n, 0.0));

  // Fill the Jacobi matrix according to the recurrence relation for Hermite polynomials
  for (std::size_t i = 1; i < n; ++i) {
    // The recurrence relation coefficients for Hermite polynomials
    if (i > 0) {
      J[i][i - 1] = std::sqrt(i);
      J[i - 1][i] = J[i][i - 1];
    }
  }

  return J;
}

template <typename T>
std::shared_ptr<Quadrature::QuadratureBase<T>> HermiteBasis<T>::getQuadrature(std::size_t                  nq,
                                                                              Quadrature::QuadratureMethod method) const
{
  return std::make_shared<Quadrature::Quadrature<T, HermiteBasis<T>>>(nq, method);
}

} // namespace Polynomials

} // namespace uq

} // namespace olb

#endif // POLYNOMIAL_HH