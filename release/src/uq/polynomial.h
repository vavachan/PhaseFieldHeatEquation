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

// polynomial.h
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "filesIO.h"
#include "quadratureBase.h"

namespace olb {

namespace uq {

namespace Polynomials {

// Abstract base class for all polynomial basis types
template <typename T>
class PolynomialBasis {
public:
  virtual ~PolynomialBasis() = default;

  // Compute polynomial coefficients of order n
  virtual std::vector<T> computeCoefficients(std::size_t n) const = 0;

  // Construct the Jacobi matrix of size n
  virtual std::vector<std::vector<T>> constructJacobiMatrix(std::size_t n) const = 0;

  // Evaluate the polynomial of order n at point x
  T evaluatePolynomial(std::size_t n, T x) const;

  // Compute the derivative of the polynomial of order n at point x
  T derivativePolynomial(std::size_t n, T x) const;

  // Dynamically create and return a quadrature object
  virtual std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(std::size_t                  nq,
                                                                       Quadrature::QuadratureMethod method) const = 0;
};

// Implement evaluatePolynomial using Horner's method for efficiency
template <typename T>
inline T PolynomialBasis<T>::evaluatePolynomial(std::size_t n, T x) const
{
  std::vector<T> coeffs = computeCoefficients(n);

  // Evaluate polynomial using Horner's method
  T result = coeffs.back();
  for (std::size_t i = coeffs.size() - 1; i-- > 0;) {
    result = result * x + coeffs[i];
  }
  return result;
}

// Implement derivativePolynomial using Horner's method
template <typename T>
inline T PolynomialBasis<T>::derivativePolynomial(std::size_t n, T x) const
{
  std::vector<T> coeffs = computeCoefficients(n);
  T              result = 0.0;
  for (std::size_t i = coeffs.size() - 1; i > 0; --i) {
    result = result * x + i * coeffs[i];
  }
  return result;
}

// LegendreBasis class that inherits from PolynomialBasis
template <typename T>
class LegendreBasis : public PolynomialBasis<T> {
public:
  // Implement all pure virtual functions
  std::vector<T>                                 computeCoefficients(std::size_t n) const override;
  std::vector<std::vector<T>>                    constructJacobiMatrix(std::size_t n) const override;
  std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(std::size_t                  nq,
                                                               Quadrature::QuadratureMethod method) const override;
};

// HermiteBasis class that inherits from PolynomialBasis
template <typename T>
class HermiteBasis : public PolynomialBasis<T> {
public:
  // Implement all pure virtual functions
  std::vector<T>                                 computeCoefficients(std::size_t n) const override;
  std::vector<std::vector<T>>                    constructJacobiMatrix(std::size_t n) const override;
  std::shared_ptr<Quadrature::QuadratureBase<T>> getQuadrature(std::size_t                  nq,
                                                               Quadrature::QuadratureMethod method) const override;
};

} // namespace Polynomials

} // namespace uq

} // namespace olb

#endif // POLYNOMIAL_H
