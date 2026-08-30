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

// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "filesIO.h"
#include "quadratureBase.h"
#include "polynomial.h"
#include "matrixOperation.h"

namespace olb {

namespace uq {

namespace Quadrature {

template <class...>
struct always_false : std::false_type {};

// Add this function to matrix_operation.h
template <typename T>
T Pythag(T a, T b)
{
  T absa = std::fabs(a);
  T absb = std::fabs(b);
  return (absa > absb ? absa * std::sqrt(1.0 + (absb / absa) * (absb / absa))
                      : (absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + (absa / absb) * (absa / absb))));
}

template <typename T, typename PolynomialBasis>
class Quadrature : public QuadratureBase<T> {
public:
  Quadrature(std::size_t nq, QuadratureMethod method = QuadratureMethod::WilkinsonShiftQR)
      : nq(nq)
      , basis(std::make_shared<PolynomialBasis>())
      , method(method)
  {
    performQRDecomposition();
  }

  const std::vector<T>& getPoints() const override { return points; }

  const std::vector<T>& getWeights() const override { return weights; }

private:
  std::size_t                      nq;
  std::shared_ptr<PolynomialBasis> basis;
  QuadratureMethod                 method;
  std::vector<T>                   points;
  std::vector<T>                   weights;

  void performQRDecomposition()
  {
    // Step 1: Construct the Jacobi matrix using the basis
    auto J = basis->constructJacobiMatrix(nq);

    // Step 2: Perform the appropriate QR decomposition
    if (method == QuadratureMethod::HouseholderQR) {
      // std::cout << "Using HouseholderQR method." << std::endl;
      points = HouseholderQRDecomposition(J); // Perform Householder QR and get eigenvalues
    }
    else if (method == QuadratureMethod::WilkinsonShiftQR) {
      // std::cerr << "Using Wilkinson's Shift QR method." << std::endl;
      points = WilkinsonShiftQRDecomposition(J); // Perform Wilkinson Shift QR
    }
    else {
      // std::cerr << "Warning: Unsupported method. Defaulting to Wilkinson's Shift QR." << std::endl;
      points = WilkinsonShiftQRDecomposition(J); // Default to Wilkinson Shift QR
    }

    // Step 3: Compute quadrature weights
    computeWeights();
  }

  // Wilkinson Shift QR decomposition to compute the eigenvalues of a tridiagonal matrix
  std::vector<T> WilkinsonShiftQRDecomposition(const std::vector<std::vector<T>>& J)
  {
    std::size_t    n = J.size();
    std::vector<T> d(n, 0.0); // Diagonal elements (eigenvalues)
    std::vector<T> e(n, 0.0); // Off-diagonal elements

    // Initialize d and e from the input matrix J
    for (std::size_t i = 0; i < n; ++i) {
      d[i] = J[i][i];
      if (i > 0)
        e[i] = J[i][i - 1];
    }

    const T eps = std::numeric_limits<T>::epsilon();
    T       g, r, p, s, c, f, b;
    int     iter, l, m, i;

    // Reduce subdiagonal elements
    for (i = 1; i < static_cast<int>(n); ++i)
      e[i - 1] = e[i];
    e[n - 1] = 0.0;

    for (l = 0; l < static_cast<int>(n); ++l) {
      iter = 0;
      do {
        // Find small subdiagonal element to isolate a diagonal block
        for (m = l; m < static_cast<int>(n) - 1; ++m) {
          T dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
          if (std::fabs(e[m]) <= eps * dd)
            break;
        }

        // If no convergence yet, apply QR with Wilkinson shift
        if (m != l) {
          if (++iter == 30) {
            std::cerr << "Too many iterations in QR decomposition!\n";
            return d; // Early return if no convergence
          }

          // Wilkinson's shift calculation
          g = (d[l + 1] - d[l]) / (2.0 * e[l]);
          r = Pythag(g, 1.0);
          g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));

          s = c = 1.0;
          p     = 0.0;
          for (i = m - 1; i >= l; --i) {
            f        = s * e[i];
            b        = c * e[i];
            e[i + 1] = (r = Pythag(f, g));
            if (r == 0.0) {
              d[i + 1] -= p;
              e[m] = 0.0;
              break;
            }
            s        = f / r;
            c        = g / r;
            g        = d[i + 1] - p;
            r        = (d[i] - g) * s + 2.0 * c * b;
            d[i + 1] = g + (p = s * r);
            g        = c * r - b;
          }

          if (r == 0.0 && i >= l)
            continue;
          d[l] -= p;
          e[l] = g;
          e[m] = 0.0;
        }
      } while (m != l);
    }

    // Return the diagonal elements as the eigenvalues
    std::sort(d.begin(), d.end());
    return d;
  }

  // Function to perform Householder QR decomposition
  void householderQR(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R)
  {
    int m = A.size();                      // Number of rows
    int n = (A.empty() ? 0 : A[0].size()); // Number of columns (assuming all rows have the same number of columns)
    //R = copyMatrix(A);
    R = A;
    Q = olb::uq::MatrixOperations::generateIdentityMatrix<T>(m);
    for (int k = 0; k < n; ++k) {
      std::vector<T> x(m - k, 0.0);
      for (int i = k; i < m; ++i) {
        x[i - k] = R[i][k];
      }

      T sign = 1.0;
      if (x[0] < 0) {
        sign = -1.0;
      }
      T norm = std::sqrt(olb::uq::MatrixOperations::dotProduct(x, x));

      std::vector<T> v = olb::uq::MatrixOperations::vectorScalarProduct(x, 1.0 / (x[0] + sign * norm));
      v[0]             = 1;
      T tau            = 2.0 / olb::uq::MatrixOperations::dotProduct(v, v);

      std::vector<std::vector<T>> H = olb::uq::MatrixOperations::generateIdentityMatrix<T>(m);

      for (int i = k; i < m; ++i) {
        for (int j = k; j < m; ++j) {
          H[i][j] -= tau * v[i - k] * v[j - k];
        }
      }

      // Update R with Householder transformation
      R = olb::uq::MatrixOperations::matrixMultiplication(H, R);

      // Update Q with Householder transformation
      Q = olb::uq::MatrixOperations::matrixMultiplication(Q, olb::uq::MatrixOperations::transposeMatrix(H));
    }
  }

  // Householder QR decomposition to compute the eigenvalues of a matrix
  std::vector<T> HouseholderQRDecomposition(const std::vector<std::vector<T>>& J)
  {
    std::size_t n = J.size();

    // Create a copy of J that we can modify
    std::vector<std::vector<T>> J_copy = olb::uq::MatrixOperations::copyMatrix(J);

    // Shift the diagonal by 1
    for (std::size_t i = 0; i < n; ++i) {
      J_copy[i][i] += 1.0;
    }

    std::vector<std::vector<T>> Q(n, std::vector<T>(n, 0.0));
    std::vector<std::vector<T>> R     = olb::uq::MatrixOperations::copyMatrix(J_copy);
    std::vector<std::vector<T>> J_new = olb::uq::MatrixOperations::copyMatrix(J_copy);
    std::vector<std::vector<T>> J_old = olb::uq::MatrixOperations::copyMatrix(Q);

    std::size_t iter     = 1;
    std::size_t max_iter = 10000;
    T           tol      = 1e-12;

    while (iter < max_iter) {
      householderQR(J_new, Q, R);
      J_new = olb::uq::MatrixOperations::matrixMultiplication(R, Q);

      if (isConverged(J_new, J_old, tol)) {
        break;
      }

      J_old.swap(J_new);
      iter++;
    }

    std::vector<T> eigenvalues(n);
    for (std::size_t i = 0; i < n; ++i) {
      eigenvalues[i] = J_old[i][i] - 1.0; // Extract diagonal as eigenvalues
    }
    std::sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
  }

  // Function to check convergence of two matrices
  bool isConverged(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B, T tol)
  {
    std::size_t n         = A.size();
    T           norm_diff = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        norm_diff += std::abs(A[i][j] - B[i][j]);
      }
    }
    // std::cout << "Norm diff: " << norm_diff << std::endl;
    return norm_diff < tol;
  }

  void computeWeights()
  {
    weights.resize(points.size());

    // If the “PolynomialBasis” is actually LegendreBasis<T>:
    if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::LegendreBasis<T>>) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        T x        = points[i];
        T Pn_prime = basis->derivativePolynomial(nq, x);
        weights[i] = 2.0 / ((1.0 - x * x) * Pn_prime * Pn_prime);
      }
    }
    // Else if it's HermiteBasis<T>:
    else if constexpr (std::is_same_v<PolynomialBasis, olb::uq::Polynomials::HermiteBasis<T>>) {
      for (std::size_t i = 0; i < points.size(); ++i) {
        T x         = points[i];
        T Hn_minus1 = basis->evaluatePolynomial(nq - 1, x);
        weights[i]  = std::pow(2, nq - 1) * std::tgamma(nq + 1) * std::sqrt(M_PI) / std::pow(nq * Hn_minus1, 2);
      }
    }
    else {
      static_assert(always_false<PolynomialBasis>::value, "Unsupported basis in Quadrature::computeWeights()");
    }

    // (optionally normalize)
    T sum = std::accumulate(weights.begin(), weights.end(), T(0));
    for (auto& w : weights) {
      w /= sum;
    }
  }
};

} // namespace Quadrature

} // namespace uq

} // namespace olb

#endif // QUADRATURE_H
