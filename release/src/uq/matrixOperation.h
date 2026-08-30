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

// matrix_operation.h
#ifndef MATRIX_OPERATION_H
#define MATRIX_OPERATION_H

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace olb {

namespace uq {

namespace MatrixOperations {

// Function declarations
template <typename T>
T dotProduct(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must be the same size for dot product.");
  }
  T result = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    result += a[i] * b[i];
  }
  return result;
}

// Function to subtract one vector from another
template <typename T>
std::vector<T> vectorSubtraction(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must be the same size for subtraction.");
  }
  std::vector<T> result(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

// Function to multiply a vector by a scalar
template <typename T>
std::vector<T> vectorScalarProduct(const std::vector<T>& v, T scalar)
{
  std::vector<T> result(v.size(), 0.0);
  for (std::size_t i = 0; i < v.size(); ++i) {
    result[i] = v[i] * scalar;
  }
  return result;
}

// Function to normalize a vector
template <typename T>
void normalize(std::vector<T>& v)
{
  T norm = std::sqrt(dotProduct(v, v));
  if (norm == 0.0) {
    throw std::runtime_error("Cannot normalize a zero vector.");
  }
  for (std::size_t i = 0; i < v.size(); ++i) {
    v[i] /= norm;
  }
}

// Function to get the column of a matrix
template <typename T>
std::vector<T> getColumn(const std::vector<std::vector<T>>& matrix, std::size_t j)
{
  if (matrix.empty() || j >= matrix[0].size()) {
    throw std::out_of_range("Invalid column index.");
  }
  std::vector<T> column(matrix.size());
  for (std::size_t i = 0; i < matrix.size(); ++i) {
    column[i] = matrix[i][j];
  }
  return column;
}

template <typename T>
std::vector<std::vector<T>> generateIdentityMatrix(int n)
{
  if (n <= 0) {
    throw std::invalid_argument("Size of identity matrix must be positive.");
  }
  std::vector<std::vector<T>> identityMatrix(n, std::vector<T>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    identityMatrix[i][i] = 1.0;
  }
  return identityMatrix;
}

// Function to add one vector to another
template <typename T>
std::vector<T> vectorAdd(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size()) {
    throw std::invalid_argument("Vectors must be the same size for addition.");
  }
  std::vector<T> result(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

// Function to multiply two matrices
template <typename T>
std::vector<std::vector<T>> matrixMultiplication(const std::vector<std::vector<T>>& A,
                                                 const std::vector<std::vector<T>>& B)
{
  if (A.empty() || B.empty()) {
    throw std::invalid_argument("Input matrices cannot be empty.");
  }
  std::size_t m = A.size();
  std::size_t n = A[0].size();
  std::size_t p = B[0].size();

  // Check dimensions
  if (n != B.size()) {
    throw std::invalid_argument("Number of columns in A must equal number of rows in B.");
  }

  // Initialize result matrix
  std::vector<std::vector<T>> C(m, std::vector<T>(p, 0.0));

  // Perform multiplication
  for (std::size_t i = 0; i < m; ++i) {
    if (A[i].size() != n) {
      throw std::invalid_argument("All rows in A must have the same number of columns.");
    }
    for (std::size_t j = 0; j < p; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

// Function to transpose a matrix
template <typename T>
std::vector<std::vector<T>> transposeMatrix(const std::vector<std::vector<T>>& A)
{
  if (A.empty()) {
    return {};
  }
  std::size_t                 m = A.size();
  std::size_t                 n = A[0].size();
  std::vector<std::vector<T>> result(n, std::vector<T>(m));
  for (std::size_t i = 0; i < m; ++i) {
    if (A[i].size() != n) {
      throw std::invalid_argument("All rows in A must have the same number of columns.");
    }
    for (std::size_t j = 0; j < n; ++j) {
      result[j][i] = A[i][j];
    }
  }
  return result;
}

// Function to compute the Frobenius norm of a matrix
template <typename T>
T frobeniusNorm(const std::vector<std::vector<T>>& matrix)
{
  T sum = 0.0;
  for (const auto& row : matrix) {
    for (T val : row) {
      sum += val * val;
    }
  }
  return std::sqrt(sum);
}

// Function to check if matrix A is converging to matrix B
template <typename T>
bool isConverging(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B, T& former_norm,
                  T tolerance)
{
  // Ensure A and B are the same size
  if (A.size() != B.size() || A.empty() || B.empty() || A[0].size() != B[0].size()) {
    throw std::invalid_argument("Matrices must be the same size and non-empty.");
  }
  std::vector<std::vector<T>> diff(A.size(), std::vector<T>(A[0].size()));
  for (std::size_t i = 0; i < A.size(); ++i) {
    for (std::size_t j = 0; j < A[0].size(); ++j) {
      diff[i][j] = A[i][j] - B[i][j];
    }
  }
  T norm = frobeniusNorm(diff);

  // Handle the first iteration
  if (former_norm == 0.0) {
    former_norm = norm;
    return false;
  }

  T relative_change = std::fabs(norm - former_norm) / former_norm;

  bool converged = relative_change <= tolerance;
  former_norm    = norm;
  return converged;
}

// Function to create a deep copy of a matrix
template <typename T>
std::vector<std::vector<T>> copyMatrix(const std::vector<std::vector<T>>& matrix)
{
  if (matrix.empty()) {
    return {};
  }

  std::vector<std::vector<T>> result;
  result.reserve(matrix.size());

  for (const auto& row : matrix) {
    result.push_back(row); // This creates a copy of each row
  }

  return result;
}

} // namespace MatrixOperations

} // namespace uq

} // namespace olb
#endif // MATRIX_OPERATION_H
