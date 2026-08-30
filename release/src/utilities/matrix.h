/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Jan Eric Marquardt, 2023 Julius Jessberger
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

#ifndef MATRIX_H
#define MATRIX_H

#include "vectorHelpers.h"

namespace olb {

/// Matrix with a defined number of ROWS and columns (COLS)
template <typename T, unsigned ROWS, unsigned COLS>
class Matrix {
private:
  Vector<Vector<T, COLS>, ROWS> _data;

public:
  /// number of rows
  static constexpr unsigned rows = ROWS;
  /// number of columns
  static constexpr unsigned cols = COLS;

  constexpr Matrix() {};
  constexpr Matrix(const T data[ROWS][COLS])
  {
    for (unsigned m = 0; m < ROWS; ++m) {
      for (unsigned n = 0; n < COLS; ++n) {
        _data[m][n] = data[m][n];
      }
    }
  }

  /// Create matrix from olb::Vector
  template <unsigned V_SIZE>
  constexpr Matrix(const Vector<T, V_SIZE>& vector)
  {
    static_assert(V_SIZE == ROWS * COLS,
                  "ERROR: Vector size must correspond to matrix size.");
    for (unsigned m = 0; m < ROWS; ++m) {
      for (unsigned n = 0; n < COLS; ++n) {
        _data[m][n] = vector[m * COLS + n];
      }
    }
  }

  constexpr Matrix(const Matrix<T, ROWS, COLS>& matrix)
  {
    _data = matrix._data;
  }
  constexpr Matrix& operator=(const Matrix& matrix)
  {
    _data = matrix._data;
    return *this;
  }
  constexpr Matrix(Matrix&& matrix) { _data = std::move(matrix._data); }
  constexpr Matrix& operator=(Matrix&& matrix)
  {
    _data = std::move(matrix._data);
    return *this;
  }

  constexpr T*       data() { return _data.data(); }
  constexpr const T* data() const { return _data.data(); }

  constexpr Vector<T, COLS>& operator[](const unsigned row)
  {
    return _data[row];
  };
  constexpr const Vector<T, COLS>& operator[](const unsigned row) const
  {
    return _data[row];
  }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix)
  {
    for (unsigned row = 0; row < matrix.rows; ++row) {
      os << matrix[row] << '\n';
    }
    return os;
  }

  template <unsigned C>
  constexpr inline Matrix operator*(const Matrix<T, COLS, C>& matrix) const
  {
    Matrix<T,ROWS,C> result;
    for (unsigned i = 0; i < ROWS; ++i) {
      for (unsigned j = 0; j < C; ++j) {
        result[i][j] = 0;
        for (unsigned k = 0; k < ROWS; ++k) {
          result[i][j] += this->operator[](i)[k] * matrix[k][j];
        }
      }
    }
    return result;
  }

  constexpr inline Vector<T, COLS>
  operator*(const Vector<T, COLS>& vector) const
  {
    Vector<T, COLS> result;
    for (unsigned row = 0; row < ROWS; ++row) {
      result[row] = _data[row] * vector;
    }
    return result;
  }

  constexpr inline Matrix<T, COLS, ROWS> transpose() const
  {
    Matrix<T,COLS,ROWS> result;
    for (unsigned m = 0; m < ROWS; ++m) {
      for (unsigned n = 0; n < COLS; ++n) {
        result[n][m] = ((*this)[m][n]);
      }
    }
    return result;
  }
};


namespace util {

template <typename T>
constexpr auto determinant(const Matrix<T,2,2>& m)
{
  return determinant<T>(m[0], m[1]);
}

template <typename T>
constexpr auto determinant(const Matrix<T,3,3>& m)
{
  return determinant<T>(m[0], m[1], m[2]);
}

/// Solve a * x = rhs
// Returns +-inf if a[0], a[1] are linearly dependent and rhs is not (no solution)
// Returns nan if a[0], a[1], rhs are linearly dependent (infinitely many solutions)
template <typename T>
auto solveLinearSystem(const Matrix<T,2,2>& a, const Vector<T,2>& rhs)
{
  return solveLinearSystem(a[0], a[1], rhs);
}

/// Solve (a0 a1 a2) x = rhs
// Returns +-inf if a0, a1, a2 are linearly dependent and rhs is not (no solution)
// Returns nan if a0, a1, a2, rhs are linearly dependent (infinitely many solutions)
template <typename T>
auto solveLinearSystem(const Matrix<T,3,3>& a, const Vector<T,3>& rhs)
{
  return solveLinearSystem(a[0], a[1], a[2], rhs);
}

}

} // namespace olb

#endif
