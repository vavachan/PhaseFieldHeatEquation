/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito, Adrian Kummerländer
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

#ifndef MATRIX_VIEW_H
#define MATRIX_VIEW_H

namespace olb {

/// Provides matrix view access to serialized vector fields
template <typename T, unsigned ROWS, unsigned COLS>
class MatrixView {
private:
  Vector<T,ROWS*COLS>& _data;

public:
  /// number of rows
  static constexpr unsigned rows = ROWS;
  /// number of columns
  static constexpr unsigned cols = COLS;

  /// expose data type
  using value_t = T;

  /// Row view to allow nested vector like access to the data
  class RowView {
  private:
    Vector<T,ROWS*COLS>& _data;
    unsigned _row;
  public:
    RowView(Vector<T,ROWS*COLS>& data, unsigned row) : _data(data), _row(row) { }

    T& operator[](unsigned col) {
      return _data[_row * COLS + col];
    }

    const T& operator[](unsigned col) const {
      return _data[_row * COLS + col];
    }
  };

  /// Create matrix from olb::Vector
  template <unsigned V_SIZE>
  explicit MatrixView(Vector<T, V_SIZE>& vector) : _data(vector)
  {
    static_assert(V_SIZE == ROWS * COLS,
                  "ERROR: Vector size must correspond to matrix size.");
  }

  explicit MatrixView(const MatrixView<T, ROWS, COLS>& matrix)
  {
    _data = matrix._data;
  }
  MatrixView& operator=(const MatrixView& matrix)
  {
    _data = matrix._data;
    return *this;
  }

  Vector<T,ROWS*COLS>& asVector() { return _data; }
  const Vector<T,ROWS*COLS>& asVector() const { return _data; }

  /// Returns row vector
  RowView operator[](unsigned row) {
    return RowView(_data, row);
  }
};

/// Provides transposed matrix view access to serialized vector fields
template <typename MATRIX_VIEW>
class TransposedMatrixView {
private:
  using T = typename MATRIX_VIEW::value_t;
  Vector<T,MATRIX_VIEW::rows*MATRIX_VIEW::cols>& _data;

public:
  /// number of rows after transposing
  static constexpr unsigned rows = MATRIX_VIEW::cols;
  /// number of columns after transposing
  static constexpr unsigned cols = MATRIX_VIEW::rows;

  /// expose data type
  using value_t = T;

  /// Column view to allow nested vector like access to the data
  class ColView {
  private:
    Vector<T,MATRIX_VIEW::rows*MATRIX_VIEW::cols>& _data;
    unsigned _col;
  public:
    ColView(Vector<T,MATRIX_VIEW::rows*MATRIX_VIEW::cols>& data, unsigned col) : _data(data), _col(col) { }

    T& operator[](unsigned row) {
      return _data[row * MATRIX_VIEW::cols + _col];
    }

    const T& operator[](unsigned row) const {
      return _data[row * MATRIX_VIEW::cols + _col];
    }
  };

  /// Create matrix from olb::Vector
  template <unsigned V_SIZE>
  explicit TransposedMatrixView(Vector<T, V_SIZE>& vector) : _data(vector)
  {
    static_assert(V_SIZE == rows * cols,
                  "ERROR: Vector size must correspond to matrix size.");
  }

  explicit TransposedMatrixView(const TransposedMatrixView<MATRIX_VIEW>& matrix)
  {
    _data = matrix._data;
  }
  TransposedMatrixView& operator=(const TransposedMatrixView& matrix)
  {
    _data = matrix._data;
    return *this;
  }

  Vector<T,rows*cols>& asVector() { return _data; }
  const Vector<T,rows*cols>& asVector() const { return _data; }

  /// Returns column vector
  ColView operator[](unsigned col) {
    return ColView(_data, col);
  }
};

template<typename T, typename DESCRIPTOR, typename LINEAR_FIELD>
using MatrixD = MatrixView<T,LINEAR_FIELD::template rows<DESCRIPTOR>(),LINEAR_FIELD::template cols<DESCRIPTOR>()>;

template<typename T, typename DESCRIPTOR, typename LINEAR_FIELD>
using TransposedMatrixD = TransposedMatrixView<MatrixView<T,LINEAR_FIELD::template rows<DESCRIPTOR>(),LINEAR_FIELD::template cols<DESCRIPTOR>()>>;

}

#endif
