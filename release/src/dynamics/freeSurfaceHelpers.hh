/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
 *                2024-2025 Danial Khazaeipoul
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

#include <atomic>

namespace olb {

namespace FreeSurface {

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T,DESCRIPTOR>& sLattice) {
  sLattice.executePostProcessors(FreeSurface::Stage0());
  sLattice.executePostProcessors(FreeSurface::Stage1());
  sLattice.executePostProcessors(FreeSurface::Stage2());
  sLattice.executePostProcessors(FreeSurface::Stage3());
  sLattice.executePostProcessors(FreeSurface::Stage4());
  sLattice.executePostProcessors(FreeSurface::Stage5());
}

template <typename CELL>
bool isCellType(CELL& cell, const FreeSurface::Type& type) {
  return cell.template getField<FreeSurface::CELL_TYPE>() == type;
}

template <typename CELL>
void setCellType(CELL& cell, const FreeSurface::Type& type) {
  cell.template setField<FreeSurface::CELL_TYPE>(type);
}

template <typename CELL>
bool hasNeighbour(CELL& cell, const FreeSurface::Type& type) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if (isCellType(nbrCell, type)) { return true; }
  }

  return false;
}

template <typename CELL>
bool hasCellFlags(CELL& cell, const FreeSurface::Flags& flags) {
  return static_cast<bool>(cell.template getField<FreeSurface::CELL_FLAGS>() & flags);
}

template <typename CELL>
void setCellFlags(CELL& cell, const FreeSurface::Flags& flags) {
  cell.template setField<FreeSurface::CELL_FLAGS>(flags);
}

template <typename CELL>
bool hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if (hasCellFlags(nbrCell, flags)) { return true; }
  }

  return false;
}

template <typename CELL>
NeighbourInfo getNeighbourInfo(CELL& cell) {
  NeighbourInfo info{};
  using DESCRIPTOR = typename CELL::descriptor_t;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    if (isCellType(nbrCell, FreeSurface::Type::Gas)) {
      info.has_gas_neighbours = true;
    }
    else if (isCellType(nbrCell, FreeSurface::Type::Fluid)) {
      info.has_fluid_neighbours = true;
    }
    else if (isCellType(nbrCell, FreeSurface::Type::Interface)) {
      ++info.interface_neighbours;
    }
  }

  return info;
}

template <typename CELL, typename V>
bool isHealthyInterface(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;

  if (!isCellType(cell, FreeSurface::Type::Interface)) { return false; }

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    if (isCellType(nbrCell, FreeSurface::Type::Gas)) {
      has_gas_neighbours = true;
      if (has_fluid_neighbours) { return true; }
    }
    else if (isCellType(nbrCell, FreeSurface::Type::Fluid)) {
      has_fluid_neighbours = true;
      if (has_gas_neighbours) { return true; }
    }
  }

  return false;
}

template <typename CELL, typename V>
void computeMassExcessWeights
(
  CELL& cell, const Vector<V, CELL::descriptor_t::d>& normal,
  const bool& enableAllInterfaces, Vector<V, CELL::descriptor_t::q>& weights
) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto direction = descriptors::c<DESCRIPTOR>(iPop);
    auto nbrCell = cell.neighbor(direction);

    Vector<V, DESCRIPTOR::d> ei{};
    if constexpr (DESCRIPTOR::d == 3) {
      ei = {V(direction[0]), V(direction[1]), V(direction[2])};
    }
    else { ei = {V(direction[0]), V(direction[1])}; }

    if (isCellType(nbrCell, FreeSurface::Type::Interface) && hasCellFlags(nbrCell, FreeSurface::Flags::None)) {
      const V nDotDirection = util::dotProduct(normal, ei);

      if (hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
        // This cell was converted from interface to fluid, so the normal vector is used as is
        weights[iPop] = nDotDirection > V(0) ? nDotDirection : V(0);
      }
      else if (hasCellFlags(cell, FreeSurface::Flags::ToGas)) {
        // This cell was converted from interface to gas, so the normal vector is inverted
        weights[iPop] = nDotDirection < V(0) ? -nDotDirection : V(0);
      }
    }
    else if (hasCellFlags(nbrCell, FreeSurface::Flags::NewInterface) && enableAllInterfaces) {
      const V nDotDirection = util::dotProduct(normal, ei);

      if (hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
        // This cell was converted from interface to fluid, so the normal vector is used as is
        weights[iPop] = nDotDirection > V(0) ? nDotDirection : V(0);
      }
      else if (hasCellFlags(cell, FreeSurface::Flags::ToGas)) {
        // This cell was converted from interface to gas, so the normal vector is inverted
        weights[iPop] = nDotDirection < V(0) ? -nDotDirection : V(0);
      }
    }
    else {
      // If the neighbour is not an interface cell, the weight is set to zero
      weights[iPop] = V(0);
    }
  }
}

template <typename CELL, typename V>
V getClampedEpsilon(const CELL& cell) {

  V epsilon = cell.template getField<FreeSurface::EPSILON>();
  return util::max(V(0), util::min(V(1), epsilon));
}

template <typename CELL, typename V>
V getSmoothEpsilon(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  // K8 kernel of Williams et al.
  // Avoid uisng expensive std::pow to boost performance
  auto kernel = [](auto radius, const auto& direction) {
    auto norm = norm_squared(direction);
    if (norm < radius * radius) {
      auto tmp = decltype(radius)(1) - norm / (radius * radius);
      return tmp * tmp * tmp * tmp;
    } else { return decltype(radius)(0); }
  };

  // Kernel support radius is fixed at 2.0
  V radius = V(2);
  V denominator = V(0);
  V tmp = V(0);

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));
    const auto direction = descriptors::c<DESCRIPTOR>(iPop);
    tmp += kernel(radius, direction) * nbrCell.template getField<FreeSurface::EPSILON>();
    if (norm(direction) < radius) { denominator += kernel(radius, direction); }
  }

  return tmp / denominator;
}

template <typename CELL, typename V>
V getClampedSmoothEpsilon(const CELL& cell) {

  V epsilon = getSmoothEpsilon(cell);
  return util::max(V(0), util::min(V(1), epsilon));
}

template <typename CELL, typename V>
Vector<V, CELL::descriptor_t::d> computeInterfaceNormal(CELL& cell) {

  // Compute interface normal using Parker-Young approximation
  Vector<V, CELL::descriptor_t::d> normal{};
  if constexpr (CELL::descriptor_t::d == 3) {
    auto eps_minusY             = getClampedEpsilon(cell.neighbor({0, -1, 0}));
    auto eps_plusY              = getClampedEpsilon(cell.neighbor({0, 1, 0}));
    auto eps_minusX             = getClampedEpsilon(cell.neighbor({-1, 0, 0}));
    auto eps_plusX              = getClampedEpsilon(cell.neighbor({1, 0, 0}));
    auto eps_minusYminusX       = getClampedEpsilon(cell.neighbor({-1, -1, 0}));
    auto eps_minusYplusX        = getClampedEpsilon(cell.neighbor({1, -1, 0}));
    auto eps_plusYminusX        = getClampedEpsilon(cell.neighbor({-1, 1, 0}));
    auto eps_plusYplusX         = getClampedEpsilon(cell.neighbor({1, 1, 0}));
    auto eps_minusZ             = getClampedEpsilon(cell.neighbor({0, 0, -1}));
    auto eps_plusZ              = getClampedEpsilon(cell.neighbor({0, 0, 1}));
    auto eps_minusZminusY       = getClampedEpsilon(cell.neighbor({0, -1, -1}));
    auto eps_minusZplusY        = getClampedEpsilon(cell.neighbor({0, 1, -1}));
    auto eps_minusZminusX       = getClampedEpsilon(cell.neighbor({-1, 0, -1}));
    auto eps_minusZplusX        = getClampedEpsilon(cell.neighbor({1, 0, -1}));
    auto eps_plusZminusY        = getClampedEpsilon(cell.neighbor({0, -1, 1}));
    auto eps_plusZplusY         = getClampedEpsilon(cell.neighbor({0, 1, 1}));
    auto eps_plusZminusX        = getClampedEpsilon(cell.neighbor({-1, 0, 1}));
    auto eps_plusZplusX         = getClampedEpsilon(cell.neighbor({1, 0, 1}));
    auto eps_minusZminusYminusX = getClampedEpsilon(cell.neighbor({-1, -1, -1}));
    auto eps_minusZplusYminusX  = getClampedEpsilon(cell.neighbor({-1, 1, -1}));
    auto eps_minusZminusYplusX  = getClampedEpsilon(cell.neighbor({1, -1, -1}));
    auto eps_minusZplusYplusX   = getClampedEpsilon(cell.neighbor({1, 1, -1}));
    auto eps_plusZminusYminusX  = getClampedEpsilon(cell.neighbor({-1, -1, 1}));
    auto eps_plusZplusYminusX   = getClampedEpsilon(cell.neighbor({-1, 1, 1}));
    auto eps_plusZminusYplusX   = getClampedEpsilon(cell.neighbor({1, -1, 1}));
    auto eps_plusZplusYplusX    = getClampedEpsilon(cell.neighbor({1, 1, 1}));

    // Rearrange arithmetic operations to minimize the number of cancellation events, i.e., positive and negative terms.
    normal[0]  = V(4) * (eps_minusX - eps_plusX);
    normal[0] += V(2) * ((eps_minusYminusX + eps_minusZminusX + eps_plusYminusX + eps_plusZminusX) - (eps_plusYplusX + eps_plusZplusX + eps_minusYplusX + eps_minusZplusX));
    normal[0] += V(1) * ((eps_minusZminusYminusX + eps_plusZminusYminusX + eps_minusZplusYminusX + eps_plusZplusYminusX) - (eps_plusZplusYplusX + eps_minusZplusYplusX + eps_plusZminusYplusX + eps_minusZminusYplusX));

    normal[1]  = V(4) * (eps_minusY - eps_plusY);
    normal[1] += V(2) * ((eps_minusYminusX + eps_minusZminusY + eps_minusYplusX + eps_plusZminusY) - (eps_plusYplusX + eps_plusZplusY + eps_plusYminusX + eps_minusZplusY));
    normal[1] += V(1) * ((eps_minusZminusYminusX + eps_plusZminusYminusX + eps_plusZminusYplusX + eps_minusZminusYplusX) - (eps_plusZplusYplusX + eps_minusZplusYplusX + eps_minusZplusYminusX + eps_plusZplusYminusX));

    normal[2]  = V(4) * (eps_minusZ - eps_plusZ);
    normal[2] += V(2) * ((eps_minusZminusX + eps_minusZminusY + eps_minusZplusX + eps_minusZplusY) - (eps_plusZplusX + eps_plusZplusY + eps_plusZminusX + eps_plusZminusY));
    normal[2] += V(1) * ((eps_minusZminusYminusX + eps_minusZplusYplusX + eps_minusZplusYminusX + eps_minusZminusYplusX) - (eps_plusZplusYplusX + eps_plusZminusYminusX + eps_plusZminusYplusX + eps_plusZplusYminusX));

    return normal;
  }
  else {
    auto eps_minusY       = getClampedEpsilon(cell.neighbor({0, -1}));
    auto eps_plusY        = getClampedEpsilon(cell.neighbor({0, 1}));
    auto eps_minusX       = getClampedEpsilon(cell.neighbor({-1, 0}));
    auto eps_plusX        = getClampedEpsilon(cell.neighbor({1, 0}));
    auto eps_minusYminusX = getClampedEpsilon(cell.neighbor({-1, -1}));
    auto eps_minusYplusX  = getClampedEpsilon(cell.neighbor({1, -1}));
    auto eps_plusYminusX  = getClampedEpsilon(cell.neighbor({-1, 1}));
    auto eps_plusYplusX   = getClampedEpsilon(cell.neighbor({1, 1}));

    // Rearrange arithmetic operations to minimize the number of cancellation events, i.e., positive and negative terms.
    normal[0] = V(2) * (eps_minusX - eps_plusX);
    normal[0] += V(1) * ((eps_minusYminusX + eps_plusYminusX) - (eps_plusYplusX + eps_minusYplusX));

    normal[1] = V(2) * (eps_minusY - eps_plusY);
    normal[1] += V(1) * ((eps_minusYminusX + eps_minusYplusX) - (eps_plusYplusX + eps_plusYminusX));

    return normal;
  }
}

template<typename T, int N>
void iterativeRefinement
(
  const std::array<T, N * N>& A, const std::array<T, N * N>& M,
  std::array<T, N>& x, const std::array<T, N>& b, const int Nsol, int maxIter
) {
  std::array<T, N> residual{};
  for (int iter = 0; iter < maxIter; ++iter) {
    // Compute the residual, i.e., r = b - A * x
    for (int i = 0; i < Nsol; ++i) {
      T sum = T(0);
      for (int j = 0; j < Nsol; ++j) { sum += A[i * N + j] * x[j]; }
      residual[i] = b[i] - sum;
    }

    // Forward substitution
    std::array<T, N> dx = residual;
    for (int i = 0; i < Nsol; ++i) {
      for (int k = 0; k < i; ++k) { dx[i] -= M[i * N + k] * dx[k]; }
    }

    // Backward substitution
    for (int i = Nsol - 1; i >= 0; --i) {
      for (int k = i + 1; k < Nsol; ++k) { dx[i] -= M[i * N + k] * dx[k]; }
      dx[i] /= M[i * N + i];
    }

    // Update the solution, i.e., x = x + dx
    for (int i = 0; i < Nsol; ++i) { x[i] += dx[i]; }
  }
}

template<typename T, int N>
void rowPivotingLUSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) {
  // Add Tikhonov Regularization to the diagonal elements of M
  for (int i = 0; i < Nsol; ++i) { M[i * N + i] += T(0); }

  // Make a copy of the original coefficient matrix.
  std::array<T, N * N> Prev_M = M;

  // Permutation vector to record row swaps.
  std::array<int, N> row_perm{};
  for (int i = 0; i < N; ++i) { row_perm[i] = i; }

  // LU decomposition method with row pivoting.
  for (int i = 0; i < Nsol; ++i) {
    int pivotRow = i;
    T maxValue = std::abs(M[i * N + i]);

    for (int j = i + 1; j < Nsol; ++j) {
      T value = std::abs(M[j * N + i]);
      if (value > maxValue) {
        maxValue = value;
        pivotRow = j;
      }
    }

    if (pivotRow != i) {
      for (int k = 0; k < Nsol; ++k) { std::swap(M[i * N + k], M[pivotRow * N + k]); }
      std::swap(row_perm[i], row_perm[pivotRow]);
      std::swap(b[i], b[pivotRow]);
    }

    for (int j = i + 1; j < Nsol; ++j) {
      M[j * N + i] /= M[i * N + i];
      for (int k = i + 1; k < Nsol; ++k) { M[j * N + k] -= M[j * N + i] * M[i * N + k]; }
    }
  }

  // Copy the right-hand side array into x array
  for (int i = 0; i < Nsol; ++i) { x[i] = b[i]; }

  // Forward substitution step, i.e., solve L * x = b in place.
  for (int i = 0; i < Nsol; ++i) {
    for (int j = 0; j < i; ++j) { x[i] -= M[i * N + j] * x[j]; }
  }

  // Backward substitution step, i.e., solve U * x = b in place.
  for (int i = Nsol - 1; i >= 0; --i) {
    for (int j = i + 1; j < Nsol; ++j) { x[i] -= M[i * N + j] * x[j]; }
    x[i] /= M[i * N + i];
  }

  // Perform iterative refinement of the solution.
  std::array<T, N * N> tmp = Prev_M;
  for (int i = 0; i < Nsol; ++i) {
    for (int j = 0; j < Nsol; ++j) { Prev_M[i * N + j] = tmp[row_perm[i] * N + j]; }
  }
  iterativeRefinement<T, N>(Prev_M, M, x, b, Nsol);
}

template<typename T, int N>
void completePivotingLUSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) {
  // Add Tikhonov Regularization to the diagonal elements of M
  for (int i = 0; i < Nsol; ++i) { M[i * N + i] += T(0); }

  // Permutation vectors to record row and columns swaps.
  std::array<int, N> row_perm{}, col_perm{};
  for (int i = 0; i < N; ++i) {
    row_perm[i] = i;
    col_perm[i] = i;
  }

  // LU decomposition method with complete pivoting.
  for (int i = 0; i < Nsol; ++i) {
    int pivotRow = i, pivotCol = i;
    T maxValue = std::abs(M[i * N + i]);

    for (int r = i; r < Nsol; ++r) {
      for (int c = i; c < Nsol; ++c) {
        T value = std::abs(M[r * N + c]);
        if (value > maxValue) {
          maxValue = value;
          pivotRow = r;
          pivotCol = c;
        }
      }
    }

    if (pivotRow != i) {
      for (int k = 0; k < Nsol; ++k) { std::swap(M[i * N + k], M[pivotRow * N + k]); }
      std::swap(row_perm[i], row_perm[pivotRow]);
      std::swap(b[i], b[pivotRow]);
    }

    if (pivotCol != i) {
      for (int k = 0; k < Nsol; ++k) { std::swap(M[k * N + i], M[k * N + pivotCol]); }
      std::swap(col_perm[i], col_perm[pivotCol]);
    }

    for (int j = i + 1; j < Nsol; ++j) {
      M[j * N + i] /= M[i * N + i];
      for (int k = i + 1; k < Nsol; ++k) { M[j * N + k] -= M[j * N + i] * M[i * N + k]; }
    }
  }

  // Copy the right-hand side array into x array
  for (int i = 0; i < Nsol; ++i) { x[i] = b[i]; }

  // Forward substitution step, i.e., solve L * x = b in place.
  for (int i = 0; i < Nsol; ++i) {
    for (int j = 0; j < i; ++j) { x[i] -= M[i * N + j] * x[j]; }
  }

  // Backward substitution step, i.e., solve U * x = b in place.
  for (int i = Nsol - 1; i >= 0; --i) {
    for (int j = i + 1; j < Nsol; ++j) { x[i] -= M[i * N + j] * x[j]; }
    x[i] /= M[i * N + i];
  }

  // Reconstruct original x by inverting the permutation, i.e.,
  // for each index i, the original index is given by col_perm[i].
  std::array<T, N> tmp = x;
  for (int i = 0; i < Nsol; ++i) { x[col_perm[i]] = tmp[i]; }
}

template<typename T, int N>
void columnPivotingQRSolver(std::array<T, N * N>& M, std::array<T, N>& x, std::array<T, N>& b, const int Nsol) {
  // Permutation array to record column swaps.
  std::array<int, N> col_perm{};
  for (int i = 0; i < N; ++i) { col_perm[i] = i; }

  // Pre-computed column norm array, i.e., L2-normalization of each column.
  std::array<T, N> col_norm{}, col_norm_tmp{};
  for (int i = 0; i < Nsol; ++i) {
    T norm = T(0);
    for (int j = 0; j < Nsol; ++j) { norm += M[j * N + i] * M[j * N + i]; }
    col_norm[i] = util::sqrt(norm);
    col_norm_tmp[i] = col_norm[i];
  }

  // QR decomposition method with column pivoting.
  for (int i = 0; i < Nsol; ++i) {
    int pivotCol = i;
    T maxValue = col_norm[i];

    for (int j = i + 1; j < Nsol; ++j) {
      if (col_norm[j] > maxValue) {
        maxValue = col_norm[j];
        pivotCol = j;
      }
    }

    if (pivotCol != i) {
      for (int k = 0; k < Nsol; ++k) { std::swap(M[k * N + i], M[k * N + pivotCol]); }
      std::swap(col_perm[i], col_perm[pivotCol]);
      std::swap(col_norm[i], col_norm[pivotCol]);
      std::swap(col_norm_tmp[i], col_norm_tmp[pivotCol]);
    }

    // Compute the L2-norm of the i-th column, i.e., range is from row (i) to (Nsol - 1),
    // and the reflection factor (r)
    T r = T(0);
    {
      T norm = T(0);
      for (int k = i; k < Nsol; ++k) { norm += M[k * N + i] * M[k * N + i]; }
      norm = util::sqrt(norm);
      r = M[i * N + i] >= T(0) ? -norm : norm;
    }

    // Prepare and normalize the Householder array for the i-th column.
    std::array<T, N> v{};
    {
      v[0] = M[i * N + i] - r;
      for (int k = i + 1; k < Nsol; ++k) { v[k - i] = M[k * N + i]; }

      T norm = T(0);
      for (int k = 0; k < Nsol - i; ++k) { norm += v[k] * v[k]; }
      norm = util::sqrt(norm);

      // Otherwise, the Householder array is already zero and no reflection is applied.
      if (norm != T(0)) {
        for (int k = 0; k < Nsol - i; ++k) { v[k] /= norm; }
      }
    }

    // Apply the Householder reflection to the remaining columns.
    for (int k = i; k < Nsol; ++k) {
      T dot_product = T(0);
      for (int j = 0; j < Nsol - i; ++j) { dot_product += v[j] * M[(j + i) * N + k]; }

      dot_product *= T(2);
      for (int j = 0; j < Nsol - i; ++j) { M[(j + i) * N + k] -= v[j] * dot_product; }
    }

    // Apply the Householder reflection to the right-hand side array.
    {
      T dot_product = T(0);
      for (int k = 0; k < Nsol - i; ++k) { dot_product += v[k] * b[k + i]; }

      dot_product *= T(2);
      for (int k = 0; k < Nsol - i; ++k) { b[k + i] -= v[k] * dot_product; }
    }

    // Update the column norm array, i.e., adapted from LAPACK’s DGEQP3 routine.
    for (int k = i + 1; k < Nsol; ++k) {
      T tmp = col_norm[k] * col_norm[k] - M[i * N + k] * M[i * N + k];

      // Recompute L2-norm from scratch, if tmp becomes negative or it drops below
      // a specified threshold, e.g., T(0.1)
      if (tmp <= T(0) || util::sqrt(tmp) < T(0.1) * col_norm_tmp[k]) {
        T norm = T(0);
        for (int j = i + 1; j < Nsol; ++j) { norm += M[j * N + k] * M[j * N + k]; }
        col_norm[k] = util::sqrt(norm);
        col_norm_tmp[k] = col_norm[k];
      }
      else {
        col_norm[k] = util::sqrt(tmp);
      }
    }

    // For the i-th column, replace diagonal elements with (r) and the non-diagonal elements with zero.
    M[i * N + i] = r;
    for (int k = i + 1; k < Nsol; ++k) { M[k * N + i] = T(0); }
  }

  // Backward substitution step, i.e., solve R*x = Qᵀ*b in place.
  for (int i = Nsol - 1; i >= 0; --i) {
    T sum = T(0);
    for (int j = i + 1; j < Nsol; ++j) { sum += M[i * N + j] * x[j]; }
    x[i] = (b[i] - sum) / M[i * N + i];
  }

  // Reconstruct x by inverting the permutation, i.e., for each
  // index i, the original index is given by col_perm[i].
  std::array<T, N> tmp = x;
  for (int i = 0; i < Nsol; ++i) { x[col_perm[i]] = tmp[i]; }
}

// Offset helper: optimized PLIC algorithm taken from Moritz Lehmann FluidX3D
template<typename T>
T plicCubeReduced(const T& volume, const Vector<T, 3>& n) {
  const T n1 = n[0], n2 = n[1], n3 = n[2];
  const T n12 = n1 + n2, n3V = n3 * volume;

  // (I) Case 5
  if (n12 <= T(2) * n3V) { return n3V + T(0.5) * n12; }

  // (II) After case 5, check n2 > 0 is true
  const T sqn1 = util::sqr(n1), n26 = T(6) * n2, v1 = sqn1 / n26;

  // (III) Case 2
  if (v1 <= n3V && n3V < v1 + T(0.5) * (n2 - n1)) {
    return T(0.5) * (n1 + util::sqrt(sqn1 + T(8) * n2 * (n3V - v1)));
  }

  // (IV) Case 1
  const T V6 = n1 * n26 * n3V;
  if (n3V < v1) { return std::cbrt(V6); }

  // (V) After case 2, check n1 > 0 is true
  const T v3 = n3 < n12 ? (util::sqr(n3) * (T(3) * n12 - n3) + sqn1 * (n1 - T(3) * n3) + util::sqr(n2) * (n2 - T(3) * n3)) / (n1 * n26) : T(0.5) * n12;
  const T sqn12 = sqn1 + util::sqr(n2), V6cbn12 = V6 - util::cube(n1) - util::cube(n2);
  const bool case34 = n3V < v3; // true: case (3), false: case (4)
  const T a = case34 ? V6cbn12 : T(0.5) * (V6cbn12 - util::cube(n3));
  const T b = case34 ? sqn12 : T(0.5) * (sqn12 + util::sqr(n3));
  const T c = case34 ? n12 : T(0.5);
  const T t = util::sqrt(util::sqr(c) - b);
  static constexpr T oneThird = T(1) / T(3);

  return c - T(2) * t * util::sin(oneThird * util::asin((util::cube(c) - T(0.5) * a - T(1.5) * b * c) / util::cube(t)));
}

// Offset computer: optimized PLIC algorithm taken from Moritz Lehmann FluidX3D
template<typename T>
T plicCube(const T& volume, const Vector<T, 3>& n) {
  // Unit cube-plane intersection: using volume and normal vector to compute the offset
  const T ax = util::abs(n[0]), ay = util::abs(n[1]), az = util::abs(n[2]);

  // Eliminate symmetry cases, and normalize n using L1 norm
  const T V = T(0.5) - util::abs(volume - T(0.5)), l = ax + ay + az;
  const T n1 = util::min(util::min(ax, ay), az) / l;
  const T n3 = util::max(util::max(ax, ay), az) / l;
  const T n2 = std::fdim(T(1), n1 + n3); // ensure n2 >= 0
  const T d = plicCubeReduced<T>(V, {n1, n2, n3});

  return l * std::copysign(T(0.5) - d, volume - T(0.5));
}

// Volume computer: optimized inverse PLIC algorithm taken from Moritz Lehmann FluidX3D
template<typename T>
T plicCubeInverse(const T& d0, const Vector<T, 3>& n) {
  // Eliminate majority of cases due to symmetry
  const T n1 = util::min(util::min(util::abs(n[0]), util::abs(n[1])), util::abs(n[2]));
  const T n3 = util::max(util::max(util::abs(n[0]), util::abs(n[1])), util::abs(n[2]));
  const T n2 = util::abs(n[0]) - n1 + util::abs(n[1]) + util::abs(n[2]) - n3;

  // Compute PLIC using reduced symmetry, shifting origin from (0.0, 0.0, 0.0) to (0.5, 0.5, 0.5)
  const T d = T(0.5) * (n1 + n2 + n3) - util::abs(d0);

  T V = T(0);
  if (util::min(n1 + n2, n3) <= d && d <= n3) {
    // (I) Case 5
    V = (d - T(0.5) * (n1 + n2)) / n3;
  }
  else if (d < n1) {
    // (II) Case 1
    V = util::cube(d) / (T(6) * n1 * n2 * n3);
  }
  else if (d <= n2) {
    // (III) Case 2
    V = (T(3) * d * (d - n1) + util::sqr(n1)) / (T(6) * n2 * n3);
  }
  else {
    // (IV) Case 3 or 4
    V = (util::cube(d) - util::cube(d - n1) - util::cube(d - n2) - util::cube(std::fdim(d, n3))) / (T(6) * n1 * n2 * n3);
  }

  return std::copysign(T(0.5) - V, d0) + V(0.5);
}

// Compute 2D interface curvature using Moritz Lehmann's FluidX3D algorithm
template <typename CELL, typename V>
V computeCurvaturePLIC2D(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  // Setup a new coordinate system where bz is perpendicular to the surface, while bx and by are
  // tangent to the surface.
  const auto normal = computeInterfaceNormal(cell);
  Vector<V, 3> by = {normal[0], normal[1], V(0)};
  if (norm_squared(by) < zeroThreshold<V>()) { return V(0); }
  by = util::normalize(by);

  const Vector<V, 3> rn{V(0), V(0), V(1)};
  const Vector<V, 3> bx = crossProduct(by, rn);

  // Number of healthy interface neighbours
  std::uint8_t nbrInterfaces = std::uint8_t(0);

  // Compute z-offset of center point using PLIC cube
  const V centerOffset = plicCube<V>(getClampedEpsilon(cell), by);

  // Number of neighbouring interface points, i.e., less than or equal to 8 (neighbours) - 1 (liquid) - 1 (gas).
  std::array<Vector<V, DESCRIPTOR::d>, 6> points;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    // Check if the neighbour is an interface cell or has a gas neighbour
    if (!isCellType(nbrCell, FreeSurface::Type::Interface) || !hasNeighbour(nbrCell, FreeSurface::Type::Gas)) {
      continue;
    }

    // Here it is assumed that neighbor normal vector is the same as center normal vector
    const auto direction = descriptors::c<DESCRIPTOR>(iPop);
    const Vector<V, 3> ei = {V(direction[0]), V(direction[1]), V(0)};
    const V offset = plicCube<V>(getClampedEpsilon(nbrCell), by) - centerOffset;

    // Transform coordinate system into (x, f(x)) and apply PLIC offset
    points[nbrInterfaces++] = {util::dotProduct(ei, bx), util::dotProduct(ei, by) + offset};
  }

  // Prepare the linear equation system, to be initialized with zeros
  std::array<V, 4> M{};
  std::array<V, 2> solution{};
  std::array<V, 2> b{};

  for (std::uint8_t i = 0; i < nbrInterfaces; ++i) {
    const V x = points[i][0];
    const V y = points[i][1];
    const V x2 = util::sqr(x);
    const V x3 = util::sqr(x) * x;

    M[0] += x2 * x2;
    M[1] += x3;
    M[3] += x2;

    b[0] += x2 * y;
    b[1] += x * y;
  }

  // Fill the lower triangle of the symmetric matrix
  M[2] = M[1];

  // Solve the linear system of equations, compile-time selection.
  if constexpr (solver == SolverType::rowPivotingLU) {
    if (nbrInterfaces >= 2) {
      rowPivotingLUSolver<V, 2>(M, solution, b, 2);
    } else {
      rowPivotingLUSolver<V, 2>(M, solution, b, util::min(2, nbrInterfaces));
    }
  }
  else if constexpr (solver == SolverType::completePivotingLU) {
    if (nbrInterfaces >= 2) {
      completePivotingLUSolver<V, 2>(M, solution, b, 2);
    } else {
      completePivotingLUSolver<V, 2>(M, solution, b, util::min(2, nbrInterfaces));
    }
  }
  else if constexpr (solver == SolverType::columnPivotingQR) {
    if (nbrInterfaces >= 2) {
      columnPivotingQRSolver<V, 2>(M, solution, b, 2);
    } else {
      columnPivotingQRSolver<V, 2>(M, solution, b, util::min(2, nbrInterfaces));
    }
  }
  else {
    // Fall back to default solver
    if (nbrInterfaces >= 2) {
      completePivotingLUSolver<V, 2>(M, solution, b, 2);
    } else {
      completePivotingLUSolver<V, 2>(M, solution, b, util::min(2, nbrInterfaces));
    }
  }

  const V A = solution[0], H = solution[1];
  const V curvature = V(2) * A * util::cube(V(1) / util::sqrt(H * H + V(1)));
  return util::max(V(-1), util::min(V(1), curvature));
}

// Compute 3D interface curvature using Moritz Lehmann's FluidX3D algorithm
template <typename CELL, typename V>
V computeCurvaturePLIC3D(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  // Setup a new coordinate system where bz is perpendicular to the surface, while bx and by are
  // tangent to the surface.
  auto bz = computeInterfaceNormal(cell);
  if (norm_squared(bz) < zeroThreshold<V>()) { return V(0); }
  bz = util::normalize(bz);

  // A random normalized vector that is just by random chance not collinear with bz
  // const Vector<V, DESCRIPTOR::d> rn{V(0.56270900), V(0.32704452), V(0.75921047)};
  const Vector<V, DESCRIPTOR::d> rn{V(0.5402151198529208), V(0.2765640977634561), V(0.7947829415070379)};

  // Normalization is necessary here because bz and rn are not perpendicular
  const Vector<V, DESCRIPTOR::d> by = util::normalize(crossProduct(bz, rn));
  const Vector<V, DESCRIPTOR::d> bx = crossProduct(by, bz);

  // Number of healthy interface neighbours
  std::uint8_t nbrInterfaces = std::uint8_t(0);

  // Compute z-offset of center point using PLIC cube
  const V centerOffset = plicCube<V>(getClampedEpsilon(cell), bz);

  // Number of neighbouring interface points, i.e., less than or equal to 26 (neighbours) - 1 (liquid) - 1 (gas).
  std::array<Vector<V, DESCRIPTOR::d>, 24> points;

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto nbrCell = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop));

    // Check if the neighbour is an interface cell or has a gas neighbour
    if (!isCellType(nbrCell, FreeSurface::Type::Interface) || !hasNeighbour(nbrCell, FreeSurface::Type::Gas)) {
      continue;
    }

    // Here it is assumed that neighbor normal vector is the same as center normal vector
    const auto direction = descriptors::c<DESCRIPTOR>(iPop);
    const Vector<V, DESCRIPTOR::d> ei = {V(direction[0]), V(direction[1]), V(direction[2])};
    const V offset = plicCube<V>(getClampedEpsilon(nbrCell), bz) - centerOffset;

    // Transform coordinate system into (x, y, f(x,y)) and apply PLIC offset
    points[nbrInterfaces++] = {util::dotProduct(ei, bx), util::dotProduct(ei, by), util::dotProduct(ei, bz) + offset};
  }

  // Prepare the linear equation system, to be initialized with zeros
  std::array<V, 25> M{};
  std::array<V, 5> solution{};
  std::array<V, 5> b{};

  for (std::uint8_t i = 0; i < nbrInterfaces; ++i) {
    const V x = points[i][0];
    const V y = points[i][1];
    const V z = points[i][2];
    const V x2 = util::sqr(x);
    const V y2 = util::sqr(y);
    const V x3 = x2 * x;
    const V y3 = y2 * y;

    M[0] += x2 * x2;
    M[1] += x2 * y2;
    M[2] += x3 * y;
    M[3] += x3;
    M[4] += x2 * y;
    M[6] += y2 * y2;
    M[7] += x * y3;
    M[8] += x * y2;
    M[9] += y3;
    M[12] += x2 * y2;
    M[13] += x2 * y;
    M[14] += x * y2;
    M[18] += x2;
    M[19] += x * y;
    M[24] += y2;

    b[0] += x2 * z;
    b[1] += y2 * z;
    b[2] += x * y * z;
    b[3] += x * z;
    b[4] += y * z;
  }

  // Fill the lower triangle of the symmetric matrix
  for (std::uint8_t i = 1; i < 5; ++i) {
    for (std::uint8_t j = 0; j < i; ++j) {
      M[i * 5 + j] = M[j * 5 + i];
    }
  }

  // Solve the linear system of equations, compile-time selection.
  if constexpr (solver == SolverType::rowPivotingLU) {
    if (nbrInterfaces >= 5) {
      rowPivotingLUSolver<V, 5>(M, solution, b, 5);
    } else {
      rowPivotingLUSolver<V, 5>(M, solution, b, util::min(5, nbrInterfaces));
    }
  }
  else if constexpr (solver == SolverType::completePivotingLU) {
    if (nbrInterfaces >= 5) {
      completePivotingLUSolver<V, 5>(M, solution, b, 5);
    } else {
      completePivotingLUSolver<V, 5>(M, solution, b, util::min(5, nbrInterfaces));
    }
  }
  else if constexpr (solver == SolverType::columnPivotingQR) {
    if (nbrInterfaces >= 5) {
      columnPivotingQRSolver<V, 5>(M, solution, b, 5);
    } else {
      columnPivotingQRSolver<V, 5>(M, solution, b, util::min(5, nbrInterfaces));
    }
  }
  else {
    // Fall back to default solver
    if (nbrInterfaces >= 5) {
      completePivotingLUSolver<V, 5>(M, solution, b, 5);
    } else {
      completePivotingLUSolver<V, 5>(M, solution, b, util::min(5, nbrInterfaces));
    }
  }

  const V A = solution[0], B = solution[1], C = solution[2], H = solution[3], I = solution[4];
  const V curvature = (A * (I * I + V(1)) + B * (H * H + V(1)) - C * H * I) * util::cube(V(1) / util::sqrt(H * H + I * I + V(1)));
  return util::max(V(-1), util::min(V(1), curvature));
}

template <typename CELL, typename V>
V computeCurvatureFDM3D(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  const auto normal = computeInterfaceNormal(cell);
  if (norm_squared(normal) < zeroThreshold<V>()) { return V(0); }

  const std::uint32_t gaussianMultipliers [DESCRIPTOR::q] =
  {
    std::uint32_t(8),
    std::uint32_t(4), std::uint32_t(4), std::uint32_t(4),
    std::uint32_t(2), std::uint32_t(2), std::uint32_t(2),
    std::uint32_t(2), std::uint32_t(2), std::uint32_t(2),
    std::uint32_t(1), std::uint32_t(1), std::uint32_t(1), std::uint32_t(1),
    std::uint32_t(4), std::uint32_t(4), std::uint32_t(4),
    std::uint32_t(2), std::uint32_t(2), std::uint32_t(2),
    std::uint32_t(2), std::uint32_t(2), std::uint32_t(2),
    std::uint32_t(1), std::uint32_t(1), std::uint32_t(1), std::uint32_t(1)
  };

  V curvature = V(0);
  V weightSum = V(0);

  for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
    auto direction = descriptors::c<DESCRIPTOR>(iPop);
    auto nbrCell = cell.neighbor(direction);

    // Check if the neighbour is an interface cell or has a gas neighbour
    if (!isCellType(nbrCell, FreeSurface::Type::Interface) || !hasNeighbour(nbrCell, FreeSurface::Type::Gas)) {
      continue;
    }

    const Vector<V, DESCRIPTOR::d> ei = {V(direction[0]), V(direction[1]), V(direction[2])};
    const auto nbrNormal = util::normalize(computeInterfaceNormal(nbrCell));
    const V weight = gaussianMultipliers[iPop];
    weightSum += weight;
    curvature += weight * util::dotProduct(ei, nbrNormal);
  }

  curvature /= weightSum;
  return util::max(V(-1), util::min(V(1), curvature));
}

template <typename CELL, typename V>
V computeCurvature(CELL& cell) {
  using DESCRIPTOR = typename CELL::descriptor_t;

  if constexpr (DESCRIPTOR::d == 2) {
    return computeCurvaturePLIC2D(cell);
  } else if (DESCRIPTOR::d == 3) {
    return computeCurvaturePLIC3D(cell);
  }

  return V(0);
}

} // namespace FreeSurface

} // namespace olb
