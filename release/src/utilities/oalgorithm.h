/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Jan E. Marquardt
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

#ifndef OLB_OALGORITHM_H
#define OLB_OALGORITHM_H

#include <algorithm>

#include "core/meta.h"
#include "core/expr.h"
#include "core/platform/platform.h"

namespace olb {

namespace util {

// Max
template <typename T>
any_platform constexpr T max(T a, meta::id_t<T> b) {
#ifdef  __CUDA_ARCH__
  return ::max(a, b);
#else
  return std::max(a, b);
#endif
}

template <typename T>
constexpr T max(std::initializer_list<T> ilist) {
  return std::max(ilist);
}

Expr max(Expr a, Expr b);

// Min
template <typename T>
any_platform constexpr T min(T a, meta::id_t<T> b) {
#ifdef  __CUDA_ARCH__
  return ::min(a, b);
#else
  return std::min(a, b);
#endif
}

template <typename T>
inline constexpr T min( std::initializer_list<T> ilist ) {
  return std::min(ilist);
}

Expr min(Expr a, Expr b);

} // namespace util

} // namespace olb

#endif
