/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CPU_SIMD_PACK_H_
#define CPU_SIMD_PACK_H_

#ifdef __AVX512F__
#include "512.h"
#else
#include "256.h"
#endif

namespace olb {

namespace cpu {

namespace simd {

namespace concepts {

template <typename T>
concept convertible_to_pack = requires(T x) {
  Pack<float>{x};
  Pack<double>{x};
};

}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator+(S lhs, Pack<T> rhs)
{
  return Pack<T>(lhs) + rhs;
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator+(Pack<T> lhs, S rhs)
{
  return lhs + Pack<T>(rhs);
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator-(S lhs, Pack<T> rhs)
{
  return Pack<T>(lhs) - rhs;
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator-(Pack<T> lhs, S rhs)
{
  return lhs - Pack<T>(rhs);
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator*(Pack<T> lhs, S rhs)
{
  return lhs * Pack<T>(rhs);
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator*(S lhs, Pack<T> rhs)
{
  return Pack<T>(lhs) * rhs;
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator/(Pack<T> lhs, S rhs)
{
  return lhs / Pack<T>(rhs);
}

template <std::floating_point T, concepts::convertible_to_pack S>
Pack<T> operator/(S lhs, Pack<T> rhs)
{
  return Pack<T>(lhs) / rhs;
}

template <std::floating_point T>
Pack<T> sqrt(Pack<T> x)
{
  return x.sqrt();
}

}

}

namespace util {

template <std::floating_point T>
cpu::simd::Pack<T> sqrt(cpu::simd::Pack<T> value)
{
  return value.sqrt();
}

template <std::floating_point T>
cpu::simd::Pack<T> fabs(cpu::simd::Pack<T> value)
{
  return cpu::simd::fabs(value);
}

template <std::floating_point T>
cpu::simd::Pack<T> pow(cpu::simd::Pack<T> base, cpu::simd::Pack<T> exp)
{
  return cpu::simd::pow(base, exp);
}

template <std::floating_point T, cpu::simd::concepts::convertible_to_pack S>
cpu::simd::Pack<T> pow(cpu::simd::Pack<T> base, S exp)
{
  return cpu::simd::pow(base, cpu::simd::Pack<T>(exp));
}

template <std::floating_point T>
cpu::simd::Pack<T> min(cpu::simd::Pack<T> rhs, cpu::simd::Pack<T> lhs)
{
  return cpu::simd::min(rhs, lhs);
}

template <std::floating_point T>
cpu::simd::Pack<T> max(cpu::simd::Pack<T> rhs, cpu::simd::Pack<T> lhs)
{
  return cpu::simd::max(rhs, lhs);
}

}

}

#endif
