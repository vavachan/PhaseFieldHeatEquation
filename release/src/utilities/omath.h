
/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Yuji Shimojima, 2020 Jan E. Marquardt
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

#ifndef OLB_OMATH_H
#define OLB_OMATH_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>

namespace olb {

namespace util {

// Pow

template <typename T, typename S>
any_platform inline auto pow(T x, S y)
    -> std::enable_if_t<std::is_floating_point_v<T> &&
                            std::is_floating_point_v<S>,
                        decltype(std::pow(x, y))>
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float> && std::is_same_v<S, float>) {
    return ::powf(x, y);
  }
  else if constexpr (std::is_same_v<T, double> && std::is_same_v<S, double>) {
    return ::pow(x, y);
  }
  else {
    return ::pow(static_cast<double>(x), static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  // __CUDA_ARCH__
  return std::pow(x, y);
#endif // __CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline auto pow(T x, S y)
    -> std::enable_if_t<std::is_integral_v<S> && std::is_floating_point_v<T>,
                        decltype(std::pow(x, y))>
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::powf(x, y);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::pow(x, y);
  }
  else {
    return ::pow(static_cast<double>(x), y);
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::pow(x, y);
#endif //__CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline auto pow(T x, S y)
    -> std::enable_if_t<std::is_integral_v<T> && std::is_floating_point_v<S>,
                        decltype(std::pow(x, y))>
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<S, float>) {
    return ::powf(x, y);
  }
  else if constexpr (std::is_same_v<S, double>) {
    return ::pow(x, y);
  }
  else {
    return ::pow(x, static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::pow(x, y);
#endif //__CUDA_ARCH__
}
template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<T> && std::is_integral_v<S>, double>
pow(T x, S y)
{
#ifdef __CUDA_ARCH__
  return ::pow(x, y);
#else  //__CUDA_ARCH__
  return std::pow(x, y);
#endif //__CUDA_ARCH__
}

any_platform inline float powf(float base, float exp)
{
#ifdef __CUDA_ARCH__
  return ::powf(base, exp);
#else  //__CUDA_ARCH__
  return std::pow(base, exp);
#endif //__CUDA_ARCH__
}

any_platform inline float powf(int base, int exp)
{
#ifdef __CUDA_ARCH__
  return ::powf(base, exp);
#else  //__CUDA_ARCH__
  return std::pow(base, exp);
#endif //__CUDA_ARCH__
}

inline long double powl(long double base, long double exp)
{
  return pow(base, exp);
  //cuda dose not support this function.
}

//fmod
template <typename T, typename S>
any_platform inline auto fmod(T x, S y)
    -> std::enable_if_t<std::is_floating_point_v<T> &&
                            std::is_floating_point_v<S>,
                        decltype(std::fmod(x, y))>
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float> && std::is_same_v<S, float>) {
    return ::fmodf(x, y);
  }
  else if constexpr (std::is_same_v<T, double> && std::is_same_v<S, double>) {
    return ::fmod(x, y);
  }
  else {
    return ::fmod(static_cast<double>(x), static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::fmod(x, y);
#endif //__CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<S> && std::is_floating_point_v<T>, T>
fmod(T x, S y)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::fmodf(x, y);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::fmod(x, y);
  }
  else {
    return ::fmod(static_cast<double>(x), y);
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::fmod(x, y);
#endif //__CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<T> && std::is_floating_point_v<S>, S>
fmod(T x, S y)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<S, float>) {
    return ::fmodf(x, y);
  }
  else if constexpr (std::is_same_v<S, double>) {
    return ::fmod(x, y);
  }
  else {
    return ::fmod(x, static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::fmod(x, y);
#endif //__CUDA_ARCH__
}

any_platform inline float fmodf(float x, float y)
{
#ifdef __CUDA_ARCH__
  return ::fmodf(x, y);
#else  //__CUDA_ARCH__
  return std::fmod(x, y);
#endif //__CUDA_ARCH__
}

any_platform inline float fmodf(int x, int y)
{
#ifdef __CUDA_ARCH__
  return ::fmodf(x, y);
#else  //__CUDA_ARCH__
  return std::fmod(x, y);
#endif //__CUDA_ARCH__
}

inline long double fmodl(long double x, long double y)
{
  return fmod(x, y);
  //cuda dose not support this function.
}

// Exp
template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> exp(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::expf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::exp(x);
  }
  else {
    return ::exp(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::exp(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> exp(T arg)
{
#ifdef __CUDA_ARCH__
  return ::exp(arg);
#else  //__CUDA_ARCH__
  return std::exp(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float expf(float arg)
{
#ifdef __CUDA_ARCH__
  return ::expf(arg);
#else  //__CUDA_ARCH__
  return std::exp(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float expf(int arg)
{
#ifdef __CUDA_ARCH__
  return ::expf(arg);
#else  //__CUDA_ARCH__
  return std::exp(arg);
#endif //__CUDA_ARCH__
}

inline long double expl(long double arg)
{
  return exp(arg);
  //cuda dose not support this function.
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> log(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::logf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::log(x);
  }
  else {
    return ::log(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::log(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> log(T arg)
{
#ifdef __CUDA_ARCH__
  return ::log(arg);
#else  //__CUDA_ARCH__
  return std::log(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float logf(float arg)
{
#ifdef __CUDA_ARCH__
  return ::logf(arg);
#else  //__CUDA_ARCH__
  return std::log(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float logf(int arg)
{
#ifdef __CUDA_ARCH__
  return ::logf(arg);
#else  //__CUDA_ARCH__
  return std::log(arg);
#endif //__CUDA_ARCH__
}

inline long double logl(long double arg)
{
  return log(arg);
  //cuda dose not support this function.
}

// Log10

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> log10(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::log10f(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::log10(x);
  }
  else {
    return ::log10(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::log10(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> log10(T arg)
{
#ifdef __CUDA_ARCH__
  return ::log10(arg);
#else  //__CUDA_ARCH__
  return std::log10(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log10f(float arg)
{
#ifdef __CUDA_ARCH__
  return ::log10f(arg);
#else  //__CUDA_ARCH__
  return std::log10(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log10f(int arg)
{
#ifdef __CUDA_ARCH__
  return ::log10f(arg);
#else  //__CUDA_ARCH__
  return std::log10(arg);
#endif //__CUDA_ARCH__
}

inline long double log10l(long double arg)
{
  return log10(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Log2

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> log2(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::log2f(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::log2(x);
  }
  else {
    return ::log2(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::log2(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> log2(T arg)
{
#ifdef __CUDA_ARCH__
  return ::log2(arg);
#else  //__CUDA_ARCH__
  return std::log2(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log2f(float arg)
{
#ifdef __CUDA_ARCH__
  return ::log2f(arg);
#else  //__CUDA_ARCH__
  return std::log2(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log2f(int arg)
{
#ifdef __CUDA_ARCH__
  return ::log2f(arg);
#else  //__CUDA_ARCH__
  return std::log2(arg);
#endif //__CUDA_ARCH__
}

inline long double log2l(long double arg)
{
  return log2(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Log1p
template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> log1p(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::log1pf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::log1p(x);
  }
  else {
    return ::log1p(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::log1p(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> log1p(T arg)
{
#ifdef __CUDA_ARCH__
  return ::log1p(arg);
#else  //__CUDA_ARCH__
  return std::log1p(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log1pf(float arg)
{
#ifdef __CUDA_ARCH__
  return ::log1pf(arg);
#else  //__CUDA_ARCH__
  return std::log1p(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float log1pf(int arg)
{
#ifdef __CUDA_ARCH__
  return ::log1pf(arg);
#else  //__CUDA_ARCH__
  return std::log1p(arg);
#endif //__CUDA_ARCH__
}

inline long double log1pl(long double arg)
{
  return log1p(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Sqrt
template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> sqrt(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::sqrtf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::sqrt(x);
  }
  else {
    return ::sqrt(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::sqrt(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> sqrt(T arg)
{
#ifdef __CUDA_ARCH__
  return ::sqrt(arg);
#else  //__CUDA_ARCH__
  return std::sqrt(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float sqrtf(float arg)
{
#ifdef __CUDA_ARCH__
  return ::sqrtf(arg);
#else  //__CUDA_ARCH__
  return std::sqrt(arg);
#endif //__CUDA_ARCH__
}

any_platform inline float sqrtf(int arg)
{
#ifdef __CUDA_ARCH__
  return ::sqrtf(arg);
#else  //__CUDA_ARCH__
  return std::sqrt(arg);
#endif //__CUDA_ARCH__
}

inline long double sqrtl(long double arg)
{
  return sqrt(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Sin

any_platform inline float sinf(float x)
{
#ifdef __CUDA_ARCH__
  return ::sinf(x);
#else  //__CUDA_ARCH__
  return std::sin(x);
#endif //__CUDA_ARCH__
}

any_platform inline float sinf(int x)
{
#ifdef __CUDA_ARCH__
  return ::sinf(x);
#else  //__CUDA_ARCH__
  return std::sin(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> sin(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::sinf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::sin(x);
  }
  else {
    return ::sin(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::sin(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> sin(T x)
{
#ifdef __CUDA_ARCH__
  return ::sin(x);
#else  //__CUDA_ARCH__
  return std::sin(x);
#endif //__CUDA_ARCH__
}

inline long double sinl(long double arg)
{
  return sin(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Sinh

any_platform inline float sinhf(float x)
{
#ifdef __CUDA_ARCH__
  return ::sinhf(x);
#else  //__CUDA_ARCH__
  return std::sinh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float sinhf(int x)
{
#ifdef __CUDA_ARCH__
  return ::sinhf(x);
#else  //__CUDA_ARCH__
  return std::sinh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> sinh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::sinhf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::sinh(x);
  }
  else {
    return ::sinh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::sinh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> sinh(T x)
{
#ifdef __CUDA_ARCH__
  return ::sinh(x);
#else  //__CUDA_ARCH__
  return std::sinh(x);
#endif //__CUDA_ARCH__
}

inline long double sinhl(long double arg)
{

  return sinh(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}
// Cos

any_platform inline float cosf(float x)
{
#ifdef __CUDA_ARCH__
  return ::cosf(x);
#else  //__CUDA_ARCH__
  return std::cos(x);
#endif //__CUDA_ARCH__
}

any_platform inline float cosf(int x)
{
#ifdef __CUDA_ARCH__
  return ::cosf(x);
#else  //__CUDA_ARCH__
  return std::cos(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> cos(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::cosf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::cos(x);
  }
  else {
    return ::cos(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::cos(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> cos(T x)
{
#ifdef __CUDA_ARCH__
  return ::cos(x);
#else  //__CUDA_ARCH__
  return std::cos(x);
#endif //__CUDA_ARCH__
}

inline long double cosl(long double arg)
{
  return cos(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Cosh

any_platform inline float coshf(float x)
{
#ifdef __CUDA_ARCH__
  return ::coshf(x);
#else  //__CUDA_ARCH__
  return std::cosh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float coshf(int x)
{
#ifdef __CUDA_ARCH__
  return ::coshf(x);
#else  //__CUDA_ARCH__
  return std::cosh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> cosh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::coshf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::cosh(x);
  }
  else {
    return ::cosh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::cosh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> cosh(T x)
{
#ifdef __CUDA_ARCH__
  return ::cosh(x);
#else  //__CUDA_ARCH__
  return std::cosh(x);
#endif //__CUDA_ARCH__
}

inline long double coshl(long double arg)
{
  return cosh(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Tan
any_platform inline float tanf(float x)
{
#ifdef __CUDA_ARCH__
  return ::tanf(x);
#else  //__CUDA_ARCH__
  return std::tan(x);
#endif //__CUDA_ARCH__
}

any_platform inline float tanf(int x)
{
#ifdef __CUDA_ARCH__
  return ::tanf(x);
#else  //__CUDA_ARCH__
  return std::tan(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> tan(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::tanf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::tan(x);
  }
  else {
    return ::tan(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::tan(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> tan(T x)
{
#ifdef __CUDA_ARCH__
  return ::tan(x);
#else  //__CUDA_ARCH__
  return std::tan(x);
#endif //__CUDA_ARCH__
}

inline long double tanl(long double arg)
{
  return tan(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

// Tanh
any_platform inline float tanhf(float x)
{
#ifdef __CUDA_ARCH__
  return ::tanhf(x);
#else  //__CUDA_ARCH__
  return std::tanh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float tanhf(int x)
{
#ifdef __CUDA_ARCH__
  return ::tanhf(x);
#else  //__CUDA_ARCH__
  return std::tanh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> tanh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::tanhf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::tanh(x);
  }
  else {
    return ::tanh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::tanh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> tanh(T x)
{
#ifdef __CUDA_ARCH__
  return ::tanh(x);
#else  //__CUDA_ARCH__
  return std::tanh(x);
#endif //__CUDA_ARCH__
}

inline long double tanhl(long double arg)
{
  return tanh(arg);
  //libcu++ of cuda dose not support long double. You will get warning in device code if you use.
}

any_platform inline float asinf(float x)
{
#ifdef __CUDA_ARCH__
  return ::asinf(x);
#else  //__CUDA_ARCH__
  return std::asin(x);
#endif //__CUDA_ARCH__
}

any_platform inline float asinf(int x)
{
#ifdef __CUDA_ARCH__
  return ::asinf(x);
#else  //__CUDA_ARCH__
  return std::asin(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> asin(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::asinf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::asin(x);
  }
  else {
    return ::asin(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::asin(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> asin(T x)
{
#ifdef __CUDA_ARCH__
  return ::asin(x);
#else  //__CUDA_ARCH__
  return std::asin(x);
#endif //__CUDA_ARCH__
}

inline long double asinl(long double arg)
{
  return asin(arg);
  //cuda dose not support this function.
}

// Asinh
any_platform inline float asinhf(float x)
{
#ifdef __CUDA_ARCH__
  return ::asinhf(x);
#else  //__CUDA_ARCH__
  return std::asinh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float asinhf(int x)
{
#ifdef __CUDA_ARCH__
  return ::asinhf(x);
#else  //__CUDA_ARCH__
  return std::asinh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> asinh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::asinhf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::asinh(x);
  }
  else {
    return ::asinh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::asinh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> asinh(T x)
{
#ifdef __CUDA_ARCH__
  return ::asinh(x);
#else  //__CUDA_ARCH__
  return std::asinh(x);
#endif //__CUDA_ARCH__
}

inline long double asinhl(long double arg)
{
  return std::asinhl(arg);
  //cuda dose not support this function.
}

// Acos
any_platform inline float acosf(float x)
{
#ifdef __CUDA_ARCH__
  return ::acosf(x);
#else  //__CUDA_ARCH__
  return std::acos(x);
#endif //__CUDA_ARCH__
}

any_platform inline float acosf(int x)
{
#ifdef __CUDA_ARCH__
  return ::acosf(x);
#else  //__CUDA_ARCH__
  return std::acos(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> acos(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::acosf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::acos(x);
  }
  else {
    return ::acos(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::acos(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> acos(T x)
{
#ifdef __CUDA_ARCH__
  return ::acos(x);
#else  //__CUDA_ARCH__
  return std::acos(x);
#endif //__CUDA_ARCH__
}

inline long double acosl(long double arg)
{
  return acos(arg);
  //cuda dose not support this function.
}

// Acosh
any_platform inline float acoshf(float x)
{
#ifdef __CUDA_ARCH__
  return ::acoshf(x);
#else  //__CUDA_ARCH__
  return std::acosh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float acoshf(int x)
{
#ifdef __CUDA_ARCH__
  return ::acoshf(x);
#else  //__CUDA_ARCH__
  return std::acosh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> acosh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::acoshf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::acosh(x);
  }
  else {
    return ::acosh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::acosh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> acosh(T x)
{
#ifdef __CUDA_ARCH__
  return ::acosh(x);
#else  //__CUDA_ARCH__
  return std::acosh(x);
#endif //__CUDA_ARCH__
}

inline long double acoshl(long double arg)
{
  return std::acoshl(arg);
  //cuda dose not support this function.
}

// Atan
any_platform inline float atanf(float x)
{
#ifdef __CUDA_ARCH__
  return ::atanf(x);
#else  //__CUDA_ARCH__
  return std::atan(x);
#endif //__CUDA_ARCH__
}

any_platform inline float atanf(int x)
{
#ifdef __CUDA_ARCH__
  return ::atanf(x);
#else  //__CUDA_ARCH__
  return std::atan(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> atan(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::atanf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::atan(x);
  }
  else {
    return ::atan(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::atan(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> atan(T x)
{
#ifdef __CUDA_ARCH__
  return ::atan(x);
#else  //__CUDA_ARCH__
  return std::atan(x);
#endif //__CUDA_ARCH__
}

inline long double atanl(long double arg)
{
  return atan(arg);
  //cuda dose not support this function.
}

// Atan2
template <typename T, typename S>
any_platform inline auto atan2(T x, S y)
    -> std::enable_if_t<std::is_floating_point_v<T> &&
                            std::is_floating_point_v<S>,
                        decltype(std::atan2(x, y))>
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float> && std::is_same_v<S, float>) {
    return ::atan2f(x, y);
  }
  else if constexpr (std::is_same_v<T, double> && std::is_same_v<S, double>) {
    return ::atan2(x, y);
  }
  else {
    return ::atan2(static_cast<double>(x), static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  // __CUDA_ARCH__
  return std::atan2(x, y);
#endif // __CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<S> && std::is_floating_point_v<T>, T>
atan2(T x, S y)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::atan2f(x, y);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::atan2(x, y);
  }
  else {
    return ::atan2(static_cast<double>(x), y);
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::atan2(x, y);
#endif //__CUDA_ARCH__
}

template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<T> && std::is_floating_point_v<S>, S>
atan2(T x, S y)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<S, float>) {
    return ::atan2f(x, y);
  }
  else if constexpr (std::is_same_v<S, double>) {
    return ::atan2(x, y);
  }
  else {
    return ::atan2(x, static_cast<double>(y));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::atan2(x, y);
#endif //__CUDA_ARCH__
}
template <typename T, typename S>
any_platform inline std::enable_if_t<
    std::is_integral_v<T> && std::is_integral_v<S>, double>
atan2(T x, S y)
{
#ifdef __CUDA_ARCH__
  return ::atan2(x, y);
#else  //__CUDA_ARCH__
  return std::atan2(x, y);
#endif //__CUDA_ARCH__
}

any_platform inline float atan2f(float base, float exp)
{
#ifdef __CUDA_ARCH__
  return ::atan2f(base, exp);
#else  //__CUDA_ARCH__
  return std::atan2(base, exp);
#endif //__CUDA_ARCH__
}

any_platform inline float atan2f(int base, int exp)
{
#ifdef __CUDA_ARCH__
  return ::atan2f(base, exp);
#else  //__CUDA_ARCH__
  return std::atan2(base, exp);
#endif //__CUDA_ARCH__
}

// Atanh
any_platform inline float atanhf(float x)
{
#ifdef __CUDA_ARCH__
  return ::atanhf(x);
#else  //__CUDA_ARCH__
  return std::atanh(x);
#endif //__CUDA_ARCH__
}

any_platform inline float atanhf(int x)
{
#ifdef __CUDA_ARCH__
  return ::atanhf(x);
#else  //__CUDA_ARCH__
  return std::atanh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> atanh(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::atanhf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::atanh(x);
  }
  else {
    return ::atanh(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::atanh(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> atanh(T x)
{
#ifdef __CUDA_ARCH__
  return ::atanh(x);
#else  //__CUDA_ARCH__
  return std::atanh(x);
#endif //__CUDA_ARCH__
}

inline long double atanhl(long double arg)
{
  return atanh(arg);
  //cuda dose not support this function.
}

// Rounding
// floor
any_platform inline float floorf(float x)
{
#ifdef __CUDA_ARCH__
  return ::floorf(x);
#else  //__CUDA_ARCH__
  return std::floor(x);
#endif //__CUDA_ARCH__
}

any_platform inline float floorf(int x)
{
#ifdef __CUDA_ARCH__
  return ::floorf(x);
#else  //__CUDA_ARCH__
  return std::floor(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> floor(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::floorf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::floor(x);
  }
  else {
    return ::floor(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::floor(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> floor(T x)
{
#ifdef __CUDA_ARCH__
  return ::floor(x);
#else  //__CUDA_ARCH__
  return std::floor(x);
#endif //__CUDA_ARCH__
}

inline long double floorl(long double arg)
{
  return floor(arg);
  //cuda dose not support this function.
}

// ceil
any_platform inline float ceilf(float x)
{
#ifdef __CUDA_ARCH__
  return ::ceilf(x);
#else  //__CUDA_ARCH__
  return std::ceil(x);
#endif //__CUDA_ARCH__
}

any_platform inline float ceilf(int x)
{
#ifdef __CUDA_ARCH__
  return ::ceilf(x);
#else  //__CUDA_ARCH__
  return std::ceil(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> ceil(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::ceilf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::ceil(x);
  }
  else {
    return ::ceil(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::ceil(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> ceil(T x)
{
#ifdef __CUDA_ARCH__
  return ::ceil(x);
#else  //__CUDA_ARCH__
  return std::ceil(x);
#endif //__CUDA_ARCH__
}

inline long double ceill(long double arg)
{
  return ceil(arg);
  //cuda dose not support this function.
}

// trunc
any_platform inline float truncf(float x)
{
#ifdef __CUDA_ARCH__
  return ::truncf(x);
#else  //__CUDA_ARCH__
  return std::trunc(x);
#endif //__CUDA_ARCH__
}

any_platform inline float truncf(int x)
{
#ifdef __CUDA_ARCH__
  return ::truncf(x);
#else  //__CUDA_ARCH__
  return std::trunc(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> trunc(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::truncf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::trunc(x);
  }
  else {
    return ::trunc(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::trunc(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> trunc(T x)
{
#ifdef __CUDA_ARCH__
  return ::trunc(x);
#else  //__CUDA_ARCH__
  return std::trunc(x);
#endif //__CUDA_ARCH__
}

inline long double truncl(long double arg)
{
  return std::truncl(arg);
  //cuda dose not support this function.
}

// round
any_platform inline float roundf(float x)
{
#ifdef __CUDA_ARCH__
  return ::roundf(x);
#else  //__CUDA_ARCH__
  return std::round(x);
#endif //__CUDA_ARCH__
}

any_platform inline float roundf(int x)
{
#ifdef __CUDA_ARCH__
  return ::roundf(x);
#else  //__CUDA_ARCH__
  return std::round(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> round(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::roundf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::round(x);
  }
  else {
    return ::round(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::round(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> round(T x)
{
#ifdef __CUDA_ARCH__
  return ::round(x);
#else  //__CUDA_ARCH__
  return std::round(x);
#endif //__CUDA_ARCH__
}

inline long double roundl(long double arg)
{
  return std::roundl(arg);
  //cuda dose not support this function.
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, long> lround(T arg)
{
  return std::lround(arg);
}

inline long lroundf(float arg) { return std::lroundf(arg); }

inline long lroundl(long double arg) { return std::lroundl(arg); }

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, long long> llround(T arg)
{
  return std::llround(arg);
}

inline long long llroundf(float arg) { return std::llroundf(arg); }

inline long long llroundl(long double arg) { return std::llroundl(arg); }

// Absolute
any_platform inline float fabsf(float x)
{
#ifdef __CUDA_ARCH__
  return ::fabsf(x);
#else  //__CUDA_ARCH__
  return std::abs(x);
#endif //__CUDA_ARCH__
}

any_platform inline float fabsf(int x)
{
#ifdef __CUDA_ARCH__
  return ::fabsf(x);
#else  //__CUDA_ARCH__
  return std::abs(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> fabs(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::fabsf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::fabs(x);
  }
  else {
    return ::fabs(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::fabs(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, double> fabs(T x)
{
#ifdef __CUDA_ARCH__
  return ::fabs(x);
#else  //__CUDA_ARCH__
  return std::fabs(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_floating_point_v<T>, T> abs(T x)
{
#ifdef __CUDA_ARCH__
  if constexpr (std::is_same_v<T, float>) {
    return ::fabsf(x);
  }
  else if constexpr (std::is_same_v<T, double>) {
    return ::fabs(x);
  }
  else {
    return ::fabs(static_cast<double>(x));
    //libcu++ of cuda supports double and float. You will get warning in device code if you use in e.g. long double.
  }
#else  //__CUDA_ARCH__
  return std::abs(x);
#endif //__CUDA_ARCH__
}

template <typename T>
any_platform inline std::enable_if_t<std::is_integral_v<T>, T> abs(T x)
{
#ifdef __CUDA_ARCH__
  return ::abs(x);
#else  //__CUDA_ARCH__
  return std::abs(x);
#endif //__CUDA_ARCH__
}

inline long double fabs(long double arg) { return std::abs(arg); }

template <std::integral T>
any_platform bool closeToZero(T x)
{
  return fabs(x) <= std::numeric_limits<T>::epsilon();
}

//isnan
template <typename T>
any_platform inline bool isnan(T arg)
{
#ifdef __CUDA_ARCH__
  return ::isnan(arg);
#else
  return std::isnan(arg);
#endif
}
} // namespace util

} // namespace olb

#endif
