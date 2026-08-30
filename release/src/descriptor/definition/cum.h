/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Louis Kronberg, Pavel Eichler, Stephan Simonis
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

#ifndef DESCRIPTOR_DEFINITION_CUM_H
#define DESCRIPTOR_DEFINITION_CUM_H

namespace olb {

namespace descriptors {

namespace tag {

struct CUM : public CATEGORY
           , public DESCRIPTOR_TAG {};

}

namespace cum_data {

using utilities::Fraction;

template <unsigned D, unsigned Q>
Fraction K[Q] = {};

template <unsigned D, unsigned Q>
int velocityIndices[Q][D] = {};

template <>
int velocityIndices<3, 27>[27][3] = {
    {10,  8, 26},
    {12, 22, 24},
    { 6,  3, 20},
    { 4,  2, 18},
    { 1,  0, 14},
    { 5, 15, 17},
    {11,  9, 25},
    { 7, 16, 19},
    {13, 21, 23},

    {10,  6, 12},
    { 4,  1,  5},
    {11,  7, 13},
    { 8,  3, 22},
    { 2,  0, 15},
    { 9, 16, 21},
    {26, 20, 24},
    {18, 14, 17},
    {25, 19, 23},

    {10,  4, 11},
    { 6,  1,  7},
    {12,  5, 13},
    { 8,  2,  9},
    { 3,  0, 16},
    {22, 15, 21},
    {26, 18, 25},
    {20, 14, 19},
    {24, 17, 23},
};

// these are the constant parameters K that are computed from the weights of the respective lattice.
template <>
Fraction K<3, 27>[27] = {
    {1 },
    0,      { 1, 3},
    0,      0,      0,      { 1, 3},
    0,      { 1, 9},
    { 1, 6},
    0,      { 1,18},
    { 2, 3},
    0,      { 2, 9},
    { 1, 6},
    0,      { 1,18},
    { 1,36},
    { 1, 9},
    { 1,36},
    { 1, 9},
    { 4, 9},
    { 1, 9},
    { 1,36},
    { 1, 9},
    { 1,36}
};

} // namespace cum_data

template <typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::CUM) any_platform
{
  return data::t<D, Q>[iPop].template as<T>();
}
template <typename T, unsigned D, unsigned Q>
constexpr T constantK(unsigned iPop) any_platform
{
  return cum_data::K<D, Q>[iPop].template as<T>();
}

} // namespace descriptors

}

#endif
