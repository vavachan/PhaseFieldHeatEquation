/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
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

/** \file
 * Implements Euclidean norm for arrays
 */

#ifndef NORM_H
#define NORM_H

namespace olb{

namespace util{

/// Squared Euclidean norm of an array
template<typename T>
T euklidN2(const T x[], int dim)
{
  T euklid = 0;
  for (int iDim=0; iDim<dim; iDim++) {
    euklid += x[iDim]*x[iDim];
  }
  return euklid;
}

/// Squared Euclidean distance between two arrays
template<typename T>
T euklidDistance2(const T x[], const T y[], int dim)
{
  T euklid = 0;
  for (int iDim=0; iDim<dim; iDim++) {
    euklid += (x[iDim]-y[iDim])*(x[iDim]-y[iDim]);
  }
  return euklid;
}


/// Euclidean norm of an array
template<typename T>
T euklidN(const T x[], int dim)
{
  return olb::util::sqrt(euklidN2(x,dim));
}

/// Euclidean distance between two arrays
template<typename T>
T euklidDistance(const T x[], const T y[], int dim)
{
  return olb::util::sqrt(euklidDistance2(x,y,dim));
}


/// Squared Euclidean norm for Container-type array
template<typename C>
auto euklidN2(const C& x)
{
  return euklidN2(x.data(), x.size());
}

/// Squared Euclidean distance for Container-type array
template<typename C>
auto euklidDistance2(const C& x, const C& y)
{
  OLB_ASSERT((x.size() == y.size()), "Arrays must have same dimension");
  return euklidDistance2(x.data(), y.data(), x.size());
}

/// Euclidean norm for Container-type array
template<typename C>
auto euklidN(const C& x)
{
  return euklidN(x.data(), x.size());
}

/// Euclidean distance for Container-type array
template<typename C>
auto euklidDistance(const C& x, const C& y)
{
  OLB_ASSERT((x.size() == y.size()), "Arrays must have same dimension");
  return euklidDistance(x.data(), y.data(), x.size());
}

}

}

#endif