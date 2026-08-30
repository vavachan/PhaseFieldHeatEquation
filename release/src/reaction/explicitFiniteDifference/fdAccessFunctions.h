/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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
 * Access functions to access the correct element of the finite-difference external field.
 *  -- header file
 */
#ifndef FD_ACCESS_FUNCTIONS_H
#define FD_ACCESS_FUNCTIONS_H

namespace olb {

namespace fd {

/*
template<typename T, typename DESCRIPTOR, typename FIELD>
T* accessOld(Cell<T,DESCRIPTOR> cell, std::size_t iT)
{
  return &cell.template getFieldPointer<FIELD>()[iT % 2];
}

template<typename T, typename DESCRIPTOR, typename FIELD>
T* accessNew(Cell<T,DESCRIPTOR> cell, std::size_t iT)
{
  return &cell.template getFieldPointer<FIELD>()[(iT+1) % 2];
}
*/

template<typename T, typename FIELD, typename CELL>
T* accessOld(CELL cell, std::size_t iT) any_platform
{
  return &cell.template getFieldPointer<FIELD>()[iT % 2];
}

template<typename T, typename FIELD, typename CELL>
T* accessNew(CELL cell, std::size_t iT) any_platform
{
  return &cell.template getFieldPointer<FIELD>()[(iT+1) % 2];
}

// Convenience function to turn turn two indices from (iExt, iD) in [0, ..., EXTENT-1][0, ..., D-1], into i in [0, ..., D*EXTENT-1]
template <unsigned EXTENT>
constexpr unsigned getArrayPos(const unsigned iExt, const unsigned iD) any_platform
{
  return iExt*EXTENT + iD;
}

} // fd

} // olb

#endif
