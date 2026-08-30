/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef CORE_DATA_HH
#define CORE_DATA_HH

#include "data.h"

namespace olb {

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteData<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  std::size_t nBlock = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    nBlock += serializable->getNblock();
  }
  return nBlock;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteData<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  std::size_t size = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    size += serializable->getSerializableSize();
  }
  return size;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* ConcreteData<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (const auto& [_, serializable] : _serializationNames) {
    this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, *serializable, loadingMode);
  }

  return dataPtr;
}

}

#endif
