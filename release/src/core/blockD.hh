/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
 *
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

#ifndef BLOCK_D_HH
#define BLOCK_D_HH

#include "blockD.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
BlockD<T,DESCRIPTOR>::BlockD(Vector<int,DESCRIPTOR::d> size, int padding, Platform platform)
  : BlockStructure<DESCRIPTOR>(size, padding),
    _platform(platform)
{ }

template<typename T, typename DESCRIPTOR>
BlockD<T,DESCRIPTOR>::~BlockD()
{ }

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::ConcreteBlockD(Vector<int,DESCRIPTOR::d> size, int padding)
  : BlockD<T,DESCRIPTOR>(size, padding, PLATFORM),
    _data(this->getNcells()),
    _descriptorFields()
{
  DESCRIPTOR::fields_t::for_each([&](auto id) {
    using field = typename decltype(id)::type;
    auto& fieldArray = _data.template get<Array<field>>();
    _descriptorFields.template set<field>(&fieldArray);
    _communicatables[typeid(field)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename field::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<field>()>
    >(fieldArray));
  });
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::resize(Vector<int,DESCRIPTOR::d> size)
{
  static_cast<BlockStructureD<DESCRIPTOR::d>*>(this)->resize(size);
  _data.resize(this->getNcells());
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
bool ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::hasData() const
{
  return _data.template provides<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
const auto& ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getData() const
{
  OLB_ASSERT(_data.template provides<FIELD_TYPE>(),
             "FIELD_TYPE must be allocated to be accessed");
  return _data.template get<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
auto& ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getData()
{
  if (_data.template provides<FIELD_TYPE>()) {
    return _data.template get<FIELD_TYPE>();
  } else {
    // TODO: Implement more generic approach to constructing arbitrary data from specific args
    auto& data = _data.template allocate<FIELD_TYPE>(this->getNcells());
    // Manage serializables and communicatables for array fields
    using concrete_data_t = typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>;
    if constexpr (std::is_base_of_v<ColumnVectorBase,concrete_data_t>) {
      using field_t = typename concrete_data_t::field_t;
      if constexpr (field_t::isSerializable()) {
        _data.template setSerialization<FIELD_TYPE>(true);
      }
      _communicatables[typeid(field_t)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
        ColumnVector<typename ImplementationOf<typename field_t::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<field_t>()>
      >(data));
    }
    return data;
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
auto& ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getField(FIELD)
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
const auto& ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getField(FIELD) const
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<const FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  return _data.getNblock();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  return _data.getSerializableSize();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;
  this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data, loadingMode);
  return dataPtr;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockD<T,DESCRIPTOR,PLATFORM>::postLoad()
{
  if constexpr (DESCRIPTOR::template provides<descriptors::POPULATION>()) {
    auto& population = getField<descriptors::POPULATION>();
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      population[iPop].postLoad();
    }
  }
}

}

#endif
