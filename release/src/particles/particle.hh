/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_HH
#define PARTICLE_HH


#include "particles/functions/particleIoFunctions.h"
#include "particles/functions/particleDynamicsFunctions.h"

namespace olb {

namespace particles{


template <typename T, typename PARTICLETYPE>
Particle<T,PARTICLETYPE>::Particle(
  DynamicFieldGroupsD<T, typename PARTICLETYPE::fields_t>& dynamicFieldGroupsD,
  std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>>& dynamicsVector,
  std::size_t iParticle )
  : _dynamicFieldGroupsD( dynamicFieldGroupsD ),
    _dynamicsVector( dynamicsVector ),
    _iParticle( iParticle)
{}

template <typename T, typename PARTICLETYPE>
void Particle<T,PARTICLETYPE>::init()
{
  particles::dynamics::initializeParticle<T, PARTICLETYPE>( _dynamicFieldGroupsD, _iParticle);
}

template <typename T, typename PARTICLETYPE>
template<bool multiOutput>
void Particle<T,PARTICLETYPE>::print(std::size_t iParticle)
{
  particles::io::printGenericParticleInfo<T, PARTICLETYPE,multiOutput>( _dynamicFieldGroupsD, iParticle);
}

template <typename T, typename PARTICLETYPE>
template<bool multiOutput>
void Particle<T,PARTICLETYPE>::print()
{
  print<multiOutput>(_iParticle);
}

template <typename T, typename PARTICLETYPE>
std::size_t Particle<T,PARTICLETYPE>::getId() const
{
  return _iParticle;
}

template <typename T, typename PARTICLETYPE>
void Particle<T,PARTICLETYPE>::setId(std::size_t iParticle)
{
  _iParticle = iParticle;
}

template <typename T, typename PARTICLETYPE>
void Particle<T,PARTICLETYPE>::advanceId()
{
  ++_iParticle;
}

template<typename T, typename PARTICLETYPE>
void Particle<T,PARTICLETYPE>::addDynamics (
  std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr )
{
  _dynamicsVector.push_back( dynamicsSPtr );
}

template<typename T, typename PARTICLETYPE>
template <typename DYNAMICS, typename ...Args>
void Particle<T,PARTICLETYPE>::defineDynamics (Args&& ...args){
  _dynamicsVector.push_back( std::make_shared<DYNAMICS>(std::forward<Args>(args)...) );
}

template<typename T, typename PARTICLETYPE>
template<bool boundsCheck>
dynamics::ParticleDynamics<T,PARTICLETYPE>* Particle<T,PARTICLETYPE>::getDynamics(
  unsigned iDyn)
{
  if constexpr(boundsCheck){
    return _dynamicsVector.at(iDyn).get();
  } else {
    return _dynamicsVector[iDyn].get();
  }
}

template<typename T, typename PARTICLETYPE>
void Particle<T,PARTICLETYPE>::process(T timeStepSize, unsigned short iDyn)
{
  if(access::isMotionComputationEnabled(*this)) {
    getDynamics(iDyn)->process(*this, timeStepSize);
  }
}

////////// Get and Set functions //TODO: implement recursively?


template <typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
std::enable_if_t<(
  PARTICLETYPE::template size<
    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
    >() > 1),
    FieldD<T,PARTICLETYPE,
    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>>
>
Particle<T,PARTICLETYPE>::getField()
{
  using GROUP_EVAL = typename PARTICLETYPE
                     ::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL
                     ::template derivedField<FIELD>;
  return _dynamicFieldGroupsD.template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>()
         .getRowPointer(_iParticle);
}


template <typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
std::enable_if_t<(
  PARTICLETYPE::template size<
    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
    >() == 1),
    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>::template value_type<T>>
Particle<T,PARTICLETYPE>::getField()
{
  using GROUP_EVAL = typename PARTICLETYPE
                     ::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL
                     ::template derivedField<FIELD>;
  return _dynamicFieldGroupsD.template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>()
         .getRowPointer(_iParticle)[0];
}

template <typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
std::enable_if_t<(PARTICLETYPE::template size<
                    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
                    >() > 1), void>
Particle<T,PARTICLETYPE>::setField(
  const FieldD<T,PARTICLETYPE,
  typename PARTICLETYPE::template derivedField<GROUP>
  ::template derivedField<FIELD>>& v )
{
  using GROUP_EVAL = typename PARTICLETYPE
                     ::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL
                     ::template derivedField<FIELD>;
  //Analogous to BlockStaticFieldsD
  auto& field = _dynamicFieldGroupsD.template get<GROUP_EVAL>()
                .template get<FIELD_EVAL>();
  for (unsigned iDim=0; iDim < PARTICLETYPE::template size<FIELD_EVAL>(); ++iDim) {
    field[iDim][_iParticle] = v[iDim];
  }
}

template <typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
std::enable_if_t<(PARTICLETYPE::template size<
                    typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
                    >() == 1), void>
Particle<T,PARTICLETYPE>::setField(
  typename PARTICLETYPE::template derivedField<GROUP>
                     ::template derivedField<FIELD>
                     ::template value_type<T> value)
{
  using GROUP_EVAL = typename PARTICLETYPE
                     ::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL
                     ::template derivedField<FIELD>;
  _dynamicFieldGroupsD.template get<GROUP_EVAL>()
  .template get<FIELD_EVAL>()
  .getRowPointer(_iParticle)[0] = value;
}

} //namespace particles

} //namespace olb

#endif
