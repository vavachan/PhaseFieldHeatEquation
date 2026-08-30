/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen
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


#ifndef PARTICLE_SYSTEM_HH
#define PARTICLE_SYSTEM_HH


namespace olb {

namespace particles{

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem()
  : _fieldGroupsContainer( Container<T,PARTICLETYPE,DATA>(0) ),
    _serialSize(DATA(0).getSerializableSize())
{}

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem( std::size_t count )
  : _fieldGroupsContainer( Container<T,PARTICLETYPE,DATA>(count) ),
    _serialSize(DATA(0).getSerializableSize())
{}

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem( std::size_t count,
    typename associatedTypes::template decompose_into<std::tuple> associatedData )
  : _associatedData( associatedData ),
    _fieldGroupsContainer( Container<T,PARTICLETYPE,DATA>(count) ),
    _serialSize(DATA(0).getSerializableSize())
{}

template<typename T, typename PARTICLETYPE>
template<typename PCONDITION>
constexpr std::size_t ParticleSystem<T, PARTICLETYPE>::size()
{
  //Either directly return size of container
  if constexpr(std::is_same_v<PCONDITION,conditions::all_particles>){
    return _fieldGroupsContainer.size();
  //Or iterate over container evaluating condition
  } else {
    std::size_t count=0;
    for (std::size_t iP=0; iP<_fieldGroupsContainer.size(); ++iP) {
      auto particle = get(iP);
      doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
        [&](Particle<T,PARTICLETYPE> particle){ ++count; });
    }
    return count;
  }
}

template<typename T, typename PARTICLETYPE>
template<typename PCONDITION>
constexpr std::size_t ParticleSystem<T, PARTICLETYPE>::size(int globiC)
{
  //Either directly return size of container
  if constexpr(std::is_same_v<PCONDITION,conditions::all_particles>){
    return _fieldGroupsContainer.size();
  //Or iterate over container evaluating condition
  } else {
    std::size_t count=0;
    for (std::size_t iP=0; iP<_fieldGroupsContainer.size(); ++iP) {
      auto particle = get(iP);
      doWhenMeetingCondition<T,PARTICLETYPE,PCONDITION>( particle,
        [&](Particle<T,PARTICLETYPE> particle){ ++count; }, globiC);
    }
    return count;
  }
}

template<typename T, typename PARTICLETYPE>
template<bool boundsCheck>
dynamics::ParticleDynamics<T,PARTICLETYPE>* ParticleSystem<T,PARTICLETYPE>::getDynamics (
  unsigned iDyn)
{
  if constexpr(boundsCheck){
    return _dynamicsVector.at(iDyn).get();
  } else {
    return _dynamicsVector[iDyn].get();
  }
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::addDynamics (
  std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr )
{
  _dynamicsVector.push_back( dynamicsSPtr );
}

template<typename T, typename PARTICLETYPE>
template <typename DYNAMICS, typename ...Args>
void ParticleSystem<T,PARTICLETYPE>::defineDynamics (Args&& ...args){
  _dynamicsVector.push_back( std::make_shared<DYNAMICS>(std::forward<Args>(args)...) );
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::extend()
{
  _fieldGroupsContainer.push_back();
//  _dynamicsMapContainer.push_back();
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::swapParticles(std::size_t iP, std::size_t jP)
{
  _fieldGroupsContainer.swapElements(iP,jP);
//  _dynamicsMapContainer.swapElements(iP,jP);
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::process (
  std::size_t iP0, std::size_t iPmax, T timeStepSize, unsigned iDyn )
{
  auto particle = get(iP0);
  for (std::size_t iP=iP0; iP<iPmax; ++iP) {
    particle.process(timeStepSize,iDyn);
    particle.advanceId();
  }
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::process(T timeStepSize, unsigned iDyn)
{
  process(0, size(), timeStepSize, iDyn);
}


template<typename T, typename PARTICLETYPE>
auto& ParticleSystem<T,PARTICLETYPE>::get()
{
  return _fieldGroupsContainer;
}

template<typename T, typename PARTICLETYPE>
std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>>& ParticleSystem<T,PARTICLETYPE>::getDynamicsVector()
{
  return _dynamicsVector;
}

template<typename T, typename PARTICLETYPE>
template<bool boundsCheck>
Particle<T,PARTICLETYPE> ParticleSystem<T, PARTICLETYPE>::get(std::size_t iParticle)
{
  if constexpr(boundsCheck){
    if (iParticle>_fieldGroupsContainer.size()){
      throw std::out_of_range("Particle does not exist in Particlesystem!");
    }
    return Particle<T,PARTICLETYPE>( _fieldGroupsContainer.data(),
              _dynamicsVector, iParticle );
  } else {
    return Particle<T,PARTICLETYPE>( _fieldGroupsContainer.data(),
              _dynamicsVector, iParticle );
  }
}

template<typename T, typename PARTICLETYPE>
Particle<T,PARTICLETYPE> ParticleSystem<T, PARTICLETYPE>::operator[](std::size_t iParticle)
{
  return get(iParticle);
}

template<typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
auto& ParticleSystem<T, PARTICLETYPE>::getFieldD()
{
  using GROUP_EVAL = typename PARTICLETYPE::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<FIELD>;
  return _fieldGroupsContainer.data().template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>();
}

template<typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
auto ParticleSystem<T, PARTICLETYPE>::getFieldPointer( std::size_t iParticle )
{
  using GROUP_EVAL = typename PARTICLETYPE::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<FIELD>;
  return _fieldGroupsContainer.data().template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>()
         .getRowPointer( iParticle );
}

template<typename T, typename PARTICLETYPE>
template<typename TYPE>
auto& ParticleSystem<T, PARTICLETYPE>::getAssociatedData()
{
  const std::size_t idx = associatedTypes::template index<TYPE>();
  return std::get<idx>(_associatedData);
}


template<typename T, typename PARTICLETYPE>
std::size_t ParticleSystem<T,PARTICLETYPE>::getSerialSize() const
{
  return _serialSize;
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::print()
{
  OstreamManager clout(std::cout, "ParticleSystem");
  clout << "-----ParticleSystem------" << std::endl;
  clout << "nop=" << size() << std::endl;
  clout << "-------------------------" << std::endl;
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::checkForErrors()
{
  //Check, whether at least one ParticleDynamics present
  if ( _dynamicsVector.size()==0 ){
    throw std::logic_error("No ParticleDynamics set!");
  }
}

} //namespace particles

} //namespace olb


#endif
