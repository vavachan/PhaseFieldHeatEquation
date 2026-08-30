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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "core/fieldArrayD.h"
#include "core/blockDynamicsMap.h"

namespace olb {

namespace particles{

//Forward declaration
namespace dynamics {
template<typename T, typename PARTICLETYPE> struct ParticleDynamics;
}


//Particle interface as pendent to cell interface
template <typename T, typename PARTICLETYPE>
class Particle {
private:
  using DATA = DynamicFieldGroupsD<T, typename PARTICLETYPE::fields_t>;
  DATA& _dynamicFieldGroupsD;
  std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>>& _dynamicsVector;
  std::size_t _iParticle;
public:
  Particle( DATA& dynamicFieldGroupsD,
            std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>>& dynamicsVector,
            std::size_t iParticle );

  /// Initialize
  void init();

  ///Print
  template<bool multiOutput=PARTICLETYPE::template providesNested<descriptors::PARALLELIZATION>()>
  void print(std::size_t iParticle);
  template<bool multiOutput=PARTICLETYPE::template providesNested<descriptors::PARALLELIZATION>()>
  void print();

  /// Return memory ID of the currently represented particle
  std::size_t getId() const;

  /// Jump to arbitrary particle memory ID
  void setId( std::size_t iParticle );

  /// Jump to next particle in linearization sequence
  void advanceId();

  /// Add dynamics
  void addDynamics(std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr );

  /// Define dynamics (factory method)
  template <typename DYNAMICS, typename ...Args>
  void defineDynamics(Args&& ...args);

  /// Get a pointer to specific dynamics
  template<bool boundsCheck=false>
  dynamics::ParticleDynamics<T,PARTICLETYPE>* getDynamics(unsigned iDyn=0);

  /// Apply processing to the particle according to dynamics at iDyn
  void process(T timeStepSize, unsigned short iDyn=0);

  ////////// Get and Set functions //TODO: implement recursively?

  /// Return copy of descriptor-declared FIELD as a vector
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() > 1),
      FieldD<T,PARTICLETYPE,
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>>>
          getField();

  /// Return copy of descriptor-declared FIELD as a scalar (more specifically as value type of FIELD)
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() == 1),
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>::template value_type<T>>
  getField();

  /// Set descriptor-declared FIELD
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() > 1),
      void>
      setField(const FieldD<T,PARTICLETYPE,
               typename PARTICLETYPE::template derivedField<GROUP>
               ::template derivedField<FIELD>>& v);

  /// Set descriptor-declared FIELD
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() == 1),
      void>
      setField( typename PARTICLETYPE::template derivedField<GROUP>
                ::template derivedField<FIELD>
                ::template value_type<T> value );


  /// Define Communicatable
  using Communicatable = typename GroupedDataCommunicatableHelper<DATA,PARTICLETYPE>::type;
  /// Get communicatable necessary for serialization
  Communicatable getCommunicatable(){
    return Communicatable(_dynamicFieldGroupsD);
  }
  /// Get serialized size
  std::size_t getSerialSize() const {
    const std::vector<unsigned int> indices{static_cast<unsigned int>(_iParticle)};
    return Communicatable(_dynamicFieldGroupsD).size(indices);
  }
  /// Serialize data
  std::size_t serialize(std::uint8_t* buffer) const {
    const std::vector<unsigned int> indices{static_cast<unsigned int>(_iParticle)};
    return Communicatable(_dynamicFieldGroupsD).serialize(indices, buffer);
  }
  /// Deserialize data
  std::size_t deserialize(const std::uint8_t* buffer){
    const std::vector<unsigned int> indices{static_cast<unsigned int>(_iParticle)};
    return Communicatable(_dynamicFieldGroupsD).deserialize(indices, buffer);
  }

};

} //namespace particles

} //namespace olb

#endif
