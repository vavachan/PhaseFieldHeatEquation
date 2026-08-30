/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

/* Particle contitions provide the possibility to perform filtering
 * operations on the particle system. Both static (evaluation at compile time)
 * and dynamic conditions can be evaluated. For now, static conditions are not
 * intended do be used with template parameters.
 */


#ifndef PARTICLE_CONDITIONS_H
#define PARTICLE_CONDITIONS_H


namespace olb {

namespace particles {

namespace conditions {

struct all_particles{
  //Properties
  static constexpr bool value=true;
  static constexpr bool dynamic=false;
};

struct valid_particles{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    return access::isValid( particle );
  }
  //Properties
  static constexpr bool dynamic=true;
};

struct invalid_particles{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    return !access::isValid( particle );
  }
  //Properties
  static constexpr bool dynamic=true;
};

struct valid_particle_centres{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle, int globiC){
    using namespace descriptors;
    bool valid = access::isValid( particle );
    bool hasCentre = true;
    if constexpr ( PARTICLETYPE::template providesNested<PARALLELIZATION,IC>() ) {
      hasCentre = (particle.template getField<PARALLELIZATION,IC>() == globiC);
    }
    return (valid && hasCentre);
  }
  //Properties
  static constexpr bool dynamic=true;
};

struct valid_particle_surfaces{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle, int globiC){
    using namespace descriptors;
    bool valid = access::isValid( particle );
    bool hasCentre = true;
    if constexpr ( PARTICLETYPE::template providesNested<PARALLELIZATION,IC>() ) {
      hasCentre = (particle.template getField<PARALLELIZATION,IC>() == globiC);
    }
    return (valid && !hasCentre);
  }
  //Properties
  static constexpr bool dynamic=true;
};

//Active particles (implies valid particle)
struct active_particles{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    using namespace descriptors;
    bool active = access::isActive( particle );
    bool valid = access::isValid( particle );
    return (active && valid);
  }
  //Properties
  static constexpr bool dynamic=true;
};

//Inactive particles (implies valid particle)
struct inactive_particles{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    using namespace descriptors;
    bool active = access::isActive( particle );
    bool valid = access::isValid( particle );
    return (!active && valid);
  }
  //Properties
  static constexpr bool dynamic=true;
};

//Active particle centres (implies valid particle)
struct active_particle_centres{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle, int globiC){
    using namespace descriptors;
    bool active = access::isActive( particle );
    bool valid = access::isValid( particle );
    bool hasCentre = true;
    if constexpr ( PARTICLETYPE::template providesNested<PARALLELIZATION,IC>() ) {
      hasCentre = (particle.template getField<PARALLELIZATION,IC>() == globiC);
    }
    return (active && valid && hasCentre);
  }
  //Properties
  static constexpr bool dynamic=true;
};

template<std::size_t selectedID>
struct particle_matching_ID{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    using namespace descriptors;
    bool match = false;
    if constexpr ( PARTICLETYPE::template providesNested<PARALLELIZATION,ID>() ) {
      match = (particle.template getField<PARALLELIZATION,ID>() == selectedID);
    } else {
      match = (particle.getId() == selectedID);
    }
    return match;
  }
  //Properties
  static constexpr bool dynamic=true;
};

template<std::size_t selectedID>
struct valid_particle_matching_ID{
  template<typename T, typename PARTICLETYPE>
  static bool value(Particle<T,PARTICLETYPE>& particle){
    using namespace descriptors;
    bool valid = access::isValid( particle );
    bool match = false;
    if constexpr ( PARTICLETYPE::template providesNested<PARALLELIZATION>() ) {
      static_assert(PARTICLETYPE::template providesNested<PARALLELIZATION,ID>(), "Field PARALLELIZATION:ID has to be provided");
      match = (particle.template getField<PARALLELIZATION,ID>() == selectedID);
    } else {
      match = (particle.getId() == selectedID);
    }
    return (valid && match);
  }
  //Properties
  static constexpr bool dynamic=true;
};

} //namespace conditions

} //namespace particles

} //namespace olb


#endif
