/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt, Mathias J. Krause
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

#ifndef PARTICLE_CONTACT_CONTAINER_HH
#define PARTICLE_CONTACT_CONTAINER_HH

#include "contactContainer.h"
#include "contactFunctions.h"

namespace olb {
namespace particles {
namespace contact {

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>::ContactContainer(
    int resizeParticleContactTo, int resizeWallContactsTo)
{
  particleContacts.resize(resizeParticleContactTo);
  wallContacts.resize(resizeWallContactsTo);
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE,
                      WALLCONTACTTYPE>::combineContacts()
{
  for (auto it = std::begin(particleContacts); it != std::end(particleContacts);
       ++it) {
    for (auto it2 = it + 1; it2 != std::end(particleContacts); ++it2) {
      // combineWith must check if the particle ids are the same
      it->combineWith(*it2);
    }
  }

  for (auto it = std::begin(wallContacts); it != std::end(wallContacts); ++it) {
    for (auto it2 = it + 1; it2 != std::end(wallContacts); ++it2) {
      // combineWith must check if the contacts are the same
      it->combineWith(*it2);
    }
  }
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE,
                      WALLCONTACTTYPE>::cleanParticleContacts()
{
  particles::contact::cleanContacts(particleContacts);
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE,
                      WALLCONTACTTYPE>::cleanWallContacts()
{
  particles::contact::cleanContacts(wallContacts);
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>::cleanContacts()
{
  cleanParticleContacts();
  cleanWallContacts();
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE,
                      WALLCONTACTTYPE>::clearParticleContacts()
{
  const std::size_t capacity = particleContacts.capacity();
  particleContacts.clear();
  particleContacts.reserve(capacity / 2);
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE,
                      WALLCONTACTTYPE>::clearWallContacts()
{
  const std::size_t capacity = wallContacts.capacity();
  wallContacts.clear();
  wallContacts.reserve(capacity / 2);
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
void ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>::clearContacts()
{
  clearParticleContacts();
  clearWallContacts();
}

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>&
ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>::operator=(
    ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& container)
{
  particleContacts.clear();
  particleContacts.reserve(container.particleContacts.size());
  for (size_t i = 0; i < container.particleContacts.size(); ++i) {
    particleContacts.push_back(container.particleContacts[i]);
  }

  wallContacts.clear();
  wallContacts.reserve(container.wallContacts.size());
  for (size_t i = 0; i < container.wallContacts.size(); ++i) {
    wallContacts.push_back(container.wallContacts[i]);
  }

  return *this;
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
