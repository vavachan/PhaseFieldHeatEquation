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

#include "particles/contact/particleContact.h"
#include "particles/contact/wallContact.h"

#ifndef CONTACT_CONTAINER_H
#define CONTACT_CONTAINER_H

namespace olb {
namespace particles {
namespace contact {

template <typename T, typename PARTICLECONTACTTYPE, typename WALLCONTACTTYPE>
struct ContactContainer {
public:
  /// Constructor
  ContactContainer(int resizeParticleContactTo = 0,
                   int resizeWallContactsTo    = 0);
  /// resizeable vector containing all particle-particle contacts
  std::vector<PARTICLECONTACTTYPE> particleContacts;
  /// resizeable vector containg all particle-wall contacts
  std::vector<WALLCONTACTTYPE> wallContacts;

  /// Combine contacts with same ids
  void combineContacts();
  /// Clean contacts - remove "empty" contacts
  void cleanParticleContacts();
  /// Clean contacts - remove "empty" contacts
  void cleanWallContacts();
  /// Clean contacts - remove "empty" contacts
  void cleanContacts();

  /// Clear contacts - remove all content
  void clearParticleContacts();
  /// Clear contacts - remove all content
  void clearWallContacts();
  /// Clear contacts - remove all content
  void clearContacts();

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& operator=(
      ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE>& container);
};

} // namespace contact
} // namespace particles
} // namespace olb
#endif
