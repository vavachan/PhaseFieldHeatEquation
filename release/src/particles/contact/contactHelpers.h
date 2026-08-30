/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Jan E. Marquardt, Mathias J. Krause
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

#ifndef CONTACT_HELPERS_H
#define CONTACT_HELPERS_H

namespace olb {

namespace particles {

namespace contact {

template <typename T, unsigned D>
void updateMinMax(PhysR<T, D>& min, PhysR<T, D>& max, const PhysR<T, D>& pos)
{
  for (unsigned iD = 0; iD < D; ++iD) {
    min[iD] = util::min(min[iD], pos[iD]);
    max[iD] = util::max(max[iD], pos[iD]);
  }
}

//TODO: Could be moved to more general scope
template <typename T, unsigned D>
unsigned short domainDimension(PhysR<T, D>& min, PhysR<T, D>& max)
{
  unsigned short dimension=0;
  for (unsigned iD = 0; iD < D; ++iD) {
    if(max[iD]>min[iD]){
      dimension+=1;
    }
  }
  return dimension;
}

} // namespace contact
} // namespace particles
} // namespace olb
#endif
