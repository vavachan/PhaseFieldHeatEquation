/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias Krause
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

#ifdef PARALLEL_MODE_OMP

#include "ompManager.h"

namespace olb {

void ompManager::init(bool verbose)
{
  setDynamic(0);
  if (verbose) {
    olb::OstreamManager clout(std::cout,"OmpManager");
    clout << "Sucessfully initialized, numThreads=" << getSize() << std::endl;
  }
}

int ompManager::getSize() const
{
  return omp_get_max_threads();
}

int ompManager::getRank() const
{
  return omp_get_thread_num();
}

void ompManager::setDynamic(int dynamicThreads)
{
  omp_set_dynamic(dynamicThreads);
}

namespace singleton {

ompManager omp() {
  return ompManager{};
}

}

}

#endif
