/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt
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

#include "olbInit.h"

#include "io/ostreamManager.h"

#include "communication/mpiManager.h"
#include "communication/ompManager.h"

#include "core/platform/platform.h"

#ifdef FEATURE_VDB
#include <openvdb/openvdb.h>
#endif

namespace olb {

namespace singleton {

ThreadPool& pool()
{
  static ThreadPool instance;
  return instance;
}

}

void initialize(int *argc, char ***argv, bool multiOutput, bool verbose)
{
  // create an OstreamManager object in order to enable multi output
  olb::OstreamManager clout(std::cout, "olbInit");
  clout.setMultiOutput(multiOutput);
  singleton::mpi().init(argc, argv, verbose);

#ifdef PARALLEL_MODE_OMP
  singleton::omp().init(verbose);
#endif

  int nThreads = 1;
  if (const char* envOlbNumThreads = std::getenv("OLB_NUM_THREADS")) {
    nThreads = std::stoi(envOlbNumThreads);
  }
  singleton::pool().init(nThreads, verbose);

  /// Verify requirements for using all enabled platforms
  #ifdef PLATFORM_GPU_CUDA
  checkPlatform<Platform::GPU_CUDA>();
  #endif

  #ifdef FEATURE_VDB
  openvdb::initialize();
  #endif
}

void initialize(int argc, char **argv, bool multiOutput, bool verbose)
{
  initialize(&argc, &argv, multiOutput, verbose);
}

}
