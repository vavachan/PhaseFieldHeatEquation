/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Fietz, Mathias J. Krause
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

#ifndef HEURISTIC_LOAD_BALANCER_HH
#define HEURISTIC_LOAD_BALANCER_HH

#include <algorithm>
#include <vector>
#include <map>
#include <math.h>
#include "core/cell.h"
#include "core/util.h"
#include "heuristicLoadBalancer.h"
#include "geometry/cuboid.h"
#include "geometry/cuboidDecomposition.h"

namespace olb {


template<typename T>
HeuristicLoadBalancer<T>::HeuristicLoadBalancer(CuboidDecomposition3D<T>& cGeometry3d,
    const double ratioFullEmpty, const double weightEmpty)
  : _ratioFullEmpty(ratioFullEmpty)
{
  reInit(cGeometry3d, ratioFullEmpty, weightEmpty);
}

template<typename T>
HeuristicLoadBalancer<T>::HeuristicLoadBalancer(CuboidDecomposition2D<T>& cGeometry2d,
    const double ratioFullEmpty, const double weightEmpty)
  : _ratioFullEmpty(ratioFullEmpty)
{
  reInit(cGeometry2d, ratioFullEmpty, weightEmpty);
}

template<typename T>
HeuristicLoadBalancer<T>::~HeuristicLoadBalancer()
{
}

template<typename T>
void HeuristicLoadBalancer<T>::swap(HeuristicLoadBalancer<T>& loadBalancer)
{
  LoadBalancer<T>::swap(loadBalancer);

#ifdef PARALLEL_MODE_MPI
  _mpiNbHelper.swap(loadBalancer._mpiNbHelper);
#endif

  std::swap(_cGeometry3d, loadBalancer._cGeometry3d);
  std::swap(_cGeometry2d, loadBalancer._cGeometry2d);
}

template<typename T>
void HeuristicLoadBalancer<T>::reInit(CuboidDecomposition3D<T>& cGeometry3d, const double ratioFullEmpty, const double weightEmpty)
{
  _ratioFullEmpty = ratioFullEmpty;
  this->_glob.clear();
  _cGeometry3d = &cGeometry3d;
  int rank = 0;
  int size = 1;
  int nC = _cGeometry3d->size();
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = util::max<int>(singleton::mpi().getSize(), 1);
#endif
  //int xN, yN, zN;
  //T globX, globY, globZ;//, delta;
  //boost::shared_array<int> tempInCN(new int[nC]);

  std::vector<int> cuboidToThread(nC);
  std::vector<int> partitionResult(nC);
  std::vector<int> vwgt(nC); // node weights
  std::vector<int> taken(nC, 0);
  std::vector<int> currentLoad(size, 0);

  if (size == 1) {
    for (int i = 0; i < nC; ++i) {
      this->_glob.push_back(i);
      this->_loc[i] = i;
      this->_rank[i] = 0;
    };
    this->_size = nC;
    return;
  }

  if (rank == 0) {
    for ( int iC = 0; iC < nC; iC++) { // assemble neighbourhood information

      int fullCells = _cGeometry3d->get(iC).getWeight();
      vwgt[iC] = int(weightEmpty*(_cGeometry3d->get(iC).getLatticeVolume() - fullCells)) + int(ratioFullEmpty * fullCells);
    }

    int maxLoad = -1;
    int maxIC = -1;
    do {
      maxLoad = -1;
      maxIC = -1;
      for ( int iC = 0 ; iC < nC; iC++) {
        if (taken[iC] == 0 && vwgt[iC] > maxLoad) {
          maxLoad = vwgt[iC];
          maxIC = iC;
        }
      }

      if (maxIC != -1) {
        double minLoad = currentLoad[0];
        int minJ = 0;
        for (int j = 1; j < size; j++) {
          if (currentLoad[j] < minLoad) {
            minLoad = currentLoad[j];
            minJ = j;
          }
        }
        taken[maxIC] = 1;
        currentLoad[minJ] += maxLoad;
        partitionResult[maxIC] = minJ;
      }
    }
    while (maxLoad != -1);
    int count = 0;
    for (int i = 0; i < nC; ++i) {
      if (partitionResult[i] == 0) {
        this->_glob.push_back(i);
        this->_loc[i] = count;
        count++;
      };
      this->_rank[i] = partitionResult[i];
      cuboidToThread[i] = partitionResult[i];
    }
    this->_size = count;
  }
  // Send all threads their number of cuboids

#ifdef PARALLEL_MODE_MPI
  if (rank == 0) {
    // Send all threads their respective cuboids
    _mpiNbHelper.free();
    _mpiNbHelper.allocate(size-1);
    for (int i = 1; i < size; i++) {
      singleton::mpi().iSend(&cuboidToThread.front(),
                             nC, i, &_mpiNbHelper.get_mpiRequest()[i-1], 0);
    }
    singleton::mpi().waitAll(_mpiNbHelper);
  }
  else {
    int *tmpCuboids = new int[nC];
    singleton::mpi().receive(tmpCuboids, nC, 0, 0);
    int count = 0;
    for (int i = 0; i < nC; ++i) {
      if (tmpCuboids[i] == rank) {
        this->_glob.push_back(i);
        this->_loc[i] = count;
        count++;
      };
      this->_rank[i] = tmpCuboids[i];
    }
    delete[] tmpCuboids;
    this->_size = count;
  }
#endif
}

template<typename T>
void HeuristicLoadBalancer<T>::reInit(CuboidDecomposition2D<T>& cGeometry2d, const double ratioFullEmpty, const double weightEmpty)
{
  _ratioFullEmpty = ratioFullEmpty;
  this->_glob.clear();
  _cGeometry2d = &cGeometry2d;
  int rank = 0;
  int size = 1;
  int nC = _cGeometry2d->size();
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = util::max<int>(singleton::mpi().getSize(), 1);
#endif
  //int xN, yN;
  //T globX, globY;//, delta;
  //boost::shared_array<int> tempInCN(new int[nC]);

  std::vector<int> cuboidToThread(nC);
  std::vector<int> partitionResult(nC);
  std::vector<int> vwgt(nC); // node weights
  std::vector<int> taken(nC, 0);
  std::vector<int> currentLoad(size, 0);

  if (size == 1) {
    for (int i = 0; i < nC; ++i) {
      this->_glob.push_back(i);
      this->_loc[i] = i;
      this->_rank[i] = 0;
    };
    this->_size = nC;
    return;
  }

  if (rank == 0) {
    for ( int iC = 0; iC < nC; iC++) { // assemble neighbourhood information

      int fullCells = _cGeometry2d->get(iC).getWeight();
      vwgt[iC] = int(weightEmpty*(_cGeometry2d->get(iC).getLatticeVolume() - fullCells)) + int(ratioFullEmpty * fullCells);

    }

    int maxLoad = -1;
    int maxIC = -1;
    do {
      maxLoad = -1;
      maxIC = -1;
      for ( int iC = 0 ; iC < nC; iC++) {
        if (taken[iC] == 0 && vwgt[iC] > maxLoad) {
          maxLoad = vwgt[iC];
          maxIC = iC;
        }
      }

      if (maxIC != -1) {
        double minLoad = currentLoad[0];
        int minJ = 0;
        for (int j = 1; j < size; j++) {
          if (currentLoad[j] < minLoad) {
            minLoad = currentLoad[j];
            minJ = j;
          }
        }
        taken[maxIC] = 1;
        currentLoad[minJ] += maxLoad;
        partitionResult[maxIC] = minJ;
      }
    }
    while (maxLoad != -1);
    int count = 0;
    for (int i = 0; i < nC; ++i) {
      if (partitionResult[i] == 0) {
        this->_glob.push_back(i);
        this->_loc[i] = count;
        count++;
      };
      this->_rank[i] = partitionResult[i];
      cuboidToThread[i] = partitionResult[i];
    }
    this->_size = count;
  }
  // Send all threads their number of cuboids

#ifdef PARALLEL_MODE_MPI
  if (rank == 0) {
    // Send all threads their respective cuboids
    _mpiNbHelper.free();
    _mpiNbHelper.allocate(size-1);
    for (int i = 1; i < size; i++) {
      singleton::mpi().iSend(&cuboidToThread.front(),
                             nC, i, &_mpiNbHelper.get_mpiRequest()[i-1], 0);
    }
    singleton::mpi().waitAll(_mpiNbHelper);
  }
  else {
    int *tmpCuboids = new int[nC];
    singleton::mpi().receive(tmpCuboids, nC, 0, 0);
    int count = 0;
    for (int i = 0; i < nC; ++i) {
      if (tmpCuboids[i] == rank) {
        this->_glob.push_back(i);
        this->_loc[i] = count;
        count++;
      };
      this->_rank[i] = tmpCuboids[i];
    }
    delete[] tmpCuboids;
    this->_size = count;
  }
#endif
}


}  // namespace olb
#endif
