/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause, Peter Weisbrod
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

/** \file
 * Representation of the 2D geometry -- generic implementation.
 */

#ifndef SUPER_GEOMETRY_HH
#define SUPER_GEOMETRY_HH

#include "utilities/omath.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

#include "geometry/cuboid.h"
#include "geometry/cuboidDecomposition.h"
#include "geometry/superGeometry.h"

#include "communication/superStructure.h"
#include "communication/loadBalancer.h"

#include "functors/analytical/indicator/indicatorF2D.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

#include "io/ostreamManager.h"
#include "io/superVtmWriter2D.h"
#include "io/superVtmWriter3D.h"

namespace olb {


template<typename T, unsigned D>
SuperGeometry<T,D>::SuperGeometry(CuboidDecomposition<T,D>& cuboidDecomposition,
                                  LoadBalancer<T>& loadBalancer,
                                  int overlap):
  SuperStructure<T,D>(cuboidDecomposition, loadBalancer, overlap),
  _communicator(new SuperCommunicator<T,SuperGeometry<T,D>>(*this)),
  _communicationNeeded(false),
  _statistics(this),
  _clout(std::cout, ("SuperGeometry" + std::to_string(D) + "D")),
  _iConstructionT{0},
  _writeIncrementalVtkEnabled{true}
{
  for (int iCloc=0; iCloc < this->getLoadBalancer().size(); iCloc++) {
    int iCglob = this->getLoadBalancer().glob(iCloc);
    _block.emplace_back(
      new BlockGeometry<T,D>(cuboidDecomposition.get(iCglob), overlap, iCglob));
  }

  _communicator->template requestField<descriptors::MATERIAL>();
  _communicator->requestOverlap(this->_overlap);
  _communicator->exchangeRequests();

  _statistics.getStatisticsStatus() = true;
  _communicationNeeded = true;
  updateStatistics(false);

  //writeIncrementalVTK();
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::writeIncrementalVTK()
{
  if (_writeIncrementalVtkEnabled) {
    SuperVTMwriter<T,D> writer("geometry");
    if (_iConstructionT == 0) {
      writer.createMasterFile();
    }
    SuperGeometryF<T,D> geometryF(*this);
    writer.addFunctor(geometryF);
    writer.write(_iConstructionT++);
  }
}

template<typename T, unsigned D>
int SuperGeometry<T,D>::get(LatticeR<D+1> latticeR) const
{
  if ( this->getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank() ) {
    return _block[this->getLoadBalancer().loc(latticeR[0])]->get(latticeR.data()+1);
  } else {
    throw std::domain_error("read only access to data which is not available locally");
  }
}

template<typename T, unsigned D>
int SuperGeometry<T,D>::getAndCommunicate(LatticeR<D+1> latticeR) const
{
  int material = 0;
  if ( this->getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank() ) {
    material = _block[this->getLoadBalancer().loc(latticeR[0])]->get(latticeR.data()+1);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(&material, 1, this->_loadBalancer.rank(latticeR[0]));
#endif
  return material;
}

template<typename T, unsigned D>
Vector<T,D> SuperGeometry<T,D>::getPhysR(LatticeR<D+1> latticeR) const
{
  return this->_cuboidDecomposition.getPhysR(latticeR);
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::getPhysR(T output[D], const int latticeR[D+1]) const
{
  auto physR = this->_cuboidDecomposition.getPhysR(latticeR);
  std::copy(physR.data(), physR.data() + D, output);
}

template<typename T, unsigned D>
BlockGeometry<T,D>& SuperGeometry<T,D>::getBlockGeometry(int locIC)
{
  _statistics.getStatisticsStatus() = true;
  return *_block[locIC];
}

template<typename T, unsigned D>
BlockGeometry<T,D> const& SuperGeometry<T,D>::getBlockGeometry(int locIC) const
{
  return *_block[locIC];
}

template<typename T, unsigned D>
template <typename BLOCK>
BLOCK& SuperGeometry<T,D>::getBlock(int locIC)
{
  _statistics.getStatisticsStatus() = true;
  return *_block[locIC];
}

template<typename T, unsigned D>
template <typename BLOCK>
const BLOCK& SuperGeometry<T,D>::getBlock(int locIC) const
{
  return *_block[locIC];
}

template<typename T, unsigned D>
SuperGeometryStatistics<T,D>& SuperGeometry<T,D>::getStatistics()
{
  if (this->_communicationNeeded) {
    this->communicate();
    getStatisticsStatus()=true;
  }
  return _statistics;
}

template<typename T, unsigned D>
const SuperGeometryStatistics<T,D>& SuperGeometry<T,D>::getStatistics() const
{
  return _statistics;
}

template<typename T, unsigned D>
bool& SuperGeometry<T,D>::getStatisticsStatus()
{
  return _statistics.getStatisticsStatus();
}

template<typename T, unsigned D>
bool const& SuperGeometry<T,D>::getStatisticsStatus() const
{
  return _statistics.getStatisticsStatus();
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::updateStatistics(bool verbose)
{
  if (this->_communicationNeeded) {
    this->communicate();
    getStatisticsStatus()=true;
  }
  _statistics.update(verbose);
  for (unsigned iC=0; iC<_block.size(); iC++) {
    _block[iC]->getStatistics().update(verbose);
  }
}

template<typename T, unsigned D>
template<typename DESCRIPTOR>
int SuperGeometry<T,D>::clean(bool verbose, std::vector<int> bulkMaterials)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_block.size(); iC++) {
    counter+=_block[iC]->template clean <DESCRIPTOR>(false, bulkMaterials);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    _clout << "cleaned "<< counter << " outer boundary voxel(s)" << std::endl;
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
  return counter;
}

template<typename T, unsigned D>
int SuperGeometry<T,D>::outerClean(bool verbose, std::vector<int> bulkMaterials)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_block.size(); iC++) {
    counter+=_block[iC]->outerClean(false, bulkMaterials);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    _clout << "cleaned "<< counter << " outer fluid voxel(s)" << std::endl;
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
  return counter;
}

template<typename T, unsigned D>
int SuperGeometry<T,D>::innerClean(bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_block.size(); iC++) {
    counter+=_block[iC]->innerClean(false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().barrier();
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    _clout << "cleaned "<< counter << " inner boundary voxel(s)" << std::endl;
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
  return counter;
}

template<typename T, unsigned D>
int SuperGeometry<T,D>::innerClean(int bcType, bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_block.size(); iC++) {
    counter+=_block[iC]->innerClean(bcType,false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    _clout << "cleaned "<< counter << " inner boundary voxel(s) of Type " << bcType << std::endl;
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
  return counter;
}

template<typename T, unsigned D>
bool SuperGeometry<T,D>::checkForErrors(bool verbose)
{
  updateStatistics(verbose);
  bool error = false;
  for (unsigned iC=0; iC<_block.size(); iC++) {
    if (_block[iC]->checkForErrors(false)) {
      error = true;
    }
  }
  if (verbose) {
    if (error) {
      this->_clout << "error!" << std::endl;
    }
    else {
      this->_clout << "the model is correct!" << std::endl;
    }
  }
  return error;
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::reset(IndicatorF<T,D>& domain)
{
  this->communicate();
  for (unsigned iC = 0; iC < _block.size(); ++iC) {
    _block[iC]->reset(domain);
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromM, int toM)
{
  this->communicate();
  for (unsigned iC=0; iC<_block.size(); iC++) {
    _block[iC]->rename(fromM,toM);
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
  writeIncrementalVTK();
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromM, int toM, FunctorPtr<IndicatorF<T,D>>&& condition)
{
  this->communicate();
  for (unsigned iC=0; iC<_block.size(); iC++) {
    _block[iC]->rename(fromM,toM,*condition);
  }
  _statistics.getStatisticsStatus() = true;
  writeIncrementalVTK();
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromM, int toM, LatticeR<D> offset)
{
  LatticeR<D> overlap (this->_overlap);
  if ( offset <= overlap ){
    _communicator->communicate();
    for (unsigned iC=0; iC<_block.size(); iC++) {
      _block[iC]->rename(fromM,toM,offset);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
    writeIncrementalVTK();
  }else{
    _clout << "error rename only implemented for offset<=overlap" << std::endl;
  }
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromM, int toM, int testM, Vector<int,D> testDirection)
{
  if ( testDirection[0]*testDirection[0]<=(this->_overlap)*(this->_overlap)
        && testDirection[1]*testDirection[1]<=(this->_overlap)*(this->_overlap)  ){
    if constexpr (D==3){
      if(testDirection[2]*testDirection[2]<=(this->_overlap)*(this->_overlap)){
        _communicator->communicate();
        for (unsigned iC=0; iC<_block.size(); iC++) {
          _block[iC]->rename(fromM,toM,testM,testDirection);
        }
        _statistics.getStatisticsStatus() = true;
        this->_communicationNeeded = true;
      }
    }else{
      _communicator->communicate();
      for (unsigned iC=0; iC<_block.size(); iC++) {
        _block[iC]->rename(fromM,toM,testM,testDirection);
      }
      _statistics.getStatisticsStatus() = true;
      this->_communicationNeeded = true;
    }
    writeIncrementalVTK();
  } else {
    _clout << "error rename only implemented for |testDirection[i]|<=overlap" << std::endl;
  }
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromBcMat, int toBcMat, int fluidMat,
                                IndicatorF<T,D>& condition)
{
  if (this->_overlap>1) {
    this->communicate();
    rename(fromBcMat, toBcMat, condition);
    Vector<int,D> testDirection = this->getStatistics().computeDiscreteNormal(toBcMat);
    this->communicate();
    for (unsigned iC=0; iC < _block.size(); iC++) {
      _block[iC]->rename(fromBcMat,toBcMat,fluidMat,condition,testDirection);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
    writeIncrementalVTK();
  }
  else {
    _clout << "error rename only implemented for overlap>=2" << std::endl;
  }
}

template<typename T, unsigned D>
void SuperGeometry<T,D>::rename(int fromBcMat, int toBcMat, int fluidMat,
                                FunctorPtr<IndicatorF<T,D>>&& condition)
{
  if (this->_overlap>1) {
    _communicator->communicate();
    rename(fromBcMat, toBcMat, *condition);
    Vector<int,D> testDirection = this->getStatistics().computeDiscreteNormal(toBcMat);
    _communicator->communicate();
    for (unsigned iC=0; iC<_block.size(); iC++) {
      _block[iC]->rename(fromBcMat,toBcMat,fluidMat,*condition,testDirection);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
    writeIncrementalVTK();
  }
  else {
    _clout << "error rename only implemented for overlap>=2" << std::endl;
  }
}


template<typename T, unsigned D>
void SuperGeometry<T,D>::print()
{
  this->_cuboidDecomposition.print();
  getStatistics().print();
}


template<typename T, unsigned D>
void SuperGeometry<T,D>::print(const T physR[D], int offset)
{
    if (offset >= this->_overlap) {
      _clout << "Warning! the offset is larger than overlap, set offset to _overlap = " << this->_overlap << std::endl;
      offset=this->_overlap;
    }
    int latticeR[D+1];

    if ( this->_cuboidDecomposition.getLatticeR(latticeR, physR) ) {
      int iCglob = latticeR[0];

      std::vector<int> loc(D, 0);
      if constexpr (D==3) {
        loc[0] = latticeR[1];
        loc[1] = latticeR[2];
        loc[2] = latticeR[3];
        _clout << "physR(center)" << "=" << "{" << latticeR[1] << "," << latticeR[2] << "," << latticeR[3] << "}" << std::endl;
      }
      else {
        loc[0] = latticeR[1];
        loc[1] = latticeR[2];
        _clout << "physR(center)" << "=" << "{" << latticeR[1] << "," << latticeR[2] << "}" << std::endl;
      }

      _block[this->getLoadBalancer().loc(iCglob)]->printNode(loc,offset);

    }
    else {
      _clout << "this point isn't inside the geometry" << std::endl;
    }
}

template<typename T, unsigned D>
std::unique_ptr<SuperIndicatorF<T,D>> SuperGeometry<T,D>::getMaterialIndicator(
  std::vector<int>&& materials)
{
  static_assert(std::is_base_of<SuperIndicatorF<T,D>, SuperIndicatorMaterial<T,D>>::value,
                "Indicator to be constructed is SuperIndicatorF implementation");

  return std::unique_ptr<SuperIndicatorF<T,D>>(
           new SuperIndicatorMaterial<T,D>(*this, std::forward<std::vector<int>>(materials))
         );
}

template<typename T, unsigned D>
std::unique_ptr<SuperIndicatorF<T,D>> SuperGeometry<T,D>::getMaterialIndicator(int material)
{
  return this->getMaterialIndicator(std::vector<int> { material });
}

template<typename T, unsigned D>
std::size_t SuperGeometry<T,D>::getNblock() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getNblock();
  });
}


template<typename T, unsigned D>
std::size_t SuperGeometry<T,D>::getSerializableSize() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getSerializableSize();
  });
}

template<typename T, unsigned D>
bool* SuperGeometry<T,D>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (std::size_t iC=0; iC < _block.size(); ++iC) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, getBlock(iC), loadingMode);
  }

  return dataPtr;
}

} // namespace olb

#endif
