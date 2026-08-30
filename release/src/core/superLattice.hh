/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007-2024 Mathias J. Krause,
 *                Adrian Kummerlaender,
 *                Dennis Teutscher
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

#ifndef SUPER_LATTICE_HH
#define SUPER_LATTICE_HH

#include "superLattice.h"

#include "communication/mpiManager.h"
#include "cell.h"
#include "io/base64.h"

#include "functors/analytical/analyticalF.h"

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

#include "io/serializerIO.h"

#include "geometry/cuboidDecomposition.h"
#include "geometry/superGeometry.hh"

#include "communication/loadBalancer.h"

#include "io/tupleParser.h"

namespace olb {

#ifdef FEATURE_EXPORT_CODE_GENERATION_TARGETS
namespace introspection {

static unsigned iLattice = 0;

}
#endif

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collectStatistics()
{
  T weight;
  T sum_weight = 0;
  T average_rho = 0;
  T average_energy = 0;
  T maxU = 0;
  T delta = 0;

  getStatistics().reset();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidDecomposition.get(this->_loadBalancer.glob(iC)).getDeltaR();
    weight = _block[iC]->getStatistics().getNumCells() * delta
             * delta * delta;
    sum_weight += weight;
    average_rho += _block[iC]->getStatistics().getAverageRho()
                   * weight;
    average_energy += _block[iC]->getStatistics().getAverageEnergy()
                      * weight;
    if (maxU < _block[iC]->getStatistics().getMaxU()) {
      maxU = _block[iC]->getStatistics().getMaxU();
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(sum_weight, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_rho, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_energy, MPI_SUM);
  singleton::mpi().reduceAndBcast(maxU, MPI_MAX);
#endif

  average_rho = average_rho / sum_weight;
  average_energy = average_energy / sum_weight;

  getStatistics().reset(average_rho, average_energy, maxU, (int) sum_weight);
  getStatistics().incrementTime();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidDecomposition.get(this->_loadBalancer.glob(iC)).getDeltaR();
    _block[iC]->getStatistics().reset(average_rho, average_energy,
        maxU, (int) sum_weight);
    _block[iC]->getStatistics().incrementTime();
  }
}

template<typename T, typename DESCRIPTOR>
SuperLattice<T,DESCRIPTOR>::SuperLattice(CuboidDecomposition<T,DESCRIPTOR::d>& decomposition,
                                         LoadBalancer<T>& loadBalancer,
                                         unsigned overlap)
  : SuperStructure<T,DESCRIPTOR::d>(decomposition,
                                    loadBalancer,
                                    overlap),
    _statistics()
{
  using namespace stage;

  auto& load = this->getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    auto& cuboid = this->_cuboidDecomposition.get(load.glob(iC));
    #ifdef PLATFORM_GPU_CUDA
    if (load.platform(iC) == Platform::GPU_CUDA) {
      if (gpu::cuda::device::getCount() == 0) {
        throw std::runtime_error("Load balancer requested GPU processing for cuboid "
                               + std::to_string(load.glob(iC))
                               + " but no CUDA device was found on rank "
                               + std::to_string(singleton::mpi().getRank()));

      }
    }
    #endif
    _block.emplace_back(constructUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      load.platform(iC), cuboid.getExtent(), this->getOverlap()));
  }

  {
    auto& communicator = getCommunicator(PostCollide());
    communicator.template requestField<descriptors::POPULATION>();
    communicator.requestOverlap(1); // Required for inter-block propagation
    communicator.exchangeRequests();
  }

  {
    auto& communicator = getCommunicator(Full());
    DESCRIPTOR::fields_t::for_each([&](auto field) {
      communicator.template requestField<typename decltype(field)::type>();
    });
    // VTK output includes overlap of 1, some implicit dependencies for overlap of 2 exist
    communicator.requestOverlap(std::min(2, this->getOverlap()));
    communicator.exchangeRequests();
  }

  _statisticsEnabled = true;
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice<T,DESCRIPTOR>::get(LatticeR<DESCRIPTOR::d+1> latticeR)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(latticeR[0])) {
    if constexpr (DESCRIPTOR::d == 3) {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
    } else {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
    }
  }
  else {
    throw std::domain_error("Cuboid iC must be locally available");
  }
#else
  if constexpr (DESCRIPTOR::d == 3) {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
  } else {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
  }
#endif
}

template<typename T, typename DESCRIPTOR>
template<typename... R>
std::enable_if_t<sizeof...(R) == DESCRIPTOR::d+1, Cell<T,DESCRIPTOR>>
SuperLattice<T,DESCRIPTOR>::get(R... latticeR)
{
  return get({latticeR...});
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::initialize()
{
  if (!_initialized) {
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->initialize();
    }
    _communicationNeeded = true;

#ifdef FEATURE_EXPORT_CODE_GENERATION_TARGETS
    writeSummary("lattice" + std::to_string(introspection::iLattice++));
#endif
    _initialized = true;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::printSummary() const
{
  _block[0]->writeDescription(this->clout);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::writeSummary(std::string fileName) const
{
  const auto& cGeometry = this->getCuboidDecomposition();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    const int globC = this->_loadBalancer.glob(iC);
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".log";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        cGeometry.get(globC).write(fout);
        fout << std::endl;
        _block[iC]->writeDescription(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".dynamics";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        _block[iC]->writeDynamicsAsCSV(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".operator";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        _block[iC]->writeOperatorAsCSV(fout);
        fout.close();
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS>();
  }
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS<T,DESCRIPTOR>>();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS<T,DESCRIPTOR>>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                DynamicsPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineDynamics(indicator->getBlockIndicatorF(iC), std::forward<DynamicsPromise<T,DESCRIPTOR>&&>(promise));
  }
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRho(indicator->getBlockIndicatorF(iC),
                         rho);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineU(indicator->getBlockIndicatorF(iC),
                        u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRhoU(indicator->getBlockIndicatorF(iC),
                          rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop.getBlockF(iC));
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}


template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  if (field->getBlockFSize() == this->_loadBalancer.size()) {
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                              field->getBlockF(iC));
    }
  }
  else {
    FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
    int coords[DESCRIPTOR::d+1];
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        coords[0] = iC;
        for (int j = 0; j < DESCRIPTOR::d; ++j) {
          coords[j+1] = loc[j];
        }
        if (indicator(coords)) {
          field(fieldTmp.data(), coords);
          _block[iC]->get(loc).template setField<FIELD>(fieldTmp);
        }
      });
    }
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                            field);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)

{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  SuperIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicator,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  SuperIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER>
void SuperLattice<T,DESCRIPTOR>::setParameter(FieldD<T,DESCRIPTOR,PARAMETER> field)
{
  if (!PARAMETER::template isValid<T,DESCRIPTOR,PARAMETER>(field)) {
    const std::string error = "Value of parameter " + fields::name<PARAMETER>() + " is invalid";
    throw std::invalid_argument(error);
  }
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template setParameter<PARAMETER>(field);
  }
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field)
{
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template getData<OperatorParameters<DYNAMICS>>().template set<PARAMETER>(
       std::forward<decltype(field)>(field));
  }
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field)
{
  setParameterOfDynamics<PARAMETER,DYNAMICS<T,DESCRIPTOR>>(std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, typename _DESCRIPTOR, typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setParameter(AbstractFieldArrayD<T,_DESCRIPTOR,FIELD>& abstractFieldArray)
{
  static_assert(DESCRIPTOR::template size<PARAMETER>() == DESCRIPTOR::template size<FIELD>(),
                "PARAMETER size must equal FIELD size");
  static_assert(std::is_same_v<typename PARAMETER::template value_type<T>,
                               typename FIELD::template value_type<T>*>,
                "PARAMETER must store pointers to FIELD components");
  // auto test = PARAMETER::template getMaxValue<T,DESCRIPTOR>();
  // std::cout << test<< std::endl;
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template setParameter<PARAMETER>(abstractFieldArray);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniEquilibrium(indicator->getBlockIndicatorF(iC), rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                SuperF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniEquilibrium(indicator->getBlockIndicatorF(iC), rho, u.getBlockF(iC));
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                SuperF<DESCRIPTOR::d,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniRegularized(indicator->getBlockIndicatorF(iC), rho, u, pi);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  iniRegularized(sGeometry.getMaterialIndicator(material), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collideAndStream()
{
  if (!_initialized) [[unlikely]] {
    initialize();
  }

  using namespace stage;

  waitForBackgroundTasks(PreCollide());
  auto& load = this->_loadBalancer;

  if (_statisticsEnabled) {
    setParameter<statistics::AVERAGE_RHO>(getStatistics().getAverageRho());
  }

  // Optional pre processing stage
  executePostProcessors(PreCollide());

  // Execute custom tasks (arbitrary callables)
  // (used for multi-stage models such as bubble model)
  executeCustomTasks(PreCollide());

  #ifdef PARALLEL_MODE_OMP
  #pragma omp taskloop
  #endif
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->collide();
  }

  // Communicate propagation overlap, optional post processing
  executePostProcessors(PostCollide());

  // Block-local propagation
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->stream();
  }

  // Communicate (default) post processor neighborhood and apply them
  executePostProcessors(PostStream());

  // Execute custom tasks (arbitrary callables)
  // (used for multi-stage models such as free surface)
  executeCustomTasks(PostStream());

  // Final communication stage (e.g. for external coupling)
  getCommunicator(PostPostProcess()).communicate();

  if (_statisticsEnabled) {
    collectStatistics();
  }
  _communicationNeeded = true;

#ifdef FEATURE_EXPORT_CODE_GENERATION_TARGETS
  if (introspection::iLattice > 1) {
    introspection::iLattice -= 1;
  } else {
    this->clout << "Terminating after export of code generation targets." << std::endl
                << "If you are confused by this message EXPORT_CODE_GENERATION_TARGETS was wrongly enabled." << std::endl;
    std::exit(0);
  }
#endif
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::executePostProcessors(STAGE stage)
{
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
  #endif

  getCommunicator(stage).communicate();

  auto& load = this->_loadBalancer;

  #ifdef PARALLEL_MODE_OMP
  #pragma omp taskloop
  #endif
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->template postProcess<STAGE>();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>& SuperLattice<T,DESCRIPTOR>::getCommunicator(STAGE stage)
{
  auto iter = _communicator.find(typeid(STAGE));
  if (iter == _communicator.end()) {
    iter = std::get<0>(_communicator.emplace(typeid(STAGE),
                                             std::make_unique<SuperCommunicator<T,SuperLattice>>(*this)));
  }
  return *std::get<1>(*iter);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::communicate()
{
  if (_communicationNeeded) {
    getCommunicator(stage::Full()).communicate();
    _communicationNeeded = false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->stripeOffDensityOffset(offset);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& SuperLattice<T,DESCRIPTOR>::getStatistics()
{
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& SuperLattice<T,DESCRIPTOR>::getStatistics() const
{
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                  PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  indicator->getBlockIndicatorF(iC),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::executeCustomTasks(STAGE stage)
{
  for (auto& f : _customTasks[typeid(STAGE)]) {
    f();
  }
}

template<typename T, typename DESCRIPTOR>
template <typename STAGE, typename F>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundTask(F&& f)
{
  _backgroundTasks[typeid(STAGE)].emplace_back(singleton::pool().schedule(f));
}

template<typename T, typename DESCRIPTOR>
template <typename F>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundOutput(F&& f)
{
  auto& load = this->getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    if (isPlatformCPU(load.platform(iC))) {
      scheduleBackgroundTask<stage::PreCollide>(std::bind(f, iC));
    }
  }
  for (int iC=0; iC < load.size(); ++iC) {
    if (load.platform(iC) == Platform::GPU_CUDA) {
      scheduleBackgroundTask<stage::PreContextSwitchTo<ProcessingContext::Evaluation>>(std::bind(f, iC));
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename CONTEXT>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundOutputVTK(CONTEXT&& vtkContext)
{
  vtkContext([](auto& writer, std::size_t iT) {
    writer.writePVD(iT);
  });
  scheduleBackgroundOutput([vtkContext](int iC) {
    vtkContext([&](auto& writer, std::size_t iT) {
      writer.writeVTI(iT, iC);
    });
  });
}

template<typename T, typename DESCRIPTOR>
template <typename STAGE>
void SuperLattice<T,DESCRIPTOR>::waitForBackgroundTasks(STAGE)
{
  if (!_backgroundTasks[typeid(STAGE)].empty()) {
    singleton::pool().waitFor(_backgroundTasks[typeid(STAGE)]);
    _backgroundTasks[typeid(STAGE)].clear();
  }
}

template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getNblock();
  });
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getSerializableSize();
  });
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (std::size_t iC=0; iC < _block.size(); ++iC) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, getBlock(iC), loadingMode);
  }

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::postLoad()
{
  for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->postLoad();
  }
}

}

#endif
