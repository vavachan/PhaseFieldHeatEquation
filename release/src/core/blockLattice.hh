/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2021 Mathias Krause
 *                2022      Adrian Kummerlaender
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

#ifndef BLOCK_LATTICE_HH
#define BLOCK_LATTICE_HH

#include "blockLattice.h"

#include "introspection.h"

#include <iterator>

namespace olb {


template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
BlockLattice<T,DESCRIPTOR>::BlockLattice(Vector<int,DESCRIPTOR::d> size, int padding, Platform platform)
  : BlockStructure<DESCRIPTOR>(size, padding),
    _platform(platform),
    _statisticsEnabled{true},
    _statistics{nullptr},
    _introspectable{true}
{ }

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
BlockLattice<T,DESCRIPTOR>::~BlockLattice()
{
  if (_statistics) {
    delete _statistics;
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::initialize()
{
  addPostProcessor(typeid(stage::PostStream),
                   meta::id<StatisticsPostProcessor>());

  // TODO: Move to StatisticsPostProcessor::setup
  _statistics = new LatticeStatistics<T>;
  _statistics->initialize();

  setProcessingContext(ProcessingContext::Simulation);
  postProcess();
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineRho(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  if (!indicator.isEmpty()) {
    Vector<T, DESCRIPTOR::d> physR;
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        rho(&rhoTmp,physR.data());
        get(loc).defineRho(rhoTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineU(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    Vector<T, DESCRIPTOR::d> physR;
    T uTmp[DESCRIPTOR::d] = { };
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp,physR.data());
        get(loc).defineU(uTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineRhoU(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
  AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    Vector<T, DESCRIPTOR::d> physR;
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);

        rho(&rhoTmp,physR.data());
        u(uTmp,physR.data());
        get(loc).defineRhoU(rhoTmp,uTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& popF)
{
  if (!indicator.isEmpty()) {
    T pop[DESCRIPTOR::q];
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        auto physR = indicator.getBlockGeometry().getPhysR(loc);
        popF(pop,physR.data());
        get(loc).definePopulations(pop);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::definePopulations(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                   BlockF<T,DESCRIPTOR::d>& popF)
{
  if (!indicator.isEmpty()) {
    T pop[DESCRIPTOR::q];
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        popF(pop, loc.data());
        get(loc).definePopulations(pop);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(
  BlockIndicatorF<T,DESCRIPTOR::d>& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  if (!indicator.isEmpty()) {
    // Don't use FieldD here as long as AnalyticalF is fixed to T
    std::array<T,DESCRIPTOR::template size<FIELD>()> fieldTmp;
    Vector<T, DESCRIPTOR::d> physR;
    getField<FIELD> ();
    this->forSpatialLocationsParallel([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);

        field(fieldTmp.data(), physR.data());
        get(loc).template setField<FIELD>(fieldTmp.data());
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                             BlockF<T,DESCRIPTOR::d>& field)
{
  if (!indicator.isEmpty()) {
    FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        field(fieldTmp.data(), loc.data());
        get(loc).template setField<FIELD>(fieldTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
template<typename FIELD>
void BlockLattice<T,DESCRIPTOR>::defineField(BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicatorF,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  BlockIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicator(indicatorF, blockGeometry);
  defineField<FIELD>(indicator, field);
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  if (!indicator.isEmpty()) {
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    Vector<T,DESCRIPTOR::d> physR;
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, physR.data());
        rho(&rhoTmp, physR.data());
        get(loc).iniEquilibrium(rhoTmp, uTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                BlockF<T,DESCRIPTOR::d>& u)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, loc.data());
        rho(&rhoTmp, physR);
        get(loc).iniEquilibrium(rhoTmp, uTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::iniRegularized(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  if (!indicator.isEmpty()) {
    T physR[DESCRIPTOR::d] = { };
    T uTmp[DESCRIPTOR::d] = { };
    T rhoTmp = T();
    T piTmp[util::TensorVal<DESCRIPTOR>::n] = { };
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        indicator.getBlockGeometry().getPhysR(physR, loc);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        pi(piTmp, physR);
        get(loc).iniRegularized(rhoTmp, uTmp, piTmp);
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
template<typename DYNAMICS>
void BlockLattice<T,DESCRIPTOR>::defineDynamics()
{
  this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
    setDynamics(this->getCellId(loc),
                DynamicsPromise(meta::id<DYNAMICS>{}));
  });
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
template<typename DYNAMICS>
void BlockLattice<T,DESCRIPTOR>::defineDynamics(BlockIndicatorF<T,DESCRIPTOR::d>& indicator)
{
  if (!indicator.isEmpty()) {
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        setDynamics(this->getCellId(loc),
                    DynamicsPromise(meta::id<DYNAMICS>{}));
      }
    });
  }
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::defineDynamics(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                                DynamicsPromise<T,DESCRIPTOR>&& promise)
{
  if (!indicator.isEmpty()) {
    this->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        defineDynamics(loc,
                       std::forward<DynamicsPromise<T,DESCRIPTOR>&&>(promise));
      }
    });
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::collide()
{
  if (_customCollisionO) {
    _customCollisionO->operator()(*this);
  } else {
    if constexpr (isPlatformCPU(PLATFORM)) {
      _dynamicsMap.collide(CollisionDispatchStrategy::Dominant);
    } else {
      _dynamicsMap.collide(CollisionDispatchStrategy::Individual);
    }
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::stream()
{
  if constexpr (PLATFORM == Platform::GPU_CUDA) {
    DESCRIPTOR::template filter<descriptors::is_propagatable_field>
              ::for_each([&](auto field) {
      auto& population = getField(field.get());
      for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        population[iPop].rotate(this->getNeighborDistance(descriptors::c<DESCRIPTOR>(iPop)));
      }
      // Required for resolving populations in (gather,scatter)_any_fields kernel
      getDataRegistry().refreshDeviceFieldArray(population);
    });
  } else {
    DESCRIPTOR::template filter<descriptors::is_propagatable_field>
              ::for_each([&](auto field) {
      auto& population = getField(field.get());
      for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        population[iPop].rotate(this->getNeighborDistance(descriptors::c<DESCRIPTOR>(iPop)));
      }
    });
  }
}

/// Operator for striping off density offset
/**
 * Used by BlockLattice::stripeOffDensityOffset
 **/
struct StripeOffDensityOffsetO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct OFFSET : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<OFFSET>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto offset = parameters.template get<OFFSET>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] -= descriptors::t<V,DESCRIPTOR>(iPop) * offset;
    }
  }
};

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void BlockLattice<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  if (!hasPostProcessor(typeid(StripeOffDensityOffsetO),
                        meta::id<StripeOffDensityOffsetO>{})) {
    addPostProcessor(typeid(StripeOffDensityOffsetO),
                     meta::id<StripeOffDensityOffsetO>{});
  }
  // Update offset and execute
  getData<OperatorParameters<StripeOffDensityOffsetO>>()
    .template set<StripeOffDensityOffsetO::OFFSET>(offset);
  postProcess(typeid(StripeOffDensityOffsetO));
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
LatticeStatistics<T>& BlockLattice<T,DESCRIPTOR>::getStatistics()
{
  return *_statistics;
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
const LatticeStatistics<T>& BlockLattice<T,DESCRIPTOR>::getStatistics() const
{
  return *_statistics;
}


template<typename T, typename DESCRIPTOR, Platform PLATFORM>
ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::ConcreteBlockLattice(Vector<int,DESCRIPTOR::d> size,
                                                                  int padding)
  : BlockLattice<T,DESCRIPTOR>(size, padding, PLATFORM),
    _data(),
    _descriptorFields(),
    _dynamicsMap(*this)
{
  DESCRIPTOR::fields_t::for_each([&](auto id) {
    using field = typename decltype(id)::type;
    using field_type = Array<field>;
    auto& fieldArray = _data.template allocate<field_type>(this->getNcells());
    _descriptorFields.template set<field>(&fieldArray);
    _communicatables[typeid(field)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
      ColumnVector<typename ImplementationOf<typename field::template column_type<T>,PLATFORM>::type,
                 DESCRIPTOR::template size<field>()>
    >(fieldArray));
    _data.template setSerialization<field_type>(true);
  });

  for (CellID iCell=0; iCell < this->getNcells(); ++iCell) {
    _dynamicsMap.set(iCell, meta::id<NoDynamics<T,DESCRIPTOR>>{});
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
bool ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::hasData() const
{
  return _data.template provides<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
const auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getData() const
{
  OLB_ASSERT(_data.template provides<FIELD_TYPE>(),
             "FIELD_TYPE must be allocated to be accessed");
  return _data.template get<FIELD_TYPE>();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD_TYPE>
auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getData()
{
  if (_data.template provides<FIELD_TYPE>()) {
    return _data.template get<FIELD_TYPE>();
  } else {
    // TODO: Implement more generic approach to constructing arbitrary data from specific args
    auto& data = _data.template allocate<FIELD_TYPE>(this->getNcells());
    // Manage serializables and communicatables for array fields
    using concrete_data_t = typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>;
    if constexpr (std::is_base_of_v<ColumnVectorBase,concrete_data_t>) {
      using field_t = typename concrete_data_t::field_t;
      if constexpr (field_t::isSerializable()) {
        _data.template setSerialization<FIELD_TYPE>(true);
      }
      _communicatables[typeid(field_t)] = std::unique_ptr<Communicatable>(new ConcreteCommunicatable<
        ColumnVector<typename ImplementationOf<typename field_t::template column_type<T>,PLATFORM>::type,
                   DESCRIPTOR::template size<field_t>()>
      >(data));
    }
    return data;
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getField(FIELD)
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
template<typename FIELD>
const auto& ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getField(FIELD) const
{
  if constexpr (DESCRIPTOR::fields_t::template contains<FIELD>()) {
    return static_cast<const FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(*_descriptorFields.template get<FIELD>());
  } else {
    return getData<Array<FIELD>>();
  }
  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::postProcess(std::type_index stage)
{
  for (auto& [_, postProcessorsOfPriority] : _postProcessors[stage]) {
    postProcessorsOfPriority.apply();
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::writeDescription(std::ostream& clout) const
{
  const auto dynamics = _dynamicsMap.getAll();
  for (const auto& promise : dynamics) {
    auto simplifiedName = introspection::getSimplifiedDynamicsName(promise.name());
    std::istringstream stream(simplifiedName);
    std::string fragment;
    std::vector<std::string> fragments;
    while (std::getline(stream, fragment, '\n')) {
      if (!fragment.empty()) {
        fragments.push_back(fragment);
      }
    }

    clout << "---Dynamics Details---" << std::endl;
    clout << " name: " << fragments[0] << std::endl;
    for (unsigned i=1; i < fragments.size(); ++i) {
    clout << "       " << fragments[i] << std::endl;
    }
    clout << " weight:           " << _dynamicsMap.getWeight(promise) << std::endl;
    if (auto optimizable = promise.isOptimizable()) {
    clout << " isOptimizable:    " << *optimizable << std::endl;
    }
    clout << " isOptimized:      " << promise.hasOptimizedVersion() << std::endl;
    if (auto count = promise.getArithmeticOperationCount()) {
    clout << " operations[FLOP]: " << *count << std::endl;
    }
    if (auto count = promise.getMemoryBandwidth()) {
    clout << " bandwidth[byte]:  " << *count << std::endl;
    }
    clout << "----------------------" << std::endl;
  }

  for (const auto& [stage, map] : _postProcessors) {
    clout << "---Stage Details---" << std::endl;
    clout << " name: " << stage.name() << std::endl;

    for (const auto& [priority, postProcessorsOfPriority] : map) {
      clout << std::endl;
      const auto operators = postProcessorsOfPriority.getAll();
      for (const auto& promise : operators) {
        clout << "---Post Processor Details---" << std::endl;
        clout << " name:        " << promise.name() << std::endl;
        clout << " priority:    " << promise.priority() << std::endl;
        clout << " scope:       " << getName(promise.scope()) << std::endl;
        clout << " weight:      " << postProcessorsOfPriority.getWeight(promise) << std::endl;
        clout << " isOptimized: " << promise.hasOptimizedVersion() << std::endl;
        clout << "---------------------------" << std::endl;
      }
    }

    clout << "-------------------" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::writeDynamicsAsCSV(std::ostream& clout) const
{
  const auto dynamics = _dynamicsMap.getAll();
  clout << "name; weight; isOptimizable; isOptimized; FLOP; complexity; bandwidth" << std::endl;
  for (const auto& promise : dynamics) {
    clout << promise.name()
          << "; " << _dynamicsMap.getWeight(promise);
    if (auto optimizable = promise.isOptimizable()) {
    clout << "; " << *optimizable;
    } else {
    clout << "; -1";
    }
    clout << "; " << promise.hasOptimizedVersion();
    if (auto count = promise.getArithmeticOperationCount()) {
    clout << "; " << *count;
    } else {
    clout << "; -1";
    }
    if (auto count = promise.getComplexity()) {
    clout << "; " << *count;
    } else {
    clout << "; -1";
    }
    if (auto count = promise.getMemoryBandwidth()) {
    clout << "; " << *count;
    } else {
    clout << "; -1";
    }
    clout << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::writeOperatorAsCSV(std::ostream& clout) const
{
  clout << "name; descriptor; weight; isOptimizable; isOptimized, complexity" << std::endl;
  for (const auto& [stage, map] : _postProcessors) {
    for (const auto& [priority, postProcessorsOfPriority] : map) {
      const auto operators = postProcessorsOfPriority.getAll();
      for (const auto& promise : operators) {
        clout << promise.name()
              << "; " << fields::name<DESCRIPTOR>()
              << "; " << postProcessorsOfPriority.getWeight(promise);
        if (auto optimizable = promise.isOptimizable()) {
          clout << "; " << *optimizable;
        } else {
          clout << "; -1";
        }
        clout << "; " << promise.hasOptimizedVersion();
        if (auto count = promise.getArithmeticOperationCount()) {
          clout << "; " <<* count;
        } else {
          clout << "; -1";
        }
        clout << std::endl;
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  return _data.getNblock();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  return _data.getSerializableSize();
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data, loadingMode);

  return dataPtr;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>::postLoad()
{
  auto& population = getField<descriptors::POPULATION>();
  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    population[iPop].postLoad();
  }
}


#ifdef PARALLEL_MODE_MPI

/// Wrapper for a non-blocking block propagation send request
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::SendTask {
private:
  const std::vector<CellID>& _cells;

  MultiConcreteCommunicatable<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> _source;

  std::unique_ptr<std::uint8_t[]> _buffer;
  MpiSendRequest _request;

public:
  SendTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block):
    _cells(cells),
    _source(block, fields),
    _buffer(new std::uint8_t[_source.size(_cells)] { }),
    _request(_buffer.get(), _source.size(_cells),
             rank, tag, comm)
  { }

  void send()
  {
    _source.serialize(_cells, _buffer.get());
    _request.start();
  }

  void wait()
  {
    _request.wait();
  }
};

/// Wrapper for a non-blocking block propagation receive request
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::RecvTask {
private:
  const int _tag;
  const int _rank;
  const std::vector<CellID>& _cells;

  MultiConcreteCommunicatable<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> _target;

  std::unique_ptr<std::uint8_t[]> _buffer;
  MpiRecvRequest _request;

public:
  /// Manual replacement for std::reference_wrapper<RecvTask>
  /**
   * Used to track pending receive requests in std::set.
   *
   * This is a workaround for problematic external definition of
   * dependently-typed comparision operators for nested classes.
   * Reconsider as soon as depending on C++17 is allowed.
   **/
  class ref {
  private:
    RecvTask& _task;
  public:
    ref(std::unique_ptr<RecvTask>& task): _task(*task) { };

    RecvTask* operator->() const
    {
      return &_task;
    }

    bool operator <(const ref& rhs) const
    {
      return _task < rhs._task;
    }
  };

  RecvTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block):
    _tag(tag),
    _rank(rank),
    _cells(cells),
    _target(block, fields),
    _buffer(new std::uint8_t[_target.size(_cells)] { }),
    _request(_buffer.get(), _target.size(_cells),
             _rank, _tag, comm)
  { }

  bool operator<(const RecvTask& rhs) const
  {
    return  _rank  < rhs._rank
        || (_rank == rhs._rank && _tag < rhs._tag);
  }

  void receive()
  {
    _request.start();
  };

  bool isDone()
  {
    return _request.isDone();
  }

  void unpack()
  {
    _target.deserialize(_cells, _buffer.get());
  }
};

#endif // PARALLEL_MODE_MPI

/// Wrapper for a local plain-copy block communication request
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::CopyTask {
  virtual ~CopyTask() { }

  virtual void copy() = 0;
  virtual void wait() = 0;
};

/// Wrapper for a local homogeneous CPU block communication request
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::HomogeneousCopyTask
  : public ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::CopyTask {
private:
  const std::vector<CellID>& _targetCells;
  const std::vector<CellID>& _sourceCells;

  MultiConcreteCommunicatable<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> _target;
  MultiConcreteCommunicatable<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>> _source;

  std::unique_ptr<std::uint8_t[]> _buffer;

public:
  HomogeneousCopyTask(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& source):
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(target, fields),
    _source(source, fields),
    _buffer(new std::uint8_t[_source.size(_sourceCells)] { })
  {
    OLB_ASSERT(_sourceCells.size() == _targetCells.size(),
               "Source cell count must match target cell count");
  }

  void copy() override
  {
    _source.serialize(_sourceCells, _buffer.get());
    _target.deserialize(_targetCells, _buffer.get());
  };

  void wait() override { }

};


template <typename T, typename DESCRIPTOR, Platform PLATFORM>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::ConcreteBlockCommunicator(
  SuperLattice<T,DESCRIPTOR>& super,
  LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T>& tagCoordinator,
  MPI_Comm comm,
#endif
  int iC,
  const BlockCommunicationNeighborhood<T,DESCRIPTOR::d>& neighborhood):
  _iC(iC)
#ifdef PARALLEL_MODE_MPI
, _mpiCommunicator(comm)
#endif
{
#ifdef PARALLEL_MODE_MPI
  neighborhood.forNeighbors([&](int remoteC) {
    if (loadBalancer.isLocal(remoteC)) {
      const Platform remotePlatform = loadBalancer.platform(loadBalancer.loc(remoteC));
      if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
        switch (remotePlatform) {
#ifdef PLATFORM_GPU_CUDA
          case Platform::GPU_CUDA:
            // Use manual copy for local GPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_CUDA,PLATFORM>(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(loadBalancer.loc(remoteC))));
            break;
#endif
          default:
            // Use manual copy for local CPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HomogeneousCopyTask(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(loadBalancer.loc(remoteC))));
            break;
        }
      }
    } else {
      if (!neighborhood.getCellsOutboundTo(remoteC).empty()) {
        _sendTasks.emplace_back(std::make_unique<SendTask>(
          _mpiCommunicator, tagCoordinator.get(loadBalancer.glob(_iC), remoteC),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsOutboundTo(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(_iC)));
      }
      if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
        _recvTasks.emplace_back(std::make_unique<RecvTask>(
          _mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsInboundFrom(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(_iC)));
      }
    }
  });

#else // not using PARALLEL_MODE_MPI
  neighborhood.forNeighbors([&](int localC) {
    if (!neighborhood.getCellsInboundFrom(localC).empty()) {
      _copyTasks.emplace_back(new HomogeneousCopyTask(
        neighborhood.getFieldsCommonWith(localC),
        neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(_iC),
        neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>(loadBalancer.loc(localC))));
    }
  });
#endif
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::~ConcreteBlockCommunicator()
{ }

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::receive()
{
  for (auto& task : _recvTasks) {
    task->receive();
  }
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::send()
{
  for (auto& task : _sendTasks) {
    task->send();
  }
  for (auto& task : _copyTasks) {
    task->copy();
  }
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::unpack()
{
  std::set<typename RecvTask::ref> pending(_recvTasks.begin(), _recvTasks.end());
  while (!pending.empty()) {
    auto task_iterator = pending.begin();
    while (task_iterator != pending.end()) {
      auto& task = *task_iterator;
      if (task->isDone()) {
        task->unpack();
        task_iterator = pending.erase(task_iterator);
      }
      else {
        ++task_iterator;
      }
    }
  }
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::wait()
{
  for (auto& task : _copyTasks) {
    task->wait();
  }
  for (auto& task : _sendTasks) {
    task->wait();
  }
}

#else // not using PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>::copy()
{
  for (auto& task : _copyTasks) {
    task->copy();
  }
}

#endif

}

#endif
