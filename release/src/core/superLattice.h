/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_H
#define SUPER_LATTICE_H

#include "utilities/aliases.h"

#include "unitConverter.h"
#include "stages.h"
#include "cellD.h"
#include "blockLattice.hh"
#include "communication/superCommunicator.h"
#include "postProcessing.hh"
#include "serializer.h"
#include "communication/superStructure.hh"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"

namespace olb {


/// Super class maintaining block lattices for a cuboid decomposition
template <typename T, typename DESCRIPTOR>
class SuperLattice final : public SuperStructure<T,DESCRIPTOR::d>
                         , public BufferSerializable {
private:
  /// Converter between SI and lattice units
  std::optional<const UnitConverter<T,DESCRIPTOR>*> _converter;
  /// Lattices with ghost cell layer of size overlap
  std::vector<std::unique_ptr<BlockLattice<T,DESCRIPTOR>>> _block;
  /// Communicators for stages of collideAndStream
  std::map<std::type_index,std::unique_ptr<SuperCommunicator<T,SuperLattice>>> _communicator;
  /// True if there are changes to be communicated using manually-triggered Full stage
  bool _communicationNeeded;
  /// Map of custom stages and associated tasks
  std::map<std::type_index,std::vector<std::function<void()>>> _customTasks;
  /// Map of custom background stages and any associated one-off tasks
  std::map<std::type_index,std::vector<std::future<void>>> _backgroundTasks;

  /// Statistics of the super structure
  LatticeStatistics<T> _statistics;
  /// Specifies if statistics are to be calculated
  /**
   * Enabled by default, needed for e.g. ConstRhoBGK dynamics.
   **/
  bool _statisticsEnabled;
  /// Aggregate global statistics
  void collectStatistics();
  /// False iff initialize was not yet called
  bool _initialized = false;

public:
  constexpr static unsigned d = DESCRIPTOR::d;

  /// Base value type of the lattice
  using value_t = T;
  /// Descriptor / discrete velocity set of the lattice
  using descriptor_t = DESCRIPTOR;

  using block_t = ConcretizableBlockLattice<T,DESCRIPTOR>;

  SuperLattice(CuboidDecomposition<T,DESCRIPTOR::d>& cuboidDecomposition,
               LoadBalancer<T>& loadBalancer,
               unsigned overlap = 3);

  SuperLattice(CuboidDecomposition<T,DESCRIPTOR::d>& cuboidDecomposition,
               LoadBalancer<T>& loadBalancer,
               unsigned overlap,
               const UnitConverter<T,DESCRIPTOR>& converter)
    : SuperLattice(cuboidDecomposition, loadBalancer, overlap)
  {
    _converter = &converter;
  }

  SuperLattice(SuperGeometry<T,DESCRIPTOR::d>& sGeometry)
    : SuperLattice(sGeometry.getCuboidDecomposition(),
                   sGeometry.getLoadBalancer(),
                   sGeometry.getOverlap())
  { }

  SuperLattice(const SuperLattice&) = delete;
  ~SuperLattice() = default;

  const UnitConverter<T,DESCRIPTOR>& getConverter() const {
    if (_converter) {
      return **_converter;
    } else {
      throw std::runtime_error("Unit converter not provided");
    }
  }

  /// Return BlockLattice with local index locC
  BlockLattice<T,DESCRIPTOR>& getBlock(int locC)
  {
    return *_block[locC];
  };
  /// Return locC-th block lattice casted to BLOCK
  /**
   * Only safe to use under external knownledge that the locC-th block is of type BLOCK
   **/
  template <typename BLOCK>
  BLOCK& getBlock(int locC)
  {
    return *dynamic_cast<BLOCK*>(_block[locC].get());
  };
  /// Return read-only BlockLattice with local index locC
  const BlockLattice<T,DESCRIPTOR>& getBlock(int locC) const
  {
    return *_block[locC];
  };
  /// Return locC-th block lattice casted to BLOCK (read-only)
  template <typename BLOCK>
  const BLOCK& getBlock(int locIC) const
  {
    return *dynamic_cast<const BLOCK*>(_block[locIC].get());
  };

  /// Apply f to every ConcreteBlockLattice of PLATFORM
  template <Platform PLATFORM, typename F>
  void forBlocksOnPlatform(F f) {
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      if (_block[iC]->getPlatform() == PLATFORM) {
        f(static_cast<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>&>(*_block[iC]));
      }
    }
  };

  /// Set processing context of block lattices
  /**
   * Used to sync data between device and host
   **/
  void setProcessingContext(ProcessingContext context)
  {
    // Communicate overlap prior to evaluation
    if (context == ProcessingContext::Evaluation) {
      communicate();
      waitForBackgroundTasks<stage::PreContextSwitchTo<ProcessingContext::Evaluation>>();
    }
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->setProcessingContext(context);
    }
  };

  /// Set processing context of FIELD_TYPE in block lattices
  /**
   * Used to sync data between device and host
   **/
  template <typename FIELD_TYPE>
  void setProcessingContext(ProcessingContext context)
  {
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      if (_block[iC]->template hasData<FIELD_TYPE>()) {
        _block[iC]->template getData<FIELD_TYPE>().setProcessingContext(context);
      }
    }
  };

  /// Return communicator for given communication stage
  template <typename STAGE>
  SuperCommunicator<T,SuperLattice>& getCommunicator(STAGE stage=STAGE());
  /// Perform full overlap communication if needed
  void communicate() override;

  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics();
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const;

  /// Get local cell interface
  Cell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d+1> latticeR);
  /// Get local cell interface
  template <typename... R>
  std::enable_if_t<sizeof...(R) == DESCRIPTOR::d+1, Cell<T,DESCRIPTOR>>
  get(R... latticeR);

  /// Initialize lattice to be ready for simulation
  /**
   * Calls default post processors (PostStream) once and
   * switches processing context to simulation.
   **/
  void initialize();

  /// Prints a summary of all dynamics and post processors
  /**
   * For brevity only of the first block of the root process
   **/
  void printSummary() const;
  /// Writes a summary of all dynamics and post processors
  void writeSummary(std::string fileName = "lattice") const;

  /// Set dynamics of all cells to DYNAMICS
  template <typename DYNAMICS>
  void defineDynamics();
  /// Set dynamics of all cells to DYNAMICS
  template <template<typename...> typename DYNAMICS>
  void defineDynamics();
  template <typename DYNAMICS>
  void defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator);
  /// Set dynamics of material cells to DYNAMICS
  template <typename DYNAMICS>
  void defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material);

  /// Set dynamics of indicated cells to DYNAMICS<T,DESCRIPTOR>
  template <template<typename...> typename DYNAMICS>
  void defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator);
  /// Set dynamics of indicated cells to DYNAMICS<T,DESCRIPTOR>
  template <template<typename...> typename DYNAMICS>
  void defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material);

  //Define the dynamics on a domain with indicator and promise
  void defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,DynamicsPromise<T,DESCRIPTOR>&& promise);

  /// Define rho on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor
   **/
  void defineRho(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&&, AnalyticalF<DESCRIPTOR::d,T,T>& rho);
  /// Define rho on a domain with a particular material number
  void defineRho(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material, AnalyticalF<DESCRIPTOR::d,T,T>& rho);

  /// Define u on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param u         Analytical functor
   **/
  void defineU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator, AnalyticalF<DESCRIPTOR::d,T,T>& u);
  /// Define u on a domain with a particular material number
  void defineU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material, AnalyticalF<DESCRIPTOR::d,T,T>& u);

  /// Define rho and u on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor
   * \param u         Analytical functor
   **/
  void defineRhoU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                  AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);
  /// Define rho and u on a domain with a particular material number
  void defineRhoU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                  AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);

  /// Define a population on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param Pop       Analytical functor, target dimension DESCRIPTOR::q
   **/
  void definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                         AnalyticalF<DESCRIPTOR::d,T,T>& Pop);
  /// Define a population on a domain with a particular material number
  void definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                         AnalyticalF<DESCRIPTOR::d,T,T>& Pop);
  /**
   * \param indicator Indicator describing the target domain
   * \param Pop       Super functor, target dimension DESCRIPTOR::q
   **/
  void definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                         SuperF<DESCRIPTOR::d,T,T>& Pop);
  /// Define a population on a domain with a particular material number
  void definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                         SuperF<DESCRIPTOR::d,T,T>& Pop);

  /// Define an external field on a domain described by an indicator
  template <typename FIELD>
  void defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                   FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field);
  /// Defines a field on a domain described by an indicator
  template <typename FIELD>
  void defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                   AnalyticalF<DESCRIPTOR::d,T,T>& field);
  /// Define an external field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                   FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field);
  /// Defines a field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                   AnalyticalF<DESCRIPTOR::d,T,T>& field);
  /// Defines a field on a indicated domain
  template <typename FIELD>
  void defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                   IndicatorF<T,DESCRIPTOR::d>& indicator,
                   AnalyticalF<DESCRIPTOR::d,T,T>& field);
  /// Defines a field on a indicated domain
  template <typename FIELD>
  void defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                   IndicatorF<T,DESCRIPTOR::d>& indicator,
                   FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field);

  /// Update PARAMETER in all dynamics and post processors
  /**
   * e.g. SuperLattice::setParameter<OMEGA>(0.6) sets the parameter
   * OMEGA to 0.6 in all operators that declare it as one of their
   * parameters.
   **/
  template <typename PARAMETER>
  void setParameter(FieldD<T,DESCRIPTOR,PARAMETER> field);

  /// Update PARAMETER in DYNAMICS
  template <typename PARAMETER, typename DYNAMICS>
  void setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field);
  /// Update PARAMETER in DYNAMICS<T,DESCRIPTOR>
  template <typename PARAMETER, template<typename...> typename DYNAMICS>
  void setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field);

  /// Set PARAMETER to column pointers of given field array
  /**
   * PARAMETER must store size(FIELD) pointers to the value type of FIELD.
   *
   * Useful for providing dynamically sized field arrays as parameters to
   * certain post processors.
   **/
  template <typename PARAMETER, typename _DESCRIPTOR, typename FIELD>
  void setParameter(AbstractFieldArrayD<T,_DESCRIPTOR,FIELD>& abstractFieldArray);

  /// Initialize by equilibrium on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  void iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);
  /// Initialize by equilibrium on a domain with a particular material number
  void iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);
  void iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, SuperF<DESCRIPTOR::d,T,T>& u);
  void iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, SuperF<DESCRIPTOR::d,T,T>& u);
  /// Initialize by non- and equilibrium on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   * \param pi        Analytical functor (global)
   **/
  void iniRegularized(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u,
                      AnalyticalF<DESCRIPTOR::d,T,T>& pi);
  /// Initialize by non- and equilibrium on a domain with a particular material number
  void iniRegularized(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u,
                      AnalyticalF<DESCRIPTOR::d,T,T>& pi);

  /// Core implementation of a single iteration of the collide and stream loop
  /**
   * 1. Pre-collision communication (optional)
   * 1. Collide interior of all local block lattices
   * 2. Post-collision communicate for inter-block propagation
   * 3. Block local streaming
   * 4. Post-stream communication for boundary conditions / post processors (optional)
   * 5. Execute default post processors on all local block lattices
   * 6. Execute manually scheduled post processing callables (optional)
   * 7. Post-post-process communication (optional)
   * 8. Reset lattice statistics and mark overlap state as unclean
   **/
  void collideAndStream();

  /// Subtract constant offset from the density
  void stripeOffDensityOffset(T offset);

  /// Switch Statistics on (default on)
  void statisticsOn()
  {
    _statisticsEnabled = true;
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->setStatisticsEnabled(true);
    }
  };
  /// Switch Statistics off (default on)
  /**
   * This speeds up the execution time
   **/
  void statisticsOff()
  {
    _statisticsEnabled = false;
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->setStatisticsEnabled(false);
    }
  };

  /// Add a non-local post-processing step
  template <typename STAGE=stage::PostStream>
  void addPostProcessor(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                        PostProcessorPromise<T,DESCRIPTOR>&& promise);
  /// Add a non-local post-processing step
  template <typename STAGE=stage::PostStream>
  void addPostProcessor(PostProcessorPromise<T,DESCRIPTOR>&& promise);

  /// Executes post processors for STAGE
  /**
   * 1. Pre-communication for STAGE
   * 2. Execute post processors for STAGE on all block lattices
   **/
  template <typename STAGE=stage::PostStream>
  void executePostProcessors(STAGE stage=STAGE());

  /// Schedules f for execution during every invocation of STAGE
  template <typename STAGE>
  void addCustomTask(std::function<void()> f) {
    _customTasks[typeid(STAGE)].emplace_back(f);
  }
  /// Executes custom tasks assigned to STAGE
  template <typename STAGE>
  void executeCustomTasks(STAGE stage=STAGE());

  /// Schedule F for one-off background execution to be finished by STAGE
  template <typename STAGE, typename F>
  void scheduleBackgroundTask(F&& f);
  /// Schedule F(iC) for one-off background output of block data
  template <typename F>
  void scheduleBackgroundOutput(F&& f);
  /// Schedule one-off background output of given VTK CONTEXT
  /**
   * CONTEXT must accept tasks as a callable task for which
   * it provides SuperVTMwriter(2,3)D and iT arguments.
   **/
  template <typename CONTEXT>
  void scheduleBackgroundOutputVTK(CONTEXT&& vtkContext);
  /// Block until all background tasks scheduled for STAGE are completed
  template <typename STAGE>
  void waitForBackgroundTasks(STAGE stage=STAGE{});

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  void postLoad() override;

};


}

#endif
