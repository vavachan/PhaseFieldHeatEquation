/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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


#ifndef OPTI_SOLVER_PARAMETERS_H
#define OPTI_SOLVER_PARAMETERS_H


#include <limits>

#include "optimization/solver/adjointLbSolver.h"
#include "solver/solverParameters.h"

namespace olb {

namespace opti {

enum class SolverMode;

}

namespace parameters {


struct OptiSimulationBase : public ParameterBase { };

struct OptiResultsBase : public ParameterBase { };


// ---------------------- Output parameters -----------------------------------


template <typename T>
struct OptiOutput : public ParameterBase {
  std::size_t                       counterOptiStep {0};
};


template<typename T, typename TAG>
struct Reader<OptiOutput<T>, TAG> : public ReaderBase<OptiOutput<T>> {
  using ReaderBase<OptiOutput<T>>::ReaderBase;

  void read(XMLreader const& xml)
  { }
};


// ---------------------- For ad and dq optimization --------------------------

template<typename T>
struct DirectOptiSimulation : public OptiSimulationBase {

  // user has to define how the control values are applied in the simulation
  virtual void applyControl(const std::vector<T>& control) = 0;
};

template<typename T>
struct DirectOptiResults : public OptiResultsBase {

  using BT = BaseType<T>;
  // Value of goal functional. Optimization seeks to minimize it
  T objective {std::numeric_limits<BT>::infinity()};
};


// ------------------- Parameters for adjoint optimization --------------------

// Parameters for distributed optimization problems, for use in combination
// with OptiCaseDual.
// Warning: only one lattice is supported so far.
template<typename T, typename LATTICES>
struct DistributedOptiSimulationBase : public OptiSimulationBase {

  using descriptor = typename LATTICES::values_t::template get<0>;

  int                                   fieldDim {descriptor::d};
  // This is set automatically by OptiCaseDual
  std::shared_ptr<SuperLatticeF<T,descriptor>> controlledField;
  // Fluid cells without Dirichlet boundaries. This is set by user in solver::prepareGeometry
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> bulkIndicator;
  // Where the control variables are active. This is set by user in solver::prepareGeometry
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> designDomain;
  // Where the objective is computed. This is set by user in solver::prepareGeometry
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> objectiveDomain;
};


template<typename T, typename LATTICES, opti::SolverMode MODE=opti::SolverMode::Reference>
struct DistributedOptiSimulation
  : public DistributedOptiSimulationBase<T,LATTICES> { };

template<typename T, typename LATTICES>
struct DistributedOptiSimulation<T,LATTICES,opti::SolverMode::Dual>
  : public DistributedOptiSimulationBase<T,LATTICES> {

  using descriptor = typename LATTICES::values_t::template get<0>;
  static constexpr unsigned                    dim = descriptor::d;

  std::shared_ptr<SuperLatticeF<T,descriptor>> fpop;          // this is set by OptiCaseDual
  std::shared_ptr<SuperF<descriptor::d,T,T>>  dObjectiveDf;  // this is set by OptiCaseDual
};


/// xml interface for DistributedOptiSimulation parameters
template<typename T, typename LATTICES, opti::SolverMode MODE, typename TAG>
struct Reader<DistributedOptiSimulation<T,LATTICES,MODE>, TAG>
      : public ReaderBase<DistributedOptiSimulation<T,LATTICES,MODE>> {
  using ReaderBase<DistributedOptiSimulation<T,LATTICES,MODE>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    xml.readOrWarn<int>("Optimization", "FieldDimension", "",
                        this->params->fieldDim, true, false, true);
  }
};


// ------------------- Results of adjoint optimization ------------------------

template<typename T, typename LATTICES>
struct DistributedOptiSimulationResults : public OptiResultsBase
{
  using descriptor = typename LATTICES::values_t::template get<0>;
  std::shared_ptr<SuperGeometry<T,descriptor::d>>   geometry;  // this is set by user in solver::prepareGeometry
  std::shared_ptr<SuperLattice<T,descriptor>>       lattice;   // this is set in AdjointLbSolver::computeResults
};


}  // namespace parameters

}  // namespace olb


#endif
