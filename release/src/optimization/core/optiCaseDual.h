/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster, Julius Jeßberger
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
* An OptiCase using Adjoint-based Differentiation
*/


#ifndef OPTI_CASE_DUAL_H
#define OPTI_CASE_DUAL_H

#include "io/xmlReader.h"

#include "optimization/solver/adjointLbSolver.h"
#include "optimization/solver/controller.h"
#include "optimization/solver/objective.h"
#include "optimization/core/optiCase.h"
#include "optimization/core/projection.h"
#include "optimization/solver/serialization.h"

namespace olb {

namespace opti {

template<typename S, unsigned dim>
class GeometrySerializer;

template <
  typename T,
  template<typename,SolverMode> typename SOLVER
>
class DistributedObjective;

/** This class implements the evaluation of the goal functional
 * and its derivatives by using adjoint LBM.
 * The adjoint equations are problem-specific and have been computed for force
 * and porosity optimization so far.
 *
 * Requirements: S is the arithmetic data type
 * SOLVER implements the (primal/ dual) simulation, inherits from AdjointLbSolver
 * An xml file is expected to provide additional information on e.g. simulation
 * and optimization parameters
 */
template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C = std::vector<S>>
class OptiCaseDual : public OptiCase<S,C> {

private:
  mutable OstreamManager                           clout {std::cout, "OptiCaseDual"};
  using descriptor = typename SOLVER<S,SolverMode::Reference>::AdjointLbSolver::DESCRIPTOR;
  static constexpr unsigned dim = descriptor::d;
  static constexpr unsigned fieldDim = CONTROLLED_FIELD::template size<descriptor>();

public:
  bool                                             _verbose {true};
  /// upper limit for the number of control variables (#voxels * field-dimension)
  std::size_t                                      _dimCtrl;

  StartValueType                                   _startValueType {Control};
  std::string                                      _projectionName;

  /// Marks, where there are active control variables
  std::shared_ptr<SuperIndicatorF<S,dim>>          _controlIndicator;

  /// Manages the array of control variables
  Controller<S>*                                   _controller {nullptr};
  UnitConverter<S,descriptor>*                     _converter {nullptr};
  std::shared_ptr<SOLVER<S,SolverMode::Primal>>    _primalSolver;
  std::shared_ptr<SOLVER<S,SolverMode::Dual>>      _dualSolver;

  std::shared_ptr<DistributedObjective<S,SOLVER>>  _objective;

  std::shared_ptr<SuperGeometry<S,dim>>            _primalGeometry;
  // this lattice is not used for simulations, but rather for functor syntax
  std::shared_ptr<SuperLattice<S,descriptor>>      _refLattice;

  std::shared_ptr<GeometrySerializer<S,dim>>       _serializer;
  std::shared_ptr<projection::Base<S>>             _projection;
  /// maps the control to a lattice functor
  std::shared_ptr<SuperLatticeF<S,descriptor>>     _projectedControl;
  /// derivative of _projectionControl
  std::shared_ptr<SuperLatticeF<S,descriptor>>     _dProjectionDcontrol;

  OptiCaseDual(XMLreader const& xml) {
    readFromXML(xml);
    initialize(xml);
    initializeFields();
  }

  ~OptiCaseDual() {
    free();
  }

  void free() {
    delete _converter;
    _converter = nullptr;
    delete _controller;
    _controller = nullptr;
  }

  /// Solve primal problem and evaluate objective
  S evaluateObjective(const C& control, unsigned optiStep=0) override;

  /// Compute derivatives via adjoint problem
  void computeDerivatives(
    const C& control, C& derivatives, unsigned optiStep=0) override;


  void setObjective(std::shared_ptr<DistributedObjective<S,SOLVER>> objective) {
    _objective = objective;
    _objective->initialize(_controller, _serializer);
  }

private:
  void readFromXML(XMLreader const& xml);

  void initialize(XMLreader const& xml);

  void initializeFields();

  void derivativesFromDualSolution(C& derivatives);

  void derivativesFromDualSolutionPointwise(C& derivatives, LatticeR<dim+1> latticeR);
};


}

}

#endif
