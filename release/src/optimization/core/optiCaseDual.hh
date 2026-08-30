/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012-2025 Mathias J. Krause, Benjamin Förster, Julius Jeßberger
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


#ifndef OPTI_CASE_DUAL_HH
#define OPTI_CASE_DUAL_HH

#include "optimization/core/optiCaseDual.h"
#include "optimization/solver/serialization.h"

namespace olb {

namespace opti {

template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::readFromXML(XMLreader const& xml)
{
  xml.readOrWarn<std::string>("Optimization", "Projection", "", _projectionName);

  std::string type ("");
  xml.readOrWarn<std::string>("Optimization", "StartValueType", "", type);
  if (type == "Porosity") {
    _startValueType = Porosity;
  } else if (type == "Permeability") {
    _startValueType = Permeability;
  } else if (type == "Control") {
    _startValueType = Control;
  } else {
    _startValueType = ProjectedControl;
  }
}

template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::initialize(XMLreader const& xml)
{
  _converter = createUnitConverter<S,descriptor>(xml);

  _primalSolver = createLbSolver <SOLVER<S,SolverMode::Primal>> (xml);
  _dualSolver = createLbSolver <SOLVER<S,SolverMode::Dual>> (xml);
  this->_postEvaluation = std::bind(&SOLVER<S,SolverMode::Primal>::postProcess, _primalSolver);

  _primalSolver->initialize();
  _controlIndicator = _primalSolver->parameters(names::Opti()).designDomain;
  _primalGeometry = _primalSolver->parameters(names::Results()).geometry;
  _refLattice = std::make_shared<SuperLattice<S,descriptor>>(*_primalGeometry);

  _serializer = std::make_shared<SimpleGeometrySerializer<S,dim>>(*_primalGeometry);
  _dimCtrl = _serializer->getNoCells() * fieldDim;

  _controller = new Controller<S>(_dimCtrl);
}

template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::initializeFields()
{
  _projection = projection::construct<S,descriptor>(*_converter, _projectionName);

  _projectedControl = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_refLattice,
    *_controller,
    fieldDim,
    _serializer,
    [this](S x) { return _projection->project(x); });
  _dProjectionDcontrol = std::make_shared<SuperLatticeSerialDataF<S,descriptor>>(
    *_refLattice,
    *_controller,
    fieldDim,
    _serializer,
    [this](S x) { return _projection->derivative(x); });
  _primalSolver->parameters(names::Opti()).controlledField = _projectedControl;
  _dualSolver->parameters(names::Opti()).controlledField = _projectedControl;
}


template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
S OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::evaluateObjective(
  const C& control, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _primalSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;
  _primalSolver->solve();

  _objective->setPrimalSolver(_primalSolver);
  return _objective->evaluate();
}


template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::computeDerivatives(
  const C& control, C& derivatives, unsigned optiStep)
{
  _controller->setControl(control, _dimCtrl);
  _dualSolver->parameters(names::OutputOpti()).counterOptiStep = optiStep;

  const auto& primalResults = _primalSolver->parameters(names::Results());
  auto& dualParams = _dualSolver->parameters(names::Opti());

  dualParams.fpop = std::make_shared<SuperLatticeFpop<S,descriptor>>(*(primalResults.lattice));
  dualParams.dObjectiveDf = _objective->derivativeByPopulations();

  _dualSolver->solve();

  derivativesFromDualSolution(derivatives);
}

template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::derivativesFromDualSolution(
  C& derivatives)
{
  const auto primalGeometry = _primalSolver->parameters(names::Results()).geometry;
  const int nC = primalGeometry->getCuboidDecomposition().size();

  for (int iC = 0; iC < nC; iC++) {

    const Vector<int,descriptor::d> extend = primalGeometry->getCuboidDecomposition().get(iC).getExtent();

    if constexpr (descriptor::d == 3) {
      for (int iX = 0; iX < extend[0]; iX++) {
        for (int iY = 0; iY < extend[1]; iY++) {
          for (int iZ = 0; iZ < extend[2]; iZ++) {
            const LatticeR<dim+1> latticeR(iC, iX, iY, iZ);
            derivativesFromDualSolutionPointwise(derivatives, latticeR);
          }
        }
      }
    }
    else if (descriptor::d == 2) {
      for (int iX = 0; iX < extend[0]; iX++) {
        for (int iY = 0; iY < extend[1]; iY++) {
          const LatticeR<dim+1> latticeR(iC, iX, iY);
          derivativesFromDualSolutionPointwise(derivatives, latticeR);
        }
      }
    }
  }
}

/// Compute the gradient dJ/d(control(latticeR)) for a certain spatial location latticeR.
/// @param derivatives the global gradient vector. Result is written in here.
/// @param latticeR the spatial location (in global lattice coordinates).
template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
void OptiCaseDual<S,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>::derivativesFromDualSolutionPointwise(
  C& derivatives,
  LatticeR<dim+1> latticeR)
{
  C result(fieldDim, 0);

  if (evaluateSuperIndicatorFglobally<dim,S>(*_controlIndicator, latticeR.data())) {

    const auto primalGeometry = _primalSolver->parameters(names::Results()).geometry;
    if (primalGeometry->getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank()) {
      S dProjectionDcontrol[fieldDim];
      (*_dProjectionDcontrol)(dProjectionDcontrol, latticeR.data());

      // compute derivative coming from dependence of objective on state -> use adjoint solution

      // construct and initialize dummy ADf cell
      const auto iCloc = primalGeometry->getLoadBalancer().loc(latticeR[0]);
      auto& block = _primalSolver->parameters(names::Results()).lattice->getBlock(iCloc);
      const LatticeR<descriptor::d> localCoords(&latticeR[1]);
      auto cell = block.get(localCoords);
      using S_AD = util::ADf<S,fieldDim>;
      const Vector<int,dim> extent(1);
      ConcreteBlockLattice<S_AD,descriptor> blockLatticeAD(extent, 0);
      blockLatticeAD.setStatisticsEnabled(false);
      const Vector<int,dim> position(0);
      blockLatticeAD.template defineDynamics<PRIMAL_DYNAMICS>(position);
      auto adFcell = blockLatticeAD.get(position);
      descriptor::fields_t::for_each([&](auto id) {
        using field = typename decltype(id)::type;
        adFcell.template setField<field>(cell.template getField<field>());
      });
      for (unsigned iCtrl=0; iCtrl<fieldDim; ++iCtrl) {
        adFcell.template getFieldPointer<CONTROLLED_FIELD>()[iCtrl].setDiffVariable(iCtrl);
      }

      // execute primal collision
      using AdParameters = ParametersOfDynamicsD<PRIMAL_DYNAMICS<S_AD,descriptor>>;
      AdParameters adParams;
      AdParameters::fields_t::for_each([&](auto id) {
        using field = typename decltype(id)::type;
        using AdFieldD = FieldD<S_AD,descriptor,field>;
        const AdFieldD adField(cell.template getField<field>());
        adParams.template set<field>(adField);
      });
      const S omega = _converter->getLatticeRelaxationFrequency();
      adParams.template set<descriptors::OMEGA>(S_AD(omega));
      PRIMAL_DYNAMICS<S,descriptor>().collide(adFcell, adParams);

      // return derivatives
      for (unsigned iCtrl=0; iCtrl<fieldDim; ++iCtrl) {
        for (unsigned jPop=0; jPop<descriptor::q; ++jPop) {
          const S phi_j = _dualSolver->parameters(names::Results()).lattice->get(latticeR)[jPop];
          result[iCtrl] -= phi_j * adFcell[jPop].d(iCtrl) * dProjectionDcontrol[iCtrl];
        }
      }

      // add (partial) derivative coming from direct dependence of objective on control
      S dObjectiveDcontrol[fieldDim];
      (*_objective->derivativeByControl())(dObjectiveDcontrol, latticeR.data());
      for (unsigned iDim=0; iDim<fieldDim; iDim++) {
        result[iDim] += dObjectiveDcontrol[iDim];  // * dProjectionDcontrol[iDim];  // todo uncomment if objective depends on projection
      }
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&result[0], fieldDim, primalGeometry->getLoadBalancer().rank(latticeR[0]));
#endif

    for (unsigned iDim=0; iDim<fieldDim; ++iDim) {
      const auto index = _serializer->getSerializedComponentIndex(latticeR, iDim, fieldDim);
      derivatives[index] = result[iDim];
    }
  }
}

} // namespace opti

} // namespace olb

#endif
