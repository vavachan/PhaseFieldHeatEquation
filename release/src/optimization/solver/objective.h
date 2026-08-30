/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Julius Jessberger
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

#ifndef OBJECTIVE_H
#define OBJECTIVE_H


#include "optimization/solver/adjointLbSolver.h"
#include "optimization/functors/dualFunctors3D.h"
#include "optimization/solver/serialization.h"
#include "utilities/aliases.h"


namespace olb {

namespace opti {

/// @brief  Objective in optimization with adjoint LBM
/// @tparam T floating point type
// Provides interface for evaluation as well as relevant partial derivatives
template <
  typename T,
  template<typename,SolverMode> typename SOLVER
>
class DistributedObjective
{
protected:
  static constexpr unsigned dim = SOLVER<T,SolverMode::Primal>::dim;
  std::shared_ptr<SOLVER<T,SolverMode::Primal>>     _primalSolver;
  std::shared_ptr<SuperIndicatorF<T,dim>>           _objectiveDomain;
  std::shared_ptr<SuperIndicatorF<T,dim>>           _designDomain;

  Controller<T>*                                    _controller {nullptr};
  std::shared_ptr<GeometrySerializer<T,dim>>        _serializer;

public:
  virtual T evaluate() = 0;

  virtual std::shared_ptr<SuperF<dim,T,T>> derivativeByPopulations() = 0;

  virtual std::shared_ptr<SuperF<dim,T,T>> derivativeByControl() = 0;

  void setPrimalSolver(std::shared_ptr<SOLVER<T,SolverMode::Primal>> primalSolver) {
    _primalSolver = primalSolver;
  }

  virtual void initialize(Controller<T>* controller,
    std::shared_ptr<GeometrySerializer<T,dim>> serializer) {
    _controller = controller;
    _serializer = serializer;
  }
};

/// @brief Objective in optimization with adjoint LBM
/// @tparam T floating point type
/// @tparam OBJECTIVE details for evaluation are implemented here (by user)
/// @tparam FieldDim dimension of controlled field
// the partial derivatives of objective w.r.t. populations and control variables
// are implemented in this class, generically with forward AD
template <
  typename T,
  template<typename,SolverMode> typename SOLVER,
  typename OBJECTIVE,
  template<typename...> typename PRIMAL_DYNAMICS,
  unsigned FieldDim
>
class GenericObjective : public DistributedObjective<T,SOLVER> {
protected:
  using DESCRIPTOR = SOLVER<T,SolverMode::Primal>::DESCRIPTOR;
  using DistributedObjective<T,SOLVER>::dim;

  // helper structure, where the details for evaluation are provided (by user)
  std::shared_ptr<OBJECTIVE>                       _objective;
  T                                                _cellVolume;

public:
  GenericObjective() = default;

  GenericObjective(std::shared_ptr<OBJECTIVE> objective)
   : _objective(objective)
  {
    this->_objectiveDomain = objective->_objectiveDomain;
    this->_designDomain = objective->_designDomain;
    _cellVolume = util::pow(
      this->_objectiveDomain->getSuperGeometry().getCuboidDecomposition().getDeltaR(),
      dim);
  }

  T evaluate() override {
    // main term that depends on state
    SuperLatticeFfromCallableF<T,DESCRIPTOR> jF(
      this->_primalSolver->lattice(),
      [this](T* output, auto cell, int iC, const int* coords) {
        _objective->j(output, cell, iC, coords);
      });
    SuperIntegral<dim,T> J(&jF, this->_objectiveDomain);

    // additional term that depends (directly) on control
    SuperLatticeFfromCallableF<T,DESCRIPTOR> rF(
      this->_primalSolver->lattice(),
      [this](T* output, int iC, const int* coords) {
        const T* control = getControlPointer(iC, coords);
        _objective->r(output, iC, coords, control);
      });
    SuperIntegral<dim,T> R(&rF, this->_designDomain);

    int dummy[dim+1]{0}; T result{0}; T result2{0};
    J(&result, dummy);
    R(&result2, dummy);
    return result + result2;
  }

  std::shared_ptr<SuperF<dim,T,T>> derivativeByPopulations() override {
    using T_AD1 = util::ADf<T,DESCRIPTOR::q>;
    std::shared_ptr<SuperF<dim,T,T>> djdf = std::make_shared<SuperLatticeFfromCallableF<T,DESCRIPTOR>> (
      this->_primalSolver->lattice(),
      [this](T* output, auto cell, int iC, const int* coords){
        const int iCglob = this->_primalSolver->lattice().getLoadBalancer().glob(iC);
        bool isInside;
        if constexpr (DESCRIPTOR::d == 3) {
          isInside = this->_objectiveDomain->operator()(iCglob,coords[0],coords[1],coords[2]);
        } else {
          isInside = this->_objectiveDomain->operator()(iCglob,coords[0],coords[1]);
        }
        if (isInside) {
          // get ADf-typed clone of cell
          const Vector<int,dim> extent(1);
          ConcreteBlockLattice<T_AD1,DESCRIPTOR> blockLatticeAD(extent, 0);
          const Vector<int,dim> position(0);
          blockLatticeAD.template defineDynamics<PRIMAL_DYNAMICS>(position);
          auto adFcell = blockLatticeAD.get(position);
          DESCRIPTOR::fields_t::for_each([&](auto id) {
            using field = typename decltype(id)::type;
            adFcell.template setField<field>(cell.template getField<field>());
          });

          for (unsigned iPop=0; iPop<DESCRIPTOR::q; ++iPop){
            adFcell.template getFieldPointer<descriptors::POPULATION>()[iPop].setDiffVariable(iPop);
          }
          T_AD1 result;
          _objective->j(&result, adFcell, iC, coords);
          for (unsigned iPop=0; iPop<DESCRIPTOR::q; ++iPop){
            output[iPop] = _cellVolume * result.d(iPop);
          }
        }
        else {
          for (unsigned iPop=0; iPop<DESCRIPTOR::q; ++iPop){
            output[iPop] = T(0);
          }
        }
    });
    return djdf;
  }

  std::shared_ptr<SuperF<dim,T,T>> derivativeByControl() override {
    using T_AD2 = util::ADf<T,FieldDim>;
    std::shared_ptr<SuperF<dim,T,T>> drdc = std::make_shared<SuperLatticeFfromCallableF<T,DESCRIPTOR>> (
      this->_primalSolver->lattice(),
      [this](T* output, int iC, const int* coords){
        const int iCglob = this->_primalSolver->lattice().getLoadBalancer().glob(iC);
        if constexpr (DESCRIPTOR::d == 3) {
          if (this->_designDomain->operator()(iCglob,coords[0],coords[1],coords[2])) {
            // get ADf-typed copy of control
            const T* control = getControlPointer(iC, coords);
            Vector<T_AD2,FieldDim> adFcontrol(control);
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              adFcontrol[iControl].setDiffVariable(iControl);
            }

            T_AD2 result;
            _objective->r(&result, iC, coords, adFcontrol.data());
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              output[iControl] = _cellVolume * result.d(iControl);
            }
          }
          else {
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              output[iControl] = T(0);
            }
          }
        }
        else {
          if (this->_designDomain->operator()(iCglob,coords[0],coords[1])) {
            // get ADf-typed copy of control
            const T* control = getControlPointer(iC, coords);
            Vector<T_AD2,FieldDim> adFcontrol(control);
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              adFcontrol[iControl].setDiffVariable(iControl);
            }

            T_AD2 result;
            _objective->r(&result, iC, coords, adFcontrol.data());
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              output[iControl] = _cellVolume * result.d(iControl);
            }
          }
          else {
            for (unsigned iControl=0; iControl<FieldDim; ++iControl){
              output[iControl] = T(0);
            }
          }
        }
    });
    return drdc;
  }

protected:
  /// @brief Access control at some position
  /// @param iC local cuboid id
  /// @param coords spatial coordinates
  /// @return pointer to control entry (the first one, if control is multi-dimensional)
  const T* getControlPointer(int iC, const int* coords) const {
    int latticeR[dim+1] = {
      this->_primalSolver->lattice().getLoadBalancer().glob(iC),
      coords[0],
      coords[1]};
    if constexpr (dim == 3) {
      latticeR[3] = coords[2];
    }
    const auto index = this->_serializer->getSerializedCellIndex(latticeR);
    return this->_controller->getControlPointer(index);
  }
};


}
}

#endif
