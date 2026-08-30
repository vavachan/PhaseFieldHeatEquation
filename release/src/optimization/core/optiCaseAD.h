/*  This file is part of the OpenLB library
*
*  Copyright (C) 2012, 2021 Mathias J. Krause, Julius Jeßberger
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
* An OptiCase using Algorithmic Differentiation (AD)
*/


#ifndef OPTI_CASE_AD_H
#define OPTI_CASE_AD_H


#include "utilities/aDiff.h"
#include "optimization/core/optiCase.h"
#include "utilities/vectorHelpers.h"


namespace olb {

namespace opti {

/**
 * @brief Derivatives are computed with automatic differentiation
 *
 * @tparam S Fundamental (floating point) datatype
 * @tparam n Number of optimized (control) parameters
 * @tparam C Container type
 */
template<typename S, std::size_t n,
  template<typename> typename C = util::StdVector>
class OptiCaseAD : public OptiCase<S,C<S>> {

protected:
  using T = util::ADf<S,n>;
  std::function<S (const C<S>&)>      _functionHelp   { [](const C<S>&){ return S{}; } };
  std::function<T (const C<T>&)>      _adFunctionHelp { [](const C<T>&){ return T{}; } };
  std::function<S (const C<S>&, unsigned)>      _function;
  std::function<T (const C<T>&, unsigned)>      _adFunction;

public:
  explicit OptiCaseAD() = default;

  OptiCaseAD(
    std::function<S (const C<S>&, unsigned)> function,
    std::function<T (const C<T>&, unsigned)> adFunction,
    std::function<void (void)> postEvaluation = [](){})
   : OptiCase<S,C<S>>(postEvaluation), _function(function), _adFunction(adFunction)
  { }

  // If function with (only) one argument is passed, use plain wrappers
  OptiCaseAD(
    std::function<S (const C<S>&)> function,
    std::function<T (const C<T>&)> adFunction,
    std::function<void (void)> postEvaluation = [](){})
   : OptiCase<S,C<S>>(postEvaluation),
     _functionHelp(function), _adFunctionHelp(adFunction),  // store original functions for memory reasons
     _function ([&](const C<S>& arg, unsigned){ return _functionHelp(arg); }),
     _adFunction ([&](const C<T>& arg, unsigned){ return _adFunctionHelp(arg); })
  { }

  S evaluateObjective(
    const C<S>& control, unsigned optiStep=0) override
  {
    return _function(control, optiStep);
  }

  void computeDerivatives(const C<S>& control,
    C<S>& derivatives, unsigned optiStep=0) override
  {
    assert((control.size() == derivatives.size()));

    auto adControl = util::iniAD<n,S,C>(control);

    const T adResult = _adFunction(adControl, optiStep);

    for(std::size_t it = 0; it < control.size(); ++it){
      derivatives[it] = adResult.d(it);
    }
  }
};


/** Interface for OptiCaseAD that performs Lattice-Boltzmann-Solver
 * construction itself (from xml file)
 */
template<typename S, std::size_t n, template<typename> typename SOLVER,
  template<typename> typename C = util::StdVector>
class OptiCaseAdForSolver : public OptiCaseAD<S,n,C> {

public:
  using typename OptiCaseAD<S,n,C>::T;
  std::shared_ptr<SOLVER<S>>                    _solver;
  std::shared_ptr<SOLVER<T>>                    _adSolver;

  OptiCaseAdForSolver(std::shared_ptr<SOLVER<S>> solver,
    std::shared_ptr<SOLVER<T>> adSolver)
   : OptiCaseAD<S,n,C>(getCallable<S>(solver), getCallable<T>(adSolver),
       std::bind(&SOLVER<S>::postProcess, solver)),
     _solver(solver), _adSolver(adSolver)
  { }

  OptiCaseAdForSolver(XMLreader const& xml)
   : OptiCaseAdForSolver(
      createLbSolver<SOLVER<S>> (xml),
      createLbSolver<SOLVER<T>> (xml))
  { }
};

template<typename S, typename T, template<typename> typename SOLVER>
OptiCaseAdForSolver(std::shared_ptr<SOLVER<S>>,std::shared_ptr<SOLVER<T>>)
  -> OptiCaseAdForSolver<S,T::dim,SOLVER,util::StdVector>;

}

}

#endif
