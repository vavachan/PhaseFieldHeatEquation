/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Mathias J. Krause, Benjamin Förster
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
 * The description of line search optimization algorithm -- header file.
 */


#ifndef OPTIMIZER_LINE_SEARCH_H
#define OPTIMIZER_LINE_SEARCH_H

#include "optimizer.h"

// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

/// Optimization algorithm: LineSearch.
/** OptimizerLineSearch optimizes an optimization problem
 * which is given in OptiCase. OptiCase provides therefore
 * methods to evaluate an object functional and compute derivatives.
 * LineSearch searches along a particular direction to find
 * a lower value of the evaluated object functional.
 *
 * This class is intended to be derived from.
 */

template<typename S, typename C>
class OptimizerLineSearch : public Optimizer<S,C> {

private:
  mutable OstreamManager clout;

protected:
  /// Lambda start value
  S _lambda;
  /// Search direction
  C _direction;
  // Upper/Lower bound
  bool _lowerBoundFlag;
  bool _upperBoundFlag;
  /// Maximal number of step attempts for conditioned line search
  int _maxStepAttempts;

  C _nextDerivative;
  bool _nextDerivFlag;

  // Step condition
  std::string _stepCondition;
  bool (OptimizerLineSearch::*_stepConditionFunction)(const S&);
  void (OptimizerLineSearch::*_stepLengthFunction)(const S&);

public:
  /// Construction of an OptimizerLineSearch
  OptimizerLineSearch(int dimCtrl, S eps, int maxIt, S lambda, int maxStepAttempts, std::string stepCondition,
                      bool verboseOn=true, const std::string fname="", const std::string logFileName="",
                      bool withUpperBound=false, S upperBound=S(),
                      bool withLowerBound=false, S lowerBound=S(), bool vectorBounds=false,
                      S controlEps=S(std::numeric_limits<double>::epsilon() ), bool failOnMaxIter = true,
                      std::vector<OptimizerLogType> gplotAnalysis = {})
    : Optimizer<S,C>(dimCtrl, eps, maxIt, verboseOn, fname, logFileName, withUpperBound, upperBound, withLowerBound,
                   lowerBound, vectorBounds, controlEps, failOnMaxIter, gplotAnalysis),
      clout(std::cout,"OptimizerLineSearch")
  {
    _lambda = lambda;
    _direction = util::ContainerCreator<C>::create(dimCtrl);
    _nextDerivative = util::ContainerCreator<C>::create(dimCtrl);
    _maxStepAttempts = maxStepAttempts;
    _stepCondition = stepCondition;
    _lowerBoundFlag = false;
    _upperBoundFlag = false;
    _nextDerivFlag = false;

    /// Define step conditions and step length calculation
    // No Condition
    if (_stepCondition == "None") {
      _stepConditionFunction = &OptimizerLineSearch::noCondition;
      _stepLengthFunction = nullptr;
      // Smaller objective function
    }
    else if (_stepCondition == "Smaller") {
      _stepConditionFunction = &OptimizerLineSearch::smallerValue;
      _stepLengthFunction = &OptimizerLineSearch::quadraticInterpolationStep;
      // Armijo, Wolfe and Strong Wolfe conditions
    }
    else if (_stepCondition == "Armijo" || _stepCondition == "Wolfe" || _stepCondition == "StrongWolfe") {
      _stepConditionFunction = &OptimizerLineSearch::armijoWolfeConditions;
      _stepLengthFunction = &OptimizerLineSearch::quadraticInterpolationStep;
    }
    else {
      clout << "No step condition chosen! Exiting..." << std::endl;
      exit(1);
    }

  };

  virtual ~OptimizerLineSearch() { };

  virtual void computeDirection() = 0;

  void checkBound()
  {
    if (this->_withUpperBound) {
      bool atUpper = true;
      for (int iDim=0; iDim<this->_dimCtrl; ++iDim) {
        const int i = (this->_vectorBounds) ? iDim : 0;
        atUpper &= (util::nearZero(this->_boundedControl[iDim] - this->_upperBound[i]));
      }
      if (atUpper) {
        if (_upperBoundFlag) {
          clout << "Control stuck at upper bound: " <<  std::string(this->_upperBound.begin(), this->_upperBound.end()) << std::endl;
          exit(1);
        }
        _upperBoundFlag = true;
      }
      else {
        _upperBoundFlag = false;
      }
    }
    if (this->_withLowerBound) {
      bool atLower = true;
      for (int iDim=0; iDim<this->_dimCtrl; ++iDim) {
        const int i = (this->_vectorBounds) ? iDim : 0;
        atLower &= (util::nearZero(this->_boundedControl[iDim] - this->_lowerBound[i]));
      }
      if (atLower) {
        if (_lowerBoundFlag) {
          clout << "Control stuck at lower bound: " << std::string(this->_lowerBound.begin(), this->_lowerBound.end()) << std::endl;
          exit(1);
        }
        _lowerBoundFlag = true;
      }
      else {
        _lowerBoundFlag = false;
      }
    }
  };

  void boundControl()
  {
    for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
      this->_boundedControl[iDim] = this->_control[iDim];
      const int i = (this->_vectorBounds) ? iDim : 0;
      if (this->_withUpperBound) {
        const S tmp = this->_upperBound[i];
        if (this->_control[iDim] > tmp ||
            (std::isnan(this->_control[iDim])
             && !std::signbit(this->_control[iDim])) ) {
          this->_boundedControl[iDim] = tmp;
        }
      }
      if (this->_withLowerBound) {
        const S tmp = this->_lowerBound[i];
        if (this->_control[iDim] < tmp ||
            (std::isnan(this->_control[iDim])
             && std::signbit(this->_control[iDim])) ) {
          this->_boundedControl[iDim] = tmp;
        }
      }
    }
  };

  bool smallerValue(const S& tempValue)
  {
    if (!(tempValue > this->_value || std::isnan(tempValue))) {
      return true;
    }
    else {
      (this->*_stepLengthFunction)(tempValue);
      return false;
    }
  }

  bool noCondition(const S& tempValue)
  {
    return true;
  }

  // Combined Armijo, Wolfe and Strong Wolfe condition
  // Depends on _stepCondition
  bool armijoWolfeConditions(const S& tempValue)
  {
    S c1 = 1e-4;
    S c2 = 0.9; // 0.9 for (quasi)Newton; else 0.1
    S dir = 0.;
    S dirNext = 0.;

    // \nabla f_k \dot p_k
    for (int i=0; i<this->_dimCtrl; i++) {
      dir += this->_direction[i] * this->_derivative[i];
    }

    // Armijo rule
    if (!(tempValue <= this->_value + c1*this->_lambda*dir)) {
      if (this->_verboseOn) {
        clout << "Armijo failed!" << std::endl;
      }
      // Decrease step length
      (this->*_stepLengthFunction)(tempValue);
      //this->_lambda *= .5;
      return false;
    }
    // (strong) Wolfe conditions (curvature condition)
    if ( _stepCondition == "Wolfe" || _stepCondition == "StrongWolfe") {
      // Compute derivative for the changed control and store it in nextDerivative
      // Set nextDerivFlag to true to prevent the optiCase to redo this process
      this->computeDerivatives(this->_control, _nextDerivative);
      _nextDerivFlag = true;
      // \nabla f_{k+1} \dot p_k
      for (int i=0; i<this->_dimCtrl; i++) {
        dirNext += this->_direction[i]*_nextDerivative[i];
      }
      bool curvature = ( dirNext <= c2*dir );
      if (_stepCondition == "StrongWolfe") {
        curvature = ( util::abs(dirNext) <= c2*util::abs(dir) );
      }

      if (!curvature) {
        if (this->_verboseOn) {
          clout << "Curvature failed!" << std::endl;
        }
        // Increase step length
        this->_lambda *= 2.1;
        return false;
      }
    }
    return true;
  }

  // Quadratic Interpolation [Nocedal; p.58]
  void quadraticInterpolationStep(const S& tempValue)
  {
    S newLambda = S();
    S dir = S();
    if (std::isnan(tempValue) ) {
      newLambda = this->_lambda/2.;
    }
    else {
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        dir += this->_derivative[iDim]*this->_direction[iDim];
      }
      newLambda = dir*this->_lambda*this->_lambda / (2.*(tempValue - this->_value + dir*this->_lambda));
      // Step size control
      if (newLambda < this->_lambda*.1 ) {
        newLambda = this->_lambda*.1;
      }
      if (newLambda > this->_lambda*.5 ) {
        newLambda = this->_lambda*.5;
      }
    }
    this->_lambda = newLambda;
  }

  void backtrackingLineSearch(S& tempValue, S lambda, bool(OptimizerLineSearch::*condition)(const S&))
  {

    // Save this->_lambda so that it is unchanged after the line search
    S startLambda = this->_lambda;
    int refinementStep = 0;
    // Do line search until the condition is fullfilled
    // the step size (this->_lambda) will be changed in the condition
    while ( !(this->*condition)(tempValue) ) {
      if ( util::abs(lambda/this->_value) < std::numeric_limits<double>::epsilon() ) {
        clout << "Excessive refinement steps (too small step size).\nProgram terminated." <<std::endl;
        exit(1);
      }
      refinementStep++;

      // Leave program if maximum number of step attempts is exceeded.
      if ( refinementStep >= _maxStepAttempts ) {
        clout << "Excessive refinement steps (maxStepAttempts exeeded).\nProgram terminated." <<std::endl;
        exit(1);
      }

      S newLambda = this->_lambda;
      bool notSensible = true;
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        // Go back (newLambda-lambda) from previous control(=control-lambda*direction),
        // therefore this produces a step of newLambda*direction
        S tmp  = (newLambda-lambda)*this->_direction[iDim];
        this->_control[iDim] -= tmp;
        // is this move larger than machine precision wrt control
        if (util::abs(this->_control[iDim]) > 0 && util::abs(tmp) > 0) {
          notSensible &= ( util::nearZero(tmp/this->_control[iDim]) );
        }
      }
      // stop excessive refinement steps when values become no longer sensible
      if (notSensible) {
        clout << "Excessive refinement steps (not sensible).\nProgram terminated." <<std::endl;
        exit(1);
      }

      if (this->_withUpperBound||this->_withLowerBound) {
        boundControl();
      }

      if (this->_verboseOn) {
        clout << "[Step " << this->_it << "][Ref " << refinementStep << "] <<<<<<<<<< lambda=" << newLambda << " <<<<<<<<<<" << std::endl;
      }

      if (this->_withUpperBound||this->_withLowerBound) {
        this->evaluateObjective(this->_boundedControl, tempValue);
      }
      else {
        this->evaluateObjective(this->_control, tempValue);
      }
      lambda = newLambda;
    }
    this->_lambda = startLambda;
  };


  /// Optimization step: line search
  virtual void optimizationStep()
  {
    // Compute search direction
    if (this->_verboseOn) {
      clout << "Computing directions..." << std::endl;
    }

    // Use given optimization algorithm to get direction:
    // Steepest Descent, LBFGS or Barzilai-Borwein
    computeDirection();

    _upperBoundFlag = false;
    _lowerBoundFlag = false;
    checkBound();

    // save initial control
    const C initialControl = this->_control;

    // Search along line to find new control
    for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
      this->_control[iDim] -= this->_lambda*_direction[iDim];
    }

    if (this->_withUpperBound||this->_withLowerBound) {
      boundControl();
      checkBound();
    }

    if (this->_verboseOn) {
      clout << "[Step " << this->_it << "] <<<<<<<<<< lambda=" << this->_lambda << " <<<<<<<<<<" << std::endl;
    }

    S tempValue = S();
    if (this->_withUpperBound||this->_withLowerBound) {
      this->evaluateObjective(this->_boundedControl, tempValue);
    }
    else {
      this->evaluateObjective(this->_control, tempValue);
    }

    // Backtracking line search with step condition
    backtrackingLineSearch(tempValue, this->_lambda, _stepConditionFunction);

    if (this->_withUpperBound||this->_withLowerBound) {
      if (this->_verboseOn) {
        clout << "Bounding control" << std::endl;
      }
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_control[iDim] = this->_boundedControl[iDim];
      }
    }

    // Update value of the objective functional
    this->_value = tempValue;
    // Update step no.
    this->_it++;

    // check if control have changed sufficiently or break from optimisation
    this->_controlsConverged = true;
    for (int iDim=0; iDim<this->_dimCtrl; ++iDim) {
      S ave = 0.5 * (initialControl[iDim] + this->_control[iDim]);
      if (util::abs(ave) > 0) {
        if (util::abs( (initialControl[iDim] - this->_control[iDim]) / ave) >= this->_controlEps) {
          this->_controlsConverged = false;
          break;
        }
      }  // else control diff is zero anyway
    }
    if (this->_verboseOn) {
      clout << "Controls converged within " << this->_controlEps << ": " << ((this->_controlsConverged) ? "true" : "false") << std::endl;
    }
  };
};

} // namespace opti

} // namespace olb

#endif
