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
 * The description of steepest decent optimization algorithm -- header file.
 */


#ifndef OPTIMIZER_STEEPEST_DECENT_H
#define OPTIMIZER_STEEPEST_DECENT_H

#include "optimizerLineSearch.h"



// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

/// Optimization algorithm: SteepestDescent
/** OptimizerSteepestDescent optimizes an optimization problem
 * which is given in OptiCase. OptiCase provides therefore
 * methods to evaluate an object functional and compute derivatives.
 * SteepestDescent searches along the direction of the derivative to find
 * a lower value of the evaluated object functional.
 *
 * This class is not intended to be derived from.
 */

template<typename S, typename C>
class OptimizerSteepestDescent : public OptimizerLineSearch<S,C> {

private:
  mutable OstreamManager clout;


public:
  OptimizerSteepestDescent(
    int dimCtrl, S eps, int maxIt, S lambda, int maxStepAttempts,
    std::string stepCondition, bool verboseOn=true, const std::string fname="",
    const std::string logFileName="", bool withUpperBound=false, S upperBound=S(),
    bool withLowerBound=false, S lowerBound=S(),
    bool vectorBounds=false, S controlEps=S(std::numeric_limits<double>::epsilon() ),
    std::vector<OptimizerLogType> gplotAnalysis = {})
    : OptimizerLineSearch<S,C>(
      dimCtrl, eps, maxIt, lambda, maxStepAttempts, stepCondition, verboseOn, fname, logFileName,
      withUpperBound, upperBound, withLowerBound, lowerBound, vectorBounds, controlEps, true,gplotAnalysis),
      clout(std::cout,"OptimizerSteepestDescent") {};

  virtual void computeDirection()
  {
    // Update _derivative
    // computeDerivative only if not already done by wolfeCondition()
    if (!this->_nextDerivFlag) {
      this->computeDerivatives(this->_control, this->_derivative);
    }
    else for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_derivative[iDim] = this->_nextDerivative[iDim];
      }

    // S normDerivative = util::euklidN(this->_derivative.data(), this->_dimCtrl);
    for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
      this->_direction[iDim] = this->_derivative[iDim]  ;  // /normDerivative;
    }
  };
};


/// Creator Function for Steepest Decent
template<typename S, typename C = std::vector<S>>
OptimizerSteepestDescent<S,C>* createOptimizerSteepestDescent(XMLreader const& params, std::size_t dimCtrl)
{
  OstreamManager clout(std::cout,"createOptimizerSteepestDescent");

  // create variables with default values
  int maxIt = 100;

  //S controlEps = S(std::numeric_limits<double>::epsilon() );
  S controlEps = S(0);
  S eps = S(1.e-10);

  S lambda = 1.;
  int maxStepAttempts = 100;
  std::string stepCondition = "Armijo";

  bool vectorBounds = false;
  bool verboseOn=true;
  std::string fname = "control.dat";
  std::string logFileName = "log.txt";

  bool withUpperBound = false;
  S upperBound = S();
  bool withLowerBound = false;
  S lowerBound = S();

  std::vector<OptimizerLogType> gplotAnalysis  = {};
  std::string gplotAnalysisString = "";

  // Read Values from XML File from area "Optimization"
  params.readOrWarn<int>("Optimization", "MaxIter", "", maxIt);

  params.readOrWarn<S>("Optimization", "Tolerance", "", eps);
  params.readOrWarn<S>("Optimization", "ControlTolerance", "", controlEps);
  params.readOrWarn<S>("Optimization", "Lambda", "", lambda);
  params.readOrWarn<int>("Optimization", "MaxStepAttempts", "", maxStepAttempts);
  params.readOrWarn<std::string>("Optimization", "StepCondition", "", stepCondition);

  params.readOrWarn<bool>("Optimization", "Verbose", "", verboseOn);
  params.readOrWarn<std::string>("Optimization", "InputFileName", "", fname);
  params.readOrWarn<std::string>("Optimization", "LogFileName", "", logFileName);

  params.readOrWarn<bool>("Optimization", "VectorBounds", "", vectorBounds);
  if ( params.readOrWarn<S>("Optimization", "UpperBound", "", upperBound, false, false) ) {
    withUpperBound = true;
  }

  if ( params.readOrWarn<S>("Optimization", "LowerBound", "", lowerBound, false, false) ) {
    withLowerBound = true;
  }

  // get the parameters for the gnuplot Analysis from the xml file from the VisualizationGnuplot area
  params.readOrWarn<std::string>("Optimization", "VisualizationGnuplot", "VisualizedParameters", "", gplotAnalysisString, false, false);
  // transform the data from the xml file to the enums needed for continueing
  getGnuplotTagsFromString(gplotAnalysisString, gplotAnalysis);

  // Create Optimizer Object
  return new OptimizerSteepestDescent<S,C>(dimCtrl, eps, maxIt, lambda, maxStepAttempts, stepCondition,
                                         verboseOn, fname, logFileName, withUpperBound, upperBound, withLowerBound, lowerBound,
                                         vectorBounds, controlEps, gplotAnalysis);
}

} // namespace opti

} // namespace olb

#endif
