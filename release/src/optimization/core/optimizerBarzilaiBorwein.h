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
 * The description of Borzilai Borwein optimization algorithm -- header file.
 */


#ifndef OPTIMIZER_BARZILAI_BORWEIN_H
#define OPTIMIZER_BARZILAI_BORWEIN_H

#include "optimizerLineSearch.h"
#include "utilities/norm.h"



// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

/// Optimization algorithm: BarzilaiBorwein.
/** OptimizerBarzilaiBorwein optimizes an optimization problem
 * which is given in OptiCase. OptiCase provides therefore
 * methods to evaluate an object functional and compute derivatives.
 * BarzilaiBorwein searches along the direction of the derivative to find
 * a lower value of the evaluated object functional. Thereby the steplength
 * is adjusted to be an approximation of the Hessian.
 * Barzilai-Borwein is therefore a "quasi" quasi Newton method.
 *
 * This class is not intended to be derived from.
 */

template<typename S, typename C>
class OptimizerBarzilaiBorwein : public OptimizerLineSearch<S,C> {

private:
  mutable OstreamManager clout;

  S _startLambda;

  C _lastControl;
  C _lastDerivative;

  C _sStore;
  C _yStore;


public:
  OptimizerBarzilaiBorwein(
    int dim, S eps, int maxIt, S lambda, int maxStepAttempts, std::string stepCondition,
    bool verboseOn=true, const std::string fname="", const std::string logFileName="",
    bool withUpperBound=false, S upperBound=S(), bool withLowerBound=false, S lowerBound=S(),
    bool vectorBounds=false, S controlEps=S(std::numeric_limits<double>::epsilon() ),
    std::vector<OptimizerLogType> gplotAnalysis = {})
    : OptimizerLineSearch<S,C>(dim, eps, maxIt, /*lambda*/ 1., maxStepAttempts, stepCondition,
      verboseOn, fname, logFileName, withUpperBound, upperBound, withLowerBound,
      lowerBound, vectorBounds, controlEps, true, gplotAnalysis),
      clout(std::cout,"OptimizerBarzilaiBorwein")
  {

    // Starting lambda for iT==0: steepest descent
    // Line Search will be initialised with the natural lambda = 1
    _startLambda = lambda;

    _lastControl = util::ContainerCreator<C>::create(this->_dimCtrl);
    _lastDerivative = util::ContainerCreator<C>::create(this->_dimCtrl);

    _sStore = util::ContainerCreator<C>::create(this->_dimCtrl);
    _yStore = util::ContainerCreator<C>::create(this->_dimCtrl);

  };

  void checkDerivativeZero()
  {
    S normDir = util::euklidN(this->_direction.data(), this->_dimCtrl);
    if (std::isnan(normDir)) {
      clout << "Warning: Derivative is null at first iteration. Check derivative calculations.\nProgram terminated" << std::endl;
      exit(1);
    }
  }

  virtual void computeDirection()
  {
    for (int i=0; i<this->_dimCtrl; i++) {
      _sStore[i] = 0.;
      _yStore[i] = 0.;
    }
    S sTy = 0.;
    S sTs = 0.;
    S yTy = 0.;
    // Alternating between lambda_1 and lambda_2, see:
    // http://www.math.ucla.edu/~wotaoyin/math273a/slides/Lec4a_Baizilai_Borwein_method_273a_2015_f.pdf, p.6
    bool alternate = false;


    // Compute Derivatives
    // computeDerivative only if not already done by wolfeCondition()
    if (!this->_nextDerivFlag) {
      this->computeDerivatives(this->_control, this->_derivative);
    }
    else {
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_derivative[iDim] = this->_nextDerivative[iDim];
      }
    }


    if (this->_it == 0) {
      // On first step do normalised steepest descent
      S normDerivative = util::euklidN(this->_derivative.data(), this->_dimCtrl);
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_direction[iDim] = _startLambda * this->_derivative[iDim] / normDerivative;
      }
      // terminate early if we recieve nan controls
      checkDerivativeZero();
      // Store control and derivative
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        _lastControl[iDim] = this->_control[iDim];
        _lastDerivative[iDim] = this->_derivative[iDim];
      }

    }
    else {
      // Calculate lambda as quasi Newton
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        _sStore[iDim] = this->_control[iDim] - _lastControl[iDim];
        _yStore[iDim] = this->_derivative[iDim] - _lastDerivative[iDim];
        sTs += _sStore[iDim]*_sStore[iDim];
        sTy += _sStore[iDim]*_yStore[iDim];
        if (alternate) {
          yTy += _yStore[iDim]*_yStore[iDim];
        }
      }
      // lambda_1 is minimizing || sD - y ||
      this->_lambda = sTs/sTy; //lambda_1
      // lambda_2 is minimizing || s - yD^-1 ||
      //this->_lambda = sTy/yTy; //lambda_2

      //Alternate lambda_1 and lambda_2
      if (alternate && this->_it%2 != 0) {
        this->_lambda = sTy/yTy;
      }

      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        _lastControl[iDim] = this->_control[iDim];
        _lastDerivative[iDim] = this->_derivative[iDim];
      }

      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_direction[iDim] = this->_derivative[iDim];
      }
    }
  };
};


/// Creator Function for Barzilai-Borwein
template<typename S, typename C = std::vector<S>>
OptimizerBarzilaiBorwein<S,C>* createOptimizerBarzilaiBorwein(XMLreader const& params, std::size_t dimCtrl)
{
  OstreamManager clout(std::cout,"createOptimizerBarzilaiBorwein");

  // create variables with default values
  int maxIt = 100;

  S controlEps = S(0);
  S eps = S(1.e-10);

  S lambda = 1.; //starting lambda for steepest descent step (it==0)
  int maxStepAttempts  = 100;
  std::string stepCondition = "None";

  bool verboseOn=true;
  std::string fname = "control.dat";
  std::string logFileName = "log.txt";

  bool vectorBounds = false;
  bool withUpperBound = false;
  S upperBound = S();
  bool withLowerBound = false;
  S lowerBound = S();

  std::vector<OptimizerLogType> gplotAnalysis  = {};
  std::string gplotAnalysisString = "";

  // Read Values from XML File
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
  return new OptimizerBarzilaiBorwein<S,C>(dimCtrl, eps, maxIt, lambda, maxStepAttempts, stepCondition,
                                         verboseOn, fname, logFileName, withUpperBound, upperBound, withLowerBound, lowerBound,
                                         vectorBounds, controlEps, gplotAnalysis);
}

} // namespace opti

} // namespace olb

#endif
