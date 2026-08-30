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
 * The description of the LBFGS optimization algorithm -- header file.
 */


#ifndef OPTIMIZER_LBFGS_H
#define OPTIMIZER_LBFGS_H

#include "optimizerLineSearch.h"
#include "utilities/norm.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {


/// Optimization algorithm: LBFGS.
/** OptimizerSteepestDescent optimizes an optimization problem
 * which is given in OptiCase. OptiCase provides therefore
 * methods to evaluate an object functional and compute derivatives.
 * LBFGS searches along the direction of the derivative to find
 * a lower value of the evaluated object functional.
 *
 * This class is not intended to be derived from.
 */

template<typename S, typename C>
class OptimizerLBFGS : public OptimizerLineSearch<S,C> {

private:
  S _startLambda;

  int _l;

  S _startCoefH;

  S** _sStore;
  S** _yStore;
  S* _rhoStore;
  S* _alpha;

  int _firstStore;
  int _lastStore;

  C _lastControl;
  C _lastDerivative;

  mutable OstreamManager clout;

public:
  OptimizerLBFGS(
    int dimCtrl, S eps, int maxIt, S lambda, int maxStepAttempts,
    std::string stepCondition, int l, S startCoefH, bool verboseOn=true,
    const std::string fname="", const std::string logFileName="",
    bool withUpperBound=false, S upperBound=S(),
    bool withLowerBound=false, S lowerBound=S(), bool vectorBounds=false,
    S controlEps=S(std::numeric_limits<double>::epsilon() ), bool failOnMaxIter = true,
    std::vector<OptimizerLogType> gplotAnalysis = {})
    : OptimizerLineSearch<S,C>(
      dimCtrl, eps, maxIt, /*lambda*/ 1., maxStepAttempts, stepCondition, verboseOn,
      fname, logFileName, withUpperBound, upperBound,
      withLowerBound, lowerBound, vectorBounds, controlEps, failOnMaxIter, gplotAnalysis),
      clout(std::cout,"OptimizerLBFGS")
  {
    // Starting lambda for iT==0: steepest descent
    // Line Search will be initialised with the natural lambda = 1
    _startLambda = lambda;

    _l=l;
    _startCoefH = startCoefH;
    _firstStore = 1;
    _lastStore = 1;

    _lastControl = util::ContainerCreator<C>::create(this->_dimCtrl);
    _lastDerivative = util::ContainerCreator<C>::create(this->_dimCtrl);

    _sStore = new S* [_l];
    _yStore = new S* [_l];

    _rhoStore = new S [_l];
    _alpha = new S [_l];

    for (int i=0; i<_l; i++) {
      _rhoStore[i] = S(0);
      _alpha[i] = S(0);
    }

    for (int i=0; i<_l; i++) {
      _sStore[i] = new S [this->_dimCtrl];
      _yStore[i] = new S [this->_dimCtrl];
      //for (int iDim=0; iDim<this->_dimCtrl; ++iDim) {
      //_sStore[i][iDim] = S(0);
      //_yStore[i][iDim] = S(0);
      //};
    }
  };

  ~OptimizerLBFGS()
  {
    for (int i=0; i<_l; i++) {
      delete[] _sStore[i];
      delete[] _yStore[i];
    }
    delete[] _sStore;
    delete[] _yStore;
    delete[] _rhoStore;
    delete[] _alpha;
  }


  void setL(int l)
  {
    _l=l;
  };

  void setStartCoefH( S startCoefH)
  {
    _startCoefH = startCoefH;
  };

  void checkDerivativeZero()
  {
    S normDir = util::euklidN(this->_direction.data(), this->_dimCtrl);
    if (std::isnan(normDir)) {
      clout << "Warning: Derivative is null at first iteration. Check derivative calculations.\nProgram terminated" << std::endl;
      exit(1);
    }
  }

  void storeHelpers()
  {

    if (this->_it!=0) {
      _rhoStore[(this->_it)%_l] = 0;
      S temp1 = S();
      S temp2 = S();
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        _sStore[(this->_it)%_l][iDim] = this->_control[iDim]-_lastControl[iDim];
        _yStore[(this->_it)%_l][iDim] = this->_derivative[iDim]-_lastDerivative[iDim];
        _rhoStore[(this->_it)%_l] += _yStore[(this->_it)%_l][iDim]*_sStore[(this->_it)%_l][iDim];
        temp1 += _sStore[(this->_it)%_l][iDim]*_yStore[(this->_it)%_l][iDim];
        temp2 += _yStore[(this->_it)%_l][iDim]*_yStore[(this->_it)%_l][iDim];
      }
      _startCoefH = temp1/temp2;

      _rhoStore[(this->_it)%_l] = 1./_rhoStore[(this->_it)%_l];
      if (this->_it>_l) {
        _firstStore = (_firstStore+1)%_l;
      }
      _lastStore = (this->_it)%_l;
    }

    for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
      _lastControl[iDim] = this->_control[iDim];
      _lastDerivative[iDim] = this->_derivative[iDim];
    }
  }

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

    if (this->_it==0) {
      // On first step do normalised steepest descent
      S normDerivative = util::euklidN(this->_derivative.data(), this->_dimCtrl);
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_direction[iDim] = _startLambda * this->_derivative[iDim] / normDerivative;
      }
      // terminate early if we recieve nan controls
      checkDerivativeZero();

      // Store Helpers!
      storeHelpers();

    }
    else {
      storeHelpers();

      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_direction[iDim] = this->_derivative[iDim];
      }

      int iMax = _l;
      if (this->_it<_l) {
        iMax = this->_it;
      }

      int iStore;
      // first for loop
      for (int i=0; i<iMax; i++) {
        iStore = (_lastStore + _l - i) % _l;

        _alpha[iStore] = 0;
        for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
          _alpha[iStore] += _sStore[iStore][iDim]*this->_direction[iDim];
        }
        _alpha[iStore] *= _rhoStore[iStore];

        for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
          this->_direction[iDim] -= _alpha[iStore]*_yStore[iStore][iDim];
        }
      } // end first for loop
      for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
        this->_direction[iDim] *= _startCoefH;
      }

      // second for loop
      for (int i=0; i<iMax; i++) {
        iStore = (_firstStore + _l + i) % _l;

        S beta = 0;
        for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
          beta += _yStore[iStore][iDim]*this->_direction[iDim];
        }
        beta *= _rhoStore[iStore];

        for (int iDim=0; iDim<this->_dimCtrl; iDim++) {
          this->_direction[iDim] += _sStore[iStore][iDim]*(_alpha[iStore]-beta);
        }
      } // end second for loop
    }
  };
};



/// Creator Function for LBFGS
template<typename S, typename C = std::vector<S>>
OptimizerLBFGS<S,C>* createOptimizerLBFGS(XMLreader const& params, std::size_t dimCtrl)
{
  OstreamManager clout(std::cout, "createOptimizerLBFGS");

  // create variables with default values
  int maxIt = 100;
  int l = 20;
  S lambda = 1.; //starting lambda for steepest descent step (it==0)
  int maxStepAttempts = 100;
  std::string stepCondition = "StrongWolfe";

  //S controlEps = S(std::numeric_limits<double>::epsilon() );
  S controlEps = S(0);
  S eps = S(1.e-10);
  S startCoefH = S(1e-4);

  bool verboseOn=true;
  std::string fname = "control.dat";
  std::string logFileName = "log.txt";
  bool vectorBounds = false;
  bool withUpperBound = false;
  S upperBound = S();
  bool withLowerBound = false;
  S lowerBound = S();
  bool failOnMaxIter = true;

  std::vector<OptimizerLogType> gplotAnalysis  = {};
  std::string gplotAnalysisString = "";

  // Read Values from XML File from area "Optimization"
  params.readOrWarn<int>("Optimization", "MaxIter", "", maxIt);
  params.readOrWarn<int>("Optimization", "L", "", l);
  params.readOrWarn<S>("Optimization", "Lambda", "", lambda);
  params.readOrWarn<int>("Optimization", "MaxStepAttempts", "", maxStepAttempts);
  params.readOrWarn<std::string>("Optimization", "StepCondition", "", stepCondition);

  params.readOrWarn<S>("Optimization", "Tolerance", "", eps);
  params.readOrWarn<S>("Optimization", "ControlTolerance", "", controlEps);
  params.readOrWarn<S>("Optimization", "StartCoefH", "", startCoefH);
  params.readOrWarn<bool>("Optimization", "FailOnMaxIter", "", failOnMaxIter);
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
  params.readOrWarn<std::string>("Optimization",
                                        "VisualizationGnuplot",
                                        "VisualizedParameters", gplotAnalysisString, false, false);
  // transform the data from the xml file to the enums needed for continueing
  getGnuplotTagsFromString(gplotAnalysisString, gplotAnalysis);

  clout << "Creating optimizer ..." << std::endl;
  // Create Optimizer Object
  return new OptimizerLBFGS<S,C>(dimCtrl, eps, maxIt, lambda, maxStepAttempts, stepCondition, l, startCoefH,
                               verboseOn, fname, logFileName, withUpperBound, upperBound, withLowerBound, lowerBound, vectorBounds,
                               controlEps, failOnMaxIter, gplotAnalysis);
}


} // namespace opti

} // namespace olb

#endif
