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
 * The description of optimization algorithms -- header file.
 */


#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "optiCase.h"
#include "communication/mpiManager.h"
#include "io/gnuplotWriter.h"
#include "io/xmlReader.h"
#include "core/singleton.h"

namespace olb {

/// \namespace opti Optimization Code.
namespace opti {


template<typename, typename> class OptiCase;

enum OptimizerLogType {control, derivative, value, error, norm_derivative};

/// Interface for the use of various optimization algorithms.
/**
 * This class is intended to be derived from.
 * S is the underlying arithmetic data type
 * C is the container type that is used to store the vectors (e.g. std::vector<S>)
 * C is expected to allow construction via ContainerCreator class
 */
template<typename S, typename C>
class Optimizer {

private:
  mutable OstreamManager clout;

protected:
  /// Number of controlled variables
  int                                 _dimCtrl;
  /// Vector of controlled variables (size _dimCtrl)
  C                                   _control;
  /// Value of the objective functional evaluated
  /// for controlled variables saved in _control
  S                                   _value;
  /// Vector of derivatives of the object functional
  /// with respect to the controlled variables
  C                                   _derivative;
  /// Current iteration no.
  int                                 _it;
  /// Maximal number of iteration
  int                                 _maxIt;
  /// Fail when max number of iteration reached
  bool                                _failOnMaxIter;
  /// Optimizer stops if |_derivatives| < _eps
  S                                   _eps;
  /// Verbose
  bool                                _verboseOn;

  /// Bounded versions
  bool                                _withUpperBound;
  bool                                _withLowerBound;
  bool                                _vectorBounds;
  C                                   _boundedControl;
  C                                   _upperBound;
  C                                   _lowerBound;

  /// For setting tolerance of controls
  bool _controlsConverged;
  S _controlEps;

  /// Provides the Optimizer with methods to evaluate the
  /// value of an object functional and compute derivatives
  OptiCase<S,C>*                      _optiCase;

  /// control vector to compare with (for numerical evaluation)
  C                                   _referenceControl;

public:
  Gnuplot<S>                          gplot;
  /// For defining what kind of gnuplot analysis is wanted, if empty vector - no analysis,
  /// value, control and derivative are the possible options
  std::vector<OptimizerLogType>       _gplotAnalysis;

  Optimizer(int dimCtrl, S eps, int maxIt,
        bool verboseOn=true, const std::string fname="",
        const std::string logFileName="",
        bool withUpperBound=false, S upperBound=S(),
        bool withLowerBound=false, S lowerBound=S(),
        bool vectorBounds=false,
        S controlEps=S(std::numeric_limits<double>::epsilon() ),
        bool failOnMaxIter = true,
        std::vector<OptimizerLogType> gplotAnalysis = {});

  virtual ~Optimizer() { };

  virtual void optimizationStep() = 0;

  void maxIterationReached();

  virtual void optimize();

  virtual void optimize(OptiCase<S,C>& optiCase) {
    _optiCase = &optiCase;
    optimize();
  }

  void simulate() {
    evaluateObjective(_control, _value);
  }

  void simulate(OptiCase<S,C>& optiCase) {
    _optiCase = &optiCase;
    simulate();
  }

  void evaluateObjective(const C& control, S& result) {
    result = _optiCase->evaluateObjective(control, _it);
  }

  void computeDerivatives(const C& control, C& derivatives) {

    _optiCase->computeDerivatives(control, derivatives, _it);
  }

  /// Prints information of the current optimization step it
  void print(int it);

  void setControl(C& control) {
    _control = std::move(control);
  }

  const C& getControl() const {
    return _control;
  }

  const C& getDerivative() const {
    return _derivative;
  }

  const S& getObjective() const {
    return _value;
  }

  int getIteration() const {
    return _it;
  }

  /// Writes the current control variables linewise into file fname
  void writeControlToFile(const std::string fname="control.dat");

  /// Reads the latest control variables from file fname
  void readControlFromFile(const std::string fname="control.dat");

  void setStartValue(S startValue);

  void setStartValue(C startValues, S factor=S(1));

  OptiCase<S,C>* getOptiCase() { return _optiCase; }

  void setOptiCase(OptiCase<S,C>* optiCase) { _optiCase = optiCase; }

  void setGnuplotData();

  /// set the reference value for the control vector (exact solution)
  void setReferenceControl(C result) {
    _referenceControl = result;
  }
};

/// the gplotAnalysisString is gained from the xml file and the function than separates and prepares it to be used in the constructor
void getGnuplotTagsFromString(std::string gplotAnalysisString, std::vector<OptimizerLogType>& gplotAnalysis );

} // namespace opti

} // namespace olb

#endif
