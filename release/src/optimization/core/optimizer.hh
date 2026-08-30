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

#ifndef OPTIMIZER_HH
#define OPTIMIZER_HH

#include "optimizer.h"
#include "utilities/norm.h"
#include <iterator>

namespace olb {

/// \namespace opti Optimization Code.
namespace opti {


template<typename, typename> class OptiCase;


template<typename S, typename C>
Optimizer<S,C>::Optimizer(int dimCtrl, S eps, int maxIt,
        bool verboseOn, const std::string fname,
        const std::string logFileName,
        bool withUpperBound, S upperBound,
        bool withLowerBound, S lowerBound,
        bool vectorBounds,
        S controlEps,
        bool failOnMaxIter,
        std::vector<OptimizerLogType> gplotAnalysis)
: clout(std::cout,"Optimizer"), _dimCtrl(dimCtrl),
    _it(0), _maxIt(maxIt), _failOnMaxIter(failOnMaxIter), _eps(eps),
    _verboseOn(verboseOn),
    _withUpperBound(withUpperBound), _withLowerBound(withLowerBound),
    _vectorBounds(vectorBounds), _controlEps(controlEps),
    gplot(logFileName, Gnuplot<S>::LINEAR, Gnuplot<S>::OFF),
    _gplotAnalysis(gplotAnalysis)
{
  _controlsConverged = false;

  _control = util::ContainerCreator<C>::create(_dimCtrl);
  _derivative = util::ContainerCreator<C>::create(_dimCtrl);

  if (_withUpperBound||_withLowerBound) {
    _boundedControl = util::ContainerCreator<C>::create(_dimCtrl);
    if (_withLowerBound) {
      _lowerBound = util::ContainerCreator<C>::create(_dimCtrl);
      if (! _vectorBounds) {
        _lowerBound[0] = lowerBound;  // only first component is used...
      }
    }
    if (_withUpperBound) {
      _upperBound = util::ContainerCreator<C>::create(_dimCtrl);
      if (! _vectorBounds) {
        _upperBound[0] = upperBound;  // only first component is used...
      }
    }
  }

  readControlFromFile(fname);
}

template<typename S, typename C>
void Optimizer<S,C>::maxIterationReached()
{
  if (_failOnMaxIter) {
    clout << "Warning: Optimization problem failed to converge within specified iteration limit of "
          << _maxIt << " iterations with tolerance of " << _eps << ".\nProgram terminated" << std::endl;
    exit(1);
  }
  else {
    clout << "Maximum iteration reached." << std::endl;
  }
}

template<typename S, typename C>
void Optimizer<S,C>::optimize()
{
  // Evaluate objective function
  evaluateObjective(_control, _value);
  _optiCase->postEvaluation();

  // Optimization step (update of _control, _value, _derivatives and _it)
  do {
    if (_it >= _maxIt) {
      maxIterationReached();
      return;
    }

    // Call optimisation algorithm: evaluate objective function and compute derivatives
    optimizationStep();

    // Print info
    _optiCase->postEvaluation();
    if (_verboseOn) {
      print(_it);
    }
    if (_gplotAnalysis.size() > 0) {
      setGnuplotData();
    }

    // Optimize as long as |derivative| > eps or a maximum of iterations is reached
  }
  while (util::euklidN2(_derivative) > (_eps*_eps) && !_controlsConverged);

  if(_gplotAnalysis.size() > 0){
    gplot.writePNG();
  }
}

template<typename S, typename C>
void Optimizer<S,C>::print(int it)
{
  clout << "=======================================" << std::endl;
  clout << ">>>>>>>>>> Optimizer: step " << it << " <<<<<<<<<<" << std::endl;
  clout << "   Value objective = " << std::setprecision(12) << _value << std::endl;
  clout << "   Norm derivative = " << std::setprecision(12) << util::euklidN(&_derivative[0], _dimCtrl) << std::endl;
  clout << "=======================================" << std::endl;
}

template<typename S, typename C>
void Optimizer<S,C>::writeControlToFile(const std::string fname)
{
  if (singleton::mpi().isMainProcessor() ) {
    std::stringstream stream;
    stream << _dimCtrl << std::endl;
    std::for_each(_control.begin(), _control.end(), [&stream](auto c){ stream << c << '\n'; });
    std::ofstream file;
    file.open(singleton::directories().getLogOutDir() + fname.c_str());
    file << stream.str();
    file.close();
  }
}

template<typename S, typename C>
void Optimizer<S,C>::readControlFromFile(const std::string fname)
{
  std::ifstream file(fname.c_str());
  if ( file.is_open() ) {
    int dimCtrl;
    file >> dimCtrl;
    if (dimCtrl!=_dimCtrl) {
      clout << "Error: dimensions do not match! dim_controller=" << _dimCtrl << "; dim_file=" << dimCtrl << std::endl;
      assert(false);
    }
    else {
      for (int i=0; i<_dimCtrl; i++) {
        double tmpVal;
        file >> tmpVal;
        _control[i]  = tmpVal;
      }
      if (_withLowerBound && _vectorBounds) {
        for (int i=0; i<_dimCtrl; i++) {
          double tmpVal;
          file >> tmpVal;
          _lowerBound[i]  = tmpVal;
        }
      }
      if (_withUpperBound && _vectorBounds) {
        for (int i=0; i<_dimCtrl; i++) {
          double tmpVal;
          file >> tmpVal;
          _upperBound[i]  = tmpVal;
        }
      }
    }
    file.close();
  }
}

template<typename S, typename C>
void Optimizer<S,C>::setStartValue(S startValue)
{
  // update control variables
  std::fill(_control.begin(), _control.end(), startValue);
}

template<typename S, typename C>
void Optimizer<S,C>::setStartValue(C startValues, S factor)
{
  std::transform(startValues.begin(), startValues.end(), _control.begin(), [factor](auto s){
    return factor*s;
  });
}

template<typename S, typename C>
void Optimizer<S,C>::setGnuplotData()
{
  std::vector<S> yValues;
  std::vector<std::string> labels;
  std::vector<char> plotTypes;

  for(unsigned i = 0; i < _gplotAnalysis.size(); i++){
    if(_gplotAnalysis[i] == value){
      yValues.push_back(_value);
      labels.push_back("value");
      plotTypes.push_back('p');
    } else if(_gplotAnalysis[i] == control){
        // get all entries of control
      for(unsigned j = 0; j< _control.size(); j++){
        yValues.push_back(_control[j]);
        labels.push_back("control[" + std::to_string(j) + "]");
        plotTypes.push_back('p');
      }
    } else if(_gplotAnalysis[i] == derivative){

      for(unsigned j = 0; j< _derivative.size(); j++){
        yValues.push_back(_derivative[j]);
        labels.push_back("derivative[" + std::to_string(j) + "]");
        plotTypes.push_back('p');
      }
    } else if(_gplotAnalysis[i] == error){
      yValues.push_back(util::euklidDistance(_control, _referenceControl));
      labels.push_back("control_error");
      plotTypes.push_back('p');
    } else if(_gplotAnalysis[i] == norm_derivative){
      yValues.push_back(util::euklidN(&_derivative[0], _dimCtrl));
      labels.push_back("derivative_norm");
      plotTypes.push_back('p');
    }
  }
  gplot.setData(S(_it), yValues, labels, "top right", plotTypes); // _it is the iterationstep
}

void getGnuplotTagsFromString(std::string gplotAnalysisString,
                              std::vector<OptimizerLogType>& gplotAnalysis ){
  // part the string to get the enum vector
  std::vector<std::string> gplotAnalysisStringVector;
  std::istringstream iss(gplotAnalysisString);
  std::copy(std::istream_iterator<std::string>(iss),
        std::istream_iterator<std::string>(),
       std::back_inserter(gplotAnalysisStringVector));

  for(unsigned i = 0; i < gplotAnalysisStringVector.size(); i++){
    if(gplotAnalysisStringVector[i] == "VALUE"){
      gplotAnalysis.push_back(OptimizerLogType::value);
    } else if(gplotAnalysisStringVector[i] == "CONTROL"){
      gplotAnalysis.push_back(OptimizerLogType::control);
    } else if(gplotAnalysisStringVector[i] == "DERIVATIVE"){
      gplotAnalysis.push_back(OptimizerLogType::derivative);
    } else if(gplotAnalysisStringVector[i] == "ERROR"){
      gplotAnalysis.push_back(OptimizerLogType::error);
    } else if(gplotAnalysisStringVector[i] == "NORM_DERIVATIVE"){
      gplotAnalysis.push_back(OptimizerLogType::norm_derivative);
    }
  }
}

} // namespace opti

} // namespace olb

#endif
