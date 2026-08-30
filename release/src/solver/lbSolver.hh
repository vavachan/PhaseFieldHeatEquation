/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster,
 *  Julius Jessberger
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


#ifndef LBSOLVER_HH
#define LBSOLVER_HH


#include "lbSolver.h"


namespace olb {

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::initialize()
{
  if (this->parameters(names::Output()).verbose) {
    clout << "Start initialization" << std::endl;
  }

  this->prepareGeometry();
  // Removes all not needed boundary voxels outside the surface
  //this->geometry().clean(this->parameters(names::Output()).verbose);
  // Removes all not needed boundary voxels inside the surface
  this->geometry().innerClean(this->parameters(names::Output()).verbose);
  geometry().checkForErrors(this->parameters(names::Output()).verbose);
  if (this->parameters(names::Output()).verbose) {
    geometry().getStatistics().print();
  }

  this->_isInitialized = true;

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished initialization" << std::endl;
  }
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::buildAndReturn()
{
  if (! this->_isInitialized) {
    this->initialize();
  }
  build();
  setInitialValues();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::build()
{
  if (this->parameters(names::Output()).printLogConverter) {
    writeLogConverter();
  }

  _timer = std::make_unique<util::Timer<BaseType<T>>>(
    converter().getLatticeTime(this->parameters(names::Simulation()).maxTime),
    _sGeometry->getStatistics().getNvoxel(),
    1); // 1 for printing a short timer summary

  if constexpr (isStationary){
    for (unsigned i = 0; i < getNumberStationaryLattices(); ++i) {
      _convergenceCheck[i] = std::make_unique<util::ValueTracer<T>>(
        converter().getLatticeTime(this->parameters(names::Stationarity()).physInterval[i]),
        this->parameters(names::Stationarity()).epsilon[i]);
    }
  }

  _itBoundaryUpdate = util::max(
    converter().getLatticeTime(this->parameters(names::Simulation()).physBoundaryValueUpdateTime),
    (unsigned long) {1});
  _itCheckStability = util::max(
    converter().getLatticeTime(this->parameters(names::Simulation()).physTimeStabilityCheck),
    (unsigned long) {1});

  renewLattices();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::prepareSimulation()
{
  if (this->parameters(names::Output()).verbose) {
    clout << "Start prepareSimulation" << std::endl;
  }

  this->_iT = 0;

  build();

  setInitialValues();

  meta::tuple_for_each(_sLattices, [](auto& lattice){
    lattice->initialize();
    lattice->getStatistics().reset();
    lattice->getStatistics().initialize();
  });


  if constexpr (outputVTK) {
    if (this->parameters(names::VisualizationVTK()).output != "off") {
      prepareVTK();
    }
  }

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished prepareSimulation" << std::endl;
    clout << "Start collide-and-stream loop" << std::endl;
  }
  _timer->start();
}


template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::timeStep(std::size_t iT)
{
  if (iT % _itBoundaryUpdate == 0) {
    setBoundaryValues(iT);
  }

  meta::tuple_for_each(_sLattices, [](auto& lattice) {
    // TODO: dispose also communication of fields
    lattice->communicate();
  });


  if constexpr (outputVTK) {
    if ( this->parameters(names::VisualizationVTK()).printOutput == olb::parameters::OutputPlot<T>::pO_intervals
      && iT % converter().getLatticeTime(this->parameters(names::VisualizationVTK()).saveTime) == 0 ) {
      writeVTK(iT);
    }
  }
  if constexpr (outputImages) {
    if ( this->parameters(names::VisualizationImages()).printOutput == olb::parameters::OutputPlot<T>::pO_intervals
      && iT % converter().getLatticeTime(this->parameters(names::VisualizationImages()).saveTime) == 0 ) {
      writeImages(iT);
    }
  }
  if constexpr (outputGnuplot) {
    if ( this->parameters(names::VisualizationGnuplot()).printOutput == olb::parameters::OutputPlot<T>::pO_intervals
      && iT % converter().getLatticeTime(this->parameters(names::VisualizationGnuplot()).saveTime) == 0 ) {
      writeGnuplot(iT);
    }
  }
  if ( this->parameters(names::Output()).verbose
    && iT % converter().getLatticeTime(this->parameters(names::Output()).logT) == 0) {
    printLog(iT);
  }
  computeResults(iT);
  getResults(iT);

  meta::tuple_for_each(_sLattices, [](auto& lattice){
    lattice->collideAndStream();
  });

  if (this->parameters(names::Simulation()).pressureFilter) {
    // TODO: allow choosing the lattices here
    std::get<0>(_sLattices)->stripeOffDensityOffset(std::get<0>(_sLattices)->getStatistics().getAverageRho() - T(1));
  }

  if constexpr (isStationary){
    unsigned counter = 0;
    using parameters_t = typename PARAMETERS::template value<names::Stationarity>;
    parameters_t::stat_lattices::for_each( [&](auto Id) {
      // using Name = typename decltype(Id)::type;
      if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::MaxLatticeVelocity){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getMaxU());
      }
      else if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::AverageEnergy){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getAverageEnergy() );
      }
      else if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::AverageRho){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getAverageRho());
      }
      else {
        throw std::invalid_argument("Convergence type is not supported.\n");
      }
      counter++;
    });
  }

  if (iT % _itCheckStability == 0) {
    checkStability(iT);
  }
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::postSimulation()
{
  meta::tuple_for_each(_sLattices, [](auto& lattice) {
    // TODO: dispose also communication of fields...
    lattice->communicate();
  });

  if ( this->parameters(names::Output()).verbose ) {
    printLog(this->_iT);
  }
  if constexpr (outputVTK) {
    if ( this->parameters(names::VisualizationVTK()).printOutput == olb::parameters::OutputPlot<T>::pO_final) {
      writeVTK(this->_iT);
    }
  }
  if constexpr (outputImages) {
    if ( this->parameters(names::VisualizationImages()).printOutput == olb::parameters::OutputPlot<T>::pO_final) {
      writeImages(this->_iT);
    }
  }
  if constexpr (outputGnuplot) {
    if ( this->parameters(names::VisualizationGnuplot()).printOutput == olb::parameters::OutputPlot<T>::pO_final) {
      writeGnuplot(this->_iT);
    }
  }
  getResults(this->_iT);

  _timer->stop();
  if (this->parameters(names::Output()).verbose) {
    clout << "Finished collide-and-stream loop" << std::endl;
    clout << "Start postprocessing" << std::endl;
  }
  computeResults();

  if (this->parameters(names::Output()).verbose) {
    _timer->printSummary();
  }

  // required to finish output tasks during use of gpu
  meta::tuple_for_each(_sLattices, [](auto& lattice) {
    lattice->setProcessingContext(ProcessingContext::Evaluation);
  });

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished postprocessing" << std::endl;
  }
}

template<typename T,  typename PARAMETERS, typename LATTICES>
bool LbSolver<T,PARAMETERS,LATTICES>::exitCondition(std::size_t iT) const
{
  if (iT > this->converter().getLatticeTime(this->parameters(names::Simulation()).maxTime)){
    return true;
  }

  if constexpr (isStationary){
    if (std::all_of(_convergenceCheck.cbegin(), _convergenceCheck.cend(), [](auto& c){
      return c->hasConverged();
      } )) {
      if (this->parameters(names::Output()).verbose) {
        clout << "Simulation converged." << std::endl;
      }
      return true;
    }
  }

  return false;
}

template<typename T,  typename PARAMETERS, typename LATTICES>
bool LbSolver<T,PARAMETERS,LATTICES>::checkStability(std::size_t iT)
{
  bool result = true;
  meta::tuple_for_each(_sLattices, [&](auto& lattice){
    if (lattice->getStatistics().getMaxU() > this->parameters(names::Simulation()).boundMaxU) {
      clout << "PROBLEM uMax=" << lattice->getStatistics().getMaxU() << std::endl;
      lattice->getStatistics().print(iT, converter().getPhysTime(iT));
      if (this->parameters(names::Simulation()).exitMaxU) {
        exit(1);
      }
      result = false;
    }
  });
  return result;
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::renewLattices()
{
  meta::tuple_for_each(_sLattices, [&](auto& lattice){
    using lattice_type = typename std::remove_reference_t<decltype(lattice)>::element_type;
    lattice = std::make_shared<lattice_type>(geometry());
  });

  prepareLattices();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::printLog(std::size_t iT) const
{
  using OutputParameters_t = typename PARAMETERS::template value<names::Output>;
  OutputParameters_t::printLatticeStatistics::for_each( [&](auto type){
    lattice(type.get()).getStatistics().print(iT, converter().getPhysTime(iT));
  } );

  _timer->print(iT, this->parameters(names::Output()).timerPrintMode);
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::writeLogConverter() const
{
  converter().print();
  converter().write(this->parameters(names::Output()).name);
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LbSolver<T,PARAMETERS,LATTICES>::prepareVTK() const
{
  // write the geometric information. Since this does not depend on the
  // lattice, we may w.l.o.g. work with the first lattice
  auto& lattice = *std::get<0>(_sLattices);
  using DESCRIPTOR = typename LATTICES::values_t::template get<0>;


  SuperVTMwriter<T,dim> vtmWriter(this->parameters(names::VisualizationVTK()).filename);

  /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
  SuperLatticeCuboid<T,DESCRIPTOR> cuboid(lattice);
  SuperLatticeRank<T,DESCRIPTOR> rank(lattice);
  vtmWriter.write( cuboid );
  vtmWriter.write( rank );
  vtmWriter.createMasterFile();
}

} // namespace olb

#endif
