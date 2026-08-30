/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger, Adrian Kummerlaender
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


#ifndef SOLVER_PARAMETERS_H
#define SOLVER_PARAMETERS_H


#include <iostream>

#include "io/ostreamManager.h"
#include "io/xmlReader.h"
#include "utilities/benchmarkUtil.h"
#include "names.h"
#include "utilities/typeMap.h"

namespace olb {

namespace parameters {

struct ParameterBase {
  virtual void initialize() { }
};

// ------------------------------ Simulation ----------------------------------

/** Base struct to keep the parameters that are necessary for the simulation.
 * Class supports reading simulation parameters from input files (e.g. XML) as
 * well as an interface for other parts of the program.
 * This class is intended to be derived from: at least, a UnitConverter has to
 * be provided additionally.
 */
template<typename T>
struct SimulationBase : public ParameterBase
{
  using BT = BaseType<T>;

  BT                                    startUpTime                 {0};
  BT                                    maxTime                     {1};
  BT                                    physBoundaryValueUpdateTime {0.1};

  BT                                    physTimeStabilityCheck      {0.1};
  bool                                  exitMaxU                    {false};
  BT                                    boundMaxU                   {1.0};

  bool                                  pressureFilter {false};

  // no of cuboids per mpi process
  int                                   noC     {1};
  int                                   overlap {3};

  virtual void initialize() { }
};

/// All the simulation parameters are read directly from an xml file
template<typename T, typename LATTICES>
struct XmlSimulation : public SimulationBase<T> {

  XMLreader     const* xml;

  using NAME = typename LATTICES::keys_t::template get<0>;
  using descriptor = typename LATTICES::template value<NAME>;
  std::shared_ptr<UnitConverter<T,descriptor> const> converter;
};


// ------------------------------ Stationarity --------------------------------

struct StationarityBase : public ParameterBase
{ };

/** All parameters that are necessary for checking whether the simulation
 * became stationary.
 * Template parameters are the arithmetic type T and (optionally) the names of
 * lattices for which stationarity is to be checked. If none are given,
 * NavierStokes is set as a default.
 */
template<typename T, typename... STAT_LATTICES>
struct Stationarity : public StationarityBase
{
  constexpr static unsigned numberOfArguments = sizeof...(STAT_LATTICES);
  // list of lattices for which stationarity is to be checked
  using stat_lattices = std::conditional_t <(numberOfArguments == 0),
    meta::list<names::NavierStokes>,
    meta::list<STAT_LATTICES...>>;
  constexpr static unsigned numberOfStationaryLattices = stat_lattices::size;

  enum ConvergenceType {MaxLatticeVelocity, AverageEnergy, AverageRho};

  std::array <ConvergenceType,numberOfStationaryLattices> convergenceType {MaxLatticeVelocity};
  std::array <T,numberOfStationaryLattices>               physInterval {1};
  std::array <T,numberOfStationaryLattices>               epsilon {1e-3};

  Stationarity() = default;

  Stationarity(ConvergenceType type, T interval_, T epsilon_)
  {
    std::fill(convergenceType.begin(), convergenceType.end(), type);
    std::fill(physInterval.begin(), physInterval.end(), interval_);
    std::fill(epsilon.begin(), epsilon.end(), epsilon_);
  }

  Stationarity(
    std::array <ConvergenceType,numberOfStationaryLattices> type,
    std::array <T,numberOfStationaryLattices> interval_,
    std::array <T,numberOfStationaryLattices> epsilon_)
   : convergenceType(type), physInterval(interval_), epsilon(epsilon_)
  { }
};


// ------------------------------ Output --------------------------------------

/** Structs to keep parameters which characterize the output.
 * So far, bash/ gnuplot/ image/ vtk are supported.
 */


template<typename T, typename LatticeLog= names::NavierStokes>
struct OutputGeneral : public ParameterBase
{
  using BT = BaseType<T>;
  std::string                              name         {"unnamed"};
  std::string                              olbDir       {"../../../"};
  std::string                              outputDir    {"./tmp/"};

  bool                                     verbose            {true};
  bool                                     printLogConverter  {true};
  BT                                       logT;

  using printLatticeStatistics =           meta::list<LatticeLog>;

  int                                      timerPrintMode     {2};

  OutputGeneral() = default;

  OutputGeneral(
    std::string name_,
    std::string olbDir_,
    std::string outputDir_,

    bool verbose_,
    bool printLogConverter_,

    BT logT_,

    int timerPrintMode_
  ) : name(name_), olbDir(olbDir_), outputDir(outputDir_), verbose(verbose_),
    printLogConverter(printLogConverter_), logT(logT_),
    timerPrintMode(timerPrintMode_)
  { }


  virtual void initialize() override
  {
    singleton::directories().setOlbDir(olbDir);
    singleton::directories().setOutputDir(outputDir);
  }
};



template <typename T>
struct OutputPlot : public ParameterBase
{
  using BT = BaseType<T>;
  // pO for print Output to avoid name clashes using off and final as keywords
  enum PrintOutput {pO_off, pO_intervals, pO_final};

  OutputPlot() = default;

  std::string output {"off"};
  std::string filename {"unnamed"};
  BT saveTime;
  PrintOutput printOutput = pO_off;

  OutputPlot(std::string out, std::string name, BT savetime) {
    if(out == "off"){
      printOutput = pO_off;
    } else if (out == "intervals"){
      printOutput = pO_intervals;
    } else if (out == "final"){
      printOutput = pO_final;
    }
    filename = name;
    saveTime = savetime;
  }
};


// ------------------------------ Results -------------------------------------

/** Struct to keep results of the simulation in order to provide communication
 * with other parts of the program.
 * This class is intended to be derived from.
 */
struct ResultsBase : public ParameterBase { };


// ------------------ Parameter reading functionality ------------------------

template<typename PARAMETERS>
struct ReaderBase
{
  mutable OstreamManager clout {std::cout, "ParameterReader"};
  std::shared_ptr<PARAMETERS> params;

  ReaderBase(std::shared_ptr<PARAMETERS> params_) : params(params_)
  { }
};

/** Base struct for reading parameters from files.
 * So far, reading from xml is supported. If you want to let your own parameter
 * struct read from xml, you may realize this as template specializations,
 * similar to the syntax in the examples below. The default (with no
 * specialization) does nothing.
 * If you want to add your own interface, you may add read(...)-methods to the
 * given structures which read the parameters from your file type.
 */
template<typename PARAMETERS, typename TAG>
struct Reader : public ReaderBase<PARAMETERS> {
  Reader(std::shared_ptr<PARAMETERS> params_) : Reader::ReaderBase(params_)
  { }

  virtual void read(XMLreader const& xml) { }
};


template<typename PARAMETERS, typename TAG>
Reader(std::shared_ptr<PARAMETERS> params_, TAG tag) -> Reader<PARAMETERS, TAG>;


template<typename PARAMETERS>
Reader(std::shared_ptr<PARAMETERS> params_) -> Reader<PARAMETERS, void>;


template<typename T, typename LatticeLog, typename TAG>
struct Reader<OutputGeneral<T,LatticeLog>, TAG> : public ReaderBase<OutputGeneral<T,LatticeLog>>
{
  using ReaderBase<OutputGeneral<T,LatticeLog>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    // Filenames
    xml.readOrWarn<std::string>("Application", "Name", "", this->params->name, true, true, true);
    xml.readOrWarn<std::string>("Application", "OlbDir", "", this->params->olbDir, true, true, true);
    xml.readOrWarn<std::string>("Output", "OutputDir", "", this->params->outputDir, true, true, true);

    // For Console Output
    xml.readOrWarn<bool>("Output", "Log", "VerboseLog", this->params->verbose, true, false, true);
    xml.readOrWarn<bool>("Output", "PrintLogConverter", "", this->params->printLogConverter, true, false, false);
    xml.readOrWarn<BaseType<T>>("Output", "Log", "SaveTime", this->params->logT, true, false, true);
    xml.readOrWarn<int>("Output", "Timer", "PrintMode", this->params->timerPrintMode, true, false, false);
  }
};



template<typename T, typename TAG>
struct Reader<OutputPlot<T>, TAG> : public ReaderBase<OutputPlot<T>>
{
  using ReaderBase<OutputPlot<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    std::string tag {TAG().name};
    xml.readOrWarn<std::string>("Output", tag, "Output", this->params->output, true, false, true);
    if (this->params->output != "off") {
      xml.readOrWarn<std::string>("Output", tag, "Filename", this->params->filename, true, false, false);
      xml.readOrWarn<BaseType<T>>("Output", tag, "SaveTime", this->params->saveTime, true, false, true);
    }
    if(this->params->output == "off"){
      this->params->printOutput = olb::parameters::OutputPlot<T>::pO_off;
    } else if (this->params->output == "intervals"){
      this->params->printOutput = olb::parameters::OutputPlot<T>::pO_intervals;
    } else if (this->params->output == "final"){
      this->params->printOutput = olb::parameters::OutputPlot<T>::pO_final;
    }
  }
};

template<typename T, typename TAG>
struct Reader<SimulationBase<T>, TAG> : public ReaderBase<SimulationBase<T>>
{
  using ReaderBase<SimulationBase<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    using BT = BaseType<T>;

    xml.readOrWarn<BT>("Application", "PhysParameters", "StartUpTime", this->params->startUpTime, false, false, false);
    xml.readOrWarn<BT>("Application", "PhysParameters", "PhysMaxTime", this->params->maxTime, false, true, true);
    xml.readOrWarn<int>("Application", "Mesh", "noCuboidsPerProcess", this->params->noC, true, false, true);
    xml.readOrWarn<bool>("Application", "PressureFilter", "", this->params->pressureFilter, true, false, false);

    this->params->physBoundaryValueUpdateTime = this->params->maxTime / BT(100); // 1% of max. time as default value
    xml.readOrWarn<BT>("Application", "PhysParameters", "BoundaryValueUpdateTime", this->params->physBoundaryValueUpdateTime, true, false, true);
    this->params->physTimeStabilityCheck = this->params->maxTime / BT(100);

    xml.readOrWarn<BT>("Application", "StabilityCheck", "Frequency", this->params->physTimeStabilityCheck, true, false, false);
    xml.readOrWarn<bool>("Application", "StabilityCheck", "ExitMaxU", this->params->exitMaxU, true, false, false);
    xml.readOrWarn<BT>("Application", "StabilityCheck", "BoundMaxU", this->params->boundMaxU, true, false, false);
  }
};

template<typename T, typename TAG>
struct Reader<Stationarity<T>, TAG> : public ReaderBase<Stationarity<T>>
{
  using ReaderBase<Stationarity<T>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    std::string type = "MaxLatticeVelocity";
    xml.readOrWarn<std::string>("Application", "ConvergenceCheck", "Type", type, true, false, true);
    this->params->convergenceType[0] = (type == "AverageEnergy") ? Stationarity<T>::AverageEnergy : Stationarity<T>::MaxLatticeVelocity;
    BaseType<T> help; // needed to prevent typecast errors with util::ADf data type
    xml.readOrWarn<BaseType<T>>("Application", "ConvergenceCheck", "Interval", help, true, false, true);
    this->params->physInterval[0] = help;
    xml.readOrWarn<BaseType<T>>("Application", "ConvergenceCheck", "Residuum", help, true, false, true);
    this->params->epsilon[0] = help;
  }
};

template<typename T, typename LATTICES, typename TAG>
struct Reader<XmlSimulation<T,LATTICES>, TAG> : public ReaderBase<XmlSimulation<T,LATTICES>>
{
  using ReaderBase<XmlSimulation<T,LATTICES>>::ReaderBase;

  void read(XMLreader const& xml_)
  {
    Reader<SimulationBase<T>, TAG>(this->params).read(xml_);
    this->params->xml = &xml_;

    using NAME = typename LATTICES::keys_t::template get<0>;
    using descriptor = typename LATTICES::template value<NAME>;
    this->params->converter = std::shared_ptr<UnitConverter<T,descriptor>>(createUnitConverter<T,descriptor>(xml_));
  }
};

} // namespace parameters


// ------------------ Creator functions for XML interface ---------------------

template <typename PARAMETERS, typename TAG>
std::shared_ptr<PARAMETERS> createParameters(XMLreader const& xml)
{
  auto result = std::make_shared<PARAMETERS>();
  parameters::Reader<PARAMETERS, TAG>(result).read(xml);
  result->initialize();
  return result;
}

template <class parameters_map>
auto createParametersTuple(XMLreader const& xml)
{
  utilities::TypeIndexedSharedPtrTuple<parameters_map> paramsTuple;
  meta::tuple_for_each(paramsTuple.tuple, [&xml](auto& element, auto index) {
    using parameter_type = typename std::remove_reference_t<decltype(element)>::element_type;
    using tag = typename parameters_map::keys_t::template get<index()>;
    element = createParameters<parameter_type,tag>(xml);
  });
  return paramsTuple;
}

} // namespace olb

#endif
