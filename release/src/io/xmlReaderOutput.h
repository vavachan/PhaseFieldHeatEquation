/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Simon Großmann
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

#ifndef IO_XML_READER_OUTPUT_H_
#define IO_XML_READER_OUTPUT_H_

#include <string>
#include <vector>
#include <iostream>

#include "io/ostreamManager.h"

namespace olb {

enum class OutputChannel {TERMINAL, ERRCHANNEL};

class XMLreaderOutput{
public:
  // printing the output of the xmlReader regarding readOrWarn
  template<typename ParameterType>
  void parameterReading(std::vector<std::string> parameters, ParameterType& var,
      bool defaultAvailable, bool exitIfMissing, bool showWarning) const;
  void loadFile(bool loadOK, std::string fName) const;
  void readValue(bool warningsOn, std::string name, std::string fName) const;
  void printWarning(std::string name, std::string typeName, std::string value,
      bool verboseOn, bool exitIfMissing) const;
  template<typename XMLreaderType>
  void print(int indent, XMLreaderType& xmlReader) const;


  XMLreaderOutput() : clout(std::cerr, "xmlReaderOutput"){}
  XMLreaderOutput(std::ostream& stream) : clout(stream, "xmlReaderOutput"){}
  XMLreaderOutput(OutputChannel outputChannel);

protected:
  mutable OstreamManager clout;
  bool _inDebugMode = false;
};

XMLreaderOutput::XMLreaderOutput(OutputChannel outputChannel) : clout(std::cerr, "xmlReaderOutput"){
  if (outputChannel == OutputChannel::TERMINAL){
    clout = OstreamManager(std::cout, "xmlReaderOutput");
  }
  else if (outputChannel == OutputChannel::ERRCHANNEL){
    clout = OstreamManager(std::cerr, "xmlReaderOutput");
  } else{
    clout = OstreamManager(std::cerr, "xmlReaderOutput");
  }
#ifdef OLB_DEBUG
  _inDebugMode = true;
#endif
}

template<typename ParameterType>
void XMLreaderOutput::parameterReading(std::vector<std::string> parameters,
                                    ParameterType& var,
                                    bool defaultAvailable,
                                    bool exitIfMissing,
                                    bool showWarning) const
{
  if (showWarning || _inDebugMode) {
    clout << "Warning: Cannot read parameter from XML File: ";
    std::for_each(parameters.begin(), parameters.end(), [this](const std::string name_parameter)
      { clout << "<" << name_parameter << ">"; });
    clout << std::endl;

    if ( exitIfMissing ) {
      clout << "Error: This program cannot continue without ";
      std::for_each(parameters.begin(), parameters.end(), [this](const std::string name_parameter)
        { clout << "<" << name_parameter << ">"; });
      clout << ". Optimization aborted." << std::endl;
      exit(1);
    }
    if (defaultAvailable) {
      clout << "\t  Setting default value: " << parameters.back() << " = "<< var << std::endl;
    }
    else {
      clout << "\t  Setting arbitrarily: " << parameters.back() << " = " << var << std::endl;
    }
  }
}


void XMLreaderOutput::loadFile(bool loadOK, std::string fName) const
{
  if (!loadOK) {
    clout << std::string("Problem processing input XML file ") << fName << std::endl;
  }
}

void XMLreaderOutput::readValue(bool warningsOn, std::string name, std::string fName) const
{
  if ( warningsOn ) {
    clout << "Warning: cannot read value from node \"" << name << "\"" << ", \"" << fName <<"\"" << std::endl;
  }
}

/// print warning if verbose mode is on and exit, if exItMissing is true
void XMLreaderOutput::printWarning(std::string name, std::string typeName,
  std::string value, bool verboseOn, bool exitIfMissing) const
{

  if ( verboseOn ) {
    clout << "Warning: Cannot read " << typeName << " value from XML element " << name << "." << std::endl;
    if ( ! value.empty() ) {
      clout << "         Setting default value = " << value << std::endl;
    }
  }
  if ( exitIfMissing ) {
    clout << "Error: This program cannot continue without \"" << name << "\". Optimization aborted." << std::endl;
    exit(1);
  }
}

/// printing the whole structure of the XMLreader
template<typename XMLreaderType>
void XMLreaderOutput::print(int indent, XMLreaderType& xmlReader) const{
  std::string indentStr(indent, ' ');
  clout << indentStr << "[" << xmlReader.getName() << "]" << std::endl;
  if (!xmlReader.getText().empty()) {
    clout << indentStr << "  " << xmlReader.getText() << std::endl;
  }
  for (unsigned int iNode=0; iNode<xmlReader._children.size(); ++iNode) {
    xmlReader._children[iNode]->_output.print(indent+2,*xmlReader._children[iNode]);
  }
}

} // namespace olb

#endif
