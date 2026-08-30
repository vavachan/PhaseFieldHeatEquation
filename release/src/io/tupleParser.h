/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Dennis Teutscher
 *
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

#include "core/superLattice.h"

#include <tinyxml2.h>  // Include TinyXML2 header
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

namespace olb {

template<typename T, typename DESCRIPTOR>
class DynamicsTupleParser {
private:
  std::string filename;
  std::vector<std::vector<std::string>> momentaLists;
  std::vector<std::string> equilibriaList;
  std::vector<std::string> collisionList;
  std::vector<int> indicators;

public:
  DynamicsTupleParser(const std::string& filename) : filename(filename) {
    readTupleFromXML();
  }

  DynamicsTupleParser() {}

  std::vector<std::string> generateDynamicsStrings() {
    std::vector<std::string> dynamicsList;
    for (size_t i = 0; i < momentaLists.size(); i++) {
      std::string dynamicsString = "dynamics::" + momentaLists[i][0] + momentaLists[i][1] + "," +
                                    momentaLists[i][2] + "," + momentaLists[i][3] + "," +
                                    momentaLists[i][4] + ">," + equilibriaList[i] + "," +
                                    collisionList[i] + ",Default>";
      dynamicsList.push_back(dynamicsString);
    }
    return dynamicsList;
  }

  std::vector<std::string> readTupleFromXML() {
    OstreamManager clout("Tupleparser");
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS) {
      std::cerr << "Failed to load file: " << filename << std::endl;
      return {};
    }

    tinyxml2::XMLElement* paramElement = doc.FirstChildElement("Param");
    if (!paramElement) {
      std::cerr << "No 'Param' element found in XML. Exiting.\n";
      return {};
    }

    tinyxml2::XMLElement* tuplesElement = paramElement->FirstChildElement("Methods");
    if (!tuplesElement) {
      std::cerr << "No 'Methods' element found in 'Param'. Exiting.\n";
      return {};
    }

    std::vector<std::string> dynamics;

    // Iterate through all Map elements
    for (tinyxml2::XMLElement* mapElement = tuplesElement->FirstChildElement("Map");
         mapElement; mapElement = mapElement->NextSiblingElement("Map")) {
      for (tinyxml2::XMLElement* dynamicElement = mapElement->FirstChildElement("Dynamic");
            dynamicElement; dynamicElement = dynamicElement->NextSiblingElement("Dynamic")) {
        std::string dynamicString;
        createDynamicString(dynamicElement, dynamicString);
        clout << dynamicString << std::endl;
        dynamics.push_back("dynamics::Tuple<" + dynamicString + ",Default>");
      }
    }
    return dynamics;
  }

  void processElement(tinyxml2::XMLElement* element, std::string& dynamicString) {
    if (!element) {
      return;
    }

    tinyxml2::XMLElement* childElement = element->FirstChildElement();
    if (childElement) {
      bool hasTextContent = false;
      while (childElement) {
        if (hasTextContent) {
          dynamicString += ",";
        }
        if (childElement->FirstChildElement()) {
          dynamicString += childElement->Value();
          dynamicString += "<";
          processElement(childElement, dynamicString);
          dynamicString += ">";
        } else {
          dynamicString += childElement->Value();
          hasTextContent = true;
        }
        childElement = childElement->NextSiblingElement();
      }
    }
  }

  void createDynamicString(tinyxml2::XMLElement* dynamicElement, std::string& dynamicString) {
    if (!dynamicElement) {
      return;
    }

    tinyxml2::XMLElement* momentaElement = dynamicElement->FirstChildElement("Momenta");
    tinyxml2::XMLElement* equilibriaElement = dynamicElement->FirstChildElement("Equilibria");
    tinyxml2::XMLElement* collisionElement = dynamicElement->FirstChildElement("Collision");

    if (momentaElement) {
      dynamicString += "Momenta<";
      processElement(momentaElement, dynamicString);
      dynamicString += ">,";
    }

    if (equilibriaElement) {
      dynamicString += equilibriaElement->FirstChild()->Value();
      dynamicString += ",";
    }

    if (collisionElement) {
      processElement(collisionElement, dynamicString);
    }

    if (!dynamicString.empty() && dynamicString.back() == ',') {
      dynamicString.pop_back();
    }
  }

  std::vector<int> readIndicatorFromXML() {
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS) {
      std::cerr << "Failed to load file: " << filename << std::endl;
      return {};
    }

    tinyxml2::XMLElement* paramElement = doc.FirstChildElement("Param");
    if (!paramElement) {
      std::cerr << "No 'Param' element found in XML. Exiting.\n";
      return {};
    }

    tinyxml2::XMLElement* tuplesElement = paramElement->FirstChildElement("Methods");
    if (!tuplesElement) {
      std::cerr << "No 'Methods' element found in 'Param'. Exiting.\n";
      return {};
    }

    // Iterate through all Map elements
    for (tinyxml2::XMLElement* mapElement = tuplesElement->FirstChildElement("Map");
          mapElement; mapElement = mapElement->NextSiblingElement("Map")) {
      tinyxml2::XMLElement* indicatorElement = mapElement->FirstChildElement("Indicator");
      if (indicatorElement) {
        const char* indicatorText = indicatorElement->GetText();
        if (indicatorText) {
          indicators.push_back(std::stoi(indicatorText));
        }
      }
    }
    return indicators;
  }

  std::map<std::string, float> readParameterFromXML() {
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS) {
      std::cerr << "Failed to load file: " << filename << std::endl;
      return {};
    }

    tinyxml2::XMLElement* paramElement = doc.FirstChildElement("Param");
    if (!paramElement) {
      std::cerr << "No 'Param' element found in XML. Exiting.\n";
      return {};
    }

    tinyxml2::XMLElement* methodsElement = paramElement->FirstChildElement("Methods");
    if (!methodsElement) {
      std::cerr << "No 'Methods' element found in 'Param'. Exiting.\n";
      return {};
    }

    std::map<std::string, float> parameterMap;

    // Iterate through all Map elements
    for (tinyxml2::XMLElement* mapElement = methodsElement->FirstChildElement("Map");
          mapElement; mapElement = mapElement->NextSiblingElement("Map")) {
      tinyxml2::XMLElement* dynamicElement = mapElement->FirstChildElement("Dynamic");
      if (!dynamicElement) {
        std::cerr << "No 'Dynamic' element found in 'Map'. Skipping.\n";
        continue;
      }

      tinyxml2::XMLElement* parameterElement = dynamicElement->FirstChildElement("Parameter");
      if (parameterElement) {
        // Iterate through all children of the Parameter element
        for (tinyxml2::XMLElement* childElement = parameterElement->FirstChildElement();
              childElement; childElement = childElement->NextSiblingElement()) {
          std::string tagName = childElement->Value();
          std::string tagValue = childElement->GetText() ? childElement->GetText() : "";
          try {
            float value = std::stof(tagValue);
            parameterMap[tagName] = value;
          } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid value for tag " << tagName << ": " << tagValue << std::endl;
          } catch (const std::out_of_range& e) {
            std::cerr << "Value out of range for tag " << tagName << ": " << tagValue << std::endl;
          }
        }
      }
    }
    return parameterMap;
  }
};

} // namespace olb
