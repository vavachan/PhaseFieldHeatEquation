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

#ifndef OLB_OSM_PARSER_HH
#define OLB_OSM_PARSER_HH

#include "osmParser.h"

using namespace std;
using namespace tinyxml2;
namespace olb{
template<typename T>
OSMParser<T>::OSMParser(const string& fileName, Vector<T, 2> lonRange, Vector<T, 2> latRange):filename(fileName)
{
  if(lonRange ==Vector<T, 2>(0) || latRange==Vector<T, 2>(0)){
    getMapBounds();
  }
  else
  {
    _lonRange = lonRange;
    _latRange = latRange;
  }
  unordered_map<string,typename OSMParser<T>::Node> nodes = parseOSMNodes();
  parseOSMWays(nodes);
}
template<typename T>
std::string OSMParser<T>::serializeFootprint(const std::vector<OSMParser<T>::Node>& footprint) {
  std::ostringstream oss;
  for (const auto& node : footprint) {
    oss << node.lat << "," << node.lon << ";";
  }
  return oss.str();
}
template<typename T>
void OSMParser<T>::getMapBounds()
{
  OstreamManager clout("test");
  XMLDocument doc;
  if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
    cerr << "Failed to load .osm file!" << endl;
  }

  XMLElement* root = doc.FirstChildElement("osm");
  if (!root) {
    cerr << "Invalid OSM file format!" << endl;
  }
  XMLElement* bounds = root->FirstChildElement("bounds");
  if (!bounds) {
    cerr << "Invalid OSM file format!" << endl;
  }
  else
  {
    _lonRange[0] = bounds->DoubleAttribute("minlon");
    _lonRange[1] = bounds->DoubleAttribute("maxlon");
    _latRange[0] = bounds->DoubleAttribute("minlat");
    _latRange[1] = bounds->DoubleAttribute("maxlat");
  }
  clout << _lonRange << std::endl;
  clout << _latRange << std::endl;
}

// Parsing Nodes
template<typename T>
unordered_map<string, typename OSMParser<T>::Node> OSMParser<T>::parseOSMNodes() {
  OstreamManager clout ("Nodes");
  unordered_map<string, Node> nodes;
  XMLDocument doc;
  if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
    cerr << "Failed to load .osm file!" << endl;
    return nodes;
  }

  XMLElement* root = doc.FirstChildElement("osm");
  if (!root) {
    cerr << "Invalid OSM file format!" << endl;
    return nodes;
  }

  XMLElement* element = root->FirstChildElement("node");
  while (element) {
    string id = element->Attribute("id");
    T lat = element->DoubleAttribute("lat");
    T lon = element->DoubleAttribute("lon");
    if(_lonRange[0]<=lon && lon <=_lonRange[1] && _latRange[0]<=lat && lat<=_latRange[1]){
      nodes[id] = {lat, lon};
    }
    element = element->NextSiblingElement("node");
  }

  return nodes;
}

// Parsing Ways, Buildings, Trees, and Roads
template<typename T>
void OSMParser<T>::parseOSMWays(unordered_map<string, typename OSMParser<T>::Node>& nodes) {
  OstreamManager clout("parseOSMWays");
  XMLDocument doc;
  if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
    cerr << "Failed to load .osm file!" << endl;
  }

  XMLElement* root = doc.FirstChildElement("osm");
  if (!root) {
    cerr << "Invalid OSM file format!" << endl;
  }

  // Parsing buildings and roads
  XMLElement* element = root->FirstChildElement("way");
  while (element) {
    bool isBuilding = false;
    bool isRoad = false;
    bool isVegitation=false;


    T height = 8.0;  // Default building height if unspecified
    vector<string> nodeRefs;

    XMLElement* ndElement = element->FirstChildElement("nd");
    while (ndElement) {
      string ref = ndElement->Attribute("ref");
      nodeRefs.push_back(ref);
      ndElement = ndElement->NextSiblingElement("nd");
    }

    XMLElement* tagElement = element->FirstChildElement("tag");
    string roadType;
    string roadName;
    string buildingName;
    T roadWidth = 5;
    while (tagElement) {
      string key = tagElement->Attribute("k");
      string value = tagElement->Attribute("v");

      // Check if the way represents a building
      if (isBuildingType(key, value) || (key=="building:part"&&value=="yes") ) {
        isBuilding = true;
      }
      if (key == "height") {
        height = stod(value);
      } else if (key == "building:levels") {
        height = stod(value) * 3.0;
      }
      if(key=="name" && isBuilding)
      {
        buildingName= value;
      }
      // Check if the way represents a road
      if (key == "highway") {
        isRoad = true;
        roadType = value;  // Store road type (e.g., "residential", "primary")
      }
      if(key == "name" && isRoad){
        roadName = value;
      }
      if(key == "width" && isRoad){
        roadWidth = stod(value);
      }
      if(isAreaVegitationType(key,value)){
        isVegitation=true;
      }

      tagElement = tagElement->NextSiblingElement("tag");
    }

  if (isBuilding) {
    Building building;
    bool inRange = false;

    for (const auto& ref : nodeRefs) {
      if (nodes.find(ref) != nodes.end()) {
        const auto& node = nodes[ref];
        if (_lonRange[0] <= node.lon && node.lon <= _lonRange[1] &&
            _latRange[0] <= node.lat && node.lat <= _latRange[1]) {
          building.footprint.push_back(node);
          inRange = true;
        }
      }
    }

    if (inRange) {
      building.height = height;
      building.name = buildingName.empty() ? "" : buildingName; // Assign the building name or leave it empty

      if (!building.name.empty()) {
        // Named building: group by name
        auto it = std::find_if(buildings.begin(), buildings.end(),
                              [&building](const Building& b) { return b.name == building.name; });

        if (it != buildings.end()) {
          // Merge footprints for the same building
          it->footprint.insert(it->footprint.end(), building.footprint.begin(), building.footprint.end());
        } else {
          buildings.push_back(building);
        }
      } else {
        // Unnamed building: group by footprint hash
        size_t footprintHash = std::hash<std::string>{}(serializeFootprint(building.footprint));

        auto it = std::find_if(buildings.begin(), buildings.end(),
                              [footprintHash, this](const Building& b) {
                                return !b.name.empty() &&
                                        std::hash<std::string>{}(this->serializeFootprint(b.footprint)) == footprintHash;
                              });


        if (it != buildings.end()) {
          // Merge footprints for the same hash
          it->footprint.insert(it->footprint.end(), building.footprint.begin(), building.footprint.end());
        } else {
          buildings.push_back(building);
        }
      }
    }
  }

    if (isVegitation) {
      AreaVegitation vegitation;
      bool inRange;
      for (const auto& ref : nodeRefs) {
        if (nodes.find(ref) != nodes.end()) {
          if(_lonRange[0]<=nodes[ref].lon && nodes[ref].lon <=_lonRange[1] && _latRange[0]<=nodes[ref].lat && nodes[ref].lat<=_latRange[1]){
            vegitation.footprint.push_back(nodes[ref]);
            inRange= true;
          }
        }
        if(inRange){
          vegitation.height = height;
          areaVegitations.push_back(vegitation);
          inRange=false;
        }
      }
    }
    // Process as a road if "isRoad" is true
    else if (isRoad) {
      Road road;
      road.type = roadType;
      road.name = roadName;
      road.width = roadWidth;
      //road.name = roadName;
      for (const auto& ref : nodeRefs) {
        if (nodes.find(ref) != nodes.end()) {
          road.path.push_back(nodes[ref]);
        }
      }
      roads.push_back(road);
    }

    element = element->NextSiblingElement("way");
  }

  // Parsing trees
  XMLElement* nodeElement = root->FirstChildElement("node");
  while (nodeElement) {
      Tree tree;
      bool isTree = false;
      // Safely parse latitude and longitude
      const char* lonStr = nodeElement->Attribute("lon");
      const char* latStr = nodeElement->Attribute("lat");
      if (lonStr && latStr) {
          tree.position.lon = stod(lonStr);
          tree.position.lat = stod(latStr);
      } else {
          // Skip the node if it doesn't have valid coordinates
          nodeElement = nodeElement->NextSiblingElement("node");
          continue;
      }

      // Process the tags
      XMLElement* tagElement = nodeElement->FirstChildElement("tag");
      while (tagElement) {
          const char* key = tagElement->Attribute("k");
          const char* value = tagElement->Attribute("v");

          if (key && value) {
              std::string keyStr = key;
              std::string valueStr = value;
              T height =0;

              if (isTreeType(keyStr, valueStr)) {
                isTree = true;

              } if (key == "leaf_type" && isTree) {
                tree.species = value;
              }
              else if(key =="height" && isTree){
                height = stod(value);
                tree.height = height;
              }
          }

          tagElement = tagElement->NextSiblingElement("tag");
      }

      // Add the tree to the list if it's valid
      if (isTree) {
        trees.push_back(tree);
      }

      nodeElement = nodeElement->NextSiblingElement("node");
  }

}


// UTM Zone Calculation
template<typename T>
int OSMParser<T>::getUTMZone(T lon) {
  return static_cast<int>(floor((lon + 180) / 6)) + 1;
}

// Convert Buildings to UTM Local Coordinates
template<typename T>
void OSMParser<T>::convertToUTMLocal() {
  if(buildings.size() > 0){
    OstreamManager clout("convertUTMLocal");
    PJ_CONTEXT *C = proj_context_create();
    int zone = getUTMZone(_lonRange[0]);
    clout << "UTM Zone: " << zone << endl;

    string projString = "EPSG:326" + to_string(zone);
    PJ *P = proj_create_crs_to_crs(C, "EPSG:4326", projString.c_str(), NULL);
    if (!P) {
      cerr << "Proj4 transformation setup failed!" << endl;
      return;
    }

    T minEasting = numeric_limits<T>::max();
    T minNorthing = numeric_limits<T>::max();
    for (auto& building : buildings) {
      for (auto& node : building.footprint) {
        PJ_COORD coord = proj_coord(node.lon, node.lat, 0, 0);
        PJ_COORD projected = proj_trans(P, PJ_FWD, coord);

        node.lon = projected.xy.x;
        node.lat = projected.xy.y;

        minEasting = min(minEasting, node.lon);
        minNorthing = min(minNorthing, node.lat);
      }
    }

    for (auto& building : buildings) {
      for (auto& node : building.footprint) {
        node.lon -= minEasting;
        node.lat -= minNorthing;
      }
    }
    for (auto& vegitaiton : areaVegitations) {
      for (auto& node : vegitaiton.footprint) {
        PJ_COORD coord = proj_coord(node.lon, node.lat, 0, 0);
        PJ_COORD projected = proj_trans(P, PJ_FWD, coord);

        node.lon = projected.xy.x;
        node.lat = projected.xy.y;

        //minEasting = min(minEasting, node.lon);
        //minNorthing = min(minNorthing, node.lat);
      }
    }
    for (auto& vegitation : areaVegitations) {
      for (auto& node : vegitation.footprint) {
        node.lon -= minEasting;
        node.lat -= minNorthing;
      }
    }

    //minEasting = numeric_limits<double>::max();
    //minNorthing = numeric_limits<double>::max();
    for (auto& road : roads) {
      for (auto& node : road.path) {
        PJ_COORD coord = proj_coord(node.lon, node.lat, 0, 0);
        PJ_COORD projected = proj_trans(P, PJ_FWD, coord);

        node.lon = projected.xy.x;
        node.lat = projected.xy.y;

        //minEasting = min(minEasting, node.lon);
        //minNorthing = min(minNorthing, node.lat);
      }
    }

    for (auto& road : roads) {
      for (auto& node : road.path) {
        node.lon -= minEasting;
        node.lat -= minNorthing;
      }
    }
    _minEasting = minEasting;
    _minNorthing = minNorthing;
    proj_destroy(P);
    proj_context_destroy(C);
  }
}

// Convert Trees to UTM Local Coordinates
template<typename T>
void OSMParser<T>::convertToUTMLocalTrees(vector<Tree>& trees, bool buildingRef) {
  if(trees.size()>0){
    PJ_CONTEXT *C = proj_context_create();
    int zone = getUTMZone(trees[0].position.lon);
    //cout << "UTM Zone: " << zone << endl;

    string projString = "EPSG:326" + to_string(zone);
    PJ *P = proj_create_crs_to_crs(C, "EPSG:4326", projString.c_str(), NULL);
    if (!P) {
      cerr << "Proj4 transformation setup failed!" << endl;
      return;
    }
    T minEasting = numeric_limits<T>::max();
    T minNorthing = numeric_limits<T>::max();
    if(buildingRef){
      minEasting=_minEasting;
      minNorthing = _minNorthing;
    }

    for (auto& tree : trees) {
      PJ_COORD coord = proj_coord(tree.position.lon, tree.position.lat, 0, 0);
      PJ_COORD projected = proj_trans(P, PJ_FWD, coord);

      tree.position.lon = projected.xy.x;
      tree.position.lat = projected.xy.y;
      if(!buildingRef){
        minEasting = min(minEasting, tree.position.lon);
        minNorthing = min(minNorthing, tree.position.lat);
      }
    }

    for (auto& tree : trees) {
      tree.position.lon -= minEasting;
      tree.position.lat -= minNorthing;
    }

    proj_destroy(P);
    proj_context_destroy(C);
  }
}
// Convert Bounding to UTM Local
template<typename T>
void OSMParser<T>::convertToUTMLocalBounding() {
    PJ_CONTEXT *C = proj_context_create();
    if (!C) {
        cerr << "Failed to create PROJ context!" << endl;
        return;
    }

    // Calculate UTM zone
    int zone = getUTMZone(_lonRange[0]);

    // Determine hemisphere
    bool isNorthernHemisphere = _latRange[0] >= 0;

    // Construct EPSG code
    string projString = isNorthernHemisphere ? "EPSG:326" + to_string(zone) : "EPSG:327" + to_string(zone);

    // Create transformation
    PJ *P = proj_create_crs_to_crs(C, "EPSG:4326", projString.c_str(), NULL);
    if (!P) {
        cerr << "Proj4 transformation setup failed!" << endl;
        proj_context_destroy(C);
        return;
    }

    // Define the bounding box corners
    PJ_COORD corners[4] = {
        proj_coord(_lonRange[0], _latRange[0], 0, 0), // Bottom-left
        proj_coord(_lonRange[1], _latRange[0], 0, 0), // Bottom-right
        proj_coord(_lonRange[1], _latRange[1], 0, 0), // Top-right
        proj_coord(_lonRange[0], _latRange[1], 0, 0)  // Top-left
    };

    // Project the corners into UTM coordinates
    for (int i = 0; i < 4; i++) {
        corners[i] = proj_trans(P, PJ_FWD, corners[i]);
    }

    // Calculate the range of the projected coordinates
    double minX = corners[0].xy.x, maxX = corners[0].xy.x;
    double minY = corners[0].xy.y, maxY = corners[0].xy.y;
    for (int i = 1; i < 4; i++) {
        if (corners[i].xy.x < minX) minX = corners[i].xy.x;
        if (corners[i].xy.x > maxX) maxX = corners[i].xy.x;
        if (corners[i].xy.y < minY) minY = corners[i].xy.y;
        if (corners[i].xy.y > maxY) maxY = corners[i].xy.y;
    }

    // Store the range
    _xRange[0] = minY;
    _xRange[1] = maxY;
    _yRange[0] = minX;
    _yRange[1] = maxX;

    // Clean up
    proj_destroy(P);
    proj_context_destroy(C);
}

template<typename T>
Vector<T,3> OSMParser<T>::getDomainSize(T maxHeight)
{
  convertToUTMLocalBounding();

  _domainSize[0] = _xRange[1]-_xRange[0];
  _domainSize[1] = _yRange[1]-_yRange[0];
  _domainSize[2] = maxHeight;
  return _domainSize;
}
template<typename T>
void OSMParser<T>::createGeometry(olb::SuperGeometry<T, 3>& superGeometry,TYPE type, int matF,int matTo,T physDeltaX, std::string objectName) {
  if(superGeometry.getWriteIncrementalVTKState())
  {
    superGeometry.setWriteIncrementalVTK(false);
    _prevWriteIncrementalState= true;
  }
  else
  {
    _prevWriteIncrementalState= false;
  }
  if(type==TYPE::BUILDING)
  {
    for (const auto& building : buildings) {
      vector<Vector<T, 4>> points;
      T height;
      for (const auto& node : building.footprint) {
        height = (building.height == 0) ? 8 : building.height;
        points.push_back({node.lat, node.lon, height, 0});
      }
      IndicatorPolygon3D<T> polygon(points);
      superGeometry.rename(matF, matTo, polygon);
    }
  }
  if(type==TYPE::TREES){
    for (const auto& tree : trees) {
      if (true){//(tree.species == "broadleaved")) {
        T height=10.;
        T radiusTrunk = 2;
        T radiusCrown = 4.5;
        if(tree.height!=0){
          height = tree.height;
        }
        //clout << "Tree height is " << height<<std::endl;
        Vector<T, 3> pos({tree.position.lat, tree.position.lon, 0});
        Vector<T, 3> pos2({tree.position.lat, tree.position.lon, height});

        Vector<T, 3> posSphere({tree.position.lat, tree.position.lon, height});

        IndicatorCylinder3D<T> trunk(pos, {0, 0, 1}, radiusTrunk, height);
        //IndicatorSphere3D<T> crown(posSphere, radiusCrown);
        IndicatorCylinder3D<T> crown(pos2,{0,0,1},radiusCrown,height);

        superGeometry.rename(matF, matTo, trunk);
        superGeometry.rename(matF, matTo, crown);
      }
    }
    if(type==TYPE::AREAVEGITATION){
      for (const auto& vegitation : areaVegitations) {
        vector<Vector<T, 4>> points;
        T height;
        for (const auto& node : vegitation.footprint) {
          height = (vegitation.height == 0) ? 8 : vegitation.height;
          points.push_back({node.lat, node.lon, height, -1.});
        }
        IndicatorPolygon3D<T> polygon(points);
        IndicatorTranslate3D<T> transPolygon({0,0,-1},polygon);
        superGeometry.rename(matF, matTo, transPolygon);
        superGeometry.rename(5, matTo, transPolygon);
      }
    }
  }
  if(type==TYPE::ROADS){
    for (const auto& road : roads) {
      if (objectName == road.name && road.path.size() >= 2||objectName=="") { // only roads with at least two points
        T roadWidth = road.width;
        vector<Vector<T, 4>> leftSide;
        vector<Vector<T, 4>> rightSide;

        for (size_t i = 0; i < road.path.size() - 1; ++i) {
          const auto& p1 = road.path[i];
          const auto& p2 = road.path[i + 1];

          // normal
          T dx = p2.lon - p1.lon;
          T dy = p2.lat - p1.lat;
          T length = sqrt(dx * dx + dy * dy);

          if (length < 1e-6) {
            //cerr << "Skipping degenerate segment: " << p1.lat << ", " << p1.lon << " to " << p2.lat << ", " << p2.lon << endl;
            continue;
          }

          T nx = -dy / length; // Normale x
          T ny = dx / length;  // Normale y

          leftSide.push_back({p1.lat + ny * roadWidth / 2, p1.lon + nx * roadWidth / 2, physDeltaX, 0});
          rightSide.push_back({p1.lat - ny * roadWidth / 2, p1.lon - nx * roadWidth / 2, physDeltaX, 0});

          if (i == road.path.size() - 2) { // add endpoints of last segments
            leftSide.push_back({p2.lat + ny * roadWidth / 2, p2.lon + nx * roadWidth / 2, physDeltaX, 0});
            rightSide.push_back({p2.lat - ny * roadWidth / 2, p2.lon - nx * roadWidth / 2, physDeltaX, 0});
          }
        }

        if (!leftSide.empty() && !rightSide.empty()) {
          //close polygons in original order
          vector<Vector<T, 4>> fullPolygon = leftSide;
          fullPolygon.insert(fullPolygon.end(), rightSide.rbegin(), rightSide.rend());

          IndicatorPolygon3D<T> polygon(fullPolygon);
          IndicatorTranslate3D<T> transPolygon({0,0,-1*physDeltaX},polygon);
          //superGeometry.rename(matF, matTo,1, transPolygon);
          superGeometry.rename(matF, matTo, transPolygon);
        }
      } else {
        //cerr << "Skipping road: " << road.name << " (insufficient points or not matching name)" << endl;
      }
    }
  }
  if(_prevWriteIncrementalState){
    superGeometry.setWriteIncrementalVTK(true);
  }
}

//DT specific
// Function to generate a 1x1x1 cube on top of the building roof
template<typename T>
void OSMParser<T>::generateBuildingExhaust(olb::SuperGeometry<T, 3>& superGeometry, T dimension,int matF, int matTo) {
  OstreamManager clout(std::cout, "GenerateBuildingExhaust");
  for (const auto& building : buildings) {
    // Skip if no footprint
    if (building.footprint.empty()) continue;
      if(building.name != ""){
        // Calculate the geometric center of the footprint
        T centerX = 0;
        T centerY = 0;
        T roofHeight = (building.height == 0) ? 8 : building.height;

        for (const auto& node : building.footprint) {
          centerX += node.lat; // Assuming latitude is x-coordinate
          centerY += node.lon; // Assuming longitude is y-coordinate
        }

        centerX /= building.footprint.size();
        centerY /= building.footprint.size();

        // Define the bounds of the 1x1x1 cube
        T cubeHalfSize = dimension/2;
        Vector<T, 3> cubeMin = {centerX - cubeHalfSize, centerY - cubeHalfSize, roofHeight};
        Vector<T, 3> cubeMax = {centerX + cubeHalfSize, centerY + cubeHalfSize, roofHeight + dimension};

        // Create the cube indicator
        IndicatorCuboid3D<T> cube(cubeMax - cubeMin, cubeMin + (cubeMax - cubeMin) / 2);

        // Assign the cube to a new material for exhaust
        superGeometry.rename(matF, matTo, cube);

        // Log the building name and exhaust placement
        clout << "Exhaust added for building: " << building.name
                  << " at (" << centerX << ", " << centerY << ", " << roofHeight + 0.5 << ")\n"<<std::endl;
      }
    }
}

template<typename T>
std::list<typename OSMParser<T>::RoadStatistics> OSMParser<T>::calculateStreetStatistics(std::map<std::string, std::pair<T, T>> nitrogenConcentrations, T physDeltaX) {
    OstreamManager clout("Road statistics");
    const T roadHeight = 3.0; // Height in meters
    const T emissionRate = 17.8; // Average car NOx emission in micrograms per kilometer (450 mg/km = 450,000 µg/km)
    const T pmEmissionPerCar = 1.2; // PM2.5 emission per car in µg/m (per meter of road length)
    std::list<RoadStatistics> roadstatistics_list;
    // Map to aggregate road statistics by street name
    std::map<std::string, std::pair<T, T>> roadStats; // {streetName: {totalLength, totalVolume}}

    // Aggregate all segments into road statistics
    for (const auto& road : roads) {
        T totalSegmentLength = 0.0;
        T roadWidth = road.width;

        // Sum up length for each segment of the road
        for (size_t i = 0; i < road.path.size() - 1; ++i) {
            const auto& p1 = road.path[i];
            const auto& p2 = road.path[i + 1];

            T dx = p2.lon - p1.lon;
            T dy = p2.lat - p1.lat;
            T length = sqrt(dx * dx + dy * dy);

            if (length > 1e-6) { // Avoid degenerate segments
                totalSegmentLength += length;
            }
        }

        // Skip segments with zero or negligible length
        if (totalSegmentLength < 1e-6) {
            continue;
        }

        // Calculate the volume of the current road segment
        T segmentVolume = totalSegmentLength * roadWidth * roadHeight;

        // Aggregate length and volume by street name
        if (roadStats.find(road.name) == roadStats.end()) {
            roadStats[road.name] = {0.0, 0.0};
        }
        roadStats[road.name].first += totalSegmentLength;
        roadStats[road.name].second += segmentVolume;
    }

    // Perform calculations for each aggregated road
    for (const auto& [streetName, stats] : roadStats) {
        T totalLength = stats.first;   // Total length in meters
        T totalVolume = stats.second; // Total volume in cubic meters

        // Fetch nitrogen concentration for this road
        auto it = nitrogenConcentrations.find(streetName);
        if (it == nitrogenConcentrations.end()) {
            //clout << "Skipping road: " << streetName << " (no nitrogen concentration specified)" << endl;
            continue;
        }
        T nitrogenConcentration = it->second.first;
        T measuredPMConcentration= it->second.second;

        // Calculate total nitrogen dioxide (NO2) from the given concentration and the total volume
        T totalNitrogenDioxide = nitrogenConcentration * totalVolume; // in micrograms (µg)

        // Estimate number of cars based on NO2 emissions
        T numberOfCars = totalNitrogenDioxide / (emissionRate * totalLength);

        // Calculate the total PM2.5 from cars (PM2.5 emission rate per car is given as 1.2 µg/m per meter of road length)
        T totalPMFromCars = numberOfCars * pmEmissionPerCar * totalLength/totalVolume;

        // Calculate the total PM2.5 in the air from the measured concentration
        T totalPM = measuredPMConcentration * totalVolume; // Total PM2.5 in micrograms (µg)

        // Calculate the amount of PM2.5 from buildings (the rest)
        T totalPMFromBuildings = measuredPMConcentration - totalPMFromCars;

        // // Output the results
        // clout << "Road: " << streetName
        //       << ", Total Length: " << totalLength
        //       << " meters, Total Volume: " << totalVolume
        //       << " m^3, Total NO2: " << totalNitrogenDioxide
        //       << " µg, Estimated Cars: " << numberOfCars
        //       << ", Total PM2.5: " << totalPM
        //       << " µg/m^3, PM2.5 from Cars: " << totalPMFromCars
        //       << " µg/m^3, PM2.5 from Buildings: " << totalPMFromBuildings
        //       << " µg" << endl;

        RoadStatistics roadStatistics;

        roadStatistics.name= streetName;
        roadStatistics.PM25_buildings=totalPMFromBuildings;
        roadStatistics.PM25_cars= totalPMFromCars;
        roadStatistics.numberCars=numberOfCars;

        roadstatistics_list.push_back(roadStatistics);

    }
  return roadstatistics_list;
}

}
#endif
