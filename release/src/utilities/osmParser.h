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

#ifndef OLB_OSM_PARSER_H
#define OLB_OSM_PARSER_H
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "tinyxml2.h"
#ifdef FEATURE_PROJ
#include "proj.h"
#else
  #error "The PROJ library is not enabled in the configuration under the FEATURE tag. To use this class, enable it by adding 'FEATURE+=PROJ'. Ensure that the Proj library is installed on your system. After installation, if necessary, add the appropriate path based on where it was installed."
#endif

using namespace std;
using namespace tinyxml2;
namespace olb{
  template<typename T>
  class OSMParser {
  public:
    // Structs for Node, Tree, and Building
    OSMParser(const std::string& fileName, Vector<T, 2> lonRange = Vector<T, 2>(0), Vector<T, 2> latRange = Vector<T, 2>(0));
    struct Node {
      T lat;
      T lon;
    };
    struct RoadStatistics
    {
      string name;
      T PM25_buildings;
      T PM25_cars;
      T NO2;
      T numberCars;

    };


    struct Tree {
      vector<Node> path;
      Node position;       // Position des Baums (x, y)
      string species;      // Art des Baums
      string leafType;
      T height=0;
    };

    struct Road{
      vector<Node> path;
      string type;
      string name;
      T width=0;
    };
    struct AreaVegitation {
      vector<Node> footprint;
      T height=0;
    };
    struct Building {
      vector<Node> footprint;
      T height=0;
      std::string name;
    };
    enum TYPE{
      BUILDING,
      TREES,
      AREAVEGITATION,
      ROADS
    };

    vector<Tree> trees; // Liste der Bäume
    vector<AreaVegitation> areaVegitations;
    vector<Building> buildings;
    vector<Road> roads;
    string filename;
    Vector<T,2> _lonRange;
    Vector<T,2> _latRange;
    Vector<T,2> _xRange;
    Vector<T,2> _yRange;
    bool _prevWriteIncrementalState;
    T _minEasting;
    T _minNorthing;

    Vector<T,3> _domainSize;
    // Parsing Funktionen
    void parseOSMWays(unordered_map<string, Node>& nodes);

    // Konvertierungsfunktionen
    void convertToUTMLocal();
    void convertToUTMLocalTrees(vector<Tree>& trees, bool buildingRef);
    void convertToUTMLocalBounding();
    Vector<T,3> getDomainSize(T maxHeight);
    // Geometrieerstellung
    void createGeometry(olb::SuperGeometry<T, 3>& superGeometry,TYPE, int matF,int matTo, T physDeltaX ,std::string objectName="");
    std::list<RoadStatistics> calculateStreetStatistics(std::map<std::string, std::pair<T, T>> nitrogenConcentration,  T physDeltaX);
    void generateBuildingExhaust(olb::SuperGeometry<T, 3>& superGeometry, T dimension,int matF, int matTo);
    std::string serializeFootprint(const std::vector<Node>& footprint);

  private:
    int getUTMZone(T lon);
    unordered_map<string, Node> parseOSMNodes();

    void getMapBounds();
    bool isBuildingType(const std::string& key, const std::string& value) {
      // Definiere die Menge der akzeptierten Gebäudetypen
      static const std::unordered_set<std::string> buildingTypes = {
        "yes", "apartments", "barracks", "bungalow", "cabin", "commercial",
        "detached", "dormitory", "farm", "ger", "hotel", "house", "houseboat",
        "industrial", "residential", "semidetached_house", "static_caravan",
        "terrace", "allotment_house", "farm_auxiliary", "barn", "bridge",
        "bunker", "carport", "cathedral", "chapel", "church", "civic",
        "college", "construction", "container", "cowshed", "digester",
        "farmhouse", "garages", "garbage_shed", "greenhouse", "hangar",
        "hospital", "hotel", "house", "hut", "kindergarten", "kiosk", "manor",
        "manufacture", "mosque", "office", "parking", "pavilion", "public",
        "religious", "residential", "retail", "roof", "ruins", "school",
        "service", "shed", "shop", "shrine", "stable", "static_caravan",
        "storage", "storage_tank", "substation", "supermarket", "synagogue",
        "temple", "train_station", "transformer_tower", "university",
        "warehouse", "works","commercial"
        };
    return key == "building" && buildingTypes.find(value) != buildingTypes.end();   };

  bool isTreeType(const std::string& key, const std::string& value) {
      // Definiere die Menge der akzeptierten Gebäudetypen
    static const std::unordered_set<std::string> naturalTagsForTrees = {
        "tree",       // Einzelner Baum
        "sapling",    // Jungbaum
    };

    return key == "natural" && naturalTagsForTrees.find(value) != naturalTagsForTrees.end();};

  bool isAreaVegitationType(const std::string& key, const std::string& value) {
      static const std::unordered_set<std::string> naturalTagsForVegitationAreas = {
        "wood",       // Waldgebiet
        "coppice",    // Niederwald (mehrstämmige Bäume oder Sträucher)
        "scrub",      // Busch- oder Strauchvegetation
        "heath",      // Heidelandschaft, oft mit kleineren Bäumen/Sträuchern
        "grassland",  // Wiese mit vereinzelten Bäumen
        "wetland"     // Feuchtgebiet mit typischer Vegetation, z. B. Mangroven
    };

    return key == "natural" && naturalTagsForVegitationAreas.find(value) != naturalTagsForVegitationAreas.end();};
};
}

#endif
