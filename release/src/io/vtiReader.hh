/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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


#ifndef VTI_READER_HH
#define VTI_READER_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cassert>

#include "geometry/cuboid.h"
#include "vtiReader.h"
#include "communication/heuristicLoadBalancer.h"

namespace olb {

/* ------------------ BaseVTIreader -----------------*/

template<typename T>
BaseVTIreader<T>::BaseVTIreader( const std::string& fileName, int dim,
                                 std::string dataName, const std::string class_name, DataType type)
  : clout(class_name), _dim(dim), _comps(0), _origin(_dim, 0), _extent(_dim, 0),
    _delta(0), _xmlReader(fileName), _nCuboids(0), _type(type)
{
  // Read WholeExtent (2 * _dim ints) from XML File and calculate _extent
  std::vector<int> wholeExtend = readExtent(&_xmlReader["ImageData"], "WholeExtent");
  _extent = getNbNodes(wholeExtend);

  // Read _delta
  std::stringstream stream_val_1(_xmlReader["ImageData"].getAttribute("Spacing"));
  T delta1;
  T delta2;
  T delta3;
  stream_val_1 >> delta1;
  stream_val_1 >> delta2;
  stream_val_1 >> delta3;
  T eps = std::numeric_limits<T>::epsilon();
  if( std::fabs(delta1 - delta2) > eps || std::fabs(delta1 - delta3) > eps ) {
    clout << " Error: The lattice has to be homogenous" << std::endl;
    exit(1);
  }
  _delta = delta1;

  // Read _origin
  std::stringstream stream_val_2(_xmlReader["ImageData"].getAttribute("Origin"));
  for (auto& origin_i : _origin) {
    stream_val_2 >> origin_i;
  }

  /// Read _comps but only from the FIRST dataName Data Tag and count cuboids (_nCuboids)
  for ( auto& pieceReader : _xmlReader["ImageData"] ) {
    if (pieceReader->getName() == "Piece") {
      if ( _nCuboids == 0 ) {
        if (_type == PointData) {
          for (auto& dataArrayReader : (*pieceReader)["PointData"]) {
            if (dataArrayReader->getAttribute("Name") == dataName && dataArrayReader->getName() == "DataArray") {
              _comps = this->getNbComps(*dataArrayReader);
            }
          }
        } else {
          for (auto& dataArrayReader : (*pieceReader)["CellData"]) {
            if (dataArrayReader->getAttribute("Name") == dataName && dataArrayReader->getName() == "DataArray") {
              _comps = this->getNbComps(*dataArrayReader);
            }
          }
        }
      }
      _nCuboids++;
    }
  }
}

template<typename T>
void BaseVTIreader<T>::printInfo()
{
  clout << "Information on VTIreader Data:" << std::endl;
  clout << "Origin: ";
  for (auto& origin_i : _origin) {
    clout << origin_i << " ";
  }
  clout << std::endl;

  clout << "Extend: ";
  for (auto& extend_i : _extent) {
    clout << extend_i << " ";
  }
  clout << std::endl;

  clout << "Spacing: " << _delta << std::endl;
}

template<typename T>
std::vector<int> BaseVTIreader<T>::readExtent(const XMLreader* reader, std::string extAttrName)
{
  // An extent is in the form of four (or six) integers "x0 x1 y0 y1 z0 z1"
  // each representing a node number
  std::stringstream extstr(reader->getAttribute(extAttrName));
  std::vector<int> extents;
  int tmp;
  for (int i = 0; i < 2 * _dim; ++i) {
    extstr >> tmp;
    extents.push_back(tmp);
  }
  return extents;
}

template<typename T>
std::vector<int> BaseVTIreader<T>::getNbNodes(std::vector<int>& extents)
{
  // Convert 4D (or 6D) extents vector into 2D (3D) extent vector
  std::vector<int> nNodes;
  for ( int i = 0; i < _dim; i++ ) {
    if (_type == PointData) nNodes.push_back(extents[2 * i + 1 ] - extents[2 * i ] + 1);
    else nNodes.push_back(extents[2 * i + 1] - extents[2 * i]);
  }
  return nNodes;
}

template<typename T>
int BaseVTIreader<T>::getNbComps(const XMLreader& tagReader)
{
  // read the NumberOfComponents-Attribute (VTI standard) and return as integer
  int comps = std::atoi((tagReader.getAttribute("NumberOfComponents")).c_str());
  if (comps == 0) {
    clout << "Warning: NumberOfComponents zero or not given!" << std::endl;
    clout << "Example: <DataArray Name='physVelocity' NumberOfComponents='3'" << std::endl;
    clout << "Setting to default of NumberOfComponents='1'" << std::endl;
    comps = 1;
  }
  return comps;
}

template<typename T>
size_t BaseVTIreader<T>::computeTotalDataValues() {
  size_t result = 1;
  for(int i = 0; i < _dim; i++) {
    result *= _extent[i];
  }
  return result * _comps;
}


/* -------------------- BaseVTIreader3D --------------------- */

template<typename T, typename BaseType>
BaseVTIreader3D<T,BaseType>::BaseVTIreader3D( const std::string& fileName, std::string dataName,
    const std::string class_name, const DataType type)
  : BaseVTIreader<T>(fileName, 3, dataName, class_name, type)
{
}

template<typename T, typename BaseType>
void BaseVTIreader3D<T,BaseType>::readCuboid(Cuboid3D<T>& cuboid, XMLreader* pieceReader)
{
  if (pieceReader->getName() == "Piece") {
    std::vector<int> extents = this->readExtent(pieceReader, "Extent");
    std::vector<int> extent = this->getNbNodes(extents);
    // int extents[i] is node number => multiply with _delta to get coordinate
    cuboid.init(extents[0] * this->_delta,
                extents[2] * this->_delta,
                extents[4] * this->_delta,
                this->_delta,
                //Why is this 0, 1, 2 and not 1, 3, 5?
                extent[0],
                extent[1],
                extent[2]);
  }
}

template<typename T, typename BaseType>
void BaseVTIreader3D<T, BaseType>::readAsciiData(std::stringstream& stream_val, BlockData<3, T, BaseType>& blockData)
{
  // Careful: respect ordering in VTI File
  for (int iz = 0; iz < blockData.getNz(); iz++) {
    for (int iy = 0; iy < blockData.getNy(); iy++) {
      for (int ix = 0; ix < blockData.getNx(); ix++) {
        for (unsigned iSize=0; iSize < blockData.getSize(); iSize++) {

          //Throws an exception, when the input stream has reached it's end, but the blockData is not completely initialized
          if(stream_val.eof()) {
            this->clout << " Error: End of input has been reached, but there is still data to be written" << std::endl;
            exit(1);
          }
          //Test if the vertex is inside of the physical extent
          if (blockData.isInside({ix, iy, iz})) {
            BaseType tmp;
            stream_val >> tmp;
            // write tmp into blockData
            blockData.get({ix, iy, iz}, iSize) = tmp;
          } else {
            this->clout << "Found point outside of the block structure" << std::endl;
            blockData.get({ix, iy, iz}, iSize) = 0;
          }
        }
      }
    }
  }
  //Test if string stream is empty
  if (!stream_val.eof()) {
    this->clout << " Error: There are still values left in the value stream but the Block Data is already filled" << std::endl;
    exit(1);
  }
}

template<typename T, typename BaseType>
void BaseVTIreader3D<T, BaseType>::readBinaryData(std::stringstream& stream_val, BlockData<3, T, BaseType>& blockData,
                                                  size_t str_length)
{
  Base64Decoder<BaseType> decoder(stream_val, str_length);
  size_t totalDataValues = this->computeTotalDataValues();
  std::vector<BaseType> values(totalDataValues);
  decoder.decode(values.data(), totalDataValues);

  size_t value_idx = 0;
  for (int iz = 0; iz < blockData.getNz(); iz++) {
    for (int iy = 0; iy < blockData.getNy(); iy++) {
      for (int ix = 0; ix < blockData.getNx(); ix++) {
        for (unsigned iSize=0; iSize < blockData.getSize(); iSize++) {
          //Throws an exception, when the input stream has reached it's end, but the blockData is not completely initialized
          if(value_idx >= totalDataValues) {
            this->clout << " Error: End of input has been reached, but there is still data to be written" << std::endl;
            exit(1);
          }
          //Test if the vertex is inside of the physical extent
          if (blockData.isInside({ix, iy, iz})) {
            // this->clout << "got value: " << values[value_idx] << std::endl;
            blockData.get({ix, iy, iz}, iSize) = values[value_idx];
          } else {
            this->clout << "Found point outside of the block structure" << std::endl;
            blockData.get({ix, iy, iz}, iSize) = 0;
          }
          value_idx++;
        }
      }
    }
  }
  //Test if string stream is empty
  if (value_idx < totalDataValues) {
    this->clout << " Error: There are still values left in the value stream but the Block Data is already filled" << std::endl;
    exit(1);
  }
}

template<typename T, typename BaseType>
bool BaseVTIreader3D<T,BaseType>::readBlockData(BlockData<3,T,BaseType>& blockData,
    const XMLreader& pieceTagReader, const std::string dataName)
{

  std::string arrayType = (this->_type == PointData) ? "PointData" : "CellData";

  // Iterate through all <DataArray> tags and take the one with the given Name attribute
  bool attributeFound = false;
  for (auto & dataArrayReader : pieceTagReader[arrayType]) {
    if (dataArrayReader->getAttribute("Name") == dataName && dataArrayReader->getName() == "DataArray") {
      attributeFound = true;
      std::string data_str;
      if (dataArrayReader->read(data_str)) {
        std::stringstream stream_val(data_str);
        if (dataArrayReader->getAttribute("format") == "binary") {
          this->clout << "Reading binary data" << std::endl;
          readBinaryData(stream_val, blockData, data_str.length());
        } else if (dataArrayReader->getAttribute("format") == "ascii") {
          this->clout << "Reading ASCII data" << std::endl;
          readAsciiData(stream_val, blockData);
        } else {
          this->clout << "Warning: The file format has to be either 'ascii' or 'binary', but it is '" << dataArrayReader->getAttribute("format") << "'."<< std::endl;
        }
      }
      return true;
    }
  }
  if(!attributeFound) {
    this->clout << "Could not find attribute with the given Name: '" << dataName << "'" << std::endl;
  }
  return false;
}

/* ---------------- BlockVTIreader3D -------------------*/

template<typename T,typename BaseType>
BlockVTIreader3D<T,BaseType>::BlockVTIreader3D(const std::string& fileName, const std::string& dataName, const DataType type)
  : BaseVTIreader3D<T,BaseType>(fileName, dataName, "BlockVTIreader3D", type),
    _cuboid(this->_origin, this->_delta, this->_extent),
    _blockData(_cuboid, 0, this->_comps)
{
  size_t pieceTagCounter = 0;
  for (XMLreader* child : this->_xmlReader["ImageData"]) {
    if(child->getName() == "Piece") {
      pieceTagCounter++;
    }
  }
  if(pieceTagCounter == 0) {
    this->clout << " Error: No piece tag could be found in " << fileName << std::endl;
    exit(1);
  }
  if(pieceTagCounter > 1) {
    this->clout << "Found several piece tags, only reads the first one" << std::endl;
  }
  // Only read the first <Piece> tag in the XML file
  this->readBlockData(_blockData, this->_xmlReader["ImageData"]["Piece"], dataName);
}


template<typename T,typename BaseType>
BlockData<3,T,BaseType>& BlockVTIreader3D<T,BaseType>::getBlockData()
{
  return _blockData;
}

template<typename T,typename BaseType>
Cuboid3D<T>& BlockVTIreader3D<T,BaseType>::getCuboid()
{
  return _cuboid;
}


/* --------------- SuperVTIreader3D ---------------------*/

template<typename T,typename BaseType>
SuperVTIreader3D<T,BaseType>::~SuperVTIreader3D()
{
  delete _loadBalancer;
  delete _cGeometry;
  delete _superData;
}

template<typename T,typename BaseType>
SuperVTIreader3D<T,BaseType>::SuperVTIreader3D(const std::string& fName, const std::string dName )
  : BaseVTIreader3D<T,BaseType>(fName, dName, "SuperVTIreader3D")
{
  this->clout << "Start reading \"" << fName << "\"... "
              << "(" << this->_nCuboids << " cuboids)" << std::endl;

  // Create CuboidDecomposition
  _cGeometry = new CuboidDecomposition3D<T> (this->_origin, this->_delta, this->_extent);

  this->clout << "* Reading Cuboid Geometry..." << std::endl;

  // Fill CuboidDecomposition
  readCuboidDecomposition();

  _cGeometry->printExtended();

  // Create LoadBalancer
  _loadBalancer = new HeuristicLoadBalancer<T> (*_cGeometry);

  // Create SuperData (allocation of the data, this->_size is already known!)
  _superData = new SuperData<3,T,BaseType> ( *_cGeometry, *_loadBalancer, 2, this->_size);

  this->clout << "* Reading BlockData..." << std::endl;

  // Fill data objects
  readSuperData(dName);

  this->clout << "VTI Reader finished." << std::endl;
}

template<typename T,typename BaseType>
void SuperVTIreader3D<T,BaseType>::readCuboidDecomposition()
{
  //_cGeometry->clearCuboids();
  for ( auto& piece : this->_xmlReader["ImageData"] ) {
    if (piece->getName() == "Piece") {
      std::vector<int> extents = this->readExtent(piece, "Extent");
      std::vector<int> extent = this->getNbNodes(extents);
      // int extent[i] is node number => multiply with _delta to get coordinate
      Cuboid3D<T> cuboid(extents[0] * this->_delta,
                         extents[2] * this->_delta,
                         extents[4] * this->_delta,
                         this->_delta,
                         extent[0],
                         extent[1],
                         extent[2]);
      _cGeometry->add(cuboid);
    }
  }
}

template<typename T,typename BaseType>
void SuperVTIreader3D<T,BaseType>::readSuperData(const std::string dName)
{
  int counter = 0;
  // Iterate over all <Piece> tags
  for (auto & piece : this->_xmlReader["ImageData"]) {
    if (piece->getName() == "Piece") {
      this->readBlockData(_superData->getBlock(counter), *piece, dName);
      counter++;
    }
  }
}

template<typename T,typename BaseType>
SuperData<3,T,BaseType>& SuperVTIreader3D<T,BaseType>::getSuperData()
{
  return *_superData;
}

template<typename T,typename BaseType>
CuboidDecomposition3D<T>& SuperVTIreader3D<T,BaseType>::getCuboidDecomposition()
{
  return *_cGeometry;
}

template<typename T,typename BaseType>
LoadBalancer<T>& SuperVTIreader3D<T,BaseType>::getLoadBalancer()
{
  return *_loadBalancer;
}


/* -------------------- BaseVTIreader2D --------------------- */

template<typename T, typename BaseType>
BaseVTIreader2D<T,BaseType>::BaseVTIreader2D( const std::string& fileName, std::string dataName,
    const std::string class_name, const DataType type)
  : BaseVTIreader<T>(fileName, 2, dataName, class_name, type)
{
}

template<typename T, typename BaseType>
void BaseVTIreader2D<T,BaseType>::readCuboid(Cuboid2D<T>& cuboid, XMLreader* pieceReader)
{
  if (pieceReader->getName() == "Piece") {
    std::vector<int> extents = this->readExtent(pieceReader, "Extent");
    std::vector<int> extent = this->getNbNodes(extents);
    // int extents[i] is node number => multiply with _delta to get coordinate
    cuboid.init(extents[0] * this->_delta,
                extents[2] * this->_delta,
                this->_delta,
                extent[0],
                extent[1]);
  }
}

template<typename T, typename BaseType>
void BaseVTIreader2D<T, BaseType>::readAsciiData(std::stringstream& stream_val, BlockData<2, T, BaseType>& blockData)
{
  // Careful: respect ordering in VTI File
  for (int iy = 0; iy < blockData.getNy(); iy++) {
    for (int ix = 0; ix < blockData.getNx(); ix++) {
      for (unsigned iSize=0; iSize < blockData.getSize(); iSize++) {

        //Throws an exception when the input stream has reached its end, but the blockData is not completely initialized
        if(stream_val.eof()) {
          this->clout << " Error: End of input has been reached, but there is still data to be written" << std::endl;
          exit(1);
        }
        //Test if the vertex is inside of the physical extent
        if (blockData.isInside({ix, iy})) {
          BaseType tmp;
          stream_val >> tmp;
          // write tmp into blockData
          blockData.get({ix, iy}, iSize) = tmp;
        } else {
          this->clout << "Found point outside of the block structure" << std::endl;
          blockData.get({ix, iy}, iSize) = 0;
        }
      }
    }
  }
  //Test if string stream is empty
  if (!stream_val.eof()) {
    this->clout << " Error: There are still values left in the value stream but the Block Data is already filled" << std::endl;
    exit(1);
  }
}

template<typename T, typename BaseType>
void BaseVTIreader2D<T, BaseType>::readBinaryData(std::stringstream& stream_val, BlockData<2, T, BaseType>& blockData, size_t str_length)
{
    // Print the contents of stream_val to the console
  // this->clout << "stream_val contents:\n" << stream_val.str() << std::endl;

  Base64Decoder<BaseType> decoder(stream_val, str_length);
  // stream_val is correct
  size_t totalDataValues = this->computeTotalDataValues();
  std::vector<BaseType> values(totalDataValues);
  decoder.decode(values.data(), totalDataValues);

  size_t value_idx = 0;
  for (int iy = 0; iy < blockData.getNy(); iy++) {
    for (int ix = 0; ix < blockData.getNx(); ix++) {
      for (unsigned iSize=0; iSize < blockData.getSize(); iSize++) {
        //Throws an exception when the input stream has reached its end, but the blockData is not completely initialized
        if(value_idx >= totalDataValues) {
          this->clout << " Error: End of input has been reached, but there is still data to be written" << std::endl;
          exit(1);
        }
        //Test if the vertex is inside of the physical extent
        if (blockData.isInside({ix, iy})) {
          // values[value_ids] are decoded wrongly
          // this->clout << values[value_idx] << std::endl;
          blockData.get({ix, iy}, iSize) = values[value_idx];
        } else {
          this->clout << "Found point outside of the block structure" << std::endl;
          blockData.get({ix, iy}, iSize) = 0;
        }
        value_idx++;
      }
    }
  }
  //Test if string stream is empty
  if (value_idx < totalDataValues) {
    this->clout << " Error: There are still values left in the value stream but the Block Data is already filled" << std::endl;
    exit(1);
  }
}

template<typename T, typename BaseType>
bool BaseVTIreader2D<T,BaseType>::readBlockData(BlockData<2,T,BaseType>& blockData, const XMLreader& pieceTagReader, const std::string dataName)
{
  std::string arrayType = (this->_type == PointData) ? "PointData" : "CellData";

  // Iterate through all <DataArray> tags and take the one with the given Name attribute
  bool attributeFound = false;
  for (auto & dataArrayReader : pieceTagReader[arrayType]) {
    if (dataArrayReader->getAttribute("Name") == dataName && dataArrayReader->getName() == "DataArray") {
      attributeFound = true;
      std::string data_str;
      if (dataArrayReader->read(data_str)) {
        std::stringstream stream_val(data_str);
        if (dataArrayReader->getAttribute("format") == "binary") {
          // this->clout << "Reading binary data" << std::endl;
          readBinaryData(stream_val, blockData, data_str.length());
        } else if (dataArrayReader->getAttribute("format") == "ascii") {
          // this->clout << "Reading ASCII data" << std::endl;
          readAsciiData(stream_val, blockData);
        } else {
          this->clout << "Warning: The file format has to be either 'ascii' or 'binary', but it is '" << dataArrayReader->getAttribute("format") << "'."<< std::endl;
        }
      }
      return true;
    }
  }
  if(!attributeFound) {
    this->clout << "Could not find attribute with the given Name: '" << dataName << "'" << std::endl;
  }
  return false;
}

/* ---------------- BlockVTIreader2D -------------------*/

template<typename T,typename BaseType>
BlockVTIreader2D<T,BaseType>::BlockVTIreader2D(const std::string& fileName, const std::string& dataName, const DataType type)
  : BaseVTIreader2D<T,BaseType>(fileName, dataName, "BlockVTIreader2D", type),
    _cuboid(this->_origin, this->_delta, this->_extent),
    _blockData(_cuboid, 0, this->_comps)
{
  size_t pieceTagCounter = 0;
  for (XMLreader* child : this->_xmlReader["ImageData"]) {
    if(child->getName() == "Piece") {
      pieceTagCounter++;
    }
  }
  if(pieceTagCounter == 0) {
    this->clout << " Error: No piece tag could be found in " << fileName << std::endl;
    exit(1);
  }
  if(pieceTagCounter > 1) {
    this->clout << "Found several piece tags, only reads the first one" << std::endl;
  }
  // Only read the first <Piece> tag in the XML file
  this->readBlockData(_blockData, this->_xmlReader["ImageData"]["Piece"], dataName);
}

template<typename T,typename BaseType>
BlockData<2,T,BaseType>& BlockVTIreader2D<T,BaseType>::getBlockData()
{
  return _blockData;
}

template<typename T,typename BaseType>
Cuboid2D<T>& BlockVTIreader2D<T,BaseType>::getCuboid()
{
  return _cuboid;
}

}
#endif
