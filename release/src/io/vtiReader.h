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

/** \file
 * The VTI reader is able to read from VTI files and create and fill
 * corresponding data structures. The reading process starts with
 * the construction of the reader object. The name of the data type
 * to be read (e.g. "physVelocity") is mandatory.
 *
 *
 * Single cuboids (BlockVTIreader) and cuboid geometries
 * (SuperVTIreader) are supported in 2D and 3D.
 *
 * The Base reader class is reading the generic information of the
 * VTI file like origin, extend and number of nodes of the surrounding
 * cuboid as well as the size (dimension) of the data vector to be read.
 *
 * In case of cuboid geometries, the reader follows these steps:
 *
 * 1. Create a CuboidDecomposition and add all Cuboids.
 *
 * 2. Create a (Heuristic) LoadBalancer out of the CuboidDecomposition.
 *
 * 3. Create a SuperData object out of the CuboidDecomposition and
 *    LoadBalancer and by that also allocate all the neccessary memory
 *    for the actual data.
 *
 * 4. Iterate through the VTI data and fill the BlockData objects of SuperData.
 *
 *
 */


#ifndef VTI_READER_H
#define VTI_READER_H

#include <string>
#include <vector>

#include "io/xmlReader.h"
#include "geometry/cuboid.h"
#include "geometry/cuboidDecomposition.h"
#include "core/superData.h"
#include "core/blockData.h"
#include "io/ostreamManager.h"
#include "communication/loadBalancer.h"

//typedef double T;
namespace olb {

template<typename T> class LoadBalancer;
template<typename T, typename BaseType> class SuperData3D;
template<typename T, typename BaseType> class BlockData3D;

template<typename T> class CuboidGeometry2D;
// template<typename T> class Cuboid2D;
template<typename T, typename BaseType> class SuperData2D;
template<typename T, typename BaseType> class BlockData2D;

enum DataType {
  PointData,
  CellData
};

template<typename T>
class BaseVTIreader {
public:
  BaseVTIreader(const std::string& fName, int dim, std::string dName,
                const std::string class_name="BaseVTIreader", const DataType type = PointData);
  virtual ~BaseVTIreader() {};
  void printInfo();
protected:
  mutable OstreamManager clout;
  /* Dimension (2D or 3D) */
  int _dim;
  /* #Data array components */
  int _comps;
  /* Origin */
  std::vector<T> _origin;
  /* #Nodes */
  std::vector<int> _extent;
  T _delta;
  XMLreader _xmlReader;
  // Number of Cuboids
  int _nCuboids;
  DataType _type;

  /// Reads Extent from extAttrName from XML Tag and returns as vector
  std::vector<int> readExtent(const XMLreader* reader, std::string extAttrName);
  /// Converts 4D (or 6D) extents vector into 2D (3D) nb_nodes vector
  std::vector<int> getNbNodes(std::vector<int>& extents);
  /// Reads size from XML tag (attribute "NumberOfComponents")
  /// Reads the number of components from the XML DataArray Tag
  int getNbComps(const XMLreader& tagReader);
  ///Computes the total amount of data values
  ///Number of components * number of cells/point
  size_t computeTotalDataValues();
};


template<typename T, typename BaseType>
class BaseVTIreader3D : public BaseVTIreader<T> {
public:
  BaseVTIreader3D(const std::string& fileName, std::string dataName,
                  const std::string class_name="BaseVTIreader3D", const DataType type = PointData);
  ~BaseVTIreader3D() override {};
protected:
  /// Reads cuboid from piece node
  void readCuboid(Cuboid3D<T>& cuboid, XMLreader* pieceReader);
  /// Reads from DataArray and fills blockData
  bool readBlockData(BlockData<3,T,BaseType>& blockData, const XMLreader& pieceTagReader,
                     const std::string dataName);
private:
  ///Reads the ASCII encoded DataArray of a vti-file
  void readAsciiData(std::stringstream& stream_val, BlockData<3,T,BaseType>& blockData);
  ///Reads the binary encoded DataArray of a vti-file
  void readBinaryData(std::stringstream& stream_val, BlockData<3,T,BaseType>& blockData, size_t str_length);
};



template<typename T, typename BaseType>
class BlockVTIreader3D : public BaseVTIreader3D<T,BaseType> {
public:
  BlockVTIreader3D(const std::string& fileName, const std::string& dataName, const DataType type = PointData);
  ~BlockVTIreader3D() override {};
  BlockData<3,T,BaseType>& getBlockData();
  Cuboid3D<T>& getCuboid();
protected:
  Cuboid3D<T> _cuboid;
  BlockData<3,T,BaseType> _blockData;
};


template<typename T, typename BaseType>
class SuperVTIreader3D : public BaseVTIreader3D<T,BaseType> {
private:
  CuboidDecomposition3D<T>* _cGeometry;
  LoadBalancer<T>* _loadBalancer;
  SuperData<3,T,BaseType>* _superData;
public:
  SuperVTIreader3D(const std::string& fName, const std::string dName);
  ~SuperVTIreader3D() override;
  SuperData<3,T,BaseType>& getSuperData();
  CuboidDecomposition3D<T>& getCuboidDecomposition();
  LoadBalancer<T>& getLoadBalancer();
private:
  /// Reads Cuboid Geometry and creates Cuboid objects
  void readCuboidDecomposition();
  /// Fills all BlockData objects of SuperData
  void readSuperData(const std::string dName);
};










/// \todo implement 2D version above
/*
template<typename T>
class VTIreader2D : public XMLreader {
public:
  VTIreader2D();
  VTIreader2D(const std::string& fName);
  ~VTIreader2D();

  void getCuboid(Cuboid2D<T>& cuboid);
  void getCuboids(std::vector<Cuboid2D<T>* >& cuboids);
  //bool getScalarData(ScalarField2D<T>* base, const std::string dName);
  //bool getVectorData(TensorField2D<T, 2>* base, const std::string dName);
  //void getScalarMultiPieceData(std::vector<const ScalarFieldBase2D<T>* >& bases, const std::string dName);
  //void getVectorMultiPieceData(std::vector<const TensorFieldBase2D<T, 2>* >& bases, const std::string dName);
  void printInfo();
private:
  int _x0, _y0, _z0;
  int _x, _y, _z;
  T _delta;
};
*/

/* ------------------ BaseVTIreader2D -----------------*/

template<typename T, typename BaseType>
class BaseVTIreader2D : public BaseVTIreader<T> {
public:
  BaseVTIreader2D(const std::string& fileName, std::string dataName,
                  const std::string class_name="BaseVTIreader2D", const DataType type = PointData);
  ~BaseVTIreader2D() override {};
protected:
  /// Reads cuboid from piece node
  void readCuboid(Cuboid2D<T>& cuboid, XMLreader* pieceReader);
  /// Reads from DataArray and fills blockData
  bool readBlockData(BlockData<2,T,BaseType>& blockData, const XMLreader& pieceTagReader,
                     const std::string dataName);
private:
  ///Reads the ASCII encoded DataArray of a vti-file
  void readAsciiData(std::stringstream& stream_val, BlockData<2,T,BaseType>& blockData);
  ///Reads the binary encoded DataArray of a vti-file
  void readBinaryData(std::stringstream& stream_val, BlockData<2,T,BaseType>& blockData, size_t str_length);
};

/* ---------------- BlockVTIreader2D -------------------*/

template<typename T, typename BaseType>
class BlockVTIreader2D : public BaseVTIreader2D<T,BaseType> {
public:
  BlockVTIreader2D(const std::string& fileName, const std::string& dataName, const DataType type = PointData);
  ~BlockVTIreader2D() override {};
  BlockData<2,T,BaseType>& getBlockData();
  Cuboid2D<T>& getCuboid();
protected:
  Cuboid2D<T> _cuboid;
  BlockData<2,T,BaseType> _blockData;
};

} // namespace olb

#endif
