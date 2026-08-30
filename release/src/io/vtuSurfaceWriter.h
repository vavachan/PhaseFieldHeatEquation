/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Christoph Gaul, Adrian Kummerlaender
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

#ifndef VTU_SURFACE_WRITER_H
#define VTU_SURFACE_WRITER_H

#include "functors/lattice/superBaseF3D.h"
#include "functors/analytical/interpolationF3D.h"
#include "io/base64.hh"
#include "io/stlReader.h"
#include "communication/mpiManager.h"

namespace olb {

template <typename T>
class VTUsurfaceWriter {
protected:
  std::string                       _name;
  std::vector<SuperF3D<T, T>*>      _functors;
  std::vector<AnalyticalF3D<T, T>*> _functorsA;
  std::vector<Vector<T, 3>>         pos;
  std::vector<std::vector<int>>     connectivity;
  CuboidDecomposition3D<T>&         _cuboidDecomposition;
  LoadBalancer<T>&                  _loadBalancer;
  bool                              _binary;

  OstreamManager clout;

  //  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  //  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName, int num);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, int i, const std::string& fullNamePVD, const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int i, const std::string& fullNamePVDMaster, const std::string& namePiece);
  void dataArray(const std::string& fullName);
  ///  interpolates and writes functors stored at functors
  void writeFunctor(const std::string& fullName, std::ofstream& fout, SuperF3D<T, T>& f);
  /// writes functors.
  void writeAnalyticalFunctor(const std::string& fullName, std::ofstream& fout, AnalyticalF3D<T, T>& f);
  /// writes coordinates of points
  void        writePosition(std::ofstream& fout);
  static bool areVectorsEqual(const std::vector<T>& v1, const std::vector<T>& v2, T epsilon = 1e-6);
  /// removes duplicates from points, needed for STL files
  std::vector<std::vector<T>>   removeDuplicates(const std::vector<std::vector<T>>& points);
  std::vector<std::vector<int>> mapTrianglesToIndices(STLreader<T>& stlReader);

public:
  VTUsurfaceWriter(std::string name, CuboidDecomposition3D<T>& cuboidDecomposition, LoadBalancer<T>& loadBalancer,
                   bool binary = false);

  void createMasterFile();
  void write(std::size_t iT);

  /// Add single location to output set
  void addPoint(Vector<T, 3>& new_position);
  /// Add multiple locations to output set
  void addPoints(std::vector<Vector<T, 3>>& new_positions);
  /// Define connectivity between previously-added output points, rendering them a surface
  void addConnectivity(std::vector<std::vector<int>>& connections);

  /// Add full STL surface as output surface
  void addSTL(STLreader<T>& stlReader);

  /// Schedule analytical functor for output
  void addFunctor(AnalyticalF3D<T, T>& f);
  /// Schedule lattice functor for output (interpolates)
  void addFunctor(SuperF3D<T, T>& f);

};

} // namespace olb

#endif
