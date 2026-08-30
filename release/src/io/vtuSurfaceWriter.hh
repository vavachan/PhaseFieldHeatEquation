/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Christoph Gaul
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

#ifndef VTU_SURFACE_WRITER_HH
#define VTU_SURFACE_WRITER_HH

#include "vtuSurfaceWriter.h"

#include <cmath>
#include <unordered_set>
#include <unordered_map>

namespace olb {

template <typename T>
VTUsurfaceWriter<T>::VTUsurfaceWriter(std::string const name, CuboidDecomposition3D<T>& cuboidDecomposition,
                                      LoadBalancer<T>& loadBalancer, bool binary)
    : _name(name)
    , _cuboidDecomposition(cuboidDecomposition)
    , _loadBalancer(loadBalancer)
    , _binary(binary)
    , clout(std::cout, "VTUsurfaceWriter3D")
{}

template <typename T>
void VTUsurfaceWriter<T>::addFunctor(AnalyticalF3D<T, T>& f)
{
  _functorsA.push_back(&f);
}

template <typename T>
void VTUsurfaceWriter<T>::addFunctor(SuperF3D<T, T>& f)
{
  _functors.push_back(&f);
}

template <typename T>
void VTUsurfaceWriter<T>::addPoints(std::vector<Vector<T, 3>>& new_positions)
{
  for (unsigned long i = 0; i < new_positions.size(); i++) {
    pos.push_back(new_positions[i]);
  }
}

template <typename T>
void VTUsurfaceWriter<T>::addPoint(Vector<T, 3>& new_position)
{
  pos.push_back(new_position);
}

template <typename T>
void VTUsurfaceWriter<T>::addConnectivity(std::vector<std::vector<int>>& connections)
{
  for (unsigned long i = 0; i < connections.size(); i++) {
    connectivity.push_back(connections[i]);
  }
}

///writes vtu files and assigns values
template <typename T>
void VTUsurfaceWriter<T>::write(std::size_t iT)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  for (SuperF3D<T, T>* f : _functors) {
    f->getSuperStructure().communicate();
  }
  std::string fullNamePVDmaster = singleton::directories().getVtkOutDir() + createFileName(_name) + ".pvd";
  std::string namePiece         = "data/" + createFileName(_name, iT, 0) + ".vtu";

  if (rank == 0) {
    // puts name of .vti piece to a .pvd file [fullNamePVD]
    // adds a namePiece to master.pvd file.
    // To do so we overwrite closePVD() and add new entry.
    dataPVDmaster(iT, 1, fullNamePVDmaster, namePiece);
  }

  //VTU data files

  std::string fullNameVTU = singleton::directories().getVtkOutDir() + "data/" + createFileName(_name, iT, 0) + ".vtu";
  if (rank == 0) {
    preambleVTU(fullNameVTU, pos.size());
  }
  //writes data arrays into vtu file. All ranks are needed for this
  this->dataArray(fullNameVTU);
  if (rank == 0) {
    closeVTU(fullNameVTU);
  }
}

template <typename T>
void VTUsurfaceWriter<T>::preambleVTU(const std::string& fullName, int num)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>" << std::endl << std::flush;
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << num << "\" NumberOfCells=\"" << connectivity.size() << "\">" << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template <typename T>
void VTUsurfaceWriter<T>::closeVTU(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template <typename T>
void VTUsurfaceWriter<T>::dataPVD(int iT, int i, const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << i << "\" " << "file=\"" << namePiece
       << "\"/>\n";
  fout.close();
}

template <typename T>
void VTUsurfaceWriter<T>::dataPVDmaster(int iT, int i, const std::string& fullNamePVDMaster,
                                        const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(), std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 from the end of file to overwrite closePVD
    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << i << "\" " << "file=\"" << namePiece
         << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  }
  else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template <typename T>
void VTUsurfaceWriter<T>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir() + createFileName(_name) + ".pvd";
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
  }
}

template <typename T>
void VTUsurfaceWriter<T>::preamblePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"1.0\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<Collection>\n";
  fout.close();
}

template <typename T>
void VTUsurfaceWriter<T>::closePVD(const std::string& fullNamePVD)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {

    std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
    if (!fout) {
      clout << "Error: could not open " << fullNamePVD << std::endl;
    }
    fout << "</Collection>\n";
    fout << "</VTKFile>\n";
    fout.close();
  }
}

template <typename T>
void VTUsurfaceWriter<T>::dataArray(const std::string& fullName)
{

  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  for (unsigned long i = 0; i < _functors.size(); i++) {
    writeFunctor(fullName, fout, *_functors[i]);
  }

  for (unsigned long i = 0; i < _functorsA.size(); i++) {
    writeAnalyticalFunctor(fullName, fout, *_functorsA[i]);
  }
  // master writes points and connectivity
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    fout << "</PointData>" << std::endl;
    fout << "<CellData /> " << std::endl;
    fout << "<Cells>" << std::endl;
    fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    // if there are no connections between the points, connectivity simply counts up
    if (connectivity.size() == 0) {
      for (int i = 0; i < 1; ++i) {
        for (unsigned long j = 0; j < pos.size(); ++j) {
          fout << j << " ";
        }
      }
    }
    // if there are connections, connectivity is written out
    else {
      for (int i = 0; i < 1; ++i) {
        for (unsigned long j = 0; j < connectivity.size(); ++j) {
          for (unsigned long k = 0; k < connectivity[j].size(); ++k) {
            fout << connectivity[j][k] << " ";
          }
          if (j % 2 == 0) {
            fout << std::endl;
          }
        }
      }
    }

    fout << "</DataArray>" << std::endl;
    fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    // if there are no connections, offsets just count up
    if (connectivity.size() == 0) {
      for (std::size_t i = 0; i < 1; ++i) {
        for (unsigned long j = 1; j <= connectivity.size(); ++j) {
          fout << j << " ";
          if (j % 6 == 0) {
            fout << std::endl;
          }
        }
      }
      // if connections are present, offset is multiplied by 3
    }
    else { //connectivity.size() != 0
      for (std::size_t i = 0; i < 1; ++i) {
        for (unsigned long j = 1; j <= connectivity.size(); ++j) {
          fout << j * 3 << " ";
          if (j % 6 == 0) {
            fout << std::endl;
          }
        }
      }
    }
    // Connected points have the type 5, unconnected points have the type 3
    fout << "</DataArray>" << std::endl;
    fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    if (connectivity.size() != 0) {
      for (unsigned long j = 0; j < connectivity.size(); ++j) {
        fout << 5 << " ";
        if (j % 6 == 0) {
          fout << std::endl;
        }
      }
    }
    else {
      for (unsigned long j = 0; j < pos.size(); ++j) {
        fout << 3 << " ";
        if (j % 6 == 0) {
          fout << std::endl;
        }
      }
    }
    // write out all coordinates of points
    fout << "</DataArray>" << std::endl;
    fout << "</Cells>" << std::endl;
    fout << "<Points>" << std::endl;
    writePosition(fout);
    fout << "</Points>" << std::endl;
    fout << "</Piece>" << std::endl;

    fout.close();
  }
}

template <typename T>
void VTUsurfaceWriter<T>::writePosition(std::ofstream& fout)
{
  fout << "<DataArray type=\"Float32\" Name=\"Point\" NumberOfComponents=\"" << 3 << "\">" << std::endl;

  int num = pos.size();
  for (int j = 0; j < num; ++j) {
    fout << pos[j][0] << " " << pos[j][1] << " " << pos[j][2] << " " << std::endl;
  }
  fout << "</DataArray>" << std::endl;
}

template <typename T>
void VTUsurfaceWriter<T>::writeFunctor(const std::string& fullName, std::ofstream& fout, SuperF3D<T, T>& f)
{
  AnalyticalFfromSuperF3D<T> interpolateF(f, false);
  writeAnalyticalFunctor(fullName, fout, interpolateF);
}

template <typename T>
void VTUsurfaceWriter<T>::writeAnalyticalFunctor(const std::string& fullName, std::ofstream& fout,
                                                 AnalyticalF3D<T, T>& f)
{
  std::size_t dimF = f.getTargetDim();

  //Header DataArray, ascii or Base64 binary
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  if (rank == 0) {
    if (!_binary) {
      fout << "<DataArray type=\"Float32\" Name=\"" << f.getName() << "\" NumberOfComponents=\"" << dimF << "\">"
           << std::endl;
    }
    else if (_binary) {
      fout << "<DataArray type=\"Float32\" Name=\"" << f.getName()
           << "\" format=\"binary\" encoding=\"base64\" NumberOfComponents=\"" << dimF << "\">" << std::endl;
    }
  }

  if (_binary) {
    fout.close();
  }
  std::size_t                 numPoints = pos.size();
  std::ofstream               ofstr(fullName.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  size_t                      binarySize = size_t(numPoints * dimF * sizeof(float));
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int                uintBinarySize = (unsigned int)binarySize;
  if (_binary) {
    sizeEncoder.encode(&uintBinarySize, 1);
  }
  Base64Encoder<float> dataEncoder(ofstr, numPoints);

  std::vector<T> vals(dimF * numPoints, 0);

  T* val = new T[dimF];

  for (std::size_t j = 0; j < numPoints; ++j) {
    if (_loadBalancer.isLocal(_cuboidDecomposition.getC(pos[j]).value())) {
      f(val, pos[j].data());

      for (std::size_t i = 0; i < dimF; ++i) {
        vals[i * numPoints + j] = val[i];
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  // reduce the values to the master, then master writes them to the file
  std::vector<T> valsCummulated(dimF * numPoints, 0);
  singleton::mpi().reduceVect(vals, valsCummulated, MPI_SUM);

  if (rank == 0) {
    for (std::size_t j = 0; j < numPoints; ++j) {
      for (std::size_t i = 0; i < dimF; ++i) {
        if (!_binary) {
          fout << valsCummulated[i * numPoints + j] << " ";
        }
        else {
          const float helper = float(valsCummulated[j]);
          dataEncoder.encode(&helper, 1);
        }
      }
    }
  }
#endif
#ifndef PARALLEL_MODE_MPI
  // write the values to the file
  for (std::size_t j = 0; j < numPoints; ++j) {
    for (std::size_t i = 0; i < dimF; ++i) {
      if (!_binary) {
        fout << vals[i * numPoints + j] << " ";
      }
      else {
        const float helper = float(vals[i * numPoints + j]);
        dataEncoder.encode(&helper, 1);
      }
    }
  }
#endif

  if (_binary) {
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  }
  if (rank == 0) {
    fout << "</DataArray>" << std::endl;
  }
}

template <typename T>
std::vector<std::vector<T>> VTUsurfaceWriter<T>::removeDuplicates(const std::vector<std::vector<T>>& points)
{
  struct VectorHash {
    std::size_t operator()(const std::vector<T>& v) const
    {
      return std::hash<T>()(v[0]) ^ std::hash<T>()(v[1]) ^ std::hash<T>()(v[2]);
    }
  };

  struct VectorEqual {
    bool operator()(const std::vector<T>& a, const std::vector<T>& b) const
    {
      return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }
  };
  std::unordered_set<std::vector<T>, VectorHash, VectorEqual> uniquePoints(points.begin(), points.end());
  return std::vector<std::vector<T>>(uniquePoints.begin(), uniquePoints.end());
}

template <typename T>
bool VTUsurfaceWriter<T>::areVectorsEqual(const std::vector<T>& v1, const std::vector<T>& v2, T epsilon)
{
  for (std::size_t i = 0; i < v1.size(); ++i) {
    if (std::fabs(v1[i] - v2[i]) > epsilon) {
      return false;
    }
  }
  return true;
}

template <typename T>
std::vector<std::vector<int>> VTUsurfaceWriter<T>::mapTrianglesToIndices(STLreader<T>& stlReader)
{
  STLmesh<T> mesh = stlReader.getMesh();

  struct VectorHash {
    std::size_t operator()(const std::vector<T>& v) const
    {
      return std::hash<T>()(v[0]) ^ std::hash<T>()(v[1]) ^ std::hash<T>()(v[2]);
    }
  };

  struct VectorEqual {
    bool operator()(const std::vector<T>& a, const std::vector<T>& b) const { return areVectorsEqual(a, b); }
  };

  std::unordered_map<std::vector<T>, int, VectorHash, VectorEqual> pointIndexMap;
  std::vector<std::vector<int>>                                    indexedTriangles;

  // Create a map from points to their indices
  for (std::size_t i = 0; i < pos.size(); ++i) {
    pointIndexMap[{pos[i][0], pos[i][1], pos[i][2]}] = i;
  }

  // Map each triangle's points to their indices
  for (const STLtriangle<T>& triangle : mesh.getTriangles()) {
    std::vector<int> indices(3);
    for (int j = 0; j < 3; ++j) {
      std::vector<T> pt = {triangle.point[j][0], triangle.point[j][1], triangle.point[j][2]};
      auto           it = pointIndexMap.find(pt);
      if (it != pointIndexMap.end()) {
        indices[j] = it->second;
      }
      else {
        // Handle the case where the point is not found (should not happen if removeDuplicates worked correctly)
        indices[j] = -1; // Invalid index
      }
    }
    if (indices[0] != -1 && indices[1] != -1 && indices[2] != -1) {
      indexedTriangles.push_back(indices);
    }
  }

  return indexedTriangles;
}

template <typename T>
void VTUsurfaceWriter<T>::addSTL(STLreader<T>& stlReader)
{
  STLmesh<T>                  mesh = stlReader.getMesh();
  std::vector<std::vector<T>> points;

  // Collect points from triangles
  for (const STLtriangle<T>& triangle : mesh.getTriangles()) {
    std::vector<T> pt = {triangle.point[0][0], triangle.point[0][1], triangle.point[0][2]};
    points.push_back(pt);
    pt = {triangle.point[1][0], triangle.point[1][1], triangle.point[1][2]};
    points.push_back(pt);
    pt = {triangle.point[2][0], triangle.point[2][1], triangle.point[2][2]};
    points.push_back(pt);
  }
  // Remove duplicate points
  points = removeDuplicates(points);

  // Add points to pos
  for (const std::vector<T>& point : points) {
    pos.push_back(Vector<T, 3>(point[0], point[1], point[2]));
  }

  std::vector<std::vector<int>> indexedTriangles = mapTrianglesToIndices(stlReader);
  for (const auto& tri : indexedTriangles) {
    connectivity.push_back(tri);
  }
}

} // end namespace olb
#endif
