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

#ifndef VTU_POINT_WRITER_H
#define VTU_POINT_WRITER_H

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/analytical/interpolationF3D.h"
#include "functors/analytical/interpolationF2D.h"
#include "io/base64.hh"
#include "core/vector.h"
namespace olb {

template <typename T, typename W, int dim>
class VTUpointWriter;

template<typename T, typename W>
class VTUpointWriter<T,W,2>{
protected:
  std::string _name;
  std::vector<SuperF2D<T,W>*> functors;
  std::vector<AnalyticalF2D<T,W>*> functorsA;
  std::vector<Vector<T,2>> pos;
  bool _binary;
  mutable OstreamManager clout;

  //  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  //  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName, int num);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, int i, const std::string& fullNamePVD,
               const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int i, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  void dataArray( const std::string& fullName );
  ///  interpolates and writes functors stored at functors
  void writeFunctor( const std::string& fullName, std::ofstream& fout ,  SuperF2D<T,W>& f );
  void writeAnalyticalFunctor( const std::string& fullName, std::ofstream& fout , AnalyticalF2D<T,W>& f );
  /// writes coordinates of points to plot
  void writePosition( std::ofstream& fout );

public:
  VTUpointWriter(  std::string const name,
                            bool binary = false);

  void createMasterFile();
  /// write function to call during runtime, also accepts additional points during call
  void write( std::size_t iT, std::vector<Vector<T,2>>& new_positions );
  void write(std::size_t iT );
  /// accepts functors to write out
  void addFunctor( AnalyticalF2D<T,W>& f, const std::string& name );
  void addFunctor( AnalyticalF2D<T,W>& f );
  void addFunctor( SuperF2D<T,W>& f, const std::string& name );
  void addFunctor( SuperF2D<T,W>& f );
  /// functions to add multiple or single points.
  void addPoint( Vector<T,2>& new_position );
  void addPoints( std::vector<olb::Vector<T,2>>& new_positions );
};

template<typename T, typename W>
class VTUpointWriter<T,W,3>{
protected:
  std::string _name;
  std::vector<SuperF3D<T,W>*> functors;
  std::vector<AnalyticalF3D<T,W>*> functorsA;
  std::vector<Vector<T,3>> pos;
  bool _binary;
  mutable OstreamManager clout;

  //  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  //  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName, int num);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, int i, const std::string& fullNamePVD,
               const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int i, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  void dataArray( const std::string& fullName );
  ///  interpolates and writes functors stored at functors
  void writeFunctor( const std::string& fullName, std::ofstream& fout ,  SuperF3D<T,W>& f );
  void writeAnalyticalFunctor( const std::string& fullName, std::ofstream& fout , AnalyticalF3D<T,W>& f );
  /// writes coordinates of points to plot
  void writePosition( std::ofstream& fout );

public:
  VTUpointWriter(  std::string const name,
                            bool binary = false);

  void createMasterFile();
  /// write function to call during runtime, also accepts additional points during call
  void write(std::size_t iT, std::vector<Vector<T,3>>& new_positions );
  void write(std::size_t iT );
  /// accepts functors to write out
  void addFunctor( AnalyticalF3D<T,W>& f, const std::string& name );
  void addFunctor( AnalyticalF3D<T,W>& f );
  void addFunctor( SuperF3D<T,W>& f, const std::string& name );
  void addFunctor( SuperF3D<T,W>& f );
  /// functions to add multiple or single points.
  void addPoint( Vector<T,3>& new_position );
  void addPoints( std::vector<Vector<T,3>>& new_positions );
};


}

#endif
