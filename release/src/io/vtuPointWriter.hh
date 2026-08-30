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


#ifndef VTU_POINT_WRITER_HH
#define VTU_POINT_WRITER_HH

#include "vtuPointWriter.h"


namespace olb{


template<typename T, typename W>
VTUpointWriter<T,W,2>::VTUpointWriter( std::string const name, bool binary )
  : _name( name ), _binary( binary ),
    clout(std::cout, "VTUpointWriter2D")
{}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addFunctor( AnalyticalF2D<T,W>& f, const std::string& name ){
  f.getName() = name;
  functorsA.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addFunctor( AnalyticalF2D<T,W>& f ){
  functorsA.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addFunctor( SuperF2D<T,W>& f, const std::string& name ){
  f.getName() = name;
  functors.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addFunctor( SuperF2D<T,W>& f ){
  functors.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addPoints( std::vector<Vector<T,2>>& new_positions ){
  for( unsigned long i=0; i<new_positions.size(); i++ ){
    pos.push_back(new_positions[i]);
  }
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::addPoint( Vector<T,2>& new_position ){
  pos.push_back(new_position);
}

///writes vtu files and assigns values
template<typename T, typename W>
void VTUpointWriter<T,W,2>::write( std::size_t iT )
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif


for (SuperF2D<T,W>* f : functors) {
    f->getSuperStructure().communicate();
}


  if (rank == 0) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + "_master.pvd";
    std::string fullNamePVD = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName(_name, iT) + ".pvd";
    preamblePVD( fullNamePVD );         // timestep
    std::string namePiece =  "data/" + createFileName(_name, iT, 0) + ".vtu";
    // puts name of .vti piece to a .pvd file [fullNamePVD]

    dataPVD(iT, 1, fullNamePVD, namePiece);
    // adds a namePiece to master.pvd file.
    // To do so we overwrite closePVD() and add new entry.
    dataPVDmaster(iT, 1, fullNamePVDmaster, namePiece);
    closePVD(fullNamePVD);            // timestep
  } // master only
  if ( rank == 0){ //master only or different systems each
    std::string fullNameVTU = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName(_name, iT, rank) + ".vtu";
    preambleVTU(fullNameVTU, pos.size());
    //writes data arrays into vtu file

    this->dataArray( fullNameVTU );
  closeVTU(fullNameVTU);
}
}

///writes vtu files and assigns values
template<typename T, typename W>
void VTUpointWriter<T,W,2>::write( std::size_t iT, std::vector<Vector<T,2>>& new_positions )
{
  addPoints(new_positions);
  write(iT);
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::createMasterFile()
{
  std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                  + createFileName(_name) + "_master.pvd";
  preamblePVD(fullNamePVDmaster);
  closePVD(fullNamePVDmaster);

}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::preamblePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::closePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::preambleVTU(
  const std::string& fullName, int num)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>" << std::endl << std::flush;
  fout
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
      << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << num
       << "\" NumberOfCells=\"" << num << "\">"
       << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::closeVTU(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::dataPVD(int iT, int i,
    const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << i
       << "\" " << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::dataPVDmaster(int iT, int i,
    const std::string& fullNamePVDMaster, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD
    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" "
         << i << "\" " << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::dataArray(
    const std::string& fullName )
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  for(unsigned long i=0; i<functors.size(); i++){
    writeFunctor( fullName, fout, *functors[i]);

  }
  for(unsigned long i=0; i<functorsA.size(); i++){
    writeAnalyticalFunctor( fullName, fout, *functorsA[i]);
  }

  fout << "</PointData>" << std::endl;
  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  for ( int i=0; i < 1; ++i){
    for ( unsigned long j=0; j<pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  for ( int i=0; i < 1; ++i){
    for ( unsigned long j=1; j <= pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
    for ( unsigned long j=0; j < pos.size(); ++j){
      fout << 1 << " ";
    }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  writePosition( fout );
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;

  fout.close();
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::writePosition( std::ofstream& fout  ){
  fout << "<DataArray type=\"Float32\" Name=\"Point\" NumberOfComponents=\""<< 3
       << "\">" << std::endl;

    int num = pos.size();
    for ( int j=0; j < num; ++j){
      fout << pos[j][0] << " "<<pos[j][1] << " 0.0 ";
    }
  fout << "</DataArray>" << std::endl;
}


template<typename T, typename W>
void VTUpointWriter<T,W,2>::writeFunctor( const std::string& fullName, std::ofstream& fout , SuperF2D<T,W>& f ){

  AnalyticalFfromSuperF2D<T> interpolateF( f, true );
  writeAnalyticalFunctor( fullName, fout , interpolateF );

}

template<typename T, typename W>
void VTUpointWriter<T,W,2>::writeAnalyticalFunctor( const std::string& fullName, std::ofstream& fout , AnalyticalF2D<T,W>& f ){
  int dimF = f.getTargetDim();
  if(!_binary){
  fout << "<DataArray type=\"Float32\" Name=\""<< f.getName()<<"\" NumberOfComponents=\""<< dimF
       << "\">" << std::endl;
  }else if(_binary){
    fout << "<DataArray type=\"Float32\" Name=\""<< f.getName()<<"\" format=\"binary\" encoding=\"base64\" NumberOfComponents=\""<< dimF
       << "\">" << std::endl;
  }

  if(_binary){
     fout.close();
  }
  size_t num = pos.size();

  std::ofstream ofstr(fullName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary);
  size_t binarySize = size_t( num*dimF * sizeof(float));
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int) binarySize;
  if(_binary){
    sizeEncoder.encode(&uintBinarySize, 1);
  }
  Base64Encoder<float> dataEncoder(ofstr, num);
  T* val = new T[dimF]();
  for ( size_t j=0; j < num; ++j){
      T point[2]= {pos[j][0],pos[j][1]};
      for (int i = 0; i < dimF; i++){
        f(val,point);
        const float helper=float(val[i]);
        if(!_binary){
          fout << helper << " ";
        }
        else{
          dataEncoder.encode(&helper, 1);
        }
      }
  }
  delete[] val;

  ofstr.close();
  if(_binary){
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  }
  fout << "</DataArray>" << std::endl;

}




//same thing as specialization for 3D
template<typename T, typename W>
VTUpointWriter<T,W,3>::VTUpointWriter( std::string const name,
                                                 bool binary )
  : _name( name ), _binary( binary ),
    clout(std::cout, "VTUpointWriter3D")
{}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addFunctor( AnalyticalF3D<T,W>& f, const std::string& name ){
  f.getName() = name;
  functorsA.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addFunctor( AnalyticalF3D<T,W>& f ){
  functorsA.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addFunctor( SuperF3D<T,W>& f, const std::string& name ){
  f.getName() = name;
  functors.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addFunctor( SuperF3D<T,W>& f ){
  functors.push_back(&f);
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addPoints( std::vector<Vector<T,3>>& new_positions ){
  for( unsigned long i=0; i<new_positions.size(); i++ ){
    pos.push_back(new_positions[i]);
  }
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::addPoint( Vector<T,3>& new_position ){
  pos.push_back(new_position);
}

///writes vtu files and assigns values
template<typename T, typename W>
void VTUpointWriter<T,W,3>::write( std::size_t iT )
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  for (SuperF3D<T,W>* f : functors) {
      f->getSuperStructure().communicate();
  }


  if (rank == 0) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + "_master.pvd";
    std::string fullNamePVD = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName(_name, iT) + ".pvd";
    preamblePVD( fullNamePVD );         // timestep
    std::string namePiece =  "data/" + createFileName(_name, iT, 0) + ".vtu";
    // puts name of .vti piece to a .pvd file [fullNamePVD]

    dataPVD(iT, 1, fullNamePVD, namePiece);
    // adds a namePiece to master.pvd file.
    // To do so we overwrite closePVD() and add new entry.
    dataPVDmaster(iT, 1, fullNamePVDmaster, namePiece);
    closePVD(fullNamePVD);            // timestep
  } // master only
  if ( rank == 0){ //master only or different systems each
    std::string fullNameVTU = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName(_name, iT, rank) + ".vtu";
    preambleVTU(fullNameVTU, pos.size());
    //writes data arrays into vtu file
    this->dataArray( fullNameVTU );

    closeVTU(fullNameVTU);
}
}

///writes vtu files and assigns values
template<typename T, typename W>
void VTUpointWriter<T,W,3>::write( std::size_t iT, std::vector<Vector<T,3>>& new_positions  )
{
  addPoints(new_positions);
  write(iT);
}


template<typename T, typename W>
void VTUpointWriter<T,W,3>::preambleVTU(
  const std::string& fullName, int num)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>" << std::endl << std::flush;
  fout
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
      << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << num
       << "\" NumberOfCells=\"" << num << "\">"
       << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::closeVTU(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}


template<typename T, typename W>
void VTUpointWriter<T,W,3>::dataPVD(int iT, int i,
    const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << i
       << "\" " << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::dataPVDmaster(int iT, int i,
    const std::string& fullNamePVDMaster, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD
    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" "
         << i << "\" " << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::createMasterFile()
{
  std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                  + createFileName(_name) + "_master.pvd";
  preamblePVD(fullNamePVDmaster);
  closePVD(fullNamePVDmaster);

}


template<typename T, typename W>
void VTUpointWriter<T,W,3>::preamblePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::closePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::dataArray(
    const std::string& fullName )
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  for(unsigned long i=0; i<functors.size(); i++){
    writeFunctor( fullName, fout, *functors[i]);
  }
  for(unsigned long i=0; i<functorsA.size(); i++){
    writeAnalyticalFunctor( fullName, fout, *functorsA[i]);
  }
  fout << "</PointData>" << std::endl;
  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  for ( size_t i=0; i < 1; ++i){
    for ( unsigned long j=0; j<pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  for ( size_t i=0; i < 1; ++i){
    for ( unsigned long j=1; j <= pos.size(); ++j){
      fout << j << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
    for ( unsigned long j=0; j < pos.size(); ++j){
      fout << 1 << " ";
    }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  writePosition( fout );
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;

  fout.close();
}

template<typename T, typename W>
void VTUpointWriter<T,W,3>::writePosition( std::ofstream& fout  ){
  fout << "<DataArray type=\"Float32\" Name=\"Point\" NumberOfComponents=\""<< 3
       << "\">" << std::endl;

    int num = pos.size();
    for ( int j=0; j < num; ++j){
      fout << pos[j][0] << " "<< pos[j][1] << " " << pos[j][2] << " ";
    }
  fout << "</DataArray>" << std::endl;
}


template<typename T, typename W>
void VTUpointWriter<T,W,3>::writeFunctor( const std::string& fullName, std::ofstream& fout , SuperF3D<T,W>& f ){
  AnalyticalFfromSuperF3D<T> interpolateF( f, true );
  writeAnalyticalFunctor( fullName, fout , interpolateF );
}


template<typename T, typename W>
void VTUpointWriter<T,W,3>::writeAnalyticalFunctor( const std::string& fullName, std::ofstream& fout , AnalyticalF3D<T,W>& f ){
  int dimF = f.getTargetDim();
  if(!_binary){
  fout << "<DataArray type=\"Float32\" Name=\""<< f.getName()<<"\" NumberOfComponents=\""<< dimF
       << "\">" << std::endl;
  }else if(_binary){
    fout << "<DataArray type=\"Float32\" Name=\""<< f.getName()<<"\" format=\"binary\" encoding=\"base64\" NumberOfComponents=\""<< dimF
       << "\">" << std::endl;
  }

  if(_binary){
     fout.close();
  }
  size_t num = pos.size();
  std::ofstream ofstr(fullName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary);
  size_t binarySize = size_t( num*dimF * sizeof(float));
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int) binarySize;
  if(_binary){
    sizeEncoder.encode(&uintBinarySize, 1);
  }

  Base64Encoder<float> dataEncoder(ofstr, num);
  T* val = new T[dimF]();
  for ( size_t j=0; j < num; ++j){
      T point[3]= {pos[j][0],pos[j][1],pos[j][2]};
      for (int i = 0; i < dimF; i++){
        f(val,point);
        if(!_binary){
          fout << val[i] << " ";
        }
        else{
          const float helper=float(val[i]);
          dataEncoder.encode(&helper, 1);
        }
      }
  }
  delete[] val;
  ofstr.close();
  if(_binary){
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  }
  fout << "</DataArray>" << std::endl;

}


}

#endif
