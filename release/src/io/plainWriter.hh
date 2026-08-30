/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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


#ifndef PLAIN_WRITER_HH
#define PLAIN_WRITER_HH

namespace olb {

//Write scalar data (single core only)
template <typename ARRAYTYPE>
void writeScalarData( std::ofstream& dataWriterOpened,
                      std::string fullFileName, std::string headLine,
                      ARRAYTYPE& dataVector, int iE, int iEinit )
{
  //Write headline if first element
  if (iE == iEinit){
    dataWriterOpened << headLine << std::endl;
  }
  //Write Data
  dataWriterOpened << dataVector[0];
  for (unsigned int i=1; i<dataVector.size(); ++i){
    dataWriterOpened << " " << dataVector[i];
  }
  dataWriterOpened << std::endl;
}

//Write scalar data
template <typename ARRAYTYPE>
void writeScalarData( std::string fullFileName, std::string headLine,
                      ARRAYTYPE& dataVector, int iE, int iEinit )
{
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()){
#endif
    //Instantiate data writer
    std::ofstream dataWriter;
    dataWriter.open( fullFileName, std::ofstream::app );
    //Write scalar data
    writeScalarData( dataWriter, fullFileName, headLine,
                     dataVector, iE, iEinit );
    //Close File
    dataWriter.close();
#ifdef PARALLEL_MODE_MPI
  }
#endif
}

//Write scalar data (including sanity check)
template <typename ARRAYTYPE>
void writeScalarData( std::string fullFileName, std::vector<std::string>& headLineVector,
                      ARRAYTYPE& dataVector, int iT, int iTinit )
{
  //Perform sanity check
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()){
#endif
    if (headLineVector.size()!=dataVector.size()){
      std::cerr << "WARNING (" << fullFileName << "): DataVector does not match provided headline!" << std::endl;
    }
#ifdef PARALLEL_MODE_MPI
  }
#endif
  //Set up headLine string
  std::string headLineStringScalar;
  for ( unsigned int iQ=0; iQ<headLineVector.size(); ++iQ ){
    if (iQ>0){ headLineStringScalar+=" "; };
    headLineStringScalar += headLineVector[iQ];
  }
  //Call write scalar data
  writeScalarData( fullFileName, headLineStringScalar, dataVector, iT, iTinit );
}




//Write array data
void writeArrayData( std::string fullFileName, std::string headLine,
                     std::vector<std::string>& dataVector )
{
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()){
#endif
    //Instantiate data writer
    std::ofstream dataWriter;
    dataWriter.open( fullFileName );
    //Write headline
    dataWriter << headLine << std::endl;
    //Write Data
    for (unsigned int i=0; i<dataVector.size(); ++i){
      dataWriter << dataVector[i] << std::endl;
    }
    //Close File
    dataWriter.close();
#ifdef PARALLEL_MODE_MPI
  }
#endif
}

//Write array data
template <typename ARRAYTYPE>
void writeArrayData( std::string fullFileName, std::string headLine,
                     std::vector<ARRAYTYPE>& dataVector )
{
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()){
#endif
    //Instantiate data writer
    std::ofstream dataWriter;
    dataWriter.open( fullFileName );
    //Write scalar data for each array entry
    for (unsigned int i=0; i<dataVector.size(); ++i){
      writeScalarData( dataWriter, fullFileName, headLine, dataVector[i], i );
    }
    //Close File
    dataWriter.close();
#ifdef PARALLEL_MODE_MPI
  }
#endif
}

//Write array data (including sanity check)
template <typename ARRAYTYPE>
void writeArrayData( std::string fullFileName, std::vector<std::string>& headLineVector,
                     std::vector<ARRAYTYPE>& dataVector )
{
  //Perform sanity check
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()){
#endif
    if (dataVector.size()==0){
      std::cerr << "WARNING: DataVector is empty!" << std::endl;
    } else {
      if (headLineVector.size()!=dataVector[0].size()){
        std::cerr << "WARNING (" << fullFileName << "): DataVector does not match provided headline!" << std::endl;
      }
    }
#ifdef PARALLEL_MODE_MPI
  }
#endif
  //Set up headLine string
  std::string headLineStringArray;
  for ( unsigned int iQ=0; iQ<headLineVector.size(); ++iQ ){
    if (iQ>0){ headLineStringArray+=" "; };
    headLineStringArray += headLineVector[iQ];
  }
  //Call write array data
  writeArrayData( fullFileName, headLineStringArray, dataVector);
}

} //namespace olb

#endif
