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


#ifndef PLAIN_WRITER_H
#define PLAIN_WRITER_H

namespace olb {

/// FileName class
class FileName{
  /// Befriend overloaded << operator
  template<class Ch, class Tr>
  friend auto& operator<<(std::basic_ostream<Ch, Tr>&, FileName&);
private:
  std::stringstream _fileName;
  std::string _suffix;
  std::string _fullname;
public:
  /// Constructor
  FileName( std::string baseName, std::string suffix=".dat" )
    : _fileName(baseName, std::ios_base::app | std::ios_base::out), _suffix(suffix)
  {}
  /// Add parameters
  template<typename T>
  std::string addParameter( std::string quantity, T value, int digits=6, int precision=2 ){
    _fileName << "_" << quantity;
    if constexpr(std::is_same_v<T,double>){
      _fileName << std::setw(digits) << std::setfill('0') << std::setprecision(precision) << value;
    } else {
      _fileName << std::setw(digits) << std::setfill('0') << value;
    }
    return _fileName.str();
  };
  //Return full name as string
  std::string str(){
    return _fileName.str()+_suffix;
  }
  //Return full name as const char* (necessary for disk operations)
  const char* c_str(){
    _fullname = str();
    return _fullname.c_str();
  }
};

//Operator << for FielName
template<class Ch, class Tr>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,FileName& fileName)
{
  return os << fileName.str();
}

//// Write functions for scalar and array data
/// - Assuming that data is only written by main processor!

/// Write scalar data (single core only)
template <typename ARRAYTYPE>
void writeScalarData( std::ofstream& dataWriterOpened,
                      std::string fullFileName, std::string headLine,
                      ARRAYTYPE& dataVector, int iE, int iE0 = 0 );
/// Write scalar data
template <typename ARRAYTYPE>
void writeScalarData( std::string fullFileName, std::string headLine,
                      ARRAYTYPE& dataVector, int iE, int iE0 = 0 );
/// Write scalar data (including sanity check)
template <typename ARRAYTYPE>
void writeScalarData( std::string fullFileName, std::vector<std::string>& headLineVector,
                      ARRAYTYPE& dataVector, int iT, int iTinit = 0 );


/// Write array data
void writeArrayData( std::string fullFileName, std::string headLine,
                     std::vector<std::string>& dataVector );
/// Write array data
template <typename ARRAYTYPE>
void writeArrayData( std::string fullFileName, std::string headLine,
                     std::vector<ARRAYTYPE>& dataVector );
/// Write array data (including sanity check)
template <typename ARRAYTYPE>
void writeArrayData( std::string fullFileName, std::vector<std::string>& headLineVector,
                     std::vector<ARRAYTYPE>& dataVector );

} //namespace olb

#endif
