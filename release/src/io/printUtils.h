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



#ifndef PRINT_UTILS_H
#define PRINT_UTILS_H

#include <list>
#include <unordered_set>

namespace olb {

/// Print std::index_sequence
template<class Ch, class Tr, std::size_t... Is>
void print_index_sequence(std::basic_ostream<Ch,Tr>& os,
                          const std::index_sequence<Is...> is)
{
  std::size_t i=0;
  ((os << (i++==0? "" : ",") << Is), ...);
}



//// operator<< for different types (mainly std::)

/// Operator << for std::index_sequence
template<class Ch, class Tr, std::size_t... Is>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 const std::index_sequence<Is...>& is)
{
  os << "<";
  print_index_sequence(os, is);
  return os << ">";
}

/// Operator << for std::vector
template<class Ch, class Tr, typename... args>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 std::vector<args...>& vec)
{
  os << "|[";
  for(auto it = vec.begin(); it != vec.end(); ++it){
    os << (it==vec.begin()? "" : ",") << *it;
  }
  return os << "]|";
}

/// Operator << for std::array
template<class Ch, class Tr, class T, std::size_t N>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 std::array<T,N>& array)
{
  os << "|[";
  for(auto it = array.begin(); it != array.end(); ++it){
    os << (it==array.begin()? "" : ",") << *it;
  }
  return os << "]|";
}

/// Operator << for std::list
template<class Ch, class Tr, typename... args>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 std::list<args...>& list)
{
  os << "|[";
  for(auto it = list.begin(); it != list.end(); ++it){
    os << (it==list.begin()? "" : ",") << *it;
  }
  return os << "]|";
}

/// Operator << for std::set
template<class Ch, class Tr, typename... args>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 std::set<args...>& set)
{
  os << "|[";
  for(auto it = set.begin(); it != set.end(); ++it){
    os << (it==set.begin()? "" : ",") << *it;
  }
  return os << "]|";
}

/// Operator << for std::unordered_set
template<class Ch, class Tr, typename... args>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 std::unordered_set<args...>& set)
{
  os << "|[";
  for(auto it = set.begin(); it != set.end(); ++it){
    os << (it==set.begin()? "" : ",") << *it;
  }
  return os << "]|";
}

//// other convenience functions

/// Create readable bool string
/// - motivated by inconsistencies during ostream operator <<
template<typename O>
std::string boolToStr( O input ){
  //If scalar
  if constexpr (std::is_arithmetic<O>::value) {
    if (input){
      return "true";
    } else {
      return "false";
    }
  //If e.g. Vector
  } else {
    std::stringstream stream;
    stream << "[";
    for(unsigned iDim=0; iDim<O::d; ++iDim){
      stream << (iDim==0 ? "" : ",");
      if (input[iDim]){
        stream << "true";
      } else {
        stream << "false";
      }
    }
    stream << "]";
    return stream.str();
  }
}

} //namespace olb


#endif
