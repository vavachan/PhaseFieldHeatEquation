/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Nicolas Hafen Mathias J. Krause
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

#ifndef PERMEABILITY_H
#define PERMEABILITY_H


namespace olb {


template<typename T, bool check=false>
T getConfinedPermeability( T latticeViscosity, T latticeRelaxationTime, T physDeltaX,
                           T charPhysLength, T physPermeability )
{
  T minPermeabiliy = physDeltaX * physDeltaX * latticeViscosity * latticeRelaxationTime;
  if constexpr (check){
    if ( physPermeability >= minPermeabiliy ){
      throw std::invalid_argument( "Error: Min permeability exceeded!" );
      return 1.0;
    } else {
      return 1. - minPermeabiliy / physPermeability;
    }
  } else {
    return 1. - minPermeabiliy / physPermeability;
  }

  __builtin_unreachable();
}

template<typename T, typename DESCRIPTOR, bool check=false>
T getConfinedPermeability( const UnitConverter<T,DESCRIPTOR>& converter, T physPermeability ){
  return getConfinedPermeability<T,check>(
    converter.getLatticeViscosity(), converter.getLatticeRelaxationTime(),
    converter.getPhysDeltaX(), converter.getCharPhysLength(),
    physPermeability );
}

template<typename T>
T getPhysPermeability( T latticeViscosity, T latticeRelaxationTime, T physDeltaX,
                           T charPhysLength, T confinedPermeability )
{
  T minPermeabiliy = physDeltaX * physDeltaX * latticeViscosity * latticeRelaxationTime;
  return minPermeabiliy / (1. - confinedPermeability);
}

template<typename T, typename DESCRIPTOR, bool check=false>
T getPhysPermeability( const UnitConverter<T,DESCRIPTOR>& converter, T confinedPermeability ){
  return getPhysPermeability<T>(
    converter.getLatticeViscosity(), converter.getLatticeRelaxationTime(),
    converter.getPhysDeltaX(), converter.getCharPhysLength(),
    confinedPermeability );
}

} // namespace olb

#endif
