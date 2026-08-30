/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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


#ifndef UNITCONVERTER_HH
#define UNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "unitConverter.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::print(std::ostream& clout) const
{
  clout << "----------------- UnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                       N=              " << getResolution() << std::endl;
  clout << "Lattice velocity:                 latticeU=       " << getCharLatticeVelocity() << std::endl;
  clout << "Lattice relaxation frequency:     omega=          " << getLatticeRelaxationFrequency(  ) << std::endl;
  clout << "Lattice relaxation time:          tau=            " << getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):       charL=          " << getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):      charU=          " << getCharPhysVelocity() << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s): charNu=         " << getPhysViscosity() << std::endl;
  clout << "Phys. density(kg/m^d):            charRho=        " << getPhysDensity() << std::endl;
  clout << "Characteristical pressure(N/m^2): charPressure=   " << getCharPhysPressure() << std::endl;
  clout << "Mach number:                      machNumber=     " << getMachNumber() << std::endl;
  clout << "Reynolds number:                  reynoldsNumber= " << getReynoldsNumber() << std::endl;
  clout << "Knudsen number:                   knudsenNumber=  " << getKnudsenNumber() << std::endl;
  clout << "Characteristical CFL number:      charCFLnumber=  " << getCharCFLnumber() << std::endl;

  clout << std::endl;
  clout << "-- Conversion factors:" << std::endl;
  clout << "Voxel length(m):                  physDeltaX=     " << getConversionFactorLength() << std::endl;
  clout << "Time step(s):                     physDeltaT=     " << getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):             physVelocity=   " << getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):           physDensity=    " << getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                  physMass=       " << getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):          physViscosity=  " << getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                  physForce=      " << getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):           physPressure=   " << getConversionFactorPressure() << std::endl;

  clout << "-------------------------------------------------------------" << std::endl;

  if ( getLatticeRelaxationTime() < T(0.55) && getCharLatticeVelocity() > (T(8)*(getLatticeRelaxationTime() - T(0.5)) + T(1.e-8)) ) {
    T tauStable = getCharLatticeVelocity() / T(8) + T(0.5);
    T timeToCell = getConversionFactorTime() / getConversionFactorLength();
    if (getCharLatticeVelocity() >= T(0.3)) {
      tauStable = T(0.15) / T(8) + T(0.5);
      timeToCell = T(0.15) / getCharPhysVelocity();
    }
    T dxNew = timeToCell * getPhysViscosity() * descriptors::invCs2<T,DESCRIPTOR>() / (tauStable - T(0.5));
    T dtNew = dxNew * timeToCell;
    clout << "WARNING:" << std::endl;
    clout << "Potentially UNSTABLE combination of relaxation time (tau=" << getLatticeRelaxationTime() << ")" << std::endl;
    clout << "and characteristical CFL number (lattice velocity) charCFLnumber=" << getCharCFLnumber() << "!" << std::endl;
    clout << "Potentially maximum characteristical CFL number (maxCharCFLnumber=" << T(8)*(getLatticeRelaxationTime() - T(0.5)) << ")" << std::endl;
    clout << "Actual characteristical CFL number (charCFLnumber=" << getCharCFLnumber() << ") > " << T(8)*(getLatticeRelaxationTime() - T(0.5)) << std::endl;
    if (getCharLatticeVelocity() >= T(0.3)) {
      clout << "Please make the CFL number smaller than 0.3!" << std::endl;
    }
    clout << "Please reduce the the cell size or the time step size!" << std::endl;
    if (getCharLatticeVelocity() >= T(0.3)) {
      clout << "We recommend to use the cell size of " << dxNew << " m and the time step size of " << dtNew << " s (CFL = 0.15)." << std::endl;
    } else {
      clout << "We recommend to use the cell size of " << dxNew << " m and the time step size of " << dtNew << " s." << std::endl;
    }
    clout << "-------------------------------------------------------------" << std::endl;
  }

}

template <typename T, class DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::print() const
{
  print(clout);
}

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
{
  std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

  if (singleton::mpi().isMainProcessor()) {
    std::ofstream fout(dataFile.c_str(), std::ios::trunc);
    if (!fout) {
      clout << "error write() function: can not open std::ofstream" << std::endl;
    }
    else {
      print( fout );
      fout.close();
    }
  }
}

template<typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createUnitConverter");
  params.setWarningsOn(false);

  T physDeltaX{};
  T physDeltaT{};

  T charPhysLength{};
  T charPhysVelocity{};
  T physViscosity{};
  T physDensity{};
  T charPhysPressure = 0;

  int resolution{};
  T latticeRelaxationTime{};
  T charLatticeVelocity{};

  int counter = 0; // counting the number of Discretization parameters and returning a warning if more than 2 are provided

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["PhysParameters"]["CharPhysLength"].read(charPhysLength);
  params["Application"]["PhysParameters"]["CharPhysVelocity"].read(charPhysVelocity);
  params["Application"]["PhysParameters"]["PhysViscosity"].read(physViscosity);
  params["Application"]["PhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["PhysParameters"]["CharPhysPressure"].read(charPhysPressure);

  std::vector<std::string> discretizationParam = {"PhysDeltaX", "Resolution",
    "CharLatticeVelocity", "PhysDeltaT", "LatticeRelaxationTime"};

  for(int i = 0; i<discretizationParam.size(); i++){
    std::string test;
    if(params["Application"]["Discretization"][discretizationParam[i]].read(test,false)){
      counter++;
    }
  }
  if(counter>2){
    clout << "WARNING: More than 2 discretization parameters provided" << std::endl;
  }

  if (!params["Application"]["Discretization"]["PhysDeltaX"].read(physDeltaX,false)) {
    if (!params["Application"]["Discretization"]["Resolution"].read<int>(resolution,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaX, resolution or charLatticeVelocity
        throw std::runtime_error("Error: Have not found PhysDeltaX, Resolution or CharLatticeVelocity in XML file.");

      }
      else {
        // found charLatticeVelocity
        if (params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
          physDeltaX = charPhysVelocity / charLatticeVelocity * physDeltaT;
        }
        else if (params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          physDeltaX = physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5);
        }
        else {
          throw std::runtime_error("Error: Only found CharLatticeVelocity, missing PhysDeltaT or LatticeRelaxationTime");
        }
      }
    }
    else {
      // found resolution
      physDeltaX = charPhysLength / resolution;

      if (params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        latticeRelaxationTime = physViscosity * charLatticeVelocity * descriptors::invCs2<T,DESCRIPTOR>() * resolution + 0.5;
      }
      else {
        if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          throw std::runtime_error("Error: Have not found LatticeRelaxationTime and was not able to derive it using CharLatticeVelocity");
        }
      }
    }
  }
  // found physDeltaX
  if (!params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
    if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaT, latticeRelaxationTime and charLatticeVelocity
        throw std::runtime_error("Error: Have not found PhysDeltaT, LatticeRelaxationTime or CharLatticeVelocity in XML file.");
      }
      else {
        // found charLatticeVelocity
        physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
      }
    }
    else {
      // found latticeRelaxationTime
      physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
    }
  }

  return new UnitConverter<T, DESCRIPTOR>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physViscosity, physDensity, charPhysPressure);
}

}  // namespace olb

#endif
