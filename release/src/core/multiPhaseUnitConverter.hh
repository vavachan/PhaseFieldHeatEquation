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


#ifndef MULTIPHASEUNITCONVERTER_HH
#define MULTIPHASEUNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void MultiPhaseUnitConverter<T, DESCRIPTOR>::print() const
{
  clout << "----------------- MultiPhaseUnitConverter information ------------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                                 N=                              " << this->getResolution() << std::endl;
  clout << "Lattice velocity:                           latticeU=                       " << this->getCharLatticeVelocity() << std::endl;
  clout << "Lattice relaxation frequency:               omega=                          " << this->getLatticeRelaxationFrequency() << std::endl;
  clout << "Lattice relaxation time:                    tau=                            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):                 charL=                          " << this->getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):                charU=                          " << this->getCharPhysVelocity() << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s):           charNu=                         " << this->getPhysViscosity() << std::endl;
  clout << "Characteristical pressure(N/m^2):           charPressure=                   " << this->getCharPhysPressure() << std::endl;
  clout << "Reynolds number:                            reynoldsNumber=                 " << this->getReynoldsNumber() << std::endl;
  clout << "Phys. Delta X(m):                           physDeltaX=                     " << this->getPhysDeltaX() << std::endl;
  clout << "Phys. Delta T(s):                           physDeltaT=                     " << this->getPhysDeltaT() << std::endl;
  clout << "Phys. Surface Tension(J/m^2):               physSurfaceTension=             " << getPhysSurfaceTension() << std::endl;
  clout << "Characteristical Phys. Temperature(K):      charPhysTemperature=            " << getCharPhysTemperature() << std::endl;

  clout << "----------------- Conversion factors:-----------------------------------" << std::endl;
  clout << "Voxel length(m):                            physDeltaX=                     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                               physDeltaT=                     " << this->getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):                       physVelocity=                   " << this->getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):                     physDensity=                    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                            physMass=                       " << this->getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):                    physViscosity=                  " << this->getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                            physForce=                      " << this->getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):                     physPressure=                   " << this->getConversionFactorPressure() << std::endl;
  clout << "Equation of state a(Jm^3/mol^2):            physEoSa=                       " << getConversionFactorEoSa() << std::endl;
  clout << "Equation of state b(m^3/mol):               physEoSb=                       " << getConversionFactorEoSb() << std::endl;
  clout << "Molar mass(kg/mol):                         physMolarMass=                  " << getConversionFactorMolarMass() <<  std::endl;
  clout << "Surface tension(J/m^2):                     physSurfaceTension=             " << getConversionFactorSurfaceTension() << std::endl;
  clout << "Temperature(K):                             physTemperature=                " << getConversionFactorTemperature() << std::endl;

  clout << "------------------------------------------------------------------------" << std::endl;

}

template <typename T, typename DESCRIPTOR>
void MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR>::print() const
{
  clout << "----------------- MultiPhaseUnitConverter information ------------------" << std::endl;
  clout << "--------------------- Parameters in lattice units: ---------------------" << std::endl;
  clout << "Resolution:                                 N=                              " << this->getResolution() << std::endl;
  clout << "Lattice velocity:                           latticeU=                       " << this->getCharLatticeVelocity() << std::endl;
  clout << "Lattice relaxation frequency:               omega=                          " << this->getLatticeRelaxationFrequency() << std::endl;
  clout << "Lattice relaxation time:                    tau=                            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "--------------------- Parameters in physical units: --------------------" << std::endl;
  clout << "Characteristical length(m):                 charL=                          " << this->getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):                charU=                          " << this->getCharPhysVelocity() << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s):           charNu=                         " << this->getPhysViscosity() << std::endl;
  clout << "Characteristical pressure(N/m^2):           charPressure=                   " << this->getCharPhysPressure() << std::endl;
  clout << "Reynolds number:                            reynoldsNumber=                 " << this->getReynoldsNumber() << std::endl;
  clout << "Phys. Delta X(m):                           physDeltaX=                     " << this->getPhysDeltaX() << std::endl;
  clout << "Phys. Delta T(s):                           physDeltaT=                     " << this->getPhysDeltaT() << std::endl;

  clout << "----------------- Conversion factors:-----------------------------------" << std::endl;
  clout << "Voxel length(m):                            physDeltaX=                     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                               physDeltaT=                     " << this->getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):                       physVelocity=                   " << this->getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):                     physDensity=                    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                            physMass=                       " << this->getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):                    physViscosity=                  " << this->getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                            physForce=                      " << this->getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):                     physPressure=                   " << this->getConversionFactorPressure() << std::endl;
  clout << "Surface tension factor(J/m^2):              physSurfaceTension=             " << getConversionFactorSurfaceTension() << std::endl;

  clout << "------------------------------------------------------------------------" << std::endl;

}

}  // namespace olb

#endif
