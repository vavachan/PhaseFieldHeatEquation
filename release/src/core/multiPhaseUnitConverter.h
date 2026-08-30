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

/** \file
 * Unit conversion handling -- header file.
 */

#ifndef MULTIPHASEUNITCONVERTER_H
#define MULTIPHASEUNITCONVERTER_H


#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

// All OpenLB code is contained in this namespace.
namespace olb {



/** Conversion between physical and lattice units, as well as discretization for multiple component lattices.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For the basic units of length L, time T, mass M, particle number N and temperature theta,
* 5 elemental conversion factors are deduced from EoS parameters a,b and R, molar mass M and surface tension sigma.
* For multiple components, the conversion factors are computed for the lightest condensable component.
*
* TODO: Extend documentation for MultiPhaseUnitConverter
*/
template <typename T, typename DESCRIPTOR>
class MultiPhaseUnitConverter : public UnitConverter<T, DESCRIPTOR> {
public:
  /** Documentation of constructor:
    * TODO: Extend constructur documentation
    */
  constexpr MultiPhaseUnitConverter(
    size_t resolution,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physEoSa,
    T latticeEoSa,
    T physEoSb,
    T physMolarMass,
    T physSurfaceTension,
    T charPhysTemperature,
    T charPhysPressure )
    : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (charPhysLength/resolution)*util::pow(physEoSa/latticeEoSa, -0.5)*util::pow(physEoSb*10.5/1., 0.5)*util::pow(physMolarMass*1, 0.5),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physMolarMass*1/(physEoSb*10.5/1.),
        charPhysPressure),
      _conversionEoSa(physEoSa/latticeEoSa),
      _conversionEoSb(physEoSb*10.5/1.),
      _conversionMolarMass(physMolarMass*1),
      _conversionSurfaceTension(this->_conversionLength*_conversionEoSa*util::pow(_conversionEoSb, -2)),
      _conversionTemperature(_conversionEoSa/_conversionEoSb/_conversionGasConstant),
      _physEoSa(physEoSa),
      _physEoSb(physEoSb),
      _physMolarMass(physMolarMass),
      _physSurfaceTension(physSurfaceTension),
      _charPhysTemperature(charPhysTemperature),
      _latticeSurfaceTension( physSurfaceTension / _conversionSurfaceTension ),
      clout(std::cout,"MultiPhaseUnitConv")
  {
  };

  /// return characteristic temperature in physical units
  constexpr T getCharPhysTemperature(  ) const
  {
    return _charPhysTemperature;
  };
  /// return equation of state parameter a in physical units
  constexpr T getPhysEoSa(  ) const
  {
    return _physEoSa;
  };
  /// return equation of state parameter b in physical units
  constexpr T getPhysEoSb(  ) const
  {
    return _physEoSb;
  };
  /// return molar mass in physical units
  constexpr T getPhysMolarMass(  ) const
  {
    return _physMolarMass;
  };
  /// return surface tension in physical units
  constexpr T getPhysSurfaceTension(  ) const
  {
    return _physSurfaceTension;
  };
  /// return characteristic temperature in physical units
  constexpr T getPhysTemperature(  ) const
  {
    return _charPhysTemperature;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorEoSa() const
  {
    return _conversionEoSa;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorEoSb() const
  {
    return _conversionEoSb;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorMolarMass() const
  {
    return _conversionMolarMass;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorGasConstant() const
  {
    return _conversionGasConstant;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorSurfaceTension() const
  {
    return _conversionSurfaceTension;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorTemperature(  ) const
  {
    return _conversionTemperature;
  };

  /// return lattice surface tension for parameter fitting
  constexpr T getLatticeSurfaceTension(  ) const
  {
    return _latticeSurfaceTension;
  };

/// nice terminal output for conversion factors, characteristical and physical data
  void print() const override;



protected:
  // conversion factors
  const T _conversionEoSa;                      // kg m^5 / s^2 mol^2
  const T _conversionEoSb;                      // m^3 / mol
  const T _conversionMolarMass;                 // kg / mol
  const T _conversionGasConstant = 8.314462618; // J / mol K = kg m^2 / s^2 mol K
  const T _conversionSurfaceTension;            // J / m^2 = kg / s^2
  const T _conversionTemperature;               // K

  // physical units, e.g characteristic or reference values
  const T _physEoSa;                            // kg m^5 / s^2 mol^2
  const T _physEoSb;                            // m^3 / mol
  const T _physMolarMass;                       // kg / mol
  const T _physSurfaceTension;                  // J / m^2 = kg / s^2
  const T _charPhysTemperature;                 // K

  // lattice units, discretization parameters
  const T _latticeSurfaceTension;               // -

private:
  mutable OstreamManager clout;
};

/** Conversion between physical and lattice units, as well as discretization for multiple component lattices.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For the basic units of length L, time T, mass M, particle number N and temperature theta,
* 5 elemental conversion factors are deduced from EoS parameters a,b and R, molar mass M and surface tension sigma.
* For multiple components, the conversion factors are computed for the lightest condensable component.
*
* TODO: Extend documentation for MultiPhaseUnitConverter
*/
template <typename T, typename DESCRIPTOR>
class MultiPhaseUnitConverterFromRelaxationTime : public UnitConverter<T, DESCRIPTOR> {
public:
  /** Documentation of constructor:
    * TODO: Extend constructur documentation
    */
  constexpr MultiPhaseUnitConverterFromRelaxationTime(
    size_t resolution,
    T latticeRelaxationTime,
    T latticeDensity,
    T charPhysLength,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0 ) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * util::pow((charPhysLength/resolution),2) / physViscosity,
        charPhysLength,
        0,
        physViscosity,
        (physDensity/latticeDensity),
        charPhysPressure),
      _conversionSurfaceTension(
        this->_conversionDensity *
        this->_conversionViscosity * this->_conversionViscosity /
        this->_conversionLength ),
      _conversionChemicalPotential( this->_conversionVelocity * this->_conversionVelocity ),
      clout(std::cout,"MultiPhaseUnitConv")
  {
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorSurfaceTension() const
  {
    return _conversionSurfaceTension;
  };

  /// access (read-only) to private member variable
  constexpr T getConversionFactorChemicalPotential() const
  {
    return _conversionChemicalPotential;
  };

  /// compute relaxation time from physical viscosity
  constexpr T computeRelaxationTimefromPhysViscosity( T userViscosity ) const
  {
    return 0.5 + descriptors::invCs2<T,DESCRIPTOR>() *
                 userViscosity / this->_conversionViscosity;
  };

  /// compute lattice surface tension from physical one
  constexpr T computeLatticeSurfaceTension( T userSurfaceTension ) const
  {
    return userSurfaceTension / this->_conversionSurfaceTension;
  };

  /// compute Reynolds from kinematic viscosity
  constexpr T computeReynolds( T userVelocity, T userLength, T userViscosity ) const
  {
    return userVelocity * userLength / userViscosity;
  };

  /// compute Weber
  constexpr T computeWeber( T userVelocity, T userLength, T userSurfaceTension ) const
  {
    return userLength * userVelocity * userVelocity / userSurfaceTension;
  };

/// nice terminal output for conversion factors, characteristical and physical data
  void print() const override;

protected:
  // conversion factors
  const T _conversionSurfaceTension;            // J / m^2 = kg / s^2
  const T _conversionChemicalPotential;         // J / kg = m^2 / s^2

private:
  mutable OstreamManager clout;
};

}  // namespace olb

#endif
