#ifndef OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_
#define OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_

#include <fstream>
#include <iostream>
#include <unistd.h>

#include "io/ostreamManager.h"
#include "io/fileName.h"

#include "core/util.h"
#include "core/unitConverter.h"
#include "core/singleton.h"

namespace olb {
template<typename T, typename DESCRIPTOR>
class AdsorptionConverter : public UnitConverter<T, DESCRIPTOR> {
 protected:
  const T _physViscosity;
  const T _conversionViscosity;
  const T _particleConcentration;
  mutable OstreamManager clout;

 public:
  /** Converter for adsorption reactions primarily for ADE lattices based on concentration.
   * Renames viscosity to diffusivity and also takes fluid parameters.
   *
   * @param physDeltaX
   * @param physDeltaT
   * @param charPhysLength
   * @param charPhysVelocity
   * @param physDiffusivity
   * @param charConcentration
   * @param particleConcentration
   * @param fluidPhysViscosity
   */
  AdsorptionConverter(
      T physDeltaX,
      T physDeltaT,
      T charPhysLength,
      T charPhysVelocity,
      T physDiffusivity,
      T charConcentration,
      T particleConcentration,
      T fluidPhysViscosity) : UnitConverter<T, DESCRIPTOR>(physDeltaX,
                                                           physDeltaT,
                                                           charPhysLength,
                                                           charPhysVelocity,
                                                           physDiffusivity,
                                                           charConcentration),
                              _physViscosity(fluidPhysViscosity),
                              _conversionViscosity(physDeltaX * physDeltaX / physDeltaT),
                              _particleConcentration(particleConcentration),
                              clout(std::cout, "AdsConverter") {};

  // rename properties because in ADE viscosity has the function of diffusivity
  constexpr T
  getPhysDiffusivity() const {
    return UnitConverter<T, DESCRIPTOR>::getPhysViscosity();
  }
  constexpr T
  getLatticeDiffusivity() const {
    return UnitConverter<T, DESCRIPTOR>::getLatticeViscosity();
  }
  constexpr T
  getConversionFactorDiffusivity() const {
    return UnitConverter<T, DESCRIPTOR>::getConversionFactorViscosity();
  }

  constexpr T
  getPhysViscosity() const {
    return _physViscosity;
  }
  constexpr T
  getConversionFactorViscosity() const {
    return _conversionViscosity;
  }
  /// conversion factor to convert particle density from lattice units to kg/m^3
  constexpr T
  getConversionFactorParticleDensity() const {
    return _particleConcentration;
  }
  T getPhysParticleConcentration(T c) {
    return c * _particleConcentration;
  }

  T getPhysConcentration(T c) {
    return c * this->getConversionFactorDensity();
  }
  T getPhysLoading(T Cq) {
    return Cq * _particleConcentration;
  }


  constexpr T
  getReynoldsNumber() const {
    return this->_charPhysVelocity * this->_charPhysLength / this->_physViscosity;
  }
  constexpr T
  getSchmidtNumber() const {
    return  this->getPhysViscosity() / this->getPhysDiffusivity();
  }
  constexpr T
  getFourierNumber() const {
    return this->getPhysDiffusivity()*this->getPhysDeltaT() / this->getPhysDeltaX()/ this->getPhysDeltaX();
  }


  void print(std::ostream &clout) const {
    {
      clout << "----------------- UnitConverter information -----------------" << std::endl;
      clout << "-- Parameters:" << std::endl;
      clout << "Resolution:                       N=              " << this->getResolution() << std::endl;
      clout << "Lattice velocity:                 latticeU=       " << this->getCharLatticeVelocity() << std::endl;
      clout << "Lattice relaxation frequency:     omega=          " << this->getLatticeRelaxationFrequency() << std::endl;
      clout << "Lattice relaxation time:          tau=            " << this->getLatticeRelaxationTime() << std::endl;
      clout << "Characteristic length(m):         charL=          " << this->getCharPhysLength() << std::endl;
      clout << "Characteristic speed(m/s):        charU=          " << this->getCharPhysVelocity() << std::endl;
      clout << "Phys. concentration(kg/m^d):      charC=          " << this->getPhysDensity() << std::endl;
      clout << "Phys. diffusivity (m^2/s):        Diffusivity=    " << this->getPhysDiffusivity() << std::endl;
      clout << "Phys. viscosity (m^2/s):          Viscosity=      " << this->getPhysViscosity() << std::endl;
      clout << "Schmidt Number:                   Sc=             " << this->getSchmidtNumber() << std::endl;
      clout << "Reynolds number:                  reynoldsNumber= " << this->getReynoldsNumber() << std::endl;
      clout << "Mach number:                      machNumber=     " << this->getMachNumber() << std::endl;
      clout << "Fourier number:                   fourierNumber=  " << this->getFourierNumber() << std::endl;
      clout << std::endl;

      clout << "-- Conversion factors:" << std::endl;
      clout << "Voxel length(m):                  physDeltaX=     " << this->getConversionFactorLength() << std::endl;
      clout << "Time step(s):                     physDeltaT=     " << this->getConversionFactorTime() << std::endl;
      clout << "Velocity factor(m/s):             physVelocity=   " << this->getConversionFactorVelocity() << std::endl;
      clout << "Density factor(kg/m^3):           physDensity=    " << this->getConversionFactorDensity() << std::endl;
      clout << "Mass factor(kg):                  physMass=       " << this->getConversionFactorMass() << std::endl;
      clout << "Viscosity factor(m^2/s):          physViscosity=  " << this->getConversionFactorViscosity() << std::endl;
      clout << "Diffusivity factor(m^2/s):        physDiffusivity=" << this->getConversionFactorDiffusivity() << std::endl;
      clout << "Force factor(N):                  physForce=      " << this->getConversionFactorForce() << std::endl;
      clout << "Pressure factor(N/m^2):           physPressure=   " << this->getConversionFactorPressure() << std::endl;
      clout << "Particle density factor:                          " << this->getConversionFactorParticleDensity() << std::endl;
      clout << "-------------------------------------------------------------" << std::endl;

    }
  }

  void print() const {
    {
      this->print(clout);
    }
  }
};

template<typename T, typename DESCRIPTOR>
class AdsorptionConverterFromSchmidtNumberAndRelaxation : public AdsorptionConverter<T, DESCRIPTOR> {
 public:
  /** Calculate parameters from nondimensionalized values
   *
   * @param schmidtNumber
   * @param reynoldsNumber
   * @param relaxationTimeADV
   * @param charPhysLength
   * @param particleLength
   * @param resolution
   * @param charPhysVelocity
   * @param physDensity
   */
  constexpr AdsorptionConverterFromSchmidtNumberAndRelaxation(
      T schmidtNumber,
      T reynoldsNumber,
      T relaxationTimeADV,
      T charPhysLength,
      T particleLength,
      T resolution,
      T charPhysVelocity,
      T physDensity,
      T particleConcentration) : AdsorptionConverter<T, DESCRIPTOR>(charPhysLength / resolution,
                                                          (relaxationTimeADV - 0.5) / descriptors::invCs2<T,DESCRIPTOR>()
                                                              * util::pow((charPhysLength / resolution), 2)
                                                              / (charPhysVelocity * particleLength
                                                              / (schmidtNumber * reynoldsNumber)),
                                                          particleLength,
                                                          charPhysVelocity,
                                                          charPhysVelocity * particleLength
                                                              / (schmidtNumber * reynoldsNumber),
                                                          physDensity,
                                                          particleConcentration,
                                                          charPhysVelocity*particleLength/(reynoldsNumber)) {};

};
}
#endif // OLB_APPS_FLORIAN_ADSORPTIONCONVERTER_H_
