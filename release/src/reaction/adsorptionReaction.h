/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

#ifndef OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_
#define OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_

#include <vector>
#include <cmath>
#include <iostream>
#include "core/adsorptionConverter.h"
#include "isotherm.h"

namespace olb {

struct SCALAR2 : public descriptors::FIELD_BASE<1,  0, 0> { };

template <typename ISOTHERM>
struct AdsorptionReaction {
  using isotherm_t = ISOTHERM;
  // K_F: film diffusion mass transfer coefficient (m/s)
  struct K_F : public descriptors::FIELD_BASE<1> { };
  // D_S: surface diffusion constant (m^2/s)
  struct D_S : public descriptors::FIELD_BASE<1> { };
  // C_0: initial concentration (mg/L)
  struct C_0 : public descriptors::FIELD_BASE<1> { };
  // R_P: particle radius (m)
  struct R_P : public descriptors::FIELD_BASE<1> { };
  // K_S: surface diffusion mass transfer coefficient (1/s)
  struct K_S : public descriptors::FIELD_BASE<1> { };
  // Q_0: equilibrium surface loading (mg/g)
  struct Q_0 : public descriptors::FIELD_BASE<1> { };
  // Conversion Factor Density
  struct CONV_DENS : public descriptors::FIELD_BASE<1> { };
  // Conversion Factor Particle Density
  struct CONV_PARC_DENS : public descriptors::FIELD_BASE<1> { };
  // External Mass Transfer Enabled
  struct EXT_MASS : public descriptors::FIELD_BASE<1> { };

  using parameters = typename ISOTHERM::parameters::template include<
    K_F,D_S,C_0,R_P,K_S,Q_0,CONV_DENS,CONV_PARC_DENS,EXT_MASS
  >;

  template <typename V, typename COUPLING, typename CONVERTER>
  static void computeParameters(COUPLING& coupling, CONVERTER& converter) {
    std::cout << "Starting computeParameters" << std::endl;
    // Compute K_F
    auto k_f = coupling.template getParameter<AdsorptionReaction::K_F>();
    auto r_p = coupling.template getParameter<AdsorptionReaction::R_P>();
    k_f[0] = k_f[0] * converter.getConversionFactorTime() * V(3.) / r_p[0];
    coupling.template setParameter<AdsorptionReaction::K_F>(k_f[0]);
    // Compute D_S
    auto D_s = coupling.template getParameter<AdsorptionReaction::D_S>();
    D_s[0] = D_s[0] * converter.getConversionFactorTime()
                    / ( converter.getConversionFactorLength() * converter.getConversionFactorLength() );
    coupling.template setParameter<AdsorptionReaction::D_S>(D_s[0]);
    // Compute C_0
    auto c_0 = coupling.template getParameter<AdsorptionReaction::C_0>();
    c_0[0] = c_0[0] / converter.getConversionFactorDensity();
    coupling.template setParameter<AdsorptionReaction::C_0>(c_0[0]);
    // Compute R_P
    r_p[0] = r_p[0] / converter.getConversionFactorLength();
    coupling.template setParameter<AdsorptionReaction::R_P>(r_p[0]);
    // Compute K_S
    auto k_s = V(15.) * D_s[0] / ( r_p[0] * r_p[0] );
    coupling.template setParameter<AdsorptionReaction::K_S>(k_s);
    // Compute Q_0
    auto q_0 = ISOTHERM().getLoadingFromCoupling(c_0,coupling) / converter.getConversionFactorParticleDensity();
    coupling.template setParameter<AdsorptionReaction::Q_0>(q_0[0]);
    // Enable external mass transfer
    coupling.template setParameter<AdsorptionReaction::CONV_DENS>(converter.getConversionFactorDensity());
    coupling.template setParameter<AdsorptionReaction::CONV_PARC_DENS>(converter.getConversionFactorParticleDensity());
    coupling.template setParameter<AdsorptionReaction::EXT_MASS>(k_f[0] > 0.);
  }

  /**
   * @brief Get surface loading without film diffusion.
   *
   * @param soluteConcentration
   * @return loading in g/m^3
   */
  template <typename V, typename PARAMETERS>
  V getSurfaceLoading(V soluteConcentration, PARAMETERS& parameters) any_platform {
    if (soluteConcentration < 0) return 0; // prevent nan
    V load = ISOTHERM().getLoading(soluteConcentration,parameters);
    if (load < 0) load = 0;
    return load;
  }

  /**
   * @brief  Get surface loading with film diffusion.
   *
   * Solves implicit boundary condition for surface concentration
   *
   * @param soluteConcentration
   * @param particleLoading
   * @param particleConcentration
   * @return T
   */
  template <typename V, typename PARAMETERS>
  V getSurfaceLoading(V soluteConcentration, V particleLoading,
                      V particleConcentration, PARAMETERS& parameters) any_platform {
    if (soluteConcentration < 0) return 0; // prevent nan
    // using D = double;  // because fsolve is hardcoded in double type

    // previous surface concentration as initial guess
    V surfaceConcentration = soluteConcentration;
    V args[] = {particleLoading, soluteConcentration, particleConcentration};
    for(int i = 0; i<100; i++){
      //ADf<T,1> y0_ad = iniAD(surfaceConcentration);
      //ADf<T,1> y1_ad = filmDiffBoundary(y0_ad, args);
      //surfaceConcentration -= filmDiffBoundary(surfaceConcentration, args) / y1_ad.d(0);

      V dCs = surfaceConcentration/1E8;
      surfaceConcentration -= 2.*dCs*filmDiffBoundary(surfaceConcentration, args,parameters)
                                    /(filmDiffBoundary(surfaceConcentration + dCs, args, parameters)
                                    - filmDiffBoundary(surfaceConcentration - dCs, args, parameters));
    }
    V load = ISOTHERM().getLoading(surfaceConcentration,parameters);
    if (load < 0) load = 0;
    return load;
  }

  /// Calculates reaction rate in lattice units
  /// \param soluteConcentration
  /// \param particleLoading
  /// \param particleConcentration
  /// \return change in solute concentration
  template <typename V, typename PARAMETERS>
  Vector<V, 2> getReactionRate(V soluteConcentration, V particleLoading,
                               V particleConcentration, PARAMETERS& params) any_platform{
    V surfaceLoad;
    V conversionFactorDensity = params.template get<CONV_DENS>();
    V conversionFactorParticleDensity = params.template get<CONV_PARC_DENS>();
    particleConcentration *= conversionFactorParticleDensity; // convert to kg/m^3
    soluteConcentration *= conversionFactorDensity;           // convert to g/m^3
    particleLoading *= conversionFactorParticleDensity;       // convert to g/m^3
    auto externalMassTransferEnabled = params.template get<EXT_MASS>();
    if (externalMassTransferEnabled) {
      surfaceLoad = getSurfaceLoading(soluteConcentration, particleLoading, particleConcentration, params);
    } else {
      surfaceLoad = getSurfaceLoading(soluteConcentration,params);
    }
    V D_s = params.template get<D_S>();
    V r_p = params.template get<R_P>();
    V k_s = params.template get<K_S>();

    V reactionRate = k_s * (surfaceLoad * particleConcentration - particleLoading);
    Vector<V, 2> reactionRates(reactionRate / conversionFactorDensity, // solute
                               reactionRate / conversionFactorParticleDensity ); // loading

    return reactionRates;
  }

  /// equations for surface concentration
  template <typename V, typename PARAMETERS>
  V filmDiffBoundary(V c_s_, V args[], PARAMETERS& params) any_platform{
    V q = args[0];
    V c = args[1];
    V particleConcentration = args[2];
    auto c_s = c_s_;
    V q_s = ISOTHERM().getLoading(c_s, params);
    // getting parameters
    V k_f = params.template get<K_F>();
    V D_s = params.template get<D_S>();
    V r_p = params.template get<R_P>();
    V k_s = params.template get<K_S>();
    V equation = k_f * (c - c_s) - k_s * (particleConcentration * q_s - q);
    return equation;
  }

  template <typename V, typename COUPLING, typename CONVERTER>
  static V getPhysFilmTransferConstant(COUPLING& coupling, CONVERTER& converter) {
    auto k_f = coupling.template getParameter<AdsorptionReaction::K_F>();
    return k_f[0] / converter.getConversionFactorTime();
  }

  template <typename V, typename COUPLING, typename CONVERTER>
  static V getPhysSurfaceTransferConstant(COUPLING& coupling, CONVERTER& converter) {
    auto k_s = coupling.template getParameter<AdsorptionReaction::K_S>();
    return k_s[0] / converter.getConversionFactorTime();
  }

  template <typename V, typename COUPLING, typename CONVERTER>
  static V getPhysSurfaceDiffusionConstant(COUPLING& coupling, CONVERTER& converter) {
    auto D_s = coupling.template getParameter<AdsorptionReaction::D_S>();
    return D_s[0] * converter.getConversionFactorLength() * converter.getConversionFactorLength()
                  / converter.getConversionFactorTime();
  }

  template <typename V, typename COUPLING, typename CONVERTER>
  static void print(std::ostream& clout, COUPLING& coupling, CONVERTER& converter) {
    // Get values
    auto r_p = coupling.template getParameter<AdsorptionReaction::R_P>();
    auto c_0 = coupling.template getParameter<AdsorptionReaction::C_0>();
    auto q_0 = coupling.template getParameter<AdsorptionReaction::Q_0>();
    clout << "----------------- Reaction information -----------------" << std::endl;
    clout << "-- Parameters:" << std::endl;
    clout << "Particle diameter(m):                   d_p=    " << r_p[0] * 2 * converter.getConversionFactorLength() << std::endl;
    clout << "Film diffusion constant(m/s):           k_f*=   " << getPhysFilmTransferConstant<V>(coupling,converter) << std::endl;
    clout << "Surface diffusion constant(m^2/s):      D_s=    " << getPhysSurfaceDiffusionConstant<V>(coupling,converter) << std::endl;
    clout << "Initial solute concentration(mg/mL):    c_0=    " << c_0[0] * converter.getConversionFactorDensity() << std::endl;
    clout << "Equilibrium surface loading(mg/g):      q_0=    " << q_0[0] * converter.getConversionFactorParticleDensity() << std::endl;
    clout << "lattice Equilibrium surface loading:    q_0=    " << q_0[0] << std::endl;
    clout << "Surface Mass transfer coefficient(m/s): k_s*=   " << getPhysSurfaceTransferConstant<V>(coupling,converter) << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;
    ISOTHERM::template print<V>(clout, coupling);
  }

};

}

#endif //OLB_APPS_FLORIAN_ADSORPTION3D_ADSORPTIONREACTION_H_
