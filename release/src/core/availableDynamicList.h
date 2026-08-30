/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Dennis Teutscher
 *
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
#ifndef AVAILABLE_DYNAMICS_H
#define AVAILABLE_DYNAMICS_H

#include "dynamics/collisionLES.h"
#include "superLattice.hh"
//#include "core3D.h"

namespace olb {
namespace momenta{
    template <typename T, typename DESCRIPTOR>
    using BulkDynamics =
        dynamics::Tuple<T, DESCRIPTOR, Porous<BulkTuple>, equilibria::SecondOrder,
                        collision::SmagorinskyEffectiveOmega<collision::BGK>>;

    template<typename T, typename DESCRIPTOR>
    using available_dynamics = meta::list<
    // "default" with momenta::bulkTuple
      NoDynamics<T,DESCRIPTOR>,
      NoDynamicsWithZero<T,DESCRIPTOR>,
      BGKdynamics<T,DESCRIPTOR>,
      ConstRhoBGKdynamics<T,DESCRIPTOR>,
      IncBGKdynamics<T,DESCRIPTOR>,
      RLBdynamics<T,DESCRIPTOR>,
      TRTdynamics<T,DESCRIPTOR>,
      //ChopardDynamics<T,DESCRIPTOR>,
      // "default" with different momenta
      PoissonDynamics<T,DESCRIPTOR>,
      //This Dynamic is not working in this setup. What is it even.
      //P1Dynamics<T,DESCRIPTOR>,
      EquilibriumBoundaryFirstOrder<T,DESCRIPTOR>,
      EquilibriumBoundarySecondOrder<T,DESCRIPTOR>,
      BounceBack<T,DESCRIPTOR>,
      BounceBackBulkDensity<T,DESCRIPTOR>,
      PartialBounceBack<T,DESCRIPTOR>,
      // NguyenLaddCorrection collision
      BounceBackVelocity<T,DESCRIPTOR>,
      // Porous momenta, SmagorinskyEffectiveOmega collision
      BulkDynamics<T,DESCRIPTOR>,
      // With Forcing "default"
      ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>,
      ForcedShanChenBGKdynamics<T,DESCRIPTOR>,
      // Guo<> forcing
      ForcedBGKdynamics<T,DESCRIPTOR>,
      ForcedIncBGKdynamics<T,DESCRIPTOR>,
      ForcedTRTdynamics<T,DESCRIPTOR>,
      // Guo<> forcing, OmegaFromCellTauEff<> collision
      ExternalTauForcedIncBGKdynamics<T,DESCRIPTOR>
      // MCGuo<> forcing
      //MultiComponentForcedBGKdynamics<T,DESCRIPTOR>
    >;
}
    template<typename T, typename DESCRIPTOR>
    using available_fields = meta::list<
        descriptors::EXTERNAL_VELOCITY,
        descriptors::VELOCITY2,
        descriptors::AVERAGE_VELOCITY,
        descriptors::SOURCE,
        descriptors::PRESSCORR,
        descriptors::FORCE,
        descriptors::EXTERNAL_FORCE,
        descriptors::NABLARHO,
        descriptors::TAU_EFF,
        descriptors::RHO,
        descriptors::GAMMA,
        descriptors::CUTOFF_KIN_ENERGY,
        descriptors::CUTOFF_HEAT_FLUX,
        descriptors::CHEM_POTENTIAL,
        descriptors::ADDEND,
        descriptors::V6,
        descriptors::V12,
        descriptors::OMEGA,
        descriptors::INTERFACE_WIDTH,
        descriptors::MAGIC,
        descriptors::G,
        descriptors::EPSILON,
        descriptors::BODY_FORCE,
        descriptors::K,
        descriptors::NU,
        descriptors::VELOCITY_NUMERATOR,
        descriptors::VELOCITY_DENOMINATOR,
        descriptors::ZETA,
        descriptors::LOCAL_DRAG,
        descriptors::VELOCITY_SOLID,
        descriptors::COORDINATE,
        descriptors::NEIGHBOR,
        descriptors::AV_SHEAR,
        descriptors::SHEAR_RATE_MAGNITUDE,
        descriptors::TAU_W,
        descriptors::SCALAR,
        descriptors::SMAGO_CONST,
        descriptors::EFFECTIVE_OMEGA,
        descriptors::VELO_GRAD,
        descriptors::FIL_RHO,
        descriptors::LOCAL_FIL_VEL_X,
        descriptors::LOCAL_FIL_VEL_Y,
        descriptors::LOCAL_FIL_VEL_Z,
        descriptors::LOCAL_AV_DISS,
        descriptors::LOCAL_AV_TKE,
        descriptors::LOCAL_SIGMA_ADM,
        descriptors::LOCAL_NU_EDDY,
        descriptors::FILTERED_VEL_GRAD,
        descriptors::ERROR_COVARIANCE,
        descriptors::VARIANCE,
        descriptors::TAU_SGS,
        descriptors::FILTERED_POPULATION,
        descriptors::INDICATE,
        descriptors::BIOGAS_INSTANT,
        descriptors::BIOGAS_CUMULATIVE,
        descriptors::METHANE_INSTANT,
        descriptors::METHANE_CUMULATIVE,
        descriptors::CO2_INSTANT,
        descriptors::CO2_CUMULATIVE,
        descriptors::TEMPERATURE,
        descriptors::SIGMA,
        descriptors::PSI,
        descriptors::NORMGRADPSI,
        descriptors::PSI0,
        descriptors::INTERPHASE_NORMAL,
        descriptors::MASS,
        descriptors::CELL_TYPE,
        descriptors::BOUNDARY,
        descriptors::SOURCE_OLD,
        descriptors::TOP,
        descriptors::BOTTOM,
        descriptors::OLD_PHIU,
        descriptors::CONTACT_DETECTION,
        descriptors::POROSITY,
        descriptors::POROSITY2,
        descriptors::EUL2LAGR,
        descriptors::LOCATION,
        descriptors::Y1,
        descriptors::Y2,
        descriptors::VORTICITY,
        collision::LES::SMAGORINSKY,
        opti::F,
        opti::DJDF,
        opti::DJDALPHA
    >;

  template <typename T, typename DESCRIPTOR>
  DynamicsPromise<T, DESCRIPTOR> mapStringToDynamics(std::string name) {
    OstreamManager clout("DynamicsFromXML");
    bool dynamicFound = false;
    DynamicsPromise<T, DESCRIPTOR> tmp(meta::id<NoDynamics<T, DESCRIPTOR>>{});

    momenta::available_dynamics<T, DESCRIPTOR>::for_each([&](auto dynamic) {
        if (name == DynamicsPromise<T, DESCRIPTOR>(dynamic).realize()->getName()) {
            tmp = DynamicsPromise<T, DESCRIPTOR>(dynamic);
            clout << DynamicsPromise<T, DESCRIPTOR>(dynamic).realize()->getName() << std::endl;
            dynamicFound = true;
        }
    });

    if (dynamicFound) {
        clout << "Returned dynamic" << std::endl;
        return DynamicsPromise<T, DESCRIPTOR>(tmp);
    } else {
        throw std::invalid_argument("Model " + name + " not available");
    }
  }

  template <typename T, typename DESCRIPTOR>
  void defineDynamicsAndSetParameterFromXml(std::string xmlFileName, SuperGeometry<T, 3>& sGeometry, SuperLattice<T, DESCRIPTOR>& sLattice) {
    defineDynamicsFromXml(xmlFileName, sLattice,sGeometry);
    setParameterFromXml(sLattice, xmlFileName);
  }

  template <typename T, typename DESCRIPTOR>
  void defineDynamicsFromXml(std::string xmlFileName, SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T, 3>& sGeometry) {
    OstreamManager clout("DynamicsFromXML");
    DynamicsTupleParser<T, DESCRIPTOR> parser(xmlFileName);
    std::vector<int> indicators = parser.readIndicatorFromXML();
    clout << "size " << indicators.size() << std::endl;

    auto dynamicsStrings = parser.readTupleFromXML();
    for (int i = 0; i < dynamicsStrings.size(); i++) {
      clout << indicators[i] << std::endl;
      sLattice.defineDynamics(sGeometry.getMaterialIndicator(indicators[i]),
                              mapStringToDynamics<T, DESCRIPTOR>(dynamicsStrings[i]));
    }
  }

  template <typename T, typename DESCRIPTOR>
  void defineDynamicsFromString(SuperLattice<T, DESCRIPTOR>& sLattice, std::unique_ptr<olb::SuperIndicatorF<T, 3U>> indicatorF, std::string modelName)
  {

    sLattice.defineDynamics(indicatorF,
                          mapStringToDynamics<T,DESCRIPTOR>(modelName));

  }

  template<typename T, typename DESCRIPTOR>
  void setParameterFromXml(SuperLattice<T, DESCRIPTOR>& sLattice, std::string xmlFileName) {
    OstreamManager clout("ParameterFromXML");

    DynamicsTupleParser<T, DESCRIPTOR> parser(xmlFileName);
    std::map<std::string, float> parameters = parser.readParameterFromXML();

    available_fields<T, DESCRIPTOR>::for_each([&](auto field) {
      auto it = parameters.find(fields::name<decltype(field.get())>()); // Find the name in the map
      if (it != parameters.end()) {
        clout << "Name: " << fields::name<decltype(field.get())>() << ", Value: " << it->second << std::endl;
        sLattice.template setParameter<decltype(field.get())>( it->second );
      }
    });
  }

} // namespace olb

#endif // AVAILABLE_DYNAMICS_H
