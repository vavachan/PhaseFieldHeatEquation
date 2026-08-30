/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Liam Sauterleute
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

#ifndef K_EPSILON_RANS_DYNAMICS_H
#define K_EPSILON_RANS_DYNAMICS_H

#include "core/operator.h"

namespace olb {


/// Coupling between Navier-Stokes and k-epsilon lattices
template<typename T>
struct RANSKE {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct VISC : public descriptors::FIELD_BASE<1> { };
  struct RNG : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<VISC, RNG>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
   {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_KE = typename CELLS::template value_t<names::TurbKineticEnergy>::descriptor_t;

    /// Velocity coupling
    auto u = cells.template get<names::TurbKineticEnergy>().template getField<descriptors::VELOCITY>();
    T rho, pi[util::TensorVal<DESCRIPTOR>::n] { };
    cells.template get<names::NavierStokes>().computeAllMomenta(rho, u.data(), pi);
    cells.template get<names::TurbKineticEnergy>().template setField<descriptors::VELOCITY>(u);
    cells.template get<names::DissipationRate>().template setField<descriptors::VELOCITY>(u);

    T C1 = 1.44;
    T C2 = 1.92;
    T Cmu = 0.09;
    T invSigmaK = 1.0;
    T invSigmaE = 0.76923;
    T eta0 = 4.38;
    T beta = 0.012;
    T kinVisc = parameters.template get<VISC>();
    T inv2dx = 1./(2.);
    T k = cells.template get<names::TurbKineticEnergy>().computeRho();
    T epsilon = cells.template get<names::DissipationRate>().computeRho();
    T rng = parameters.template get<RNG>();

    if(epsilon != 0 && k != 0) {
      T turbVisc = Cmu*k*k/epsilon;
      T diffK = turbVisc*invSigmaK;
      T diffE = turbVisc*invSigmaE;
      if(rng != 0){
        diffK += kinVisc;
        diffE += kinVisc;
      }
      T tau_turb_K = diffK * descriptors::invCs2<T,DESCRIPTOR_KE>() + 0.5;
      T tau_turb_E = diffE * descriptors::invCs2<T,DESCRIPTOR_KE>() + 0.5;

      /// Finite Difference for Derivatives
      auto& cellK = cells.template get<names::TurbKineticEnergy>();
      auto& cellE = cells.template get<names::DissipationRate>();

      auto u_pXp = cellK.neighbor({1,0,0}).template getField<descriptors::VELOCITY>();
      auto u_pXm = cellK.neighbor({-1,0,0}).template getField<descriptors::VELOCITY>();
      auto u_pYp = cellK.neighbor({0,1,0}).template getField<descriptors::VELOCITY>();
      auto u_pYm = cellK.neighbor({0,-1,0}).template getField<descriptors::VELOCITY>();
      auto u_pZp = cellK.neighbor({0,0,1}).template getField<descriptors::VELOCITY>();
      auto u_pZm = cellK.neighbor({0,0,-1}).template getField<descriptors::VELOCITY>();

      T kXp = cellK.neighbor({1,0,0}).template getField<descriptors::K>();
      T kXm = cellK.neighbor({-1,0,0}).template getField<descriptors::K>();
      T kYp = cellK.neighbor({0,1,0}).template getField<descriptors::K>();
      T kYm = cellK.neighbor({0,-1,0}).template getField<descriptors::K>();
      T kZp = cellK.neighbor({0,0,1}).template getField<descriptors::K>();
      T kZm = cellK.neighbor({0,0,-1}).template getField<descriptors::K>();

      T eXp = cellE.neighbor({1,0,0}).template getField<descriptors::EPSILON>();
      T eXm = cellE.neighbor({-1,0,0}).template getField<descriptors::EPSILON>();
      T eYp = cellE.neighbor({0,1,0}).template getField<descriptors::EPSILON>();
      T eYm = cellE.neighbor({0,-1,0}).template getField<descriptors::EPSILON>();
      T eZp = cellE.neighbor({0,0,1}).template getField<descriptors::EPSILON>();
      T eZm = cellE.neighbor({0,0,-1}).template getField<descriptors::EPSILON>();

      T dxK = 0.5*(kXp - kXm);
      T dyK = 0.5*(kYp - kYm);
      T dzK = 0.5*(kZp - kZm);

      T dxE = 0.5*(eXp - eXm);
      T dyE = 0.5*(eYp - eYm);
      T dzE = 0.5*(eZp - eZm);

      T dxTurbVisc = 0.5*Cmu*(kXp*kXp/eXp - kXm*kXm/eXm);
      T dyTurbVisc = 0.5*Cmu*(kYp*kYp/eYp - kYm*kYm/eYm);
      T dzTurbVisc = 0.5*Cmu*(kZp*kZp/eZp - kZm*kZm/eZm);

      // Velocity Derivatives
      T UxDx = (u_pXp[0] - u_pXm[0])*inv2dx;
      T UxDy = (u_pYp[0] - u_pYm[0])*inv2dx;
      T UxDz = (u_pZp[0] - u_pZm[0])*inv2dx;
      T UyDx = (u_pXp[1] - u_pXm[1])*inv2dx;
      T UyDy = (u_pYp[1] - u_pYm[1])*inv2dx;
      T UyDz = (u_pZp[1] - u_pZm[1])*inv2dx;
      T UzDx = (u_pXp[2] - u_pXm[2])*inv2dx;
      T UzDy = (u_pYp[2] - u_pYm[2])*inv2dx;
      T UzDz = (u_pZp[2] - u_pZm[2])*inv2dx;

      T turbKinEnergyProduction = turbVisc*( (UxDx+UxDx)*UxDx + (UxDy+UyDx)*UxDy + (UxDz+UzDx)*UxDz
                                           + (UyDx+UxDy)*UyDx + (UyDy+UyDy)*UyDy + (UyDz+UzDy)*UyDz
                                           + (UzDx+UxDz)*UzDx + (UzDy+UyDz)*UzDy + (UzDz+UzDz)*UzDz
                                           - 2./3. * (UxDx + UyDy + UzDz) * (UxDx + UyDy + UzDz) );

      /// Modification of C2 for RNG model
      if(rng != 0){
        std::vector<T> meanStrain(6, T());
        meanStrain[0] = UxDx;
        meanStrain[1] = 0.5*( UxDy + UyDx );
        meanStrain[2] = 0.5*( UxDz + UzDx );
        meanStrain[3] = UyDy;
        meanStrain[4] = 0.5*( UyDz + UzDy );
        meanStrain[5] = UzDz;

        T meanStrainNormSqr = meanStrain[0]*meanStrain[0] + 2.0*meanStrain[1]*meanStrain[1] + 2.0*meanStrain[2]*meanStrain[2]
                          + meanStrain[3]*meanStrain[3] + 2.0*meanStrain[4]*meanStrain[4] + meanStrain[5]*meanStrain[5];
        T weightedMeanStrainNorm = util::sqrt(2.0*meanStrainNormSqr);
        T eta = weightedMeanStrainNorm * k / epsilon;
        T C2add = Cmu*eta*eta*eta*(1. - eta/eta0)/(1. + beta*eta*eta*eta);
        C2 += C2add;
      }

      /// Source Terms
      T sourceK = turbKinEnergyProduction - epsilon - 2./3.*k*(UxDx + UyDy + UzDz) + (dxK*dxTurbVisc*invSigmaK + dyK*dyTurbVisc*invSigmaK + dzK*dzTurbVisc*invSigmaK);
      T sourceE = C1*epsilon/k*turbKinEnergyProduction - C2*epsilon*epsilon/k + (dxE*dxTurbVisc*invSigmaE + dyE*dyTurbVisc*invSigmaE + dzE*dzTurbVisc*invSigmaE);
      cells.template get<names::TurbKineticEnergy>().template setField<descriptors::SOURCE>(sourceK);
      cells.template get<names::DissipationRate>().template setField<descriptors::SOURCE>(sourceE);

      T tau_turb_RANS = (turbVisc + kinVisc) * descriptors::invCs2<T,DESCRIPTOR>() + 0.5;

      /// Modification of the relaxation times
      cells.template get<names::NavierStokes>().template setField<descriptors::OMEGA>(T{1} / (tau_turb_RANS));
      cells.template get<names::TurbKineticEnergy>().template setField<descriptors::OMEGA>(T{1} / (tau_turb_K));
      cells.template get<names::DissipationRate>().template setField<descriptors::OMEGA>(T{1} / (tau_turb_E));
    }
  }

};

}

#endif