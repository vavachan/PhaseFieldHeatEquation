/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

#ifndef PHASE_FIELD_COUPLING_H
#define PHASE_FIELD_COUPLING_H

#include "core/operator.h"

namespace olb {

namespace stage {
  struct ChemPotCalc { };
  struct InitOutlet { };
  struct PhiLimiter { };
}

// =========================================================================//
// ==================Allen-Cahn + signed distance function==================//
// =========================================================================//

struct initialPsi {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();
    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    V psi = 0;
    if (phi > (0.99)) {
      psi += -4.595120;
    }
    else if (phi < (0.01)) {
      psi += 4.595120;
    }
    else {
      psi += util::log(1/phi-1);
    }
    psi *= -w/4.;
    cell.template setField<descriptors::SCALAR>(psi/util::sqrt(psi*psi+1.));
    cell.template setField<descriptors::PSI>(psi);
    cell.template setField<descriptors::PSI0>(psi);
  }

};

struct normGradPsi {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::EPSILON>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    auto psi = cell.template getField<descriptors::PSI>();
    auto psi0 = cell.template getField<descriptors::PSI0>();
    auto epsilon = parameters.template get<descriptors::EPSILON>();
    Vector<V,DESCRIPTOR::d> gradPsi{};
    V normGradPsi = 1.;
    if (util::fabs(psi) <= epsilon) {
      V a = psi - cell.neighbor({-1,0}).template getField<descriptors::PSI>();
      V b = cell.neighbor({1,0}).template getField<descriptors::PSI>() - psi;
      V c = psi - cell.neighbor({0,-1}).template getField<descriptors::PSI>();
      V d = cell.neighbor({0,1}).template getField<descriptors::PSI>() - psi;
      V gradX, gradY;
      if (psi0 > 0) {
        if (a<0) a = 0;
        if (b>0) b = 0;
        if (c<0) c = 0;
        if (d>0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else if (psi0 < 0) {
        if (a>0) a = 0;
        if (b<0) b = 0;
        if (c>0) c = 0;
        if (d<0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else {}
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    } else {
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    }
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct normGradPsiBoundary2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::EPSILON>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto psi = cell.template getField<descriptors::PSI>();
    auto psi0 = cell.template getField<descriptors::PSI0>();
    auto epsilon = parameters.template get<descriptors::EPSILON>();
    Vector<T,DESCRIPTOR::d> gradPsi{};
    T normGradPsi = 1.;
    if (util::fabs(psi) <= epsilon) {
      T a = psi - cell.neighbor({-1,0}).template getField<descriptors::PSI>();
      T b = cell.neighbor({1,0}).template getField<descriptors::PSI>() - psi;
      T c = psi - cell.neighbor({0,-1}).template getField<descriptors::PSI>();
      T d = cell.neighbor({0,1}).template getField<descriptors::PSI>() - psi;
      if (xNormal == 1) b = a;
      else if (xNormal == -1) a = b;
      if (yNormal == 1) d = c;
      else if (yNormal == -1) c = d;
      T gradX, gradY;
      if (psi0 > 0) {
        if (a<0) a = 0;
        if (b>0) b = 0;
        if (c<0) c = 0;
        if (d>0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else if (psi0 < 0) {
        if (a>0) a = 0;
        if (b<0) b = 0;
        if (c>0) c = 0;
        if (d<0) d = 0;
        gradX = a;
        if (b*b>a*a) gradX = b;
        gradY = c;
        if (d*d>c*c) gradY = d;
        normGradPsi = util::sqrt(gradX*gradX+gradY*gradY);
      }
      else {}
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    } else {
      cell.template setField<descriptors::NORMGRADPSI>(normGradPsi);
    }
  }
};

struct psiEvolve {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct DELTAT    : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<DELTAT,descriptors::EPSILON>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    auto dt = parameters.template get<DELTAT>();
    auto epsilon = parameters.template get<descriptors::EPSILON>();
    auto s = cell.template getField<descriptors::SCALAR>();
    auto psi = cell.template getField<descriptors::PSI>();
    V psi_new = psi;
    if (util::fabs(psi) >= epsilon) cell.template setField<descriptors::PSI>(psi_new);
    else {
      auto normGradPsi = cell.template getField<descriptors::NORMGRADPSI>();
      psi_new = psi-dt*s*(normGradPsi-1.);
      cell.template setField<descriptors::PSI>(psi_new);
    }
  }
};

struct dispersionLimiter {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    V w = parameters.template get<descriptors::INTERFACE_WIDTH>();
    V psi = cell.template getField<descriptors::PSI>();
    auto phi = cell.template getField<descriptors::STATISTIC>();
    if (fabs(psi) >= 1.5*w) {
      if (phi[0] > 0.5 && phi[0] < 0.995) {
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] += descriptors::t<V,DESCRIPTOR>(iPop) * (1-phi[0]);
        }
      } else if (phi[0] < 0.5 && phi[0] > 0.005) {
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] += descriptors::t<V,DESCRIPTOR>(iPop) * (0-phi[0]);
        }
      }
      phi[0] = cell.computeRho();
      cell.template setField<descriptors::STATISTIC>(phi);
    }
  }
};

struct AllenCahnNonLocalHelper {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::EPSILON>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    auto epsilon = parameters.template get<descriptors::EPSILON>();

    V psi = cell.template getField<descriptors::PSI>();
    V interfaceIndicator = 0;
    if (util::fabs(psi) < epsilon) interfaceIndicator = 1;
    V top = (1.-interfaceIndicator)*phi*(phi-1.)*(phi-0.5);
    V bottom = (1.-interfaceIndicator)*(1.-phi)*phi;

    cell.template setField<descriptors::TOP>(top);
    cell.template setField<descriptors::BOTTOM>(bottom);
  }
};

struct AllenCahnPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct TAUS        : public descriptors::FIELD_BASE<3> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  struct NONLOCALITY : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<SIGMA,W,TAUS,RHOS,NONLOCALITY>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    auto scale = cellNS.template getField<descriptors::SCALAR>();
    auto sigma = parameters.template get<SIGMA>()*(0.99*scale+0.01);
    auto w = parameters.template get<W>();
    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phi;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = (tau_v + (tau_l-tau_v)*phi);//*(-9.*scale+10.);
    auto gamma = parameters.template get<NONLOCALITY>();
    auto tau_mobil = parameters.template get<TAUS>()[2];
    V M = (tau_mobil-0.5)/descriptors::invCs2<V,DESCRIPTOR>();

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    //cellNS.computeU(u);
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    forceNS += mu*gradPhi;
    auto externalForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::FORCE>(externalForce + forceNS/rho);
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);
    cellNS.computeU(u);

    V psi = cellAC.template getField<descriptors::PSI>();
    V interfaceIndicator = 0.;
    if (util::fabs(psi) < 3.*w) interfaceIndicator = 1;
    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n{};
    Vector<V,DESCRIPTOR::d> phiU{};
    V gradPhiSqr = 0;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      gradPhiSqr += gradPhi[iD]*gradPhi[iD];
      phiU[iD] = phi*u[iD];
    }
    if (gradPhiSqr >= 1e-28) {
      n += gradPhi/util::sqrt(gradPhiSqr);
    }
    old_phiU = cellAC.template getField<descriptors::OLD_PHIU>();
    cellAC.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    V D_N = 4*beta/k*(interfaceIndicator-1.)*(phi*(phi-1.)*(phi-0.5)-gamma*(1.-phi)*phi);
    V source_old = cellAC.template getField<descriptors::SOURCE_OLD>();
    cellAC.template setField<descriptors::SOURCE_OLD>(M*D_N);
    forceAC += (phiU - old_phiU) + interfaceIndicator*lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cellAC.template setField<descriptors::FORCE>(forceAC);
    V source = 1.5*M*D_N-0.5*source_old;
    cellAC.template setField<descriptors::SOURCE>(source);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};

struct VelocityCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cells) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    // Computation and storage of forces
    V u[DESCRIPTOR::d] {};
    cellNS.computeU(u);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};

struct AllenCahnOutletCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct NORMAL      : public descriptors::FIELD_BASE<2> { };
  using parameters = meta::list<SIGMA,W,NORMAL>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    auto w = parameters.template get<W>();
    auto sigma = parameters.template get<SIGMA>();
    Vector<int, DESCRIPTOR::d> normal = parameters.template get<NORMAL>();
    auto rho = cellNS.template getField<descriptors::RHO>();
    auto nablaRho = cellNS.template getField<descriptors::NABLARHO>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    auto cellAC_n1 = cellAC.neighbor({-normal[0],-normal[1]});
    auto cellNS_n1 = cellNS.neighbor({-normal[0],-normal[1]});
    V phi_n1 = cellAC_n1.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }
    V laplacePhi_n1 = 0;
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i_n1 = cellAC_n1.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      laplacePhi_n1 += 2*(phi_i_n1 - phi_n1) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    // Computation of thermodynamic pressure and update populations
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    V mu_n1 = 4*beta*phi_n1*(phi_n1-1.)*(phi_n1-0.5)-k*laplacePhi_n1;

    V u[DESCRIPTOR::d] {};
    //cellAC_n1.computeU(u);
    V rhoU = u[0]*nablaRho[0]+u[1]*nablaRho[1];
    V p_n1 = cellNS_n1.computeRho();
    V p_n2 = cellNS.neighbor({-2*normal[0],-2*normal[1]}).computeRho();
    V p_old = cellNS.computeRho();
    V p_out = (p_n1 + (mu+mu_n1)/2*(phi-phi_n1));//0.8*(0.25*sigma*curv*(2*phi-1)*((2*phi-1)*(2*phi-1)-3)+0.5*sigma*curv);
    //V p_star = p_out;
    if (phi >= 0.95) {
      p_out = 0;
    }
    //V p_star = p_out-3./10.*rhoU;
    cellNS.defineRho(p_out);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cellNS[iPop] += descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>() * (p_out-p_old);
    }
  }
};

struct LiangSinglePostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n{};
    Vector<V,DESCRIPTOR::d> phiU{};
    V u[DESCRIPTOR::d] {};
    cell.computeU(u);
    V gradPhiSqr = 0;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      gradPhiSqr += gradPhi[iD]*gradPhi[iD];
      phiU[iD] = phi*u[iD];
    }
    if (gradPhiSqr >= 1e-28) {
      n += gradPhi/util::sqrt(gradPhiSqr);
    }

    old_phiU = cell.template getField<descriptors::OLD_PHIU>();
    cell.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    forceAC += (phiU - old_phiU) + lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cell.template setField<descriptors::FORCE>(forceAC);
  }
};

struct LiangPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct SIGMA       : public descriptors::FIELD_BASE<1> { };
  struct W           : public descriptors::FIELD_BASE<1> { };
  struct TAUS        : public descriptors::FIELD_BASE<2> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  struct SWITCH      : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::list<SIGMA,W,TAUS,RHOS,SWITCH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAC = cells.template get<names::Component1>();

    V phi = cellAC.template getFieldComponent<descriptors::STATISTIC>(0);
    Vector<V,DESCRIPTOR::d> gradPhi{};
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellAC.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      V tcs2 = descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * tcs2;
      laplacePhi += 2*(phi_i - phi) * tcs2;
    }

    auto scale = cellNS.template getField<descriptors::SCALAR>();
    auto sigma = parameters.template get<SIGMA>()*scale;
    auto w = parameters.template get<W>();
    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phi;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = (tau_v + (tau_l-tau_v)*phi)*scale + (1-scale)*5.;

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)-k*laplacePhi;
    forceNS += mu*gradPhi;
    auto externalForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::FORCE>(externalForce + forceNS/rho);
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);
    cellNS.computeU(u);

    Vector<V,DESCRIPTOR::d> forceAC{};
    Vector<V,DESCRIPTOR::d> old_phiU{};
    Vector<V,DESCRIPTOR::d> n{};
    Vector<V,DESCRIPTOR::d> phiU{};
    auto velo_switch = parameters.template get<SWITCH>();
    V gradPhiSqr = 0;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      gradPhiSqr += gradPhi[iD]*gradPhi[iD];
      u[iD] *= velo_switch;
      phiU[iD] = phi*u[iD];
    }
    if (gradPhiSqr >= 1e-28) {
      n += gradPhi/util::sqrt(gradPhiSqr);
    }
    old_phiU = cellAC.template getField<descriptors::OLD_PHIU>();
    cellAC.template setField<descriptors::OLD_PHIU>(phiU);
    V lambda = 4*phi*(1.-phi)/w;
    forceAC += (phiU - old_phiU) + lambda*n/descriptors::invCs2<V,DESCRIPTOR>();
    cellAC.template setField<descriptors::FORCE>(forceAC);
    cellAC.template setField<descriptors::VELOCITY>(u);
  }
};

struct RhoWettingStatistics  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    auto statistic = cell.template getField<descriptors::STATISTIC>();
    statistic[0] = cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(statistic);
    cell.template setField<descriptors::PHIWETTING>(statistic[0]);
  }
};

struct WellBalancedCahnHilliardPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct TAUS        : public descriptors::FIELD_BASE<2> { };
  struct RHOS        : public descriptors::FIELD_BASE<2> { };
  using parameters = meta::list<TAUS,RHOS>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellCH = cells.template get<names::Component1>();

    V phi = cellCH.template getFieldComponent<descriptors::STATISTIC>(0);
    V mu = cellCH.template getField<descriptors::CHEM_POTENTIAL>();
    Vector<V,DESCRIPTOR::d> gradPhi{};
    Vector<V,DESCRIPTOR::d> gradMu{};
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cellCH.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0);
      const V mu_i = cellCH.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::CHEM_POTENTIAL>();
      //const V phiMu_i = phi_i*mu_i;
      //gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      //gradMu += mu_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      if(int(cellCH.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::BOUNDARY>())) gradMu += mu * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      else gradMu += mu_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      if(int(cellCH.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::BOUNDARY>())) gradPhi += phi * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      else gradPhi += phi_i * descriptors::c<DESCRIPTOR>(iPop) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    auto rho_v = parameters.template get<RHOS>()[0];
    auto rho_l = parameters.template get<RHOS>()[1];
    V rho = rho_v + (rho_l-rho_v)*phi;
    auto tau_v = parameters.template get<TAUS>()[0];
    auto tau_l = parameters.template get<TAUS>()[1];
    V tau = tau_v + (tau_l-tau_v)*phi;

    // Computation and storage of forces
    Vector<V,DESCRIPTOR::d> forceNS{};
    V u[DESCRIPTOR::d] {};
    forceNS += -phi*gradMu;
    //V mu = cellCH.template getField<descriptors::CHEM_POTENTIAL>();
    //forceNS += mu*gradPhi;
    auto externalBlockForce = cellNS.template getField<descriptors::EXTERNAL_FORCE>();
    cellNS.template setField<descriptors::FORCE>(externalBlockForce + forceNS/rho);
    cellNS.template setField<descriptors::NABLARHO>((rho_l-rho_v)*gradPhi);
    cellNS.template setField<descriptors::TAU_EFF>(tau);
    cellNS.template setField<descriptors::RHO>(rho);

    V source{};
    cellNS.computeU(u);
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      source += gradPhi[iD]*u[iD];
    }
    V source_old = cellCH.template getField<descriptors::SOURCE_OLD>();
    cellCH.template setField<descriptors::SOURCE>(1.5*source-0.5*source_old);
    cellCH.template setField<descriptors::SOURCE_OLD>(source);
    cellCH.template setField<descriptors::VELOCITY>(u);
  }
};

struct ChemPotentialPhaseFieldProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::SCALAR,descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    V phi = cell.template getFieldComponent<descriptors::STATISTIC>(0);
    V laplacePhi = 0;
    //TODO: use lattice gradient schemes here from Cahn-Hilliard implementation
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      const V phi_i = cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::PHIWETTING>();
      laplacePhi += 2*(phi_i - phi) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    }

    auto sigma = parameters.template get<descriptors::SCALAR>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    // Computation and storage of chemical potential
    V k = 1.5*sigma*w;
    V beta = 12*sigma/w;
    V mu = 4*beta*phi*(phi-1.)*(phi-0.5)+0.25*(phi<0)*phi-k*laplacePhi;
    cell.template setField<descriptors::CHEM_POTENTIAL>(mu);
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct WellBalancedWallProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::THETA,descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto theta = parameters.template get<descriptors::THETA>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    auto phi_n = cell.template getField<descriptors::STATISTIC>();
    phi_n[0] = cell.neighbor({-xNormal,-yNormal}).computeRho();
    //auto phi_w = phi_n[0] + 5.*4./w * cos(theta) *phi_n[0]*(1-phi_n[0])*phi_n[0]*(1-phi_n[0]);
    auto phi_w = phi_n[0] + 4./w * cos(theta) *phi_n[0]*(1-phi_n[0]);
    //auto phi_w = phi_n[0] + 0.55662613315165*5.*4./w * cos(theta) *((0.2+phi_n[0])*(1.2-phi_n[0])*(0.2+phi_n[0])*(1.2-phi_n[0])*(0.2+phi_n[0])*(1.2-phi_n[0])-0.0138824);
    cell.template setField<descriptors::PHIWETTING>(phi_w);

    auto cp_n = cell.template getField<descriptors::CHEM_POTENTIAL>();
    cp_n = cell.neighbor({-xNormal,-yNormal}).template getField<descriptors::CHEM_POTENTIAL>();
    cell.template setField<descriptors::STATISTIC>(phi_n);
    cell.template setField<descriptors::CHEM_POTENTIAL>(cp_n);
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
struct WellBalancedWallProcessor3D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::THETA,descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto theta = parameters.template get<descriptors::THETA>();
    auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();

    auto phi_n = cell.template getField<descriptors::STATISTIC>();
    phi_n[0] = cell.neighbor({-xNormal,-yNormal,-zNormal}).computeRho();
    //auto phi_w = phi_n[0] + 5.*4./w * cos(theta) *phi_n[0]*(1-phi_n[0])*phi_n[0]*(1-phi_n[0]);
    auto phi_w = phi_n[0] + 4./w * cos(theta) *phi_n[0]*(1-phi_n[0]);
    //auto phi_w = phi_n[0] + 0.55662613315165*5.*4./w * cos(theta) *((0.2+phi_n[0])*(1.2-phi_n[0])*(0.2+phi_n[0])*(1.2-phi_n[0])*(0.2+phi_n[0])*(1.2-phi_n[0])-0.0138824);
    cell.template setField<descriptors::PHIWETTING>(phi_w);

    auto cp_n = cell.template getField<descriptors::CHEM_POTENTIAL>();
    cp_n = cell.neighbor({-xNormal,-yNormal,-zNormal}).template getField<descriptors::CHEM_POTENTIAL>();
    cell.template setField<descriptors::STATISTIC>(phi_n);
    cell.template setField<descriptors::CHEM_POTENTIAL>(cp_n);
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct GeometricPhaseFieldWallProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::THETA>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto theta = parameters.template get<descriptors::THETA>();
    Vector<int,DESCRIPTOR::d> tangent{yNormal*(-1),xNormal*1};
    Vector<int,DESCRIPTOR::d> opp_tang{tangent[0]*(-1),tangent[1]*(-1)};

    auto phi_1 = cell.neighbor({-xNormal,-yNormal}).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_1r = cell.neighbor({-xNormal,-yNormal}).neighbor(tangent).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_1l = cell.neighbor({-xNormal,-yNormal}).neighbor(opp_tang).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_2r = cell.neighbor({-2*xNormal,-2*yNormal}).neighbor(tangent).template getFieldComponent<descriptors::STATISTIC>(0);
    auto phi_2l = cell.neighbor({-2*xNormal,-2*yNormal}).neighbor(opp_tang).template getFieldComponent<descriptors::STATISTIC>(0);

    T dphi_1 = ( phi_1r - phi_1l ) / 2.;
    T dphi_2 = ( phi_2r - phi_2l ) / 2.;
    T tau_x_dphi = 1.5*dphi_1 - 0.5*dphi_2;

    auto phi = cell.template getField<descriptors::STATISTIC>();
    phi[0] = phi_1 + tan( M_PI/2. - theta ) * abs(tau_x_dphi);
    cell.template setField<descriptors::STATISTIC>(phi);
  }
};

template<typename T, typename DESCRIPTOR>
struct GeometricPhaseFieldCurvedWallProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::THETA>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto phi = cell.template getField<descriptors::STATISTIC>();
    phi[0] = 0;
    T weightSum = 0;
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      if(cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::BOUNDARY>() == 1.) {
        phi[0] += cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getFieldComponent<descriptors::STATISTIC>(0) * descriptors::t<T,DESCRIPTOR>(iPop);
        weightSum += descriptors::t<T,DESCRIPTOR>(iPop);
      }
    }
    phi[0] = phi[0]/weightSum;
    cell.template setField<descriptors::STATISTIC>(phi);
  }
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal>
struct IsoPhaseFieldCurvedWallProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::INTERFACE_WIDTH>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto phi_s = cell.template getField<descriptors::STATISTIC>();
    auto phi_b = cell.neighbor({-xNormal,-yNormal}).template getFieldComponent<descriptors::STATISTIC>(0);
    if (phi_b >= 0.995) {
      phi_s[0] = 1.;
    }
    else if (phi_b <= 0.005) {
      phi_s[0] = 0.;
    }
    else {
      auto theta = cell.template getField<descriptors::THETA>();
      auto w = parameters.template get<descriptors::INTERFACE_WIDTH>();
      T z_z0 = w/2.*util::atanh(2.*phi_b-1.);
      T cx = util::cos(-theta)*xNormal-util::sin(-theta)*yNormal;
      T cy = util::sin(-theta)*xNormal+util::cos(-theta)*yNormal;
      Vector<T,DESCRIPTOR::d> c{cx,cy};
      Vector<T,DESCRIPTOR::d> n{T(xNormal),T(yNormal)};
      c = util::normalize(c);
      T d = c[0]*n[0]+c[1]*n[1];
      phi_s[0] = 0.5*(1+util::tanh(2.*(z_z0+d)/w));
    }
    cell.template setField<descriptors::STATISTIC>(phi_s);
  }
};

template<int xNormal, int yNormal>
struct SetOutletCells {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    using T = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    T phi, u[DESCRIPTOR::d];
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();
    cell.computeRhoU(phi, u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      outlet_cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, phi, u);
    }
    cell.template setField<descriptors::CONV_POPS>(outlet_cell);
    auto phiGhost = cell.neighbor({xNormal,yNormal}).template getField<descriptors::STATISTIC>();
    phiGhost[0] = phi;
    cell.neighbor({xNormal,yNormal}).template setField<descriptors::STATISTIC>(phiGhost);
  }
};

struct SetIncOutletCells {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    //T p, u[DESCRIPTOR::d];
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();
    //cell.computeRhoU(p, u);
    //auto rho = cell.template getField<descriptors::RHO>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      outlet_cell[iPop] = cell[iPop];
    }
    cell.template setField<descriptors::CONV_POPS>(outlet_cell);
  }
};

};

#endif
