/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Czelusniak,
 *                     Adrian Kummerlaender
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

#ifndef EQUATIONS_OF_STATE_H
#define EQUATIONS_OF_STATE_H


namespace olb {

namespace EOS {

// Landau EOS
struct Landau {
  struct RHOV           : public descriptors::FIELD_BASE<1> { };
  struct RHOL           : public descriptors::FIELD_BASE<1> { };
  struct THICKNESS      : public descriptors::FIELD_BASE<1> { };
  struct SURFTENSION    : public descriptors::FIELD_BASE<1> { };

  struct RHOC           : public descriptors::FIELD_BASE<1> { };
  struct PC             : public descriptors::FIELD_BASE<1> { };
  struct BETA           : public descriptors::FIELD_BASE<1> { };
  struct KAPPA          : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<RHOV,RHOL,THICKNESS,SURFTENSION,RHOC,PC,BETA,KAPPA>;

  template <typename V, typename PARAMETERS>
  static void computeParameters(PARAMETERS& params) any_platform {
    // Compute and set rhoc (critical density)
    auto rhov = params.template getParameter<EOS::Landau::RHOV>()[0];
    auto rhol = params.template getParameter<EOS::Landau::RHOL>()[0];
    const V rhoc = 0.5 * ( rhov + rhol );
    params.template setParameter<EOS::Landau::RHOC>(rhoc);
    // Compute and set beta (density ratio parameter)
    const V beta = ( V(0.5)*(rhol-rhov)/rhoc )*( V(0.5)*(rhol-rhov)/rhoc );
    params.template setParameter<EOS::Landau::BETA>(beta);
    // Compute p_c
    auto thickness = params.template getParameter<EOS::Landau::THICKNESS>()[0];
    auto surfTension = params.template getParameter<EOS::Landau::SURFTENSION>()[0];
    const V pc = surfTension/thickness*V(3./8.)/sqrt(2.)/beta/beta;
    params.template setParameter<EOS::Landau::PC>(pc);
    // Compute kappap
    const V kappa = surfTension*thickness*V(3./4.)*sqrt(2.)/beta/rhoc/rhoc;
    params.template setParameter<EOS::Landau::KAPPA>(kappa);
  }

  template <typename V, typename PARAMETERS>
  V computeMU( V rho, PARAMETERS& params ) any_platform {
    auto rho_c = params.template get<RHOC>(); // critical density
    auto p_c = params.template get<PC>(); // critical pressure
    auto beta = params.template get<BETA>(); // additional parameter to control density ratio
    V nu_rho = (rho-rho_c)/rho_c;
    V mu = V(4.)*p_c/rho_c*nu_rho*( nu_rho*nu_rho - beta );
    return mu;
  };
};

} // end namespace EOS

} // end namespace olb

#endif

