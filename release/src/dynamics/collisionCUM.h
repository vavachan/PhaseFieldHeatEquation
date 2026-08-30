/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Louis Kronberg, Pavel Eichler, Stephan Simonis
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

#ifndef DYNAMICS_COLLISION_CUM_H
#define DYNAMICS_COLLISION_CUM_H

#include "descriptor/definition/cum.h"

namespace olb {

template <typename DESCRIPTOR>
struct cum {
  static_assert(std::is_same<typename DESCRIPTOR::category_tag,
                             descriptors::tag::CUM>::value,
                "DESCRIPTOR is tagged as CUM");

  /// Uses the Chimera Transformation to compute central moments from the population distribution.
  ///
  /// FIXME: Somebody should check whether this method in the current form translates well to descriptors
  /// other than D3Q27.
  /// FIXME: should the chimera transformation even be used over comuputing the central moments directly?
  template <typename MOMENTA, typename CELL, typename U,
            typename V = typename CELL::value_t>
  static void computeMomenta(MOMENTA& momenta, CELL& cell, U& u) any_platform
  {
    // initialize momenta with the populations
    for (int i = 0; i < DESCRIPTOR::q; i++) {
      momenta[i] = cell[i];
    }
    constexpr auto passes = DESCRIPTOR::q / DESCRIPTOR::d;
    for (int i = DESCRIPTOR::d - 1; i >= 0; i--) {
      for (int j = 0; j < passes; j++) {
        auto [a, b, c] = descriptors::cum_data::velocityIndices<
            DESCRIPTOR::d, DESCRIPTOR::q>[passes * i + j];

        V k = descriptors::constantK<V, DESCRIPTOR::d, DESCRIPTOR::q>(
            passes * i + j);
        V sum        = momenta[a] + momenta[c];
        V difference = momenta[c] - momenta[a];
        momenta[a]   = momenta[a] + momenta[b] + momenta[c];
        momenta[b]   = difference - (momenta[a] + k) * u[i];
        momenta[c] =
            sum - V(2) * difference * u[i] + u[i] * u[i] * (momenta[a] + k);
      }
    }
  }

  /// Backwards Chimera Transformation to compute population distribution from the central moments
  /// FIXME: Somebody should check whether this method in the current form translates well to descriptors
  /// other than D3Q27
  template <typename MOMENTA, typename CELL, typename U,
            typename V = typename CELL::value_t>
  static void computePopulations(MOMENTA& momenta, CELL& cell,
                                 U& u) any_platform
  {

    constexpr auto passes = DESCRIPTOR::q / DESCRIPTOR::d;
    for (int i = 0; i < DESCRIPTOR::d; i++) {
      for (int j = 0; j < passes; j++) {
        auto [a, b, c] = descriptors::cum_data::velocityIndices<
            DESCRIPTOR::d, DESCRIPTOR::q>[passes * i + j];
        V k = descriptors::constantK<V, DESCRIPTOR::d, DESCRIPTOR::q>(
            passes * i + j);

        V ma = ((momenta[c] - momenta[b]) * V(0.5) + momenta[b] * u[i] +
                (momenta[a] + k) * (u[i] * u[i] - u[i]) * V(0.5));
        V mb = (momenta[a] - momenta[c]) - V(2) * momenta[b] * u[i] -
               (momenta[a] + k) * u[i] * u[i];
        V mc = ((momenta[c] + momenta[b]) * V(0.5) + momenta[b] * u[i] +
                (momenta[a] + k) * (u[i] * u[i] + u[i]) * V(0.5));

        momenta[a] = ma;
        momenta[b] = mb;
        momenta[c] = mc;
      }
    }

    // initialize momenta with the populations
    for (int i = 0; i < DESCRIPTOR::q; i++) {
      cell[i] = momenta[i];
    }
  }

  template <typename CELL, typename U, typename V = typename CELL::value_t>
  static V cumCollision(CELL& cell, const V& omega, V& rho, U& u) any_platform
  {
    V drho        = rho - 1;
    V inverse_rho = 1. / rho;
    V uSqr        = util::normSqr<V, DESCRIPTOR::d>(u);

    // compute central moments from the populations and save them in `moments`
    V moments[DESCRIPTOR::q];
    computeMomenta(moments, cell, u);
    auto [mbbb, mabb, mbab, mbba, maab, macb, maba, mabc, mbaa, mbac, maaa,
          maac, maca, macc, mcbb, mbcb, mbbc, mccb, mcab, mcbc, mcba, mbcc,
          mbca, mccc, mcca, mcac, mcaa] = moments;

    const V omega2 =
        1; // If modified, constants A and B must be modified too!!!
    const V omega3  = 1;
    const V omega4  = 1;
    const V omega5  = 1;
    const V omega6  = 1;
    const V omega7  = 1;
    const V omega10 = 1; // Collision of the sixth order cumulant

    // COMPUTE cumulant moments
    V CUMcbb = mcbb - ((mcaa + 1. / 3) * mabb + 2 * mbba * mbab) * inverse_rho;
    V CUMbcb = mbcb - ((maca + 1. / 3) * mbab + 2 * mbba * mabb) * inverse_rho;
    V CUMbbc = mbbc - ((maac + 1. / 3) * mbba + 2 * mbab * mabb) * inverse_rho;

    V CUMcca =
        mcca - (((mcaa * maca + 2 * mbba * mbba) + 1. / 3 * (mcaa + maca)) *
                    inverse_rho -
                1. / 9 * (drho * inverse_rho));
    V CUMcac =
        mcac - (((mcaa * maac + 2 * mbab * mbab) + 1. / 3 * (mcaa + maac)) *
                    inverse_rho -
                1. / 9 * (drho * inverse_rho));
    V CUMacc =
        macc - (((maac * maca + 2 * mabb * mabb) + 1. / 3 * (maac + maca)) *
                    inverse_rho -
                1. / 9 * (drho * inverse_rho));

    // fifth order cumulant moments
    V CUMbcc = mbcc - ((maac * mbca + maca * mbac + 4 * mabb * mbbb +
                        2 * (mbab * macb + mbba * mabc)) +
                       1. / 3 * (mbca + mbac)) *
                          inverse_rho;
    V CUMcbc = mcbc - ((maac * mcba + mcaa * mabc + 4 * mbab * mbbb +
                        2 * (mabb * mcab + mbba * mbac)) +
                       1. / 3 * (mcba + mabc)) *
                          inverse_rho;
    V CUMccb = mccb - ((mcaa * macb + maca * mcab + 4 * mbba * mbbb +
                        2 * (mbab * mbca + mabb * mcba)) +
                       1. / 3 * (macb + mcab)) *
                          inverse_rho;

    // sixth order cumulant moments
    V CUMccc =
        mccc +
        ((-4 * mbbb * mbbb - (mcaa * macc + maca * mcac + maac * mcca) -
          4 * (mabb * mcbb + mbab * mbcb + mbba * mbbc) -
          2 * (mbca * mbac + mcba * mabc + mcab * macb)) *
             inverse_rho +
         (4 * (mbab * mbab * maca + mabb * mabb * mcaa + mbba * mbba * maac) +
          2 * (mcaa * maca * maac) + 16 * mbba * mbab * mabb) *
             inverse_rho * inverse_rho -
         1. / 3 * (macc + mcac + mcca) * inverse_rho -
         1. / 9 * (mcaa + maca + maac) * inverse_rho +
         (2 * (mbab * mbab + mabb * mabb + mbba * mbba) +
          (maac * maca + maac * mcaa + maca * mcaa) +
          1. / 3 * (maac + maca + mcaa)) *
             inverse_rho * inverse_rho * 2. / 3 +
         1. / 27 * ((drho * drho - drho) * inverse_rho * inverse_rho));

    // Moment renaming
    V mxxPyyPzz = mcaa + maca + maac;
    V mxxMyy    = mcaa - maca;
    V mxxMzz    = mcaa - maac;

    V mxxyPyzz = mcba + mabc;
    V mxxyMyzz = mcba - mabc;

    V mxxzPyyz = mcab + macb;
    V mxxzMyyz = mcab - macb;

    V mxyyPxzz = mbca + mbac;
    V mxyyMxzz = mbca - mbac;

    // Relaxation of second order cumulants with no correction terms
    mxxPyyPzz += omega2 * (maaa - mxxPyyPzz);
    mxxMyy += -(-omega) * (-mxxMyy);
    mxxMzz += -(-omega) * (-mxxMzz);

    mabb += omega * (-mabb);
    mbab += omega * (-mbab);
    mbba += omega * (-mbba);

    // Relaxation of third order cumulants
    //  no limiter
    mbbb += omega4 * (-mbbb);
    mxxyPyzz += omega3 * (-mxxyPyzz);
    mxxyMyzz += omega4 * (-mxxyMyzz);
    mxxzPyyz += omega3 * (-mxxzPyyz);
    mxxzMyyz += omega4 * (-mxxzMyyz);
    mxyyPxzz += omega3 * (-mxyyPxzz);
    mxyyMxzz += omega4 * (-mxyyMxzz);

    // Compute inverse linear combinations of second and third order cumulants
    mcaa = 1. / 3 * (mxxMyy + mxxMzz + mxxPyyPzz);
    maca = 1. / 3 * (-2 * mxxMyy + mxxMzz + mxxPyyPzz);
    maac = 1. / 3 * (mxxMyy - 2 * mxxMzz + mxxPyyPzz);

    mcba = (mxxyMyzz + mxxyPyzz) * 0.5;
    mabc = (-mxxyMyzz + mxxyPyzz) * 0.5;
    mcab = (mxxzMyyz + mxxzPyyz) * 0.5;
    macb = (-mxxzMyyz + mxxzPyyz) * 0.5;
    mbca = (mxyyMxzz + mxyyPxzz) * 0.5;
    mbac = (-mxyyMxzz + mxyyPxzz) * 0.5;

    CUMacc = (1 - omega6) * (CUMacc);
    CUMcac = (1 - omega6) * (CUMcac);
    CUMcca = (1 - omega6) * (CUMcca);
    CUMbbc = (1 - omega6) * (CUMbbc);
    CUMbcb = (1 - omega6) * (CUMbcb);
    CUMcbb = (1 - omega6) * (CUMcbb);

    CUMbcc += omega7 * (-CUMbcc);
    CUMcbc += omega7 * (-CUMcbc);
    CUMccb += omega7 * (-CUMccb);

    CUMccc += omega10 * (-CUMccc);

    // Compute central moments from post collision cumulants
    mcbb = CUMcbb +
           1. / 3 * ((3 * mcaa + 1) * mabb + 6 * mbba * mbab) * inverse_rho;
    mbcb = CUMbcb +
           1. / 3 * ((3 * maca + 1) * mbab + 6 * mbba * mabb) * inverse_rho;
    mbbc = CUMbbc +
           1. / 3 * ((3 * maac + 1) * mbba + 6 * mbab * mabb) * inverse_rho;

    mcca = CUMcca + (((mcaa * maca + 2 * mbba * mbba) * 9 + 3 * (mcaa + maca)) *
                         inverse_rho -
                     (drho * inverse_rho)) *
                        1. / 9;
    mcac = CUMcac + (((mcaa * maac + 2 * mbab * mbab) * 9 + 3 * (mcaa + maac)) *
                         inverse_rho -
                     (drho * inverse_rho)) *
                        1. / 9;
    macc = CUMacc + (((maac * maca + 2 * mabb * mabb) * 9 + 3 * (maac + maca)) *
                         inverse_rho -
                     (drho * inverse_rho)) *
                        1. / 9;

    mbcc = CUMbcc + 1. / 3 *
                        (3 * (maac * mbca + maca * mbac + 4 * mabb * mbbb +
                              2 * (mbab * macb + mbba * mabc)) +
                         (mbca + mbac)) *
                        inverse_rho;
    mcbc = CUMcbc + 1. / 3 *
                        (3 * (maac * mcba + mcaa * mabc + 4 * mbab * mbbb +
                              2 * (mabb * mcab + mbba * mbac)) +
                         (mcba + mabc)) *
                        inverse_rho;
    mccb = CUMccb + 1. / 3 *
                        (3 * (mcaa * macb + maca * mcab + 4 * mbba * mbbb +
                              2 * (mbab * mbca + mabb * mcba)) +
                         (macb + mcab)) *
                        inverse_rho;

    mccc =
        CUMccc -
        ((-4 * mbbb * mbbb - (mcaa * macc + maca * mcac + maac * mcca) -
          4 * (mabb * mcbb + mbab * mbcb + mbba * mbbc) -
          2 * (mbca * mbac + mcba * mabc + mcab * macb)) *
             inverse_rho +
         (4 * (mbab * mbab * maca + mabb * mabb * mcaa + mbba * mbba * maac) +
          2 * (mcaa * maca * maac) + 16 * mbba * mbab * mabb) *
             inverse_rho * inverse_rho -
         1. / 9 * (macc + mcac + mcca) * inverse_rho -
         1. / 9 * (mcaa + maca + maac) * inverse_rho +
         (2 * (mbab * mbab + mabb * mabb + mbba * mbba) +
          (maac * maca + maac * mcaa + maca * mcaa) +
          1. / 3 * (maac + maca + mcaa)) *
             inverse_rho * inverse_rho * 2. / 3 +
         1. / 27 * ((drho * drho - drho) * inverse_rho * inverse_rho));

    // Add acceleration (body force) to first order cumulants
    mbaa = -mbaa;
    maba = -maba;
    maab = -maab;

    // Chimera transform from central moments to well conditioned distributions
    computePopulations(moments, cell, u);

    return uSqr;
  }
};


namespace collision {
  struct CUM {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "CUM";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {

    static_assert(DESCRIPTOR::d==3 && DESCRIPTOR::q == 27, "Cumulant Dynamics only implemented in D3Q27");
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      const V omega = parameters.template get<descriptors::OMEGA>();
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V uSqr = cum<DESCRIPTOR>::cumCollision(cell, omega, rho, u);
      return {rho, uSqr};
    };
  };
};
}
}
#endif
