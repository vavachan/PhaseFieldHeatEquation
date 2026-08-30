/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Robin Trunk, Sam Avis
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

#ifndef FREE_ENERGY_COUPLING_3D_H
#define FREE_ENERGY_COUPLING_3D_H

namespace olb {

struct ChemicalPotentialCoupling3D {

  // Possible scope options: PerCell, PerBlock, PerCellWithParameters
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct ALPHA: public descriptors::FIELD_BASE<1> { };
  struct KAPPA1: public descriptors::FIELD_BASE<1> { };
  struct KAPPA2: public descriptors::FIELD_BASE<1> { };
  struct KAPPA3: public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<ALPHA, KAPPA1, KAPPA2, KAPPA3>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform {

    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    // Get the cell of the first lattice
    auto& cellA = cells.template get<names::A>();

    // Get the cell of the second lattice
    auto& cellB = cells.template get<names::B>();

    V rhoA = cellA.computeRho();
    V rhoB = cellB.computeRho();

    V densitySum = rhoA + rhoB;
    V densityDifference = rhoA - rhoB;

    V kappa1 = parameters.template get<KAPPA1>();
    V kappa2 = parameters.template get<KAPPA2>();

    V term1 = 0.125 * kappa1 * (densitySum)
                  * (densitySum-1.) * (densitySum-2.);

    V term2 = 0.125 * kappa2 * (densityDifference)
                  * (densityDifference-1.) * (densityDifference-2.);

    V laplaceRho1 = -24.0 * rhoA;
    V laplaceRho2 = -24.0 * rhoB;
    V laplaceRho3 = 0;

    for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {

      auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
      auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));


      V nextRhoA = nextCellA.computeRho();
      V nextRhoB = nextCellB.computeRho();

      if(iPop == 1 || iPop == 2 || iPop == 3 || iPop == 10 || iPop == 11 || iPop == 12) {
        laplaceRho1 += (2 * nextRhoA);
        laplaceRho2 += (2 * nextRhoB);
      } else {
        laplaceRho1 += nextRhoA;
        laplaceRho2 += nextRhoB;
      }
    }

    laplaceRho1 *= 1.0 / 6.0;
    laplaceRho2 *= 1.0 / 6.0;

    auto alpha = parameters.template get<ALPHA>();

    cellA.template setField<descriptors::CHEM_POTENTIAL>(term1 + term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * laplaceRho2
                                 +(kappa2 + kappa1) * (laplaceRho3 - laplaceRho1) ));

    cellB.template setField<descriptors::CHEM_POTENTIAL>(term1 - term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * (laplaceRho1 - laplaceRho3)
                                 -(kappa2 + kappa1) * laplaceRho2 ));

  }

};




struct ForceCoupling3D {

  // Possible scope options: PerCell, PerBlock, PerCellWithParameters
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  template <typename CELLS>
  void apply(CELLS& cells) any_platform {

    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    // Get the cell of the first lattice
    auto& cellA = cells.template get<names::A>();

    // Get the cell of the second lattice
    auto& cellB = cells.template get<names::B>();

    V phi = cellA.computeRho();
    V rho = cellB.computeRho();

    // Calculatations for lattice A.
    auto cellA1 = cellA.neighbor(descriptors::c<DESCRIPTOR>(1));
    auto cellA2 = cellA.neighbor(descriptors::c<DESCRIPTOR>(2));
    auto cellA3 = cellA.neighbor(descriptors::c<DESCRIPTOR>(3));
    auto cellA4 = cellA.neighbor(descriptors::c<DESCRIPTOR>(4));
    auto cellA5 = cellA.neighbor(descriptors::c<DESCRIPTOR>(5));
    auto cellA6 = cellA.neighbor(descriptors::c<DESCRIPTOR>(6));
    auto cellA7 = cellA.neighbor(descriptors::c<DESCRIPTOR>(7));
    auto cellA8 = cellA.neighbor(descriptors::c<DESCRIPTOR>(8));
    auto cellA9 = cellA.neighbor(descriptors::c<DESCRIPTOR>(9));
    auto cellA10 = cellA.neighbor(descriptors::c<DESCRIPTOR>(10));
    auto cellA11 = cellA.neighbor(descriptors::c<DESCRIPTOR>(11));
    auto cellA12 = cellA.neighbor(descriptors::c<DESCRIPTOR>(12));
    auto cellA13 = cellA.neighbor(descriptors::c<DESCRIPTOR>(13));
    auto cellA14 = cellA.neighbor(descriptors::c<DESCRIPTOR>(14));
    auto cellA15 = cellA.neighbor(descriptors::c<DESCRIPTOR>(15));
    auto cellA16 = cellA.neighbor(descriptors::c<DESCRIPTOR>(16));
    auto cellA17 = cellA.neighbor(descriptors::c<DESCRIPTOR>(17));
    auto cellA18 = cellA.neighbor(descriptors::c<DESCRIPTOR>(18));

    V gradMuPhiX = 1./12. * (
    -     cellA6.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA7.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA4.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellA1.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA5.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA14.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellA10.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA13.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA16.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA15.template getField<descriptors::CHEM_POTENTIAL>());

    V gradMuPhiY = 1./12. * (
    -     cellA4.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA14.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA8.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellA2.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA9.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA18.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellA11.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA17.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA5.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA13.template getField<descriptors::CHEM_POTENTIAL>());

    V gradMuPhiZ = 1./12. * (
    -     cellA8.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA18.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA6.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellA3.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellA16.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA7.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellA12.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA15.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA9.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellA17.template getField<descriptors::CHEM_POTENTIAL>());


    // Calculatations for lattice B.
    auto cellB1 = cellB.neighbor(descriptors::c<DESCRIPTOR>(1));
    auto cellB2 = cellB.neighbor(descriptors::c<DESCRIPTOR>(2));
    auto cellB3 = cellB.neighbor(descriptors::c<DESCRIPTOR>(3));
    auto cellB4 = cellB.neighbor(descriptors::c<DESCRIPTOR>(4));
    auto cellB5 = cellB.neighbor(descriptors::c<DESCRIPTOR>(5));
    auto cellB6 = cellB.neighbor(descriptors::c<DESCRIPTOR>(6));
    auto cellB7 = cellB.neighbor(descriptors::c<DESCRIPTOR>(7));
    auto cellB8 = cellB.neighbor(descriptors::c<DESCRIPTOR>(8));
    auto cellB9 = cellB.neighbor(descriptors::c<DESCRIPTOR>(9));
    auto cellB10 = cellB.neighbor(descriptors::c<DESCRIPTOR>(10));
    auto cellB11 = cellB.neighbor(descriptors::c<DESCRIPTOR>(11));
    auto cellB12 = cellB.neighbor(descriptors::c<DESCRIPTOR>(12));
    auto cellB13 = cellB.neighbor(descriptors::c<DESCRIPTOR>(13));
    auto cellB14 = cellB.neighbor(descriptors::c<DESCRIPTOR>(14));
    auto cellB15 = cellB.neighbor(descriptors::c<DESCRIPTOR>(15));
    auto cellB16 = cellB.neighbor(descriptors::c<DESCRIPTOR>(16));
    auto cellB17 = cellB.neighbor(descriptors::c<DESCRIPTOR>(17));
    auto cellB18 = cellB.neighbor(descriptors::c<DESCRIPTOR>(18));

    V gradMuRhoX = 1./12. * (
    -     cellB6.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB7.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB4.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellB1.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB5.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB14.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellB10.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB13.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB16.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB15.template getField<descriptors::CHEM_POTENTIAL>());

    V gradMuRhoY = 1./12. * (
    -     cellB4.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB14.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB8.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellB2.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB9.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB18.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellB11.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB17.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB5.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB13.template getField<descriptors::CHEM_POTENTIAL>());

    V gradMuRhoZ = 1./12. * (
    -     cellB8.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB18.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB6.template getField<descriptors::CHEM_POTENTIAL>()
    - 2 * cellB3.template getField<descriptors::CHEM_POTENTIAL>()
    -     cellB16.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB7.template getField<descriptors::CHEM_POTENTIAL>()
    + 2 * cellB12.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB15.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB9.template getField<descriptors::CHEM_POTENTIAL>()
    +     cellB17.template getField<descriptors::CHEM_POTENTIAL>());

    V psi = 0.;
    V gradMuPsiX = 0.;
    V gradMuPsiY = 0.;
    V gradMuPsiZ = 0.;

    V forceX = - rho*gradMuRhoX - phi*gradMuPhiX - psi*gradMuPsiX;
    V forceY = - rho*gradMuRhoY - phi*gradMuPhiY - psi*gradMuPsiY;
    V forceZ = - rho*gradMuRhoZ - phi*gradMuPhiZ - psi*gradMuPsiZ;

    cellB.template setField<descriptors::FORCE>({forceX, forceY, forceZ});

    V u[3];
    cellB.computeU(u);

    cellA.template setField<descriptors::FORCE>(u);

  }
};


struct InletOutletCoupling3D {

  // Possible scope options: PerCell, PerBlock, PerCellWithParameters
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  template <typename CELLS>
  void apply(CELLS& cells) any_platform {

    using V = typename CELLS::template value_t<names::A>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::A>::descriptor_t;

    // Get the cell of the first lattice
    auto& cellA = cells.template get<names::A>();

    // Get the cell of the second lattice
    auto& cellB = cells.template get<names::B>();

    V u[DESCRIPTOR::d];
    cellB.computeU(u);
    cellA.defineU(u);
  }
};

}  // namespace olb

#endif
