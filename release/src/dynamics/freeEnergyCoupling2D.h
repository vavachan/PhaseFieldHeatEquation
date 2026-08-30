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

#ifndef FREE_ENERGY_COUPLING_2D_H
#define FREE_ENERGY_COUPLING_2D_H

namespace olb {

struct ChemicalPotentialCoupling2D {

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

    V laplaceRho1 = -12.0 * rhoA;
    V laplaceRho2 = -12.0 * rhoB;
    V laplaceRho3 = 0.; // Only used in case of three fluids / lattices

    auto alpha = parameters.template get<ALPHA>();

    // Three fluids / lattices involved: A, B, and C
    if constexpr (CELLS::map_t::keys_t::template contains<names::C>()){
      // Cell of the third lattice:
      auto& cellC = cells.template get<names::C>();

      V rhoC = cellC.computeRho();

      densitySum -= rhoC;
      densityDifference -= rhoC;

      V kappa3 = parameters.template get<KAPPA3>();
      V term3 = kappa3 * rhoC * (rhoC - 1.) * (2. * rhoC - 1.);

      laplaceRho3 = -12.0 * rhoC;

      /*
      * The neighbours of the currently looked at cell
      * (i.e. cell 0) are ordered as following:
      *
      *    1 - 8 - 7
      *    |   |   |
      *    2 - 0 - 6
      *    |   |   |
      *    3 - 4 - 5
      */
      for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {

        auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
        auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));
        auto nextCellC = cellC.neighbor(descriptors::c<DESCRIPTOR>(iPop));

        V rhoA = nextCellA.computeRho();
        V rhoB = nextCellB.computeRho();
        V rhoC = nextCellC.computeRho();

        if(iPop % 2 == 0) {
          laplaceRho1 += (2. * rhoA);
          laplaceRho2 += (2. * rhoB);
          laplaceRho3 += (2. * rhoC);
        } else {
          laplaceRho1 += rhoA;
          laplaceRho2 += rhoB;
          laplaceRho3 += rhoC;
        }
      }

      laplaceRho1 *= 0.25;
      laplaceRho2 *= 0.25;
      laplaceRho3 *= 0.25;

      cellC.template setField<descriptors::CHEM_POTENTIAL>( - term1 - term2 + term3
            + 0.25*alpha*alpha*( (kappa2 + kappa1) * laplaceRho1
                                 -(kappa2 - kappa1) * laplaceRho2
                                 -(kappa2 + kappa1 + 4.*kappa3) * laplaceRho3 ));


    } else { // Only 2 fluids / lattices involved: A and B

      for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {

        auto nextCellA = cellA.neighbor(descriptors::c<DESCRIPTOR>(iPop));
        auto nextCellB = cellB.neighbor(descriptors::c<DESCRIPTOR>(iPop));

        V rhoA = nextCellA.computeRho();
        V rhoB = nextCellB.computeRho();

        if(iPop % 2 == 0) {
          laplaceRho1 += (2. * rhoA);
          laplaceRho2 += (2. * rhoB);
        } else {
          laplaceRho1 += rhoA;
          laplaceRho2 += rhoB;
        }
      }

      laplaceRho1 *= 0.25;
      laplaceRho2 *= 0.25;

    }

    cellA.template setField<descriptors::CHEM_POTENTIAL>(term1 + term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * laplaceRho2
                                 +(kappa2 + kappa1) * (laplaceRho3 - laplaceRho1) ));

    cellB.template setField<descriptors::CHEM_POTENTIAL>(term1 - term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * (laplaceRho1 - laplaceRho3)
                                 -(kappa2 + kappa1) * laplaceRho2 ));

  }

};




struct ForceCoupling2D {

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

    V u[2] { };
    V phi = cellA.computeRho();
    V rho = cellB.computeRho();

    /*
    * Calculatations for lattice A. The neighbours of the currently
    * looked at cell (i.e. cell 0) are ordered as following:
    *
    *    1 - 8 - 7
    *    |   |   |
    *    2 - 0 - 6
    *    |   |   |
    *    3 - 4 - 5
    *
    */
    auto cellA1 = cellA.neighbor(descriptors::c<DESCRIPTOR>(1));
    auto cellA2 = cellA.neighbor(descriptors::c<DESCRIPTOR>(2));
    auto cellA3 = cellA.neighbor(descriptors::c<DESCRIPTOR>(3));

    auto cellA7 = cellA.neighbor(descriptors::c<DESCRIPTOR>(7));
    auto cellA6 = cellA.neighbor(descriptors::c<DESCRIPTOR>(6));
    auto cellA5 = cellA.neighbor(descriptors::c<DESCRIPTOR>(5));

    // 1/12 * ( -A1 -4*A2 -A3   +A7 +4*A6 +A5)
    V gradMuPhiX = 1./12. * (
    -      cellA1.template getField<descriptors::CHEM_POTENTIAL>()
    - 4. * cellA2.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellA3.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellA7.template getField<descriptors::CHEM_POTENTIAL>()
    + 4. * cellA6.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellA5.template getField<descriptors::CHEM_POTENTIAL>());

    auto cellA8 = cellA.neighbor(descriptors::c<DESCRIPTOR>(8));
    auto cellA4 = cellA.neighbor(descriptors::c<DESCRIPTOR>(4));

    // 1/12 * ( +A1 +4*A8 +A7    -A3 -4*A4 -A5)
    V gradMuPhiY = 1./12. * (
    +      cellA1.template getField<descriptors::CHEM_POTENTIAL>()
    + 4. * cellA8.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellA7.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellA3.template getField<descriptors::CHEM_POTENTIAL>()
    - 4. * cellA4.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellA5.template getField<descriptors::CHEM_POTENTIAL>());

    /*
    * Calculatations for lattice B. The neighbours of the currently
    * looked at cell (i.e. cell 0) are ordered as following:
    *
    *    1 - 8 - 7
    *    |   |   |
    *    2 - 0 - 6
    *    |   |   |
    *    3 - 4 - 5
    */
    auto cellB1 = cellB.neighbor(descriptors::c<DESCRIPTOR>(1));
    auto cellB2 = cellB.neighbor(descriptors::c<DESCRIPTOR>(2));
    auto cellB3 = cellB.neighbor(descriptors::c<DESCRIPTOR>(3));

    auto cellB7 = cellB.neighbor(descriptors::c<DESCRIPTOR>(7));
    auto cellB6 = cellB.neighbor(descriptors::c<DESCRIPTOR>(6));
    auto cellB5 = cellB.neighbor(descriptors::c<DESCRIPTOR>(5));

    // 1/12 * ( -B1 -4*B2 -B3   +B7 +4*B6 +B5)
    V gradMuRhoX = 1./12. * (
    -      cellB1.template getField<descriptors::CHEM_POTENTIAL>()
    - 4. * cellB2.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellB3.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellB7.template getField<descriptors::CHEM_POTENTIAL>()
    + 4. * cellB6.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellB5.template getField<descriptors::CHEM_POTENTIAL>());

    auto cellB8 = cellB.neighbor(descriptors::c<DESCRIPTOR>(8));
    auto cellB4 = cellB.neighbor(descriptors::c<DESCRIPTOR>(4));

    // 1/12 * ( +B1 +4*B8 +B7    -B3 -4*B4 -B5)
    V gradMuRhoY = 1./12. * (
    +      cellB1.template getField<descriptors::CHEM_POTENTIAL>()
    + 4. * cellB8.template getField<descriptors::CHEM_POTENTIAL>()
    +      cellB7.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellB3.template getField<descriptors::CHEM_POTENTIAL>()
    - 4. * cellB4.template getField<descriptors::CHEM_POTENTIAL>()
    -      cellB5.template getField<descriptors::CHEM_POTENTIAL>());



    // Get the cell of the third lattice, if it is provided
    if constexpr (CELLS::map_t::keys_t::template contains<names::C>()){
      auto& cellC = cells.template get<names::C>();
      V psi = cellC.computeRho();

      auto cellC1 = cellC.neighbor(descriptors::c<DESCRIPTOR>(1));
      auto cellC2 = cellC.neighbor(descriptors::c<DESCRIPTOR>(2));
      auto cellC3 = cellC.neighbor(descriptors::c<DESCRIPTOR>(3));

      auto cellC7 = cellC.neighbor(descriptors::c<DESCRIPTOR>(7));
      auto cellC6 = cellC.neighbor(descriptors::c<DESCRIPTOR>(6));
      auto cellC5 = cellC.neighbor(descriptors::c<DESCRIPTOR>(5));


      // 1/12 * ( -C1 -4*C2 -C3   +C7 +4*C6 +C5)
      V gradMuPsiX = 1./12. * (
      -     cellC1.template getField<descriptors::CHEM_POTENTIAL>()
      - 4 * cellC2.template getField<descriptors::CHEM_POTENTIAL>()
      -     cellC3.template getField<descriptors::CHEM_POTENTIAL>()
      +     cellC7.template getField<descriptors::CHEM_POTENTIAL>()
      + 4 * cellC6.template getField<descriptors::CHEM_POTENTIAL>()
      +     cellC5.template getField<descriptors::CHEM_POTENTIAL>());

      auto cellC8 = cellC.neighbor(descriptors::c<DESCRIPTOR>(8));
      auto cellC4 = cellC.neighbor(descriptors::c<DESCRIPTOR>(4));

      // 1/12 * ( +C1 +4*C8 +C7    -C3 -4*C4 -C5)
      V gradMuPsiY = 1./12. * (
      +     cellC1.template getField<descriptors::CHEM_POTENTIAL>()
      + 4 * cellC8.template getField<descriptors::CHEM_POTENTIAL>()
      +     cellC7.template getField<descriptors::CHEM_POTENTIAL>()
      -     cellC3.template getField<descriptors::CHEM_POTENTIAL>()
      - 4 * cellC4.template getField<descriptors::CHEM_POTENTIAL>()
      -     cellC5.template getField<descriptors::CHEM_POTENTIAL>());

      cellB.template setField<descriptors::FORCE>({
        - rho*gradMuRhoX - phi*gradMuPhiX - psi*gradMuPsiX,
        - rho*gradMuRhoY - phi*gradMuPhiY - psi*gradMuPsiY
      });

      cellB.computeU(u);
      cellA.template setField<descriptors::FORCE>(u);
      cellC.template setField<descriptors::FORCE>(u);

    } else {

      cellB.template setField<descriptors::FORCE>({
        - rho*gradMuRhoX - phi*gradMuPhiX,
        - rho*gradMuRhoY - phi*gradMuPhiY
      });

      cellB.computeU(u);
      cellA.template setField<descriptors::FORCE>(u);
    }

  }
};




struct InletOutletCoupling2D {

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

    // If relevant, get the cell of the third lattice
    if constexpr (CELLS::map_t::keys_t::template contains<names::C>()){
      auto& cellC = cells.template get<names::C>();
      cellC.defineU(u);
    }
  }
};



struct DensityOutletCoupling2D {
  // Possible scope options: PerCell, PerBlock, PerCellWithParameters
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct RHO: public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<RHO>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform {

    using V = typename CELLS::template value_t<names::A>::value_t;

    auto outletRho = parameters.template get<RHO>();

    // Get the cell of the first lattice
    auto& cellA = cells.template get<names::A>();

    // Get the cell of the second lattice
    auto& cellB = cells.template get<names::B>();

    V rho0, phi;
    rho0 = cellA.computeRho();
    phi = cellB.computeRho();

    cellA.defineRho(outletRho);

    V temp = phi * outletRho / rho0;
    cellB.defineRho(temp);


    // Get the cell of the third lattice, in case it exists
    if constexpr (CELLS::map_t::keys_t::template contains<names::C>()){
      auto& cellC = cells.template get<names::C>();

      V psi = cellC.computeRho();
      temp = psi * outletRho / rho0;
      cellC.defineRho(temp);
    }

  }
};

}  // namespace olb

#endif
