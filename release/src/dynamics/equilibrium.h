/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_EQUILIBRIUM_H
#define DYNAMICS_EQUILIBRIUM_H

#include "lbm.h"
#include "descriptor/fields.h"

namespace olb {

template<typename T> struct CellStatistic;

namespace equilibria {

struct None {
  using parameters = meta::list<>;

  static std::string getName() {
    return "None";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = V{0};
      }
      return {0, 0};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct ZerothOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "ZerothOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = descriptors::t<V,DESCRIPTOR>(iPop) * rho - descriptors::t<V,DESCRIPTOR>(iPop);
      }
      return {rho, 0};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct CahnHilliardZerothOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "CahnHilliardZerothOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      const auto mu = cell.template getField<descriptors::CHEM_POTENTIAL>();
      fEq[0] = rho - (1-descriptors::t<V,DESCRIPTOR>(0))*mu-descriptors::t<V,DESCRIPTOR>(0);
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = descriptors::t<V,DESCRIPTOR>(iPop)*mu-descriptors::t<V,DESCRIPTOR>(iPop);
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    }

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct FirstOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "FirstOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct SecondOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "SecondOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct ThirdOrder {
  using parameters = meta::list<>;

  static std::string getName() {
    return "ThirdOrder";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq) any_platform {
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::thirdOrder(iPop, rho, u, uSqr);
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      return compute(cell, rho, u, fEq);
    };
  };
};

struct Incompressible {
  using parameters = meta::list<>;

  static std::string getName() {
    return "Incompressible";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename J, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, J& j, FEQ& fEq) any_platform {
      const V pressure = rho / descriptors::invCs2<V,DESCRIPTOR>();
      const V jSqr = util::normSqr<V,DESCRIPTOR::d>(j);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::incompressible(iPop, j, jSqr, pressure);
      }
      return {rho, jSqr / (rho*rho)};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, j[DESCRIPTOR::d];
      MomentaF().computeRhoJ(cell, rho, j);
      return compute(cell, rho, j, fEq);
    };
  };
};

struct MPIncompressible {
  using parameters = meta::list<>;

  static std::string getName() {
    return "Multi-Phase-Incompressible";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename P, typename U, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, P& p, U& u, FEQ& fEq) any_platform {
      const auto rho = cell.template getField<descriptors::RHO>();
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::mpincompressible(iPop, rho, u, uSqr, p);
      }
      return {rho, uSqr};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V p, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, p, u);
      return compute(cell, p, u, fEq);
    };
  };
};

struct P1 {
  using parameters = meta::list<>;

  static std::string getName() {
    return "P1";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename J, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, J& j, FEQ& fEq) any_platform {
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::P1(iPop, rho, j);
      }
      return {rho, 0};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, j[DESCRIPTOR::d];
      MomentaF().computeRhoJ(cell, rho, j);
      return compute(cell, rho, j, fEq);
    };
  };
};

struct Chopard {
  struct SPEED_OF_SOUND : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<SPEED_OF_SOUND>;

  static std::string getName() {
    return "Chopard";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename RHO, typename U, typename FEQ, typename VS2, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, RHO& rho, U& u, FEQ& fEq, VS2 vs2 = V{1} / descriptors::invCs2<V,DESCRIPTOR>()) any_platform {
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        if (iPop==0) {
          return rho * (V{1} - vs2 * descriptors::invCs2<V,DESCRIPTOR>()*(V{1}-descriptors::t<V,DESCRIPTOR>(0))
                              - descriptors::t<V,DESCRIPTOR>(0)/V{2}*descriptors::invCs2<V,DESCRIPTOR>()*uSqr)
                  - descriptors::t<V,DESCRIPTOR>(0);
        }
        else {
          V c_u{};
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
          }
          return rho * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>()
                      * (vs2 + c_u
                        + descriptors::invCs2<V,DESCRIPTOR>() / V{2} * c_u*c_u
                        - uSqr / V{2})
            - descriptors::t<V,DESCRIPTOR>(iPop);
        }
      }
      return {rho, uSqr};
    };

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) any_platform {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      const V vs2 = parameters.template get<SPEED_OF_SOUND>();
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        if (iPop==0) {
          return rho * (V{1} - vs2 * descriptors::invCs2<V,DESCRIPTOR>()*(V{1}-descriptors::t<V,DESCRIPTOR>(0))
                              - descriptors::t<V,DESCRIPTOR>(0)/V{2}*descriptors::invCs2<V,DESCRIPTOR>()*uSqr)
                  - descriptors::t<V,DESCRIPTOR>(0);
        }
        else {
          V c_u{};
          for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
            c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
          }
          return rho * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>()
                      * (vs2 + c_u
                        + descriptors::invCs2<V,DESCRIPTOR>() / V{2} * c_u*c_u
                        - uSqr / V{2})
            - descriptors::t<V,DESCRIPTOR>(iPop);
        }
      }
      return {rho, uSqr};
    };
  };
};

}

}

#endif
