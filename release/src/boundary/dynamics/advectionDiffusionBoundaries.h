/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *                2022 Nando Suntoyo, Adrian Kummerlaender, Shota Ito
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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H

#include "descriptor/descriptor.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "dynamics/dynamics.h"

namespace olb {

//===================================================================================
//================= AdvectionDiffusionDynamcis on Flat Boundaries =========
//===================================================================================
// D2Q5/D3Q7:   moment-based BC (see ref. doi:10.1504/PCFD.2016.077296)
// D2Q9/D3Q19:  non-equilibrium bounce back method (see ref. doi:10.1007/978-3-319-44649-3)
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
struct AdvectionDiffusionBoundariesDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::template exchange_momenta<MOMENTA>::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = AdvectionDiffusionBoundariesDynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA,
    direction, orientation
  >;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,DYNAMICS,M,direction,orientation>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionBoundariesDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<AdvectionDiffusionBoundariesDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V dirichletTemperature = MomentaF().computeRho(cell);

    // Placeholder - the index calculation of this section can and should be done at compile time
    constexpr auto unknownIndices = util::subIndexOutgoing<DESCRIPTOR, direction, orientation>();
    constexpr auto knownIndices = util::subIndexOutgoingRemaining<DESCRIPTOR, direction, orientation>();

    if constexpr ((DESCRIPTOR::d == 3 && DESCRIPTOR::q == 7)||(DESCRIPTOR::d == 2 && DESCRIPTOR::q == 5)) {
      static_assert(unknownIndices.size() == 1 && knownIndices.size() == DESCRIPTOR::q - 1,
                    "More than one population missing");
      V sum = V{0};
      for (unsigned iPop : knownIndices) {
        sum += cell[iPop];
      }

      // on cell there are non-shifted values -> temperature has to be changed
      V difference = dirichletTemperature - V{1} - sum;
      cell[unknownIndices[0]] = difference;
      return typename DYNAMICS::template exchange_momenta<MOMENTA>::CollisionO().apply(cell, parameters);
    }
    else {
      auto u = cell.template getField<descriptors::VELOCITY>();
      // part for q=19 copied from AdvectionDiffusionEdgesDynamics.collide()
      // but here just all missing directions, even at border of inlet area
      // has to be checked!
      for (unsigned iPop : unknownIndices) {
        cell[iPop] =
          equilibrium<DESCRIPTOR>::firstOrder(iPop, dirichletTemperature, u)
          - (cell[descriptors::opposite<DESCRIPTOR>(iPop)]
             - equilibrium<DESCRIPTOR>::firstOrder(
               descriptors::opposite<DESCRIPTOR>(iPop),
               dirichletTemperature, u));
      }
      return {-1,-1};
    }
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "AdvectionDiffusionBoundariesDynamics";
  };

};

template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int...ARGS>
using AdvectionLocalDiffusionBoundariesDynamics = dynamics::ParameterFromCell<descriptors::OMEGA,
AdvectionDiffusionBoundariesDynamics< T, DESCRIPTOR, DYNAMICS, MOMENTA, ARGS...>
>;

//===================================================================================
//================= AdvectionDiffusionDynamcis On Edges =========
//===================================================================================
// non-equilibrium bounce back method (see ref. doi:10.1007/978-3-319-44649-3)
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgesDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::template exchange_momenta<MOMENTA>::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = AdvectionDiffusionEdgesDynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA,
    plane, normal1, normal2
  >;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,DYNAMICS,M,plane,normal1,normal2>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionEdgesDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<AdvectionDiffusionEdgesDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    constexpr auto unknownIndexes = util::subIndexOutgoing3DonEdges<DESCRIPTOR,plane,normal1,normal2>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop : unknownIndexes) {
      cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, temperature, u)
                 - (  cell[descriptors::opposite<DESCRIPTOR>(iPop)]
                    - equilibrium<DESCRIPTOR>::firstOrder(descriptors::opposite<DESCRIPTOR>(iPop),
                                                                   temperature,
                                                                   u));
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::template exchange_momenta<MOMENTA>::CollisionO().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "AdvectionDiffusionEdgesDynamics";
  };
};

template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int...ARGS>
using AdvectionLocalDiffusionEdgesDynamics = dynamics::ParameterFromCell<descriptors::OMEGA,
AdvectionDiffusionEdgesDynamics< T, DESCRIPTOR, DYNAMICS, MOMENTA, ARGS...>
>;

//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 2D Boundaries =========
//===================================================================================
// non-equilibrium bounce back method (see ref. doi:10.1007/978-3-319-44649-3)
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int xNormal, int yNormal>
struct AdvectionDiffusionCornerDynamics2D final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::template exchange_momenta<MOMENTA>::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = AdvectionDiffusionCornerDynamics2D<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA,
    xNormal, yNormal
  >;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,DYNAMICS,M,xNormal,yNormal>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionCornerDynamics2D);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<AdvectionDiffusionCornerDynamics2D>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    constexpr auto unknownIndexes = util::subIndexOutgoing2DonCorners<DESCRIPTOR,xNormal,yNormal>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop : unknownIndexes) {
      cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, temperature, u)
                 - (  cell[descriptors::opposite<DESCRIPTOR>(iPop)]
                    - equilibrium<DESCRIPTOR>::firstOrder(descriptors::opposite<DESCRIPTOR>(iPop),
                                                                   temperature,
                                                                   u));
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::CollisionO().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "AdvectionDiffusionCornerDynamics2D";
  };

};

template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int... NORMAL>
using AdvectionLocalDiffusionCornerDynamics2D = dynamics::ParameterFromCell<descriptors::OMEGA,
 AdvectionDiffusionCornerDynamics2D<T, DESCRIPTOR, DYNAMICS, MOMENTA, NORMAL... >
>;

//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 3D Boundaries =========
//===================================================================================
// non-equilibrium bounce back method (see ref. doi:10.1007/978-3-319-44649-3)
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int xNormal, int yNormal, int zNormal>
struct AdvectionDiffusionCornerDynamics3D final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::template exchange_momenta<MOMENTA>::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = AdvectionDiffusionCornerDynamics3D<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA,
    xNormal, yNormal, zNormal
  >;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,DYNAMICS,M,xNormal,yNormal,zNormal>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionCornerDynamics3D);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<AdvectionDiffusionCornerDynamics3D>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    constexpr auto unknownIndexes = util::subIndexOutgoing3DonCorners<DESCRIPTOR,xNormal,yNormal,zNormal>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop : unknownIndexes) {
      cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, temperature, u)
                 - (  cell[descriptors::opposite<DESCRIPTOR>(iPop)]
                    - equilibrium<DESCRIPTOR>::firstOrder(descriptors::opposite<DESCRIPTOR>(iPop),
                                                                   temperature,
                                                                   u));
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::template exchange_momenta<MOMENTA>::CollisionO().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "AdvectionDiffusionCornerDynamics3D";
  };

};

template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int... NORMAL>
using AdvectionLocalDiffusionCornerDynamics3D = dynamics::ParameterFromCell<descriptors::OMEGA,
 AdvectionDiffusionCornerDynamics3D<T, DESCRIPTOR, DYNAMICS, MOMENTA, NORMAL... >
>;


}  // namespace olb

#endif
