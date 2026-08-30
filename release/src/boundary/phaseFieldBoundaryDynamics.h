/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

#ifndef PHASE_FIELD_BOUNDARY_DYNAMICS_H
#define PHASE_FIELD_BOUNDARY_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of Dirichlet boundary condition for the order parameter.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
class PhaseFieldInletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename CORRECTED_DYNAMICS::EquilibriumF;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template<typename M>
  using exchange_momenta = PhaseFieldInletDynamics<T,DESCRIPTOR,DYNAMICS,M,direction,orientation>;

  std::type_index id() override {
    return typeid(PhaseFieldInletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldInletDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    // Along all the commented parts of this code there will be an example based
    // on the situation where the wall's normal vector if (0,1) and the
    // numerotation of the velocites are done according to the D2Q9
    // lattice of the OpenLB library.

    // Find all the missing populations
    // (directions 3,4,5)
    constexpr auto missingIndices = util::subIndexOutgoing<DESCRIPTOR,direction,orientation>();

    V phi = MomentaF().computeRho(cell);
    V missingPhi = phi - 1.;
    V missingWeightSum = 0;
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {

      bool contains = false;
      for(unsigned i=0; i < missingIndices.size(); i++){
        if(missingIndices[i] == (unsigned)iPop){
          contains = true;
        }
      }

      if(contains){
        missingWeightSum += descriptors::t<V,DESCRIPTOR>(iPop);
      } else {
        missingPhi -= cell[iPop];
      }
    }

    for (unsigned iPop=0; iPop < missingIndices.size(); ++iPop) {
      cell[missingIndices[iPop]] = missingPhi * descriptors::t<V,DESCRIPTOR>(missingIndices[iPop]) / missingWeightSum;
    }

    return typename CORRECTED_DYNAMICS::CollisionO().apply(cell, parameters);
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "PhaseFieldInletDynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

/**
* Implementation of convective boundary condition for the order parameter.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
struct PhaseFieldConvectiveOutletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  //using parameters = typename DYNAMICS::parameters;
  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  std::type_index id() override {
    return typeid(PhaseFieldConvectiveOutletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldConvectiveOutletDynamics>>();
  }

  constexpr static bool is_vectorizable = false;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V phi = MomentaF().computeRho(cell);
    V u[DESCRIPTOR::d];
    MomentaF().computeU(cell, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    return {phi, uSqr};
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "PhaseFieldConvectiveOutletDynamics<" + DYNAMICS().getName() + ">";
  };

};

template <typename T, typename DESCRIPTOR>
using NoCollideDynamicsExternalVelocity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::ExternalVelocityTuple,
  equilibria::FirstOrder,
  collision::None
>;
///Initialising the setConvectivePhaseFieldBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setConvectivePhaseFieldBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material)
{
  setConvectivePhaseFieldBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setConvectivePhaseFieldBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setConvectivePhaseFieldBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setConvectivePhaseFieldBoundary");

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); iCloc++) {
    setConvectivePhaseFieldBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc));
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  int _overlap = 1;
  {
    auto& communicator = sLattice.getCommunicator(stage::PostCollide());
    communicator.template requestField<descriptors::POPULATION>();

    SuperIndicatorBoundaryNeighbor<T,DESCRIPTOR::d> neighborIndicator(std::forward<decltype(indicator)>(indicator), _overlap);
    communicator.requestOverlap(_overlap, neighborIndicator);
    communicator.exchangeRequests();
  }
}


/// Set ConvectivePhaseFieldBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setConvectivePhaseFieldBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator)
{
  using namespace boundaryhelper;
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = 1;
  std::vector<int> discreteNormal(3,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin && indicator(iX, iY)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY);
      if ((abs(discreteNormal[1]) + abs(discreteNormal[2])) == 1) {
        block.addPostProcessor(
          typeid(stage::PostCollide), {iX,iY},
          promisePostProcessorForNormal<T,DESCRIPTOR,FlatConvectivePhaseFieldPostProcessorA2D>(
            Vector<int,2>(discreteNormal.data() + 1)));
        block.addPostProcessor(
          typeid(stage::PostStream), {iX,iY},
          promisePostProcessorForNormal<T,DESCRIPTOR,FlatConvectivePhaseFieldPostProcessorB2D>(
            Vector<int,2>(discreteNormal.data() + 1)));
        block.template defineDynamics<NoCollideDynamicsExternalVelocity>({iX, iY});
      } else {
        throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
      }
    }
  });
}

}

#endif
