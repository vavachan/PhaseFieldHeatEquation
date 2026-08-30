/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef SET_BOUNDARY_H
#define SET_BOUNDARY_H

#include <optional>

#include "setBoundary2D.h"
#include "setBoundary3D.h"

#include "normalDynamicsContructors.h"

#include "core/concepts.h"

#include "geometry/discreteNormals.h"

namespace olb {

namespace concepts {

/// Concept of a discrete-normal-dependent boundary condition
/**
 * This is a very common case for LBM boundary conditions.
 *
 * In OpenLB such conditions are modeled as the combination
 * of Dynamics for the local treatment and a post processor
 * for the non-local part. Individual BCs may define only
 * one of these options depending on the discrete normal.
 *
 * Types fullfilling this concept can be applied using the
 * `boundary::set` setter function.
 **/
template <typename BC>
concept BoundaryCondition = requires(BC bc,
                                     DiscreteNormalType type,
                                     DiscreteNormal<typename BC::descriptor_t> n) {
  requires BaseType<typename BC::value_t>;
  requires LatticeDescriptor<typename BC::descriptor_t>;

  /// Optionally maps discrete normal to dynamics
  { bc.getDynamics(type, n) } -> std::same_as<
    std::optional<DynamicsPromise<typename BC::value_t,
                                  typename BC::descriptor_t>>
  >;
  /// Optionally maps discrete normal to post processor
  { bc.getPostProcessor(type, n) } -> std::same_as<
    std::optional<PostProcessorPromise<typename BC::value_t,
                                       typename BC::descriptor_t>>
  >;
  /// Defines required neighborhood radius for communication (of post processor)
  { bc.getNeighborhoodRadius() } -> std::same_as<CellDistance>;
};

}

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, concepts::BoundaryCondition BC>
void set(BlockLattice<T,DESCRIPTOR>& block,
         BlockIndicatorF<T,DESCRIPTOR::d>& boundaryI,
         BlockIndicatorF<T,DESCRIPTOR::d>& fluidI,
         BlockIndicatorF<T,DESCRIPTOR::d>& outsideI)
{
  OstreamManager clout(std::cout, "boundary::set");
  auto& blockGeometryStructure = boundaryI.getBlockGeometry();
  blockGeometryStructure.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
    if (   blockGeometryStructure.getNeighborhoodRadius(latticeR) >= util::max(BC().getNeighborhoodRadius(), 1)
        && boundaryI(latticeR)) {
      const auto [normalType, normal] = computeBoundaryTypeAndNormal(fluidI, outsideI, latticeR);
      if (auto dynamics = BC().getDynamics(normalType, normal)) {
        block.defineDynamics(latticeR, std::move(*dynamics));
      }
      if (auto op = BC().getPostProcessor(normalType, normal)) {
        block.addPostProcessor(typeid(stage::PostStream),
                               latticeR,
                               std::move(*op));
      }
    }
  });
}

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, concepts::BoundaryCondition BC>
void set(SuperLattice<T,DESCRIPTOR>& sLattice,
         FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& boundaryI,
         FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& fluidI,
         FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& outsideI)
{
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    set<T,DESCRIPTOR,BC>(sLattice.getBlock(iCloc),
                         boundaryI->getBlockIndicatorF(iCloc),
                         fluidI->getBlockIndicatorF(iCloc),
                         outsideI->getBlockIndicatorF(iCloc));
  }
  addPoints2CommBC(sLattice,
                   std::forward<decltype(boundaryI)>(boundaryI),
                   BC().getNeighborhoodRadius());
}

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, concepts::BoundaryCondition BC>
void set(SuperLattice<T,DESCRIPTOR>& sLattice,
         FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  auto fluidI = indicator->getSuperGeometry().getMaterialIndicator(1);
  auto outsideI = indicator->getSuperGeometry().getMaterialIndicator(0);
  set<T,DESCRIPTOR,BC>(sLattice,
                       std::forward<decltype(indicator)>(indicator),
                       fluidI,
                       outsideI);
}

template <template<typename...> typename BC, concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void set(SuperLattice<T,DESCRIPTOR>& sLattice,
         FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  set<T,DESCRIPTOR,BC<T,DESCRIPTOR>>(sLattice,
                                     std::forward<decltype(indicator)>(indicator));
}

template <concepts::BoundaryCondition BC>
void set(SuperLattice<typename BC::value_t,typename BC::descriptor_t>& sLattice,
         FunctorPtr<SuperIndicatorF<typename BC::value_t,BC::descriptor_t::d>>&& indicator)
{
  set<typename BC::value_t,
      typename BC::descriptor_t,
      BC>(sLattice,
          std::forward<decltype(indicator)>(indicator));
}

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, concepts::BoundaryCondition BC>
void set(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
         int material)
{
  set<T,DESCRIPTOR,BC>(sLattice,
                       sGeometry.getMaterialIndicator(material));
}

template <template<typename...> typename BC, concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void set(SuperLattice<T,DESCRIPTOR>& sLattice,
         SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
         int material)
{
  set<T,DESCRIPTOR,BC<T,DESCRIPTOR>>(sLattice,
                                     sGeometry.getMaterialIndicator(material));
}

template <concepts::BoundaryCondition BC>
void set(SuperLattice<typename BC::value_t,typename BC::descriptor_t>& sLattice,
         SuperGeometry<typename BC::value_t,BC::descriptor_t::d>& sGeometry,
         int material)
{
  set<typename BC::value_t,
      typename BC::descriptor_t,
      BC>(sLattice,
          sGeometry.getMaterialIndicator(material));
}

}

}

#endif
