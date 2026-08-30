/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito, Adrian Kummerländer
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

#ifndef IMMERSED_BOUNDARY_METHOD_3D_H
#define IMMERSED_BOUNDARY_METHOD_3D_H

#include "stencil.h"

namespace olb {

namespace ibm {

template <int WIDTH>
struct InterpolateVelocityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::converter::PHYS_DELTA_X
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;

    Vector<V,3> u{}; // interpolated vel. of particle with physR
    auto particle = cells.template get<names::Points>(); // current particle assoc. with cells
    auto f_cell = cells.template get<names::NavierStokes>(); // closest lattice cell to particle
    auto f_physR = f_cell.template getField<descriptors::LOCATION>(); // physR of fluid cell
    auto p_physR = particle.template getField<fields::PHYS_R>(); // physR of particle
    V dx = parameters.template get<fields::converter::PHYS_DELTA_X>(); // delta x of fluid cell
    Vector<int,3> rel_lat_min;
    for (int d = 0; d < 3; ++d) {
      rel_lat_min[d] = static_cast<int>(util::floor((p_physR[d] - f_physR[d]) / dx - WIDTH / 2.) + 1);
    }
    for (int x = rel_lat_min[0]; x < WIDTH + rel_lat_min[0]; ++x) {
      for (int y = rel_lat_min[1]; y < WIDTH + rel_lat_min[1]; ++y) {
        for (int z = rel_lat_min[2]; z < WIDTH + rel_lat_min[2]; ++z) {
          auto n_physR = f_cell.neighbor({x, y, z}).template getField<descriptors::LOCATION>(); // physR of cell inside stencil
          auto dist = (n_physR - p_physR) / dx; // lat. dist.
          Vector<V,3> u_tmp{};
          f_cell.neighbor({x, y, z}).computeU(u_tmp.data()); // get fluid vel.
          u_tmp *= Stencil<V,WIDTH>::weight(dist[0]) *
                   Stencil<V,WIDTH>::weight(dist[1]) *
                   Stencil<V,WIDTH>::weight(dist[2]); // mult. interp. weights
          u += u_tmp; // add contribution
        }
      }
    }
    particle.template setField<fields::membrane::VELOCITY>(u);
  }
};

template <int WIDTH>
struct SpreadForceO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<
    fields::converter::PHYS_DELTA_X
  >;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;

    auto particle = cells.template get<names::Points>(); // current particle assoc. with cells
    auto f_cell = cells.template get<names::NavierStokes>(); // closest lattice cell to particle
    auto f_physR = f_cell.template getField<descriptors::LOCATION>(); // physR of fluid cell
    auto p_physR = particle.template getField<fields::PHYS_R>(); // physR of particle
    V dx = parameters.template get<fields::converter::PHYS_DELTA_X>(); // delta x of fluid cell
    Vector<int,3> rel_lat_min;
    for (int d = 0; d < 3; ++d) {
      rel_lat_min[d] = static_cast<int>(util::floor((p_physR[d] - f_physR[d]) / dx - WIDTH / 2.) + 1);
    }

    for (int x = rel_lat_min[0]; x < WIDTH + rel_lat_min[0]; ++x) {
      for (int y = rel_lat_min[1]; y < WIDTH + rel_lat_min[1]; ++y) {
        for (int z = rel_lat_min[2]; z < WIDTH + rel_lat_min[2]; ++z) {
          auto n_physR = f_cell.neighbor({x, y, z}).template getField<descriptors::LOCATION>(); // physR of cell inside stencil
          auto dist = (n_physR - p_physR) / dx; // lat. dist.
          auto f_F = f_cell.neighbor({x,y,z}).template getField<descriptors::FORCE>(); // old fluid force
          auto p_F = particle.template getField<fields::membrane::FORCE>(); // get solid-node force
          p_F *= Stencil<V,WIDTH>::weight(dist[0]) *
                 Stencil<V,WIDTH>::weight(dist[1]) *
                 Stencil<V,WIDTH>::weight(dist[2]); // mult. interp. weights
          f_cell.neighbor({x,y,z}).template setField<descriptors::FORCE>(f_F + p_F); // update by adding force contribution
        }
      }
    }
  }
};

}

/// Wrapper class to perform the immersed boundary method
template<int WIDTH, typename COUPLEES>
class SuperImmersedBoundaryCoupling3D {
private:
  template <typename VALUED_DESCRIPTOR>
  using ptr_to_lattice = std::conditional_t<
    VALUED_DESCRIPTOR::descriptor_t::template provides<descriptors::POPULATION>(),
    SuperLattice<typename VALUED_DESCRIPTOR::value_t,
                 typename VALUED_DESCRIPTOR::descriptor_t>,
    SuperD<typename VALUED_DESCRIPTOR::value_t,
           typename VALUED_DESCRIPTOR::descriptor_t>
  >*;

  utilities::TypeIndexedTuple<typename COUPLEES::template map_values<
    ptr_to_lattice
  >> _lattices;

  std::unique_ptr<SuperLatticePointCoupling<ibm::InterpolateVelocityO<WIDTH>,COUPLEES>> _interpolateO;
  std::unique_ptr<SuperLatticePointCoupling<ibm::SpreadForceO<WIDTH>,COUPLEES>> _spreadForceO;
public:
  template<typename... MAP>
  SuperImmersedBoundaryCoupling3D(std::integral_constant<int,WIDTH>, MAP&&... args)
  {
    auto map = std::make_tuple(&args...);
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      constexpr unsigned idx = COUPLEES::keys_t::template index<name_t>();
      _lattices.template set<name_t>(std::get<2*idx+1>(map));
    });

    _interpolateO = std::make_unique<SuperLatticePointCoupling<ibm::InterpolateVelocityO<WIDTH>,COUPLEES>>(
      ibm::InterpolateVelocityO<WIDTH>{},
      names::NavierStokes{}, *(_lattices.get(meta::id<names::NavierStokes>{})),
      names::Points{}, *(_lattices.get(meta::id<names::Points>{})));
    _spreadForceO = std::make_unique<SuperLatticePointCoupling<ibm::SpreadForceO<WIDTH>,COUPLEES>>(
      ibm::SpreadForceO<WIDTH>{},
      names::NavierStokes{}, *(_lattices.get(meta::id<names::NavierStokes>{})),
      names::Points{}, *(_lattices.get(meta::id<names::Points>{})));

    // set the LOCATION FIELD
    auto sLattice = _lattices.get(meta::id<names::NavierStokes>{});
    auto& load = sLattice->getLoadBalancer();
    for (int iC=0; iC < load.size(); ++iC) {
      sLattice->getBlock(iC).forSpatialLocations([&](LatticeR<3> latticeR) {
        int latR[4] {iC, latticeR[0], latticeR[1], latticeR[2]};
        auto physR = sLattice->getCuboidDecomposition().getPhysR(latR);
        sLattice->getBlock(iC).get(latticeR).template setField<descriptors::LOCATION>(physR);
      });
    }
  }

  template <typename T>
  void setPhysDeltaX(T dx)
  {
    _interpolateO->template setParameter<fields::converter::PHYS_DELTA_X>(dx);
    _spreadForceO->template setParameter<fields::converter::PHYS_DELTA_X>(dx);
  }

  void execute()
  {
    // purge force field
    auto sLattice = _lattices.get(meta::id<names::NavierStokes>{});
    auto& load = sLattice->getLoadBalancer();
    using V = typename std::remove_reference<decltype(*sLattice)>::type::value_t;
    for (int iC=0; iC < load.size(); ++iC) {
      sLattice->getBlock(iC).forSpatialLocations([&](LatticeR<3> latticeR) {
        int latR[4] {0, latticeR[0], latticeR[1], latticeR[2]};
        sLattice->getBlock(iC).get(latR[1],latR[2],latR[3]).template setField<descriptors::FORCE>(V{0.0});
      });
    }

    // spread force from membrane to fluid
    _spreadForceO->execute();

    // interpolate fluid velocity on membrane surface
    _interpolateO->execute();
  }
};

template <int WIDTH, typename... MAP>
SuperImmersedBoundaryCoupling3D(std::integral_constant<int,WIDTH>, MAP&&...)
  -> SuperImmersedBoundaryCoupling3D<
    WIDTH,
    typename meta::map<MAP...>::template map_values<descriptors::extract_valued_descriptor_t>
>;

}

#endif
