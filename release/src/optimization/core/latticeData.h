/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef LATTICE_DATA_H
#define LATTICE_DATA_H

namespace olb {

template<typename T, typename DESCRIPTOR> class LatticeResults;

/// Encapsulates all necessary classes to run a simulation
template<typename T, typename DESCRIPTOR>
class LatticeData {
private:
  std::unique_ptr<UnitConverter<T,DESCRIPTOR>> _converter;
  std::unique_ptr<CuboidDecomposition<T,DESCRIPTOR::d>> _decomposition;
  std::unique_ptr<LoadBalancer<T>> _balancer;
  std::unique_ptr<SuperGeometry<T,DESCRIPTOR::d>> _geometry;
  std::unique_ptr<SuperLattice<T,DESCRIPTOR>> _lattice;
  std::unique_ptr<LatticeResults<T,DESCRIPTOR>> _results;
  std::string _name;

public:
  LatticeData (std::string name) : _name(name)
    { };

  /// Create instances and return unique pointer access
  template<typename MEMBER, typename... ARGS>
  auto& create(ARGS&&... args) {
    if constexpr (std::is_base_of_v<UnitConverter<T,DESCRIPTOR>,MEMBER>) {
      _converter = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_converter;
    }
    else if constexpr (std::is_same_v<CuboidDecomposition<T,DESCRIPTOR::d>,MEMBER>) {
      _decomposition = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_decomposition;
    }
    else if constexpr (std::is_base_of_v<LoadBalancer<T>,MEMBER>) {
      _balancer = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_balancer;
    }
    else if constexpr (std::is_same_v<SuperGeometry<T,DESCRIPTOR::d>,MEMBER>) {
      _geometry = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_geometry;
    }
    else if constexpr (std::is_same_v<SuperLattice<T,DESCRIPTOR>,MEMBER>) {
      _lattice = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_lattice;
    }
    else if constexpr (std::is_same_v<LatticeResults<T,DESCRIPTOR>,MEMBER>) {
      _results = std::make_unique<MEMBER>(std::forward<ARGS>(args)...);
      return *_results;
    } else {
      throw std::runtime_error("Unsupported type to create.");
    }
  }

  auto& getUnitConverter() {
    if (!_converter) {
      throw std::runtime_error("LatticeData: Missing UnitConverter.");
    }
    return *_converter;
  }

  auto& getCuboidDecomposition() {
    if (!_decomposition) {
      throw std::runtime_error("LatticeData: Missing CuboidDecomposition.");
    }
    return *_decomposition;
  }

  auto& getLoadBalancer() {
    if (!_balancer) {
      throw std::runtime_error("LatticeData: Missing LoadBalancer.");
    }
    return *_balancer;
  }

  auto& getSuperGeometry() {
    if (!_geometry) {
      throw std::runtime_error("LatticeData: Missing SuperGeometry.");
    }
    return *_geometry;
  }

  auto& getSuperLattice() {
    if (!_lattice) {
      throw std::runtime_error("LatticeData: Missing SuperLattice.");
    }
    return *_lattice;
  }

  auto& getLatticeResults() {
    if (!_results) {
      throw std::runtime_error("LatticeData: Missing LatticeResults.");
    }
    return *_results;
  }

  /// Used to reset simulations
  void resetLattice() {
    _lattice.reset(new SuperLattice<T,DESCRIPTOR>(*_geometry));
  }

  /// Used for differentiating different latticeData instances
  std::string getName() {
    return _name;
  }
};

}

#endif
