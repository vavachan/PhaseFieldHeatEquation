/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
 *                2025 Adrian Kummerlaender
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

#ifndef CUBOID_H
#define CUBOID_H

#include "core/vector.h"
#include "core/serializer.h"

#include "utilities/vectorHelpers.h"

#include <optional>

namespace olb {

template <typename T, unsigned D>
class Cuboid final {
private:
  /// Global position of the lower left corner in SI meters
  Vector<T,D> _origin;
  /// Extent in number of nodes
  Vector<int,D> _extent;
  /// Spacing of nodes in SI meters
  T _delta;
  /// Weight for load balancing purposes
  std::size_t _weight;

  // Don't want to use concept here as it would have to be repeated for all method implementations
  static_assert(D == 2 || D == 3, "Only 2D/3D supported");

public:
  Cuboid() = default;
  Cuboid(Vector<T,D> origin, T delta, Vector<int,D> extent):
    _origin(origin),
    _extent(extent),
    _delta(delta),
    _weight(getLatticeVolume()) { }
  Cuboid(const Cuboid& rhs):
    _origin(rhs._origin),
    _extent(rhs._extent),
    _delta(rhs._delta),
    _weight(rhs._weight) { }
  Cuboid(IndicatorF<T,D>& indicatorF, T delta):
    _origin(indicatorF.getMin()),
    _extent((indicatorF.getMax() - indicatorF.getMin()) / delta + 1.5),
    _delta(delta),
    _weight(getLatticeVolume()) { }

  Cuboid& operator=(const Cuboid& rhs);
  bool operator==(const Cuboid& rhs) const;

  /// Returns lower left corner coordinates
  Vector<T,D> getOrigin() const { return _origin; }
  /// Returns spacing of cuboid nodes
  T getDeltaR() const { return _delta; }
  /// Returns extent in number of voxels
  Vector<int,D> getExtent() const { return _extent; }

  int getNx() const { return _extent[0]; }
  int getNy() const { return _extent[1]; }
  int getNz() const { return _extent[2]; }

  /// Returns the volume of the cuboid
  T getPhysVolume() const;

  std::size_t getWeight() const { return _weight; }
  void setWeight(std::size_t weight) { _weight = weight; }

  Vector<T,D> getPhysR(LatticeR<D> latticeR) const {
    return _origin + latticeR*_delta;
  }

  LatticeR<D> getLatticeR(Vector<T,D> physR) const {
    return util::floor((physR - _origin) / _delta + 0.5);
  }

  LatticeR<D> getFloorLatticeR(Vector<T,D> physR) const {
    return util::floor((physR - _origin) / _delta);
  }

  /// Returns closest latticeR within eps of physR if it exists
  std::optional<LatticeR<D>> getCloseLatticeR(Vector<T,D> physR, T eps=1e-5) const;

  /// Returns the number of Nodes in the volume
  std::size_t getLatticeVolume() const;
  /// Returns the perimeter of the cuboid
  T getPhysPerimeter() const requires (D == 2);
  /// Returns the perimeter of the cuboid
  T getPhysPerimeter() const requires (D == 3);
  /// Returns the number of Nodes at the perimeter
  std::size_t getLatticePerimeter() const requires (D == 2);
  /// Returns the number of Nodes at the perimeter
  std::size_t getLatticePerimeter() const requires (D == 3);

  /// Checks whether pos is contained in the cuboid extended with an layer of size overlap*delta
  bool isInside(Vector<T,D> pos, int overlap = 0) const;

  /// Checks whether there is an intersection with the cuboid extended by a layer of size overlap*delta
  bool intersects(Vector<T,D> globMin, Vector<T,D> globMax, int overlap = 0) const;
  /// Returns true iff self intersects cuboid
  bool intersects(const Cuboid<T,D>& cuboid) const;

  /// Divides the cuboid in p*q*r cuboids of equal volume and add them to the given vector
  void divide(Vector<int,D> division, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 3);
  /// Divides the cuboid in p*q cuboids of equal volume and add them to the given vector
  void divide(Vector<int,D> division, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 2);

  /// Divides the cuboid in p cuboids and add them to the given vector
  void divideP(int p, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 2);
  /// Divides the cuboid in p cuboids and add them to the given vector
  void divideP(int p, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 3);

  /// Divides the cuboid into fractions along the iDth dimension
  void divideFractional(int iD, std::vector<T> fractions, std::vector<Cuboid<T,D>>& childrenC) const;

  void resize(Vector<int,D> offset, Vector<int,D> extent);

  void refine(int factor);

  void write(std::ostream& cout) const;
  void print() const;

  void writeAsXML(std::ostream&) const;

};

template <typename T>
using Cuboid2D = Cuboid<T,2>;
template <typename T>
using Cuboid3D = Cuboid<T,3>;

}

#endif
