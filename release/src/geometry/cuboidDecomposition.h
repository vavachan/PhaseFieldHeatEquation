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

#ifndef CUBOID_DECOMPOSITION_H
#define CUBOID_DECOMPOSITION_H

#include <vector>
#include <fstream>

#include "core/singleton.h"
#include "core/vector.h"

#include "geometry/cuboid.h"

#include "io/ostreamManager.h"
#include "io/xmlReader.h"

#include "cuboid.h"
#include "cuboidDecompositionMinimizer.h"

namespace olb {

template <typename T> class LoadBalancer;

/// Decomposition of a physical volume into a set of disjoint cuboids
/**
 * The union of all cuboids is a superset of the decomposed volume.
 * Cuboids are dimensioned s.t. they exactly represent a regular
 * lattice with fixed spacing.
 * Two cuboids are neighbors if the distance between them is less
 * than the lattice spacing.
 */
template <typename T, unsigned D>
class CuboidDecomposition final {
private:
  /// Cuboid which contains all other cuboids
  Cuboid<T,D> _motherCuboid;
  /// Vector of the cuboids
  std::vector<Cuboid<T,D>> _cuboids;
  /// Periodicity flag
  Vector<bool,D> _periodicityOn;

public:
  /// Constructs cuboid decomposition of cuboid with origin and extent
  CuboidDecomposition(Vector<T,D> origin, T deltaR, Vector<int,D> extent, int nC=1);
  /// Construction from an given mother cuboid
  CuboidDecomposition(const Cuboid<T,D>& motherCuboid, int nC);
  /// Constructs a cuboid decomposition with uniform spacing of voxelSize which consists of nC cuboids
  CuboidDecomposition(IndicatorF<T,D>& indicatorF, T voxelSize, int nC=1);
  /// Constructs a cuboid decomposition with uniform spacing of voxelSize which consists of nC cuboids
  /// shrunken using the minimizeBy strategy.
  CuboidDecomposition(IndicatorF<T,D>& indicatorF, T voxelSize, int nC, std::string minimizeBy);

  /// Returns number of cuboids in decomposition
  int size() const;

  /// Returns ID of cuboid containing physR within padding
  std::optional<int> getC(Vector<T,D> physR, int padding = 0) const;
  /// Returns physical position of lattice position
  Vector<T,D> getPhysR(LatticeR<D+1> latticeR) const;
  /// Returns lattice position for given physical position if it exists
  std::optional<LatticeR<D+1>> getLatticeR(Vector<T,D> physR) const;
  /// Returns floor lattice position for given physical position if it exists
  std::optional<LatticeR<D+1>> getFloorLatticeR(Vector<T,D> physR) const;

  /// Returns spacing between lattice points in physical units
  T getDeltaR() const;

  /// Returns true iff physR is covered by the decomposition
  bool isInside(Vector<T,D> physR) const;

  /// Read access to a single cuboid
  const Cuboid<T,D>& get(int iC) const;
  /// Read and write access to a single cuboid
  Cuboid<T,D>& get(int iC);
  /// Returns the smallest cuboid that includes all cuboids of the structure
  const Cuboid<T,D>& getMotherCuboid() const;

  /// Returns the smallest cuboid that includes all cuboids of the structure
  Cuboid<T,D>& getMotherCuboid();

  /// Set flag to enable/disable periodicity depending of direction. Be aware that not all directions are true to ensure boundary conditions like for velocity are not disturbed.
  void setPeriodicity(Vector<bool,D> periodicity);

  std::vector<Cuboid<T,D>>& cuboids() {
    return _cuboids;
  }

  /// Returns set of neighbors to cuboid iCglob within overlap
  std::set<int> getNeighborhood(int iCglob, int overlap = 0) const;

  /// Sets the number of full cells of each cuboid
  void setWeights(IndicatorF<T,D>& indicatorF) requires (D == 3);

  /// Splits cuboid iC, removes it and adds p cuboids of same volume
  void split(int iC, int p);
  /// Splits cuboid iC, removes it, adds approx. width^3 sized new cuboids
  void splitRegular(int iC, int width);
  /// Splits cuboid iC, removes it and adds p cuboids of same weight
  void splitByWeight(int iC, int p, IndicatorF<T,D>& indicatorF) requires (D == 3);
  /// Splits cuboid iC along dimension iD into cuboids of fractions
  void splitFractional(int iC, int iD, std::vector<T> fractions);

  /// Removes the cuboid iC
  void remove(int iC);
  /// Removes all cuboids where indicatorF = 0
  void remove(IndicatorF<T,D>& indicatorF);
  /// Removes all cuboids where weight = 0
  void removeByWeight();

  /// Shrink cuboid iC so that no empty planes are left
  void shrink(int iC, IndicatorF<T,D>& indicatorF);
  /// Shrink all cuboids so that no empty planes are left
  void shrink(IndicatorF<T,D>& indicatorF);

  /// Refines mesh by splitting each cell into factor^3 cells
  void refine(int factor);
  /// Tries to refine mesh to given deltaR
  bool tryRefineTo(T deltaR);

  /// Prints cuboid geometry details
  void print() const;
  /// Prints cuboid geometry details plus details of all cuboids
  void printExtended();

  /// Save CuboidDecomposition into an existing XML File
  void writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer);
  /// Save CuboidDecomposition into XML File
  void writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer);

  /// Returns the minimum of the ratio nX/nY/nZ in the structure
  T getMinRatio() const;
  /// Returns the maximum of the ratio nX/nY/nZ in the structure
  T getMaxRatio() const;

  /// Returns the minimum coordinate in the structure
  Vector<T,D> getMinPhysR() const;
  /// Returns the maximum coordinate in the structure
  Vector<T,D> getMaxPhysR() const;

  /// Returns the minimum volume in the structure
  T getMinPhysVolume() const;
  /// Returns the maximum volume in the structure
  T getMaxPhysVolume() const;

  /// Returns the minimum number of nodes in the structure
  std::size_t getMinLatticeVolume() const;
  /// Returns the maximum number of nodes in the structure
  std::size_t getMaxLatticeVolume() const;

  /// Returns the minimum number of nodes in the structure inside the indicator
  std::size_t getMinLatticeWeight() const;
  /// Returns the maximum number of nodes in the structure inside the indicator
  std::size_t getMaxLatticeWeight() const;

  /// Returns the total number cells (without overlap)
  std::size_t getNumNodes() const;
};

template <typename T>
using CuboidDecomposition2D = CuboidDecomposition<T,2>;
template <typename T>
using CuboidDecomposition3D = CuboidDecomposition<T,3>;

/// Load CuboidDecomposition from XML File
template<typename T, unsigned D>
std::unique_ptr<CuboidDecomposition<T,D>> createCuboidDecomposition(std::string fileName);

}

#endif
