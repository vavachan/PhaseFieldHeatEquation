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

#ifndef SERIALIZE_H
#define SERIALIZE_H

namespace olb {

/// Returns serialized vector of field data inside indicated region
template <typename FIELD, typename T, typename DESCRIPTOR>
std::vector<T> getSerializedFromField(SuperLattice<T,DESCRIPTOR>& lattice,
                                      FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  auto& loadBalancer = lattice.getLoadBalancer();
  auto& cDecomposition = indicator->getSuperStructure().getCuboidDecomposition();
  auto& superGeometry = indicator->getSuperGeometry();
  const std::size_t serializedSize = getSerializedFieldSize<FIELD>(lattice,
    std::forward<FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>>(indicator));

  std::vector<T> serial;
  serial.resize(serializedSize, 0);
  std::size_t offset = 0;
  for (int iC=0; iC<cDecomposition.size(); ++iC) {
    std::size_t numNodes = 0;
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockGeometry = superGeometry.getBlockGeometry(loc);
      auto& blockLattice = lattice.getBlock(loc);
      auto& blockIndicator = indicator->getBlockIndicatorF(loc);
      blockGeometry.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> pos) {
        if (blockIndicator(pos)) {
          std::size_t index = offset + numNodes;
          for (std::size_t iD=0; iD<DESCRIPTOR::template size<FIELD>(); ++iD) {
            serial.at(index + iD) = blockLattice.get(pos).template getFieldComponent<FIELD>(iD);
            ++numNodes;
          }
        }
      });
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(numNodes, MPI_SUM);
#endif
    offset += numNodes;
  }
#ifdef PARALLEL_MODE_MPI
  std::vector<T> tmp;
  tmp.resize(serializedSize, 0);
  singleton::mpi().allreduce(serial.data(), tmp.data(), serializedSize, MPI_SUM);
  serial = tmp;
#endif
  return serial;
}

/// Set field content from a serialized data vector inside indicated region
template <typename FIELD, typename T, typename DESCRIPTOR>
void setFieldFromSerialized(const std::vector<T>& serial, SuperLattice<T,DESCRIPTOR>& lattice,
                            FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  auto& loadBalancer = lattice.getLoadBalancer();
  auto& cDecomposition = indicator->getSuperStructure().getCuboidDecomposition();
  auto& superGeometry = indicator->getSuperGeometry();

  std::size_t offset = 0;
  for (int iC=0; iC<cDecomposition.size(); ++iC) {
    std::size_t numNodes = 0;
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockGeometry = superGeometry.getBlockGeometry(loc);
      auto& blockLattice = lattice.getBlock(loc);
      auto& blockIndicator = indicator->getBlockIndicatorF(loc);
      blockGeometry.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> pos) {
        if (blockIndicator(pos)) {
          std::size_t index = offset + numNodes;
          FieldD<T,DESCRIPTOR,FIELD> field;
          for (std::size_t iD=0; iD<DESCRIPTOR::template size<FIELD>(); ++iD) {
            field[iD] = serial.at(index + iD);
            ++numNodes;
          }
          blockLattice.get(pos).template setField<FIELD>(field);
        }
      });
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(numNodes, MPI_SUM);
#endif
    offset += numNodes;
  }
}

/// Returns size of serialized data vector of a field
template <typename FIELD, typename T, typename DESCRIPTOR>
std::size_t getSerializedFieldSize(SuperLattice<T,DESCRIPTOR>& lattice) {
  const std::size_t numNodes = lattice.getCuboidDecomposition().getNumNodes();
  const std::size_t fieldDim = DESCRIPTOR::template size<FIELD>();
  return numNodes * fieldDim;
}

/// Returns size of serialized data vector of a field
template <typename FIELD, typename T, typename DESCRIPTOR>
std::size_t getSerializedFieldSize(SuperLattice<T,DESCRIPTOR>& lattice,
                                   FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator) {
  auto& loadBalancer = lattice.getLoadBalancer();
  auto& cDecomposition = indicator->getSuperStructure().getCuboidDecomposition();
  auto& superGeometry = indicator->getSuperGeometry();
  const std::size_t fieldDim = DESCRIPTOR::template size<FIELD>();

  std::size_t numNodes = 0;
  for (int iC=0; iC<cDecomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockGeometry = superGeometry.getBlockGeometry(loc);
      auto& blockIndicator = indicator->getBlockIndicatorF(loc);
      blockGeometry.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> pos) {
        if (blockIndicator(pos)) {
          ++numNodes;
        }
      });
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numNodes, MPI_SUM);
#endif
  return numNodes*fieldDim;
}

}

#endif
