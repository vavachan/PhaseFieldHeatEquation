/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022-2023 Nicolas Hafen, Jan E. Marquardt
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

#ifndef SUPER_PARTICLE_SYSTEM_H
#define SUPER_PARTICLE_SYSTEM_H

namespace olb {

namespace particles {

template <typename T, typename PARTICLETYPE>
class SuperParticleSystem {
private:
  std::vector<ParticleSystem<T, PARTICLETYPE>*> _blockParticleSystems;
  SuperStructure<T, PARTICLETYPE::d>            _superStructure;
  /// Cache for all ranks responsible for neighboring cuboids
  std::set<int>                       _neighbourRanks;
  /// Cache for all ranks responsible for neighboring cuboids and their neighboring cuboids
  std::set<int>                       _extendedNeighbourRanks;
  /// Cached cuboid neighborhood: all cuboid neighbors of all cuboids
  std::vector<std::set<int>>          _cuboidNeighborhood;
  std::size_t                                   _currentGlobalID = 0;
  const std::size_t                             _serialSize;
  T                                             _maximalCircumRadius = T{0.};

  unsigned                                      getOffset();
  void                                          updateCuboidNeighborhood();

public:
  SuperParticleSystem(SuperStructure<T, PARTICLETYPE::d>& superStructure, T maximalCircumRadius = T{0.});
  /// Create SuperParticleSystem from SuperGeometry and optional maximal circumference radius (if already known)
  SuperParticleSystem(SuperGeometry<T, PARTICLETYPE::d>& superGeometry, T maximalCircumRadius = T{0.});
  void print();
  void addDynamics(std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr );
  template <typename DYNAMICS, typename ...Args>
  void defineDynamics(Args&& ...args);
  std::vector<ParticleSystem<T, PARTICLETYPE>*>& getBlockParticleSystems();
  SuperStructure<T, PARTICLETYPE::d>&            getSuperStructure();
  const std::set<int>&                 getNeighbourRanks();
  const std::set<int>&                 getExtendedNeighbourRanks();
  const std::vector<std::set<int>>&    getCuboidNeighborhood();
  std::size_t                                    getGlobID();
  std::size_t                                    getSerialSize() const;
  Particle<T, PARTICLETYPE>                      get(std::size_t globalParticleID);
  template<typename PCONDITION=conditions::all_particles>
  constexpr std::size_t                          size(); //Returns total number of particles
  void                                           checkForErrors();
  void                                           updateOffsetFromCircumRadius(T circumRadius);
};

//Define type that can be both ParticleSystem and SuperParticleSystem
//- depending on particle type
template <typename T, typename PARTICLETYPE>
using XParticleSystem = std::conditional_t<
    PARTICLETYPE::template providesNested<descriptors::PARALLELIZATION>(),
    SuperParticleSystem<T, PARTICLETYPE>, ParticleSystem<T, PARTICLETYPE>>;


//Simple stuctur that stores information to avoid an additional particle search
//-intended to be returnd during particle creation
//-synchronize(): ensures the same localID on all cores (most likely never necessary,
//  as localID is only valid in the respective cuboid )
struct ParallelParticleLocator{
  //Stored quantities
  int globiC;
  std::size_t globalID;
  std::size_t localID;
  //Constructor
  ParallelParticleLocator()
    : globiC(0), globalID(0), localID(0)
  {}
  ParallelParticleLocator(int globiC_, std::size_t globalID_, std::size_t localID_)
    : globiC(globiC_), globalID(globalID_), localID(localID_)
  {}

  //Synchronize
  void synchronize(){
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(globiC, MPI_MAX);
    singleton::mpi().reduceAndBcast(globalID, MPI_MAX);
    singleton::mpi().reduceAndBcast(localID, MPI_MAX);
#endif
  }
  //Print
  void print(){
    std::cout  << "ParallelParticleLocator: ("
               << "globiC="   << globiC
               << ", globID=" << globalID
               << ", locID="  << localID << ")"
               << std::endl;
  }
};


} //namespace particles

} //namespace olb

#endif
