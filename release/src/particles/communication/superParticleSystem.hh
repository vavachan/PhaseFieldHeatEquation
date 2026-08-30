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

#ifndef SUPER_PARTICLE_SYSTEM_HH
#define SUPER_PARTICLE_SYSTEM_HH

namespace olb {

namespace particles {

template <typename T, typename PARTICLETYPE>
SuperParticleSystem<T, PARTICLETYPE>::SuperParticleSystem(
    SuperStructure<T, PARTICLETYPE::d>& superStructure, T maximalCircumRadius)
    : _superStructure(superStructure)
    , _serialSize(ParticleSystem<T, PARTICLETYPE>(1).getSerialSize())
    , _maximalCircumRadius(maximalCircumRadius)
{
  static_assert(access::providesValid<PARTICLETYPE>(), "Field GENERAL:VALID has to be provided");
  auto loadBalancer = superStructure.getLoadBalancer();
  // always create at least one BlockParticleSystem to counteract possible problems
  // (for example in creator functions, where the first blockParticleSystem is hard coded)
  _blockParticleSystems.push_back(new ParticleSystem<T, PARTICLETYPE>);
  for (int iC = 1; iC < loadBalancer.size(); ++iC) {
    _blockParticleSystems.push_back(new ParticleSystem<T, PARTICLETYPE>);
  }

  _cuboidNeighborhood.resize(_superStructure.getCuboidDecomposition().size());
  updateCuboidNeighborhood();
}

template <typename T, typename PARTICLETYPE>
SuperParticleSystem<T, PARTICLETYPE>::SuperParticleSystem(
    SuperGeometry<T, PARTICLETYPE::d>& superGeometry, T maximalCircumRadius)
    : SuperParticleSystem<T, PARTICLETYPE>(
          static_cast<SuperStructure<T, PARTICLETYPE::d>&>(superGeometry), maximalCircumRadius)
{}

template <typename T, typename PARTICLETYPE>
void SuperParticleSystem<T, PARTICLETYPE>::print()
{
  OstreamManager clout(std::cout, "SuperParticleSystem");
  auto           loadBalancer = _superStructure.getLoadBalancer();
  auto           rank         = singleton::mpi().getRank();
  clout << "---SuperParticleSystem---" << std::endl;
  clout.setMultiOutput(true);
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().synchronizeIO();
#endif
  for (int iC = 0; iC < loadBalancer.size(); ++iC) {
    int    globiC          = loadBalancer.glob(iC);
    auto&  bParticleSystem = *_blockParticleSystems[iC];
    size_t size            = bParticleSystem.size();
    clout << "BlockParticleSystem " << globiC << " (" << rank << "," << iC
          << ")"
          << ": nop=" << size << std::endl;
  }
  clout.setMultiOutput(false);
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().synchronizeIO();
#endif
  clout << "-------------------------" << std::endl;
}

template <typename T, typename PARTICLETYPE>
void SuperParticleSystem<T, PARTICLETYPE>::addDynamics(
  std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr )
{
  auto loadBalancer = _superStructure.getLoadBalancer();
  for (int iC = 0; iC < loadBalancer.size(); ++iC) {
    auto& bParticleSystem = *_blockParticleSystems[iC];
    bParticleSystem.addDynamics( dynamicsSPtr );
  }
}

template<typename T, typename PARTICLETYPE>
template <typename DYNAMICS, typename ...Args>
void SuperParticleSystem<T,PARTICLETYPE>::defineDynamics (Args&& ...args){
  auto loadBalancer = _superStructure.getLoadBalancer();
  for (int iC = 0; iC < loadBalancer.size(); ++iC) {
    auto& bParticleSystem = *_blockParticleSystems[iC];
    bParticleSystem.template defineDynamics<DYNAMICS>( std::forward<Args>(args)... );
  }
}

template <typename T, typename PARTICLETYPE>
const std::set<int>& SuperParticleSystem<T, PARTICLETYPE>::getNeighbourRanks()
{
  return _neighbourRanks;
}

template <typename T, typename PARTICLETYPE>
const std::set<int>& SuperParticleSystem<T, PARTICLETYPE>::getExtendedNeighbourRanks()
{
  return _extendedNeighbourRanks;
}

template <typename T, typename PARTICLETYPE>
const std::vector<std::set<int>>& SuperParticleSystem<T, PARTICLETYPE>::getCuboidNeighborhood()
{
  return _cuboidNeighborhood;
}

template <typename T, typename PARTICLETYPE>
std::vector<ParticleSystem<T, PARTICLETYPE>*>&
SuperParticleSystem<T, PARTICLETYPE>::getBlockParticleSystems()
{
  return _blockParticleSystems;
}

template <typename T, typename PARTICLETYPE>
SuperStructure<T, PARTICLETYPE::d>&
SuperParticleSystem<T, PARTICLETYPE>::getSuperStructure()
{
  return _superStructure;
}

template <typename T, typename PARTICLETYPE>
std::size_t SuperParticleSystem<T, PARTICLETYPE>::getGlobID()
{
  return _currentGlobalID++;
}

template <typename T, typename PARTICLETYPE>
std::size_t SuperParticleSystem<T, PARTICLETYPE>::getSerialSize() const
{
  return _serialSize;
}

//TODO: consider using searchParticle() in here instead
template <typename T, typename PARTICLETYPE>
Particle<T, PARTICLETYPE>
SuperParticleSystem<T, PARTICLETYPE>::get(std::size_t globalParticleID)
{
  using namespace descriptors;

  for (ParticleSystem<T, PARTICLETYPE>* blockParticleSystem :
       _blockParticleSystems) {
    for (std::size_t localiP = 0; localiP < blockParticleSystem->size();
         ++localiP) {
      auto              particle = blockParticleSystem->get(localiP);
      const std::size_t currentGlobalParticleID =
          particle.template getField<PARALLELIZATION, ID>();
      if (globalParticleID == currentGlobalParticleID) {
        return particle;
      }
    }
  }

  std::cerr << "ERROR: Rank " << singleton::mpi().getRank()
            << " cannot find particle with global id " << globalParticleID
            << std::endl;
  assert(false);
  return _blockParticleSystems[0]->get(0);
}

template <typename T, typename PARTICLETYPE>
template <typename PCONDITION>
constexpr std::size_t SuperParticleSystem<T, PARTICLETYPE>::size()
{
  std::size_t size = 0;
  communication::forSystemsInSuperParticleSystem<T,PARTICLETYPE>(
    *this,
    [&](ParticleSystem<T, PARTICLETYPE>& particleSystem, int iC, int globiC) {
      size += particleSystem.template size<PCONDITION>(globiC);
    });
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(size, MPI_SUM);
#endif
  return size;
}

template <typename T, typename PARTICLETYPE>
void SuperParticleSystem<T, PARTICLETYPE>::checkForErrors()
{
  //Check, whether at least one ParticleDynamics present
#ifdef PARALLEL_MODE_MPI
  if(singleton::mpi().isMainProcessor()){
#endif
    _blockParticleSystems[0]->checkForErrors();
#ifdef PARALLEL_MODE_MPI
  }
#endif
}

template <typename T, typename PARTICLETYPE>
void SuperParticleSystem<T, PARTICLETYPE>::updateOffsetFromCircumRadius(
    T circumRadius)
{
  if(_maximalCircumRadius < circumRadius) {
    _maximalCircumRadius = circumRadius;
    updateCuboidNeighborhood();
  }
}

template <typename T, typename PARTICLETYPE>
unsigned SuperParticleSystem<T, PARTICLETYPE>::getOffset()
{
  auto& cuboidDecomposition = _superStructure.getCuboidDecomposition();
  const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR();
  unsigned offset = util::max( util::ceil(_maximalCircumRadius/physDeltaX),
      _superStructure.getOverlap() );
  return offset;
}


template <typename T, typename PARTICLETYPE>
void SuperParticleSystem<T, PARTICLETYPE>::updateCuboidNeighborhood()
{
  auto& cuboidDecomposition = _superStructure.getCuboidDecomposition();
  int offset = getOffset();

  std::map<int,std::vector<int>> neighbourhood;
  evaluateCuboidDecompositionNeighbourhood(cuboidDecomposition, neighbourhood, offset);

  constexpr bool correctFaultyNeighbourhood = true;
  checkCuboidNeighbourhoodConsistency(neighbourhood, correctFaultyNeighbourhood);

  for (auto [globiC, neighbors] : neighbourhood) {
    _cuboidNeighborhood[globiC] = std::set<int>(neighbors.begin(), neighbors.end());
  }

  _neighbourRanks = communication::getNeighbourRanksFromCuboidNeighbourhood<T, PARTICLETYPE::d, false>(
    _superStructure, _cuboidNeighborhood);

  _extendedNeighbourRanks = _neighbourRanks;
  for(int rank : _neighbourRanks) {
    const std::set<int> neighborRanks = communication::getNeighbourRanksFromCuboidNeighbourhood<
      T, PARTICLETYPE::d, false>(_superStructure, rank, _cuboidNeighborhood);
    _extendedNeighbourRanks.insert(neighborRanks.begin(), neighborRanks.end());
  }
}


} //namespace particles

} //namespace olb

#endif
