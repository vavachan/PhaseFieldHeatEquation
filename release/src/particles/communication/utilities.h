/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022-2023 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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


#ifndef PARTICLE_COMMUNICATION_UTILITIES_H
#define PARTICLE_COMMUNICATION_UTILITIES_H

#include "particles/contact/contactFunctions.h"

namespace olb {

///Evaluate complete neighbourhood and store in std::map
template<typename T, unsigned D>
void evaluateCuboidDecompositionNeighbourhood(
  CuboidDecomposition<T,D>& cuboidDecomposition,
  std::map<int,std::vector<int>>& neighbourhood,
  int offset )
{
  for (int iC=0; iC<cuboidDecomposition.size(); ++iC){
    neighbourhood[iC].clear();
    auto neighbours = cuboidDecomposition.getNeighborhood(iC, offset);
    for (int jC : neighbours) {
      neighbourhood[iC].emplace_back(jC);
    }
  }
}

/// Consistency check for neighbour retrieval
/// - workaround for issue #319
bool checkCuboidNeighbourhoodConsistency(
  std::map<int,std::vector<int>>& neighbourhood,
  bool correct = false,
  bool verbose = false)
{
  bool consistent = true;
  for ( auto cuboidNeighbourPair : neighbourhood){
    //Retrieve iC and neighbours
    int iC = cuboidNeighbourPair.first;
    std::vector<int>& neighbours = cuboidNeighbourPair.second;
    //Loop over neighbours
    for (int iCN : neighbours){
      //Retreive neighbour's neighbours
      std::vector<int>& neighboursNeighbours = neighbourhood[iCN];
      bool iCfound = false;
      //Loop over neighbour's neighbours
      for (int iCNN : neighboursNeighbours ){
        if (iCNN == iC){
          iCfound=true;
          break;
        }
      }
      if (!iCfound){
        //Set consistency boolean to false
        consistent = false;
        //Output, if desired
        if (verbose){
          std::cout << "iC " << iC << " not found in list of neighbour iC "
                             << iCN << std::endl;
        }
        //Correct, if desired
        if (correct){
          neighbourhood[iCN].push_back(iC);
        }
      }
    }
  }
  //Return whether consistent
  return consistent;
}




namespace particles {

namespace communication {


/// Get neighbour cuboids
template<typename T, unsigned D, bool verbose=false>
std::vector<int> getNeighborCuboids(
  CuboidDecomposition<T,D>& cuboidDecomposition, unsigned offset, int globiC )
{
  auto neighbors = cuboidDecomposition.getNeighborhood(globiC, offset);

  //Print, if desired
  if constexpr(verbose){
    std::cout << globiC << ": " << neighbors << std::endl;
  }

  std::vector<int> ns;
  for (int iC : neighbors) {
    ns.emplace_back(iC);
  }
  return ns;
}

/// Get neighbour ranks
template<typename T, unsigned D, bool verbose=false>
std::set<int> getNeighbourRanksFromCuboidNeighbourhood(
  SuperStructure<T,D>& superStructure, int rank, std::map<int,std::vector<int>> neighbourhood )
{
  std::set<int> rankNeighbours;

#ifdef PARALLEL_MODE_MPI
  auto& loadBalancer = superStructure.getLoadBalancer();

  for ( auto cuboidNeighbourPair : neighbourhood){
    int globiC = cuboidNeighbourPair.first;
    if (loadBalancer.rank(globiC) == rank) {
      std::vector<int>& neighbours = cuboidNeighbourPair.second;
      for (int globiCN : neighbours ){
        rankNeighbours.insert( loadBalancer.rank(globiCN));
      }
    }
  }

  if constexpr(verbose){
    std::cout << rank << ": " << rankNeighbours << std::endl;
  }
#endif

  return rankNeighbours;
}


/// Get neighbour ranks
template<typename T, unsigned D, bool verbose=false>
std::set<int> getNeighbourRanksFromCuboidNeighbourhood(
  SuperStructure<T,D>& superStructure, int rank, const std::vector<std::set<int>>& neighbourhood )
{
  std::map<int,std::vector<int>> neighbourhoodMap;

  for(unsigned i=0; i<neighbourhood.size(); ++i) {
    neighbourhoodMap[i] = std::vector<int>();
    for(int iC : neighbourhood[i]) {
      neighbourhoodMap[i].push_back(iC);
    }
  }

  return getNeighbourRanksFromCuboidNeighbourhood<T,D,verbose>(
    superStructure, rank, neighbourhoodMap );
}

/// Get neighbour ranks
template<typename T, unsigned D, bool verbose=false>
std::set<int> getNeighbourRanksFromCuboidNeighbourhood(
  SuperStructure<T,D>& superStructure, const std::vector<std::set<int>>& neighbourhood )
{
  const int rank = singleton::mpi().getRank();
  return getNeighbourRanksFromCuboidNeighbourhood<T,D,verbose>(
    superStructure, rank, neighbourhood);
}


/// Get neighbour ranks
template<typename T, unsigned D, bool verbose=false>
std::set<int> getNeighbourRanks(
  SuperStructure<T,D>& superStructure, unsigned offset, int rank )
{
  std::set<int> rankNeighbours;

#ifdef PARALLEL_MODE_MPI
  auto& loadBalancer = superStructure.getLoadBalancer();
  auto& cuboidDecomposition = superStructure.getCuboidDecomposition();

  ///HOTFIX VERSION (according to #319)
  // 1. sample complete neighbourhood
  std::map<int,std::vector<int>> neighbourhood;
  evaluateCuboidDecompositionNeighbourhood(cuboidDecomposition, neighbourhood, offset);

  //2. check consistency
  constexpr bool correctFaultyNeighbourhood = true;
  checkCuboidNeighbourhoodConsistency( neighbourhood, correctFaultyNeighbourhood );

  //3. evaluate relevant ranks
  rankNeighbours = getNeighbourRanksFromCuboidNeighbourhood<T,D,verbose>(
    superStructure, rank, neighbourhood );
#endif

  return rankNeighbours;
}

/// Get neighbour ranks
template<typename T, unsigned D, bool verbose=false>
std::set<int> getNeighbourRanks(
  SuperStructure<T,D>& superStructure, unsigned offset)
{
  //Retrive parallization infos
  const int rank = singleton::mpi().getRank();
  return getNeighbourRanks<T,D,false>(superStructure, offset, rank);
}


template<typename T, typename PARTICLETYPE>
std::size_t serializeIcDest( int globiCdest,
  std::uint8_t* bufferiCandParticle )
{
  //Define iC field as FieldArrayD with single entry
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,descriptors::IC> fieldiC(1);
  fieldiC.setField(0, globiCdest);
  const std::vector<unsigned int> index{0};
  auto communicatable = ConcreteCommunicatable(fieldiC);
  //Serialize iC
  communicatable.serialize(index, bufferiCandParticle);
  //Return communicatable size
  return communicatable.size(index);
}

template<typename T, typename PARTICLETYPE>
void serializeIcDestAndParticle( int globiCdest,
  Particle<T,PARTICLETYPE>& particle,
  std::uint8_t* bufferiCandParticle )
{
  //Serialize iC and retrieve particle buffer
  std::size_t serialSizeField =
    serializeIcDest<T,PARTICLETYPE>( globiCdest, bufferiCandParticle );
  std::uint8_t* bufferRawParticle = &bufferiCandParticle[serialSizeField];
  //Serialize particle
  particle.serialize(bufferRawParticle);
}

template<typename T, typename PARTICLETYPE>
std::size_t deserializeIcDest( int& globiCdest,
  std::uint8_t* bufferRaw )
{
  //Define iC field as FieldArrayD with single entry
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,descriptors::IC> fieldiC(1);
  const std::vector<unsigned int> index{0};
  auto communicatable = ConcreteCommunicatable(fieldiC);
  //Deserialice iC
  communicatable.deserialize(index, bufferRaw);
  globiCdest = fieldiC.getField(0);
  //Return communicatable size
  return communicatable.size(index);
}


template<typename T, typename PARTICLETYPE>
void deserializeIcDestAndParticle( int& globiCdest,
  Particle<T,PARTICLETYPE>& particle,
  std::uint8_t* bufferiCandParticle )
{
  //Deserialize iC and retrieve particle buffer
  std::size_t serialSizeField =
    deserializeIcDest<T,PARTICLETYPE>(globiCdest,bufferiCandParticle);
  std::uint8_t* bufferRawParticle = &bufferiCandParticle[serialSizeField];
  //Deserialize particle
  particle.deserialize(bufferRawParticle);
}


/*
 * The function provides a check if the cuboids fit the constraints from the
 * implementation of the particle parallization
 */
template<typename T, typename PARTICLETYPE>
void checkCuboidSizes( SuperParticleSystem<T, PARTICLETYPE>& sParticleSystem )
{
  constexpr unsigned D = PARTICLETYPE::d;
  auto& cuboidDecomposition = sParticleSystem.getSuperStructure().getCuboidDecomposition();

  const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR();

  T maxCircumRadius{0};
  communication::forParticlesInSuperParticleSystem(
      sParticleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
          maxCircumRadius = util::max(maxCircumRadius,
              particles::contact::evalCircumRadius(particle, physDeltaX));
      });
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(maxCircumRadius, MPI_MAX,
                                    singleton::mpi().bossId());
#endif

  int minCuboidExtent{std::numeric_limits<int>::max()};
  for (int iC=0; iC<cuboidDecomposition.size(); ++iC){
    for(unsigned iD=0; iD<D; ++iD){
      minCuboidExtent = util::min(minCuboidExtent,
          cuboidDecomposition.get(iC).getExtent()[iD]);
    }
  }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(minCuboidExtent, MPI_MIN,
                                    singleton::mpi().bossId());
#endif

  // It's enough that one rank checks for an error
  if (singleton::mpi().getRank() == 0) {
    OLB_ASSERT(minCuboidExtent*physDeltaX > maxCircumRadius,
        "At least one cuboid too small so that the currently implemented "
        "particle parallelization strategy does not work.");
  }
}

template<typename T>
T movePositionToStart(const T position, const T max, const T min)
{
  return min + (position - max);
}

template<typename T>
T movePositionToEnd(const T position, const T max, const T min)
{
  return max - (min - position);
}

template<typename T>
T applyPeriodicityToPosition(const bool isPeriodic, T position, const T max, const T min)
{
  if(isPeriodic) {
    if(position < min) {
      position = movePositionToEnd(position, max, min);
    }
    else if(position > max) {
      position = movePositionToStart(position, max, min);
    }
  }
  return position;
}

// The following functions use an offset of 0.5 * physDeltaX, because there is
// a distance of 1 * physDeltaX between the maximal coordinate and minimal coordinate

/// Returns minimal coordinate of domain for periodic particle boundaries
template<typename T, unsigned D>
Vector<T,D> getCuboidMin(const CuboidDecomposition<T,D>& cuboidDecomposition)
{
  const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR();
  return (cuboidDecomposition.getMotherCuboid().getOrigin() - 0.5 * physDeltaX);

}

/// Returns maximal coordinate of domain for periodic particle boundaries
template<typename T, unsigned D>
Vector<T,D> getCuboidMax(const CuboidDecomposition<T,D>& cuboidDecomposition,
    const PhysR<T,D>& min)
{
  const T physDeltaX = cuboidDecomposition.getMotherCuboid().getDeltaR();
  return min + physDeltaX * cuboidDecomposition.getMotherCuboid().getExtent();
}

/// Updates a position if out of bounds and periodic setup is used
template<typename T, unsigned D>
Vector<T,D> applyPeriodicityToPosition(
    const CuboidDecomposition<T,D>& cuboidDecomposition,
    const Vector<bool,D>& periodicity,
    PhysR<T,D> position)
{
  const PhysR<T,D> min = getCuboidMin<T,D>(cuboidDecomposition);
  const PhysR<T,D> max = getCuboidMax<T,D>(cuboidDecomposition, min);

  for(unsigned iD=0; iD<D; ++iD) {
    position[iD] = applyPeriodicityToPosition(periodicity[iD],
        position[iD], max[iD], min[iD]);
  }

  return position;
}

/// Function returns true if cuboid was found and gives iC
template<typename T, unsigned D>
bool getCuboid(const CuboidDecomposition<T,D>& cuboidDecomposition,
    const Vector<bool,D>& periodicity,
    const PhysR<T,D>& position, int& iC)
{
  PhysR<T,D> newPosition(
      applyPeriodicityToPosition(cuboidDecomposition, periodicity, position));
  if (auto iC_ = cuboidDecomposition.getC(newPosition)) {
    iC = *iC_;
    return true;
  }
  return false;
}

} //namespace communication

} //namespace particles

} //namespace olb


#endif
