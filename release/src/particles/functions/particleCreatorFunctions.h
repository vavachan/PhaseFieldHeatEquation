/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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

#ifndef PARTICLE_CREATOR_FUNCTIONS_H
#define PARTICLE_CREATOR_FUNCTIONS_H

#include <sstream>
#include <string>
#include <vector>

#include "particles/communication/particleParallelCreatorFunctions3D.h"
#include "particles/functions/particleContactForceFunctions.h"
#include "particles/functions/particleCreatorFunctions2D.h"
#include "particles/functions/particleCreatorFunctions3D.h"
#include "particles/particles.h"

namespace olb {

namespace particles {

namespace creators {

template <typename T, unsigned D>
struct SpawnData {
  SpawnData(const PhysR<T, D>& _position,
            const Vector<T, utilities::dimensions::convert<D>::rotation>&
                _angleInDegree)
      : position(_position)
      , angleInDegree(_angleInDegree)
  {}

  PhysR<T, D>                                            position;
  Vector<T, utilities::dimensions::convert<D>::rotation> angleInDegree;
};

template <typename T, unsigned D>
T calcParticleVolume(
    const std::vector<SpawnData<T, D>>&         spawnData,
    const std::function<T(const std::size_t&)>& getParticleVolume)
{
  T           particleVolume(0);
  std::size_t pID(0);

  // Use with C++20:
  //for(std::size_t pID = 0; [[maybe_unused]] const SpawnData<T,D>& entry : spawnData) {
  for ([[maybe_unused]] const SpawnData<T, D>& entry : spawnData) {
    particleVolume += getParticleVolume(pID++);
  }
  return particleVolume;
}

template <typename T, unsigned D>
T calcParticleVolumeFraction(
    const std::vector<SpawnData<T, D>>& spawnData, const T fluidVolume,
    const std::function<T(const std::size_t&)>& getParticleVolume)
{
  return calcParticleVolume(spawnData, getParticleVolume) / fluidVolume;
}

template <typename T, unsigned D>
void extendSpawnData(
    std::vector<SpawnData<T, D>>& spawnData,
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T                                     fluidDomainVolume,
    const std::function<T(const std::size_t&)>& getCircumRadius,
    const std::function<T(const std::size_t&)>& getParticleVolume)
{
  constexpr unsigned maxTries = 1e6;

  OstreamManager clout(std::cout, "ParticleSeeding");
  //Create randomizer
  util::Randomizer<T> randomizer;

  std::size_t pID            = spawnData.size();
  unsigned    count          = 0;
  T    currVolumeOfParticles = calcParticleVolume(spawnData, getParticleVolume);
  bool isPositionValid;

  const Vector<T, D> extent = fluidDomain.getMax() - fluidDomain.getMin();

  while (currVolumeOfParticles / fluidDomainVolume <
             wantedParticleVolumeFraction &&
         count < maxTries) {
    const T circumRadius = getCircumRadius(pID);
    isPositionValid      = true;

    // Randomize position
    PhysR<T, D> randomPosition = randomizer.template generate<PhysR<T, D>>();
    randomPosition *= extent;
    randomPosition += fluidDomain.getMin();

    // check overlap with other particles
    std::size_t i = 0;
    for (const SpawnData<T, D>& entry : spawnData) {
      const T distanceCenterOfMass = norm(entry.position - randomPosition);
      if (distanceCenterOfMass <= circumRadius + getCircumRadius(i++)) {
        ++count;
        isPositionValid = false;
        break;
      }
    }
    if (!isPositionValid) {
      continue;
    }

    // check overlap with boundaries (should work for convex geometries)
    // for better results an iteration over the whole volume could be added
    auto pointsOnHull = particles::discrete_points_on_hull::calculate(
        randomPosition, circumRadius);
    for (const PhysR<T, D>& pointOnHull : pointsOnHull) {
      bool isInside;
      fluidDomain(&isInside, pointOnHull.data());
      if (!isInside) {
        isPositionValid = false;
        break;
      }
    }

    if (isPositionValid) {
      spawnData.push_back(SpawnData<T, D>(
          randomPosition,
          Vector<T, utilities::dimensions::convert<D>::rotation>(T {0})));
      count = 0;
      currVolumeOfParticles += getParticleVolume(pID++);
    }
    else {
      ++count;
    }
  }
  if (count >= maxTries) {
    clout << "WARNING: Could not set a particle volume fraction of "
          << wantedParticleVolumeFraction
          << ", only a particle volume fraction of "
          << currVolumeOfParticles / fluidDomainVolume << " was possible."
          << std::endl;
  }

  return;
}

template <typename T, unsigned D>
std::vector<SpawnData<T, D>> setParticles(
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T                                     fluidDomainVolume,
    const std::function<T(const std::size_t&)>& getCircumRadius,
    const std::function<T(const std::size_t&)>& getParticleVolume,
    const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
        createParticle)
{
  std::vector<SpawnData<T, D>> spawnData;
  extendSpawnData<T, D>(spawnData, wantedParticleVolumeFraction, fluidDomain,
                        fluidDomainVolume, getCircumRadius, getParticleVolume);

  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    createParticle(entry, pID++);
  }

  return spawnData;
}

std::vector<std::vector<std::string>>
readParticlePositions(const std::string& filename)
{
  std::vector<std::vector<std::string>> content;
  std::vector<std::string>              row;
  std::string                           line, word;

  std::fstream file(filename, std::ios::in);
  if (file.is_open()) {
    while (getline(file, line)) {
      row.clear();

      std::stringstream str(line);

      while (getline(str, word, ';')) {
        row.push_back(word);
      }
      content.push_back(row);
    }
  }
  else {
    std::cerr << "Could not open the file " << filename << std::endl;
  }

  return content;
}

template <typename T, unsigned D>
void saveParticlePositions(
    const std::string& filename, const std::vector<SpawnData<T, D>>& spawnData,
    const std::function<std::string(const std::size_t&)>& evalIdentifier)
{
  std::ofstream file;
  file.open(filename.c_str(), std::ios::trunc);
  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    for (unsigned iD = 0; iD < D; ++iD) {
      file << std::setprecision(16) << entry.position[iD] << ';';
    }
    for (unsigned iD = 0; iD < utilities::dimensions::convert<D>::rotation;
         ++iD) {
      file << std::setprecision(16) << entry.angleInDegree[iD] << ';';
    }
    file << evalIdentifier(pID++) << std::endl;
  }
  file.close();
}

template <typename T, unsigned D>
std::vector<SpawnData<T, D>> setParticles(
    const std::string& filename, const T wantedParticleVolumeFraction,
    IndicatorF<T, D>& fluidDomain, const T fluidDomainVolume,
    const std::function<T(const std::size_t&)>& getParticleVolume,
    const std::function<void(const SpawnData<T, D>&, const std::string&)>
        createParticle)
{
  std::vector<SpawnData<T, D>>          spawnData;
  std::vector<std::vector<std::string>> filecontent =
      particles::creators::readParticlePositions(filename);
  for (const std::vector<std::string>& line : filecontent) {
    PhysR<T, D>                                            particlePosition;
    Vector<T, utilities::dimensions::convert<D>::rotation> particleAngle;
    unsigned                                               iD = 0;
    for (; iD < D; ++iD) {
      particlePosition[iD] = std::stod(line[iD]);
    }
    for (; iD < (D + utilities::dimensions::convert<D>::rotation); ++iD) {
      particleAngle[iD - D] = std::stod(line[iD]);
    }
    SpawnData<T, D> tmpSpawnData(particlePosition, particleAngle);
    spawnData.push_back(tmpSpawnData);
    createParticle(tmpSpawnData, line[iD]);
  }

  return spawnData;
}

// Should only work if the particles have a similar size
// otherwise a bigger bounding box might overlap a smaller one without detection
// (could be improved by using a finer lattice / using much more points to check)
template <typename T, unsigned D>
void extendSpawnData(
    std::vector<SpawnData<T, D>>& spawnData,
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T fluidDomainVolume,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
        getMin,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
                                                getMax,
    const std::function<T(const std::size_t&)>& getParticleVolume)
{
  constexpr unsigned maxTries  = 1e6;
  const auto         getPoints = [&](const PhysR<T, D>& position,
                             const PhysR<T, D>& min, const PhysR<T, D>& max) {
    std::vector<PhysR<T, D>> points;
    if constexpr (D == 3) {
      points.reserve(9);
      points.push_back(PhysR<T, D>(min[0], min[1], min[2]));
      points.push_back(PhysR<T, D>(min[0], min[1], max[2]));
      points.push_back(PhysR<T, D>(min[0], max[1], min[2]));
      points.push_back(PhysR<T, D>(min[0], max[1], max[2]));
      points.push_back(PhysR<T, D>(max[0], min[1], max[2]));
      points.push_back(PhysR<T, D>(max[0], min[1], min[2]));
      points.push_back(PhysR<T, D>(max[0], max[1], min[2]));
      points.push_back(PhysR<T, D>(max[0], max[1], max[2]));
    }
    else {
      points.reserve(5);
      points.push_back(PhysR<T, D>(min[0], min[1]));
      points.push_back(PhysR<T, D>(min[0], max[1]));
      points.push_back(PhysR<T, D>(max[0], min[1]));
      points.push_back(PhysR<T, D>(max[0], max[1]));
    }
    points.push_back(position);
    return points;
  };

  OstreamManager clout(std::cout, "ParticleSeeding");
  //Create randomizer
  util::Randomizer<T> randomizer;

  std::size_t pID            = spawnData.size();
  unsigned    count          = 0;
  T    currVolumeOfParticles = calcParticleVolume(spawnData, getParticleVolume);
  bool isPositionValid;

  const Vector<T, D> extent = fluidDomain.getMax() - fluidDomain.getMin();

  while (currVolumeOfParticles / fluidDomainVolume <
             wantedParticleVolumeFraction &&
         count < maxTries) {
    isPositionValid = true;

    // Randomize position
    PhysR<T, D> randomPosition = randomizer.template generate<PhysR<T, D>>();
    randomPosition *= extent;
    randomPosition += fluidDomain.getMin();
    const std::vector<PhysR<T, D>> pointsToCheck =
        getPoints(randomPosition, getMin(pID, randomPosition),
                  getMax(pID, randomPosition));

    // check overlap with other particles
    for (std::size_t iP = 0; iP < spawnData.size(); ++iP) {
      const PhysR<T, D> min = getMin(iP, spawnData[iP].position);
      const PhysR<T, D> max = getMax(iP, spawnData[iP].position);

      for (const PhysR<T, D>& point : pointsToCheck) {
        /*
        const int count = std::count(cornerPoints.begin(), cornerPoints.end(), corner);
        if (count > 1) {
          clout<< "WARNING: Same point was checked twice during the particle seeding";
        }
        */

        bool overlap = true;
        for (unsigned iD = 0; iD < D; ++iD) {
          overlap = overlap && min[iD] <= point[iD] && max[iD] >= point[iD];
        }
        if (overlap) {
          isPositionValid = false;
          break;
        }
      }

      if (!isPositionValid) {
        break;
      }
    }
    if (!isPositionValid) {
      ++count;
      continue;
    }

    // check overlap with boundaries (should work for convex geometries)
    // for better results an iteration over the whole volume could be added
    for (const PhysR<T, D>& point : pointsToCheck) {
      bool isInside;
      fluidDomain(&isInside, point.data());
      if (!isInside) {
        isPositionValid = false;
        break;
      }
    }

    if (isPositionValid) {
      spawnData.push_back(SpawnData<T, D>(
          randomPosition,
          Vector<T, utilities::dimensions::convert<D>::rotation>(T {0.})));
      count = 0;
      currVolumeOfParticles += getParticleVolume(pID++);
    }
    else {
      ++count;
    }
  }
  if (count >= maxTries) {
    clout << "WARNING: Could not set a particle volume fraction of "
          << wantedParticleVolumeFraction
          << ", only a particle volume fraction of "
          << currVolumeOfParticles / fluidDomainVolume << " was possible."
          << std::endl;
  }

  return;
}

template <typename T, unsigned D>
std::vector<SpawnData<T, D>> setParticles(
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T fluidDomainVolume,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
        getMin,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
                                                getMax,
    const std::function<T(const std::size_t&)>& getParticleVolume,
    const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
        createParticle)
{
  std::vector<SpawnData<T, D>> spawnData;
  extendSpawnData<T, D>(spawnData, wantedParticleVolumeFraction, fluidDomain,
                        fluidDomainVolume, getMin, getMax, getParticleVolume);

  for (std::size_t pID = 0; pID < spawnData.size(); ++pID) {
    createParticle(spawnData[pID], pID);
  }
  return spawnData;
}

template <typename T, unsigned D>
void extendSpawnData(
    std::vector<SpawnData<T, D>>& spawnData,
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T fluidDomainVolume, const T deltaX, const T contactDetectionDistance,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
        getMin,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
                                                getMax,
    const std::function<T(const std::size_t&, const SpawnData<T, D>&,
                          const PhysR<T, D>&)>& signedDistanceToParticle,
    const std::function<T(const std::size_t&)>& getParticleVolume)
{
  OstreamManager     clout(std::cout, "ParticleSeeding");
  constexpr unsigned maxTries = 1e6;

  //Create randomizer
  util::Randomizer<T> randomizer;

  std::size_t pID         = spawnData.size();
  unsigned    count       = 0;
  T currVolumeOfParticles = calcParticleVolume(spawnData, getParticleVolume);

  const Vector<T, D> extent = fluidDomain.getMax() - fluidDomain.getMin();

  while (currVolumeOfParticles / fluidDomainVolume <
             wantedParticleVolumeFraction &&
         count < maxTries) {
    bool isPositionValid = true;

    // Randomize position
    PhysR<T, D> randomPosition = randomizer.template generate<PhysR<T, D>>();
    randomPosition *= extent;
    randomPosition += fluidDomain.getMin();

    const PhysR<T, D> min = getMin(pID, randomPosition);
    const PhysR<T, D> max = getMax(pID, randomPosition);
    Vector<int, D>    start;
    Vector<int, D>    end;
    for (unsigned iD = 0; iD < D; ++iD) {
      start[iD] = util::floor(min[iD] / deltaX);
      end[iD]   = util::ceil(max[iD] / deltaX);
    }

    const auto evalCurrentPosition = [&](const Vector<int, D>& pos,
                                         bool&                 breakLoop) {
      const PhysR<T, D> physPos = deltaX * pos;
      if (signedDistanceToParticle(
              pID,
              SpawnData<T, D>(
                  randomPosition,
                  Vector<T, utilities::dimensions::convert<D>::rotation>(
                      T {0.})),
              physPos) < contactDetectionDistance) {

        // Check wall intersection
        if (fluidDomain.signedDistance(physPos) > -contactDetectionDistance) {
          breakLoop       = true;
          isPositionValid = false;
          return;
        }

        // Check intersection with particles
        for (std::size_t i = 0; i < spawnData.size(); ++i) {
          if (signedDistanceToParticle(i, spawnData[i], physPos) <
              contactDetectionDistance) {
            breakLoop       = true;
            isPositionValid = false;
            return;
          }
        }
      }
    };

    particles::contact::forEachPositionWithBreak(start, end,
                                                 evalCurrentPosition);

    if (!isPositionValid) {
      ++count;
    }
    else {
      spawnData.push_back(SpawnData<T, D>(
          randomPosition,
          Vector<T, utilities::dimensions::convert<D>::rotation>(T {0})));
      currVolumeOfParticles += getParticleVolume(pID++);
      count = 0;
    }
  }

  if (count >= maxTries) {
    clout << "WARNING: Could not set a particle volume fraction of "
          << wantedParticleVolumeFraction
          << ", only a particle volume fraction of "
          << currVolumeOfParticles / fluidDomainVolume << " was possible."
          << std::endl;
  }

  return;
}

template <typename T, unsigned D>
std::vector<SpawnData<T, D>> setParticles(
    const T wantedParticleVolumeFraction, IndicatorF<T, D>& fluidDomain,
    const T fluidDomainVolume, const T deltaX, const T contactDetectionDistance,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
        getMin,
    const std::function<PhysR<T, D>(const std::size_t&, const PhysR<T, D>&)>&
                                                getMax,
    const std::function<T(const std::size_t&, const SpawnData<T, D>&,
                          const PhysR<T, D>&)>& signedDistanceToParticle,
    const std::function<T(const std::size_t&)>& getParticleVolume,
    const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
        createParticle)

{
  std::vector<SpawnData<T, D>> spawnData;
  extendSpawnData<T, D>(spawnData, wantedParticleVolumeFraction, fluidDomain,
                        fluidDomainVolume, deltaX, contactDetectionDistance,
                        getMin, getMax, signedDistanceToParticle,
                        getParticleVolume);

  for (std::size_t pID = 0; pID < spawnData.size(); ++pID) {
    createParticle(spawnData[pID], pID);
  }
  return spawnData;
}

/// Updates particle positions so that they can be easily written to a txt file with the function above
// TODO: Needs testing
template <typename T, typename PARTICLETYPE>
std::vector<SpawnData<T, PARTICLETYPE::d>>
updateParticlePositions(std::vector<SpawnData<T, PARTICLETYPE::d>> spawnData,
                        XParticleSystem<T, PARTICLETYPE>& particleSystem
#ifdef PARALLEL_MODE_MPI
                        ,
                        MPI_Comm particleCreatorComm = MPI_COMM_WORLD
#endif
)
{
  using namespace descriptors;

#ifdef PARALLEL_MODE_MPI
  // Prepare a set of destination ranks from all other ranks
  // -> every rank will know the position
  auto& cuboidDecomposition = particleSystem.getSuperStructure().getCuboidDecomposition();
  auto& loadBalancer   = particleSystem.getSuperStructure().getLoadBalancer();
  std::set<int> destRanksSet;
  for (int iC = 0; iC < cuboidDecomposition.size(); ++iC) {
    int rank = loadBalancer.rank(iC);
    if (rank != singleton::mpi().getRank()) {
      destRanksSet.insert(rank);
    }
  }

  // Preparations
  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;
  using GENERAL_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::GENERAL>;
  using POSITION_EVAL =
      typename GENERAL_EVAL ::template derivedField<descriptors::POSITION>;
  using SURFACE_EVAL =
      typename PARTICLETYPE ::template derivedField<descriptors::SURFACE>;
  using ANGLE_EVAL =
      typename SURFACE_EVAL ::template derivedField<descriptors::ANGLE>;

  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, descriptors::ID> fieldID(1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, POSITION_EVAL> fieldPosition(
      1);
  FieldArrayD<T, PARTICLETYPE, Platform::CPU_SISD, ANGLE_EVAL> fieldAngle(1);

  auto communicatableID       = ConcreteCommunicatable(fieldID);
  auto communicatablePosition = ConcreteCommunicatable(fieldPosition);
  auto communicatableAngle    = ConcreteCommunicatable(fieldAngle);
  const std::vector<unsigned int> indices {0};
  const std::size_t               serialSize = communicatableID.size(indices) +
                                 communicatablePosition.size(indices) +
                                 communicatableAngle.size(indices);

  // Obtain data for communication
  communication::forParticlesInSuperParticleSystem<
      T, PARTICLETYPE,
      conditions::valid_particle_centres //only consider center for resolved
      >(particleSystem, [&](Particle<T, PARTICLETYPE>&       particle,
                            ParticleSystem<T, PARTICLETYPE>& bParticleSystem,
                            int                              globiC) {
    fieldID.setField(0, particle.template getField<PARALLELIZATION, ID>());
    fieldPosition.setField(0, access::getPosition(particle));
    fieldAngle.setField(0, access::getAngle(particle));
    spawnData[fieldID.getField(0)].position = fieldPosition.getField(0);
    spawnData[fieldID.getField(0)].angleInDegree =
        util::radianToDegree(fieldAngle.getField(0));
    for (const int destRank : destRanksSet) {
      std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSize] {});
      std::uint8_t*                   bufferRaw = buffer.get();
      std::size_t serialIdx = communicatableID.serialize(indices, bufferRaw);
      serialIdx +=
          communicatablePosition.serialize(indices, &bufferRaw[serialIdx]);
      serialIdx +=
          communicatableAngle.serialize(indices, &bufferRaw[serialIdx]);
      dataMap.insert(std::make_pair(destRank, std::move(buffer)));
    }
  });

  //Create non blocking mpi helper
  singleton::MpiNonBlockingHelper mpiNbHelper;

  std::map<int, std::vector<std::uint8_t>> rankDataMapSorted;
  communication::fillSendBuffer(dataMap, rankDataMapSorted, serialSize);

  // Send mapped data
  communication::sendMappedData(rankDataMapSorted, destRanksSet, serialSize,
                                particleCreatorComm, mpiNbHelper);

  // Receive positions from all other ranks
  communication::receiveAndExecuteForData(
      destRanksSet, serialSize, particleCreatorComm, mpiNbHelper,
      [&](int rankOrig, std::uint8_t* buffer) {
        std::size_t serialIdx = communicatableID.deserialize(indices, buffer);
        serialIdx +=
            communicatablePosition.deserialize(indices, &buffer[serialIdx]);
        serialIdx +=
            communicatableAngle.deserialize(indices, &buffer[serialIdx]);

        spawnData[fieldID.getField(0)].position = fieldPosition.getField(0);
        spawnData[fieldID.getField(0)].angleInDegree =
            util::radianToDegree(fieldAngle.getField(0));
      });

#else
  for (int i = 0; i < particleSystem.size(); ++i) {
    spawnData[i].position = access::getPosition(particleSystem.get(i));
    spawnData[i].angleInDegree =
        util::radianToDegree(access::getAngle(particleSystem.get(i)));
  }
#endif

  return spawnData;
}

template <typename T, typename PARTICLETYPE>
std::vector<SpawnData<T, PARTICLETYPE::d>> saveUpdatedParticlePositions(
    const std::string&                                    filename,
    const std::vector<SpawnData<T, PARTICLETYPE::d>>&     originalSpawnData,
    const std::function<std::string(const std::size_t&)>& evalIdentifier,
    XParticleSystem<T, PARTICLETYPE>&                     particleSystem
#ifdef PARALLEL_MODE_MPI
    ,
    MPI_Comm particleCreatorComm = MPI_COMM_WORLD
#endif
)
{
  std::vector<SpawnData<T, PARTICLETYPE::d>> spawnData =
      updateParticlePositions<T, PARTICLETYPE>(originalSpawnData, particleSystem
#ifdef PARALLEL_MODE_MPI
                                               ,
                                               particleCreatorComm
#endif
      );
  saveParticlePositions(filename, spawnData, evalIdentifier);
  return spawnData;
}

} //namespace creators

} //namespace particles

} //namespace olb

#endif
