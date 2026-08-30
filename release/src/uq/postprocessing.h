/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
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

#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include "uncertaintyQuantification.h"

namespace olb {

namespace uq {

// Internal helper struct for grid dimensions and spacing
template <unsigned D, typename T>
struct FieldGeometryInfo {
  std::array<int, D> dims {};
  unsigned           size {};
  unsigned           nodeCount {};
  Vector<T, D>       spacing {};
};

/**
    * \brief Read data from multiple samples' .vti files into a 4D vector
*/
template <typename T, typename DESCRIPTOR>
void readData(int samples_number, const std::string& uqFolder, const std::string& name, const std::string& dataName,
              const std::vector<std::size_t>& iT, std::size_t rank, std::vector<std::vector<std::vector<T>>>& data,
              FieldGeometryInfo<DESCRIPTOR::d, T>& info)
{
  // Read the first sample's VTI file to get geometry information
  std::string basePath = uqFolder + "0/tmp/";
  singleton::directories().setOutputDir(basePath);

  std::string fileName = singleton::directories().getVtkOutDir() + "data/" + createFileName(name, iT[0], rank) + ".vti";

  BlockVTIreader<T, T, DESCRIPTOR::d> reader0(fileName, dataName);
  BlockData<DESCRIPTOR::d, T, T>&     blockDataSample0 = reader0.getBlockData();

  // Read spacing from VTI file via XML
  XMLreader         xmlConfig(fileName);
  std::stringstream spacingStream(xmlConfig["ImageData"].getAttribute("Spacing"));

  for (int i = 0; i < DESCRIPTOR::d; ++i) {
    spacingStream >> info.spacing[i];
  }

  info.dims[0] = blockDataSample0.getNx();
  info.dims[1] = blockDataSample0.getNy();
  if constexpr (DESCRIPTOR::d == 3) {
    info.dims[2] = blockDataSample0.getNz();
  }
  info.nodeCount = std::accumulate(info.dims.begin(), info.dims.end(), 1, std::multiplies<> {});

  info.size = blockDataSample0.getSize();

  // Initialize data vector: data[nx*ny][size][samples_number] or data[nx*ny*nz][size][samples_number]
  // data[nodeIndex][iSize][iSample]
  data.assign(info.nodeCount, std::vector<std::vector<T>>(info.size, std::vector<T>(samples_number)));

  // Read data from each sample's .vti
  for (int iSample = 0; iSample < samples_number; ++iSample) {
    std::string samplePath = uqFolder + std::to_string(iSample) + "/tmp/";
    singleton::directories().setOutputDir(samplePath);

    std::string sampleFileName =
        singleton::directories().getVtkOutDir() + "data/" + createFileName(name, iT[iSample], rank) + ".vti";

    BlockVTIreader<T, T, DESCRIPTOR::d> reader(sampleFileName, dataName);
    BlockData<DESCRIPTOR::d, T, T>&     blockDataSample = reader.getBlockData();

    // Extract data for each lattice point
    if constexpr (DESCRIPTOR::d == 2) {
      for (int iy = 0; iy < info.dims[1]; ++iy) {
        for (int ix = 0; ix < info.dims[0]; ++ix) {
          int nodeIndex = ix + info.dims[0] * iy;
          // std::array<int, 2> coords = {ix, iy};

          for (unsigned iSize = 0; iSize < info.size; ++iSize) {
            data[nodeIndex][iSize][iSample] = blockDataSample.get({ix, iy}, iSize);
          }
        }
      }
    }
    else if constexpr (DESCRIPTOR::d == 3) {
      for (int iz = 0; iz < info.dims[2]; ++iz) {
        for (int iy = 0; iy < info.dims[1]; ++iy) {
          for (int ix = 0; ix < info.dims[0]; ++ix) {
            int nodeIndex = ix + info.dims[0] * (iy + info.dims[1] * iz);
            // std::array<int, 3> coords = {ix, iy, iz};

            for (unsigned iSize = 0; iSize < info.size; ++iSize) {
              data[nodeIndex][iSize][iSample] = blockDataSample.get({ix, iy, iz}, iSize);
            }
          }
        }
      }
    }
  }
}

/**
  * \brief Compute mean and std dev for each lattice node across samples
  */
template <typename T, unsigned D>
void computeMeanAndStd(UncertaintyQuantification<T>&                   uq,
                       const std::vector<std::vector<std::vector<T>>>& data, // [node][size][sample]
                       const FieldGeometryInfo<D, T>& info, BlockData<D, T, T>& blockDataMean,
                       BlockData<D, T, T>& blockDataStd)
{
  if constexpr (D == 2) {
    for (int iy = 0; iy < info.dims[1]; ++iy) {
      for (int ix = 0; ix < info.dims[0]; ++ix) {
        int nodeIndex = ix + info.dims[0] * iy;
        // std::array<int, 2> coords = {ix, iy};

        for (unsigned iSize = 0; iSize < info.size; ++iSize) {
          const auto& samples                = data[nodeIndex][iSize];
          blockDataMean.get({ix, iy}, iSize) = uq.mean(samples);
          blockDataStd.get({ix, iy}, iSize)  = uq.std(samples);
        }
      }
    }
  }
  else if constexpr (D == 3) {
    for (int iz = 0; iz < info.dims[2]; ++iz) {
      for (int iy = 0; iy < info.dims[1]; ++iy) {
        for (int ix = 0; ix < info.dims[0]; ++ix) {
          int nodeIndex = ix + info.dims[0] * (iy + info.dims[1] * iz);
          // std::array<int, 3> coords = {ix, iy, iz};

          for (unsigned iSize = 0; iSize < info.size; ++iSize) {
            const auto& samples                    = data[nodeIndex][iSize];
            blockDataMean.get({ix, iy, iz}, iSize) = uq.mean(samples);
            blockDataStd.get({ix, iy, iz}, iSize)  = uq.std(samples);
          }
        }
      }
    }
  }
}

/**
  * \brief Compute mean, std dev, and write them out to VTI for postprocessing
  */
template <typename T, typename DESCRIPTOR>
void computeMeanAndStdAndWriteVTI(UncertaintyQuantification<T>& uq, const std::string& uqFolder,
                                  const std::string& name, const std::string& dataName,
                                  CuboidDecomposition<T, DESCRIPTOR::d>& cuboidDecomposition,
                                  SuperGeometry<T, DESCRIPTOR::d>&       sGeometry)
{
  constexpr unsigned D = DESCRIPTOR::d;

  OstreamManager clout(std::cout, "computeMeanAndStdAndWriteVTI");

  // Build a super lattice for reading/writing the final data
  SuperLattice<T, DESCRIPTOR> sLattice(sGeometry);
  sLattice.template defineDynamics<NoDynamics>(sGeometry, /*material=*/0);
  sLattice.initialize();

  // Load iteration logs for each sample
  std::vector<std::vector<size_t>> iTList(uq.getSamplesNumber());
  size_t                           numIterations = std::numeric_limits<size_t>::max();

  for (size_t n = 0; n < uq.getSamplesNumber(); ++n) {
    std::string   filePath = uqFolder + std::to_string(n) + "/tmp/iteration_log.txt";
    std::ifstream inFile(filePath);
    if (inFile.is_open()) {
      size_t iVal;
      while (inFile >> iVal) {
        iTList[n].push_back(iVal);
      }
      inFile.close();
      numIterations = std::min(numIterations, iTList[n].size());
    }
    else {
      clout << "Could not open file: " << filePath << std::endl;
    }
  }

  SuperVTMwriter<T, D> vtmWriter(name);

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  // Rename materials for visualization
  for (int rank = 0; rank < noOfCuboids; ++rank) {
    Cuboid<T, D>& cuboid = cuboidDecomposition.get(rank);
    Vector<T, D>  origin = cuboid.getOrigin();
    origin[0] -= cuboid.getDeltaR();
    origin[1] -= cuboid.getDeltaR();
    if constexpr (DESCRIPTOR::d == 3) {
      origin[2] -= cuboid.getDeltaR();
    }

    Vector<T, D> extent = cuboid.getExtent() * cuboid.getDeltaR();
    extent[0] += cuboid.getDeltaR();
    extent[1] += cuboid.getDeltaR();
    if constexpr (DESCRIPTOR::d == 3) {
      extent[2] += cuboid.getDeltaR();
    }

    IndicatorCuboid<T, D> readCuboid(extent, origin);
    sGeometry.rename(/*fromMaterial=*/0, /*toMaterial=*/10 + rank, readCuboid);
  }

  // Process each iteration
  singleton::directories().setOutputDir("./tmp/");
  SuperLatticeCuboid<T, DESCRIPTOR> cuboid(sLattice);
  SuperLatticeRank<T, DESCRIPTOR>   rankFunctor(sLattice);

  vtmWriter.write(cuboid);
  vtmWriter.write(rankFunctor);
  vtmWriter.createMasterFile();

  for (size_t iter = 0; iter < numIterations; ++iter) {
    // Collect the iteration index from each sample
    std::vector<size_t> iT;
    iT.reserve(iTList.size());
    for (const auto& list : iTList) {
      iT.push_back(list[iter]);
    }

    clout << "Processing iteration " << iT[0] << std::endl;

    // Process data for each cuboid/rank
    for (int rank = 0; rank < noOfCuboids; ++rank) {
      std::string basePath = uqFolder + "0/tmp/";
      singleton::directories().setOutputDir(basePath);

      std::string fileName =
          singleton::directories().getVtkOutDir() + "data/" + createFileName(name, iT[0], rank) + ".vti";

      // We'll read the block data from the 0th sample's file,
      // then overwrite it with the mean/std from all samples.
      BlockVTIreader<T, T, D> readerMean(fileName, dataName);
      BlockVTIreader<T, T, D> readerStd(fileName, dataName);

      BlockData<D, T, T>& blockDataMean = readerMean.getBlockData();
      BlockData<D, T, T>& blockDataStd  = readerStd.getBlockData();

      // Gather all samples' data
      std::vector<std::vector<std::vector<T>>> data;
      FieldGeometryInfo<D, T>                  info;

      readData<T, DESCRIPTOR>(uq.getSamplesNumber(), uqFolder, name, dataName, iT, rank, data, info);

      // Overwrite blockDataMean, blockDataStd with computed mean, std
      computeMeanAndStd<T, D>(uq, data, info, blockDataMean, blockDataStd);

      // Create AnalyticalF2D wrappers to map blockData -> a lattice function
      BlockDataF<T, T, D> blockDataFMean(blockDataMean);
      BlockDataF<T, T, D> blockDataFStd(blockDataStd);

      SpecialAnalyticalFfromBlockF<T, T, D> meanField(blockDataFMean, cuboidDecomposition.get(rank), info.spacing,
                                                      /*scale=*/1.0);
      SpecialAnalyticalFfromBlockF<T, T, D> stdField(blockDataFStd, cuboidDecomposition.get(rank), info.spacing,
                                                     /*scale=*/1.0);

      // Define the newly computed fields on the lattice
      auto materialIndicator = sGeometry.getMaterialIndicator({10 + rank});
      for (int iBalancer = 0; iBalancer < sGeometry.getLoadBalancer().size(); ++iBalancer) {
        sLattice.getBlock(iBalancer).template defineField<descriptors::VELOCITY>(
            materialIndicator->getBlockIndicatorF(iBalancer), meanField);
        sLattice.getBlock(iBalancer).template defineField<descriptors::VELOCITY2>(
            materialIndicator->getBlockIndicatorF(iBalancer), stdField);
      }
    }

    // Write data to VTM file
    SuperLatticePhysField<T, DESCRIPTOR, descriptors::VELOCITY> meanPhysField(sLattice,
                                                                              /*scale=*/1.0, dataName + "Mean");
    vtmWriter.addFunctor(meanPhysField);

    SuperLatticePhysField<T, DESCRIPTOR, descriptors::VELOCITY2> stdPhysField(sLattice,
                                                                              /*scale=*/1.0, dataName + "Std");
    vtmWriter.addFunctor(stdPhysField);

    // Finally, write the iteration data
    singleton::directories().setOutputDir("./tmp/");
    vtmWriter.write(iT[0]);
  } // end for (iter)
}

} // namespace uq

} // namespace olb

#endif // POSTPROCESSING_H