/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * Interface for post-processing steps -- header file.
 */
#ifndef LATTICE_STATISTICS_H
#define LATTICE_STATISTICS_H

#include <vector>
#include "io/ostreamManager.h"

namespace olb {

/////////////////// Statistics Postprocessing ////////////////////////////

namespace statistics {

struct AVERAGE_RHO : public descriptors::FIELD_BASE<1> { };

}

template<typename T>
class LatticeStatistics {
public:
  enum { avRho=0, avEnergy=1 } AverageT;
  enum { maxU=0 } MaxT;

  struct Aggregatable {
    std::size_t nCells = 0;
    T avRho { };
    T avEnergy { };
    T maxU { };

    void increment(T rho, T uSqr) {
      nCells   += 1;
      avRho    += rho;
      avEnergy += uSqr;
      if constexpr (!std::is_same_v<Expr,T>) {
        maxU = std::max(uSqr, maxU);
      }
    }

    Aggregatable& operator+=(const Aggregatable& rhs) {
      nCells += rhs.nCells;
      avRho += rhs.avRho;
      avEnergy += rhs.avEnergy;
      maxU = std::max(maxU, rhs.maxU);
      return *this;
    }
  };

public:
  LatticeStatistics();
  ~LatticeStatistics() = default;
  void reset();
  void reset(T average_rho_, T average_energy_, T maxU_, std::size_t numCells_);

  int subscribeAverage();
  int subscribeSum();
  int subscribeMin();
  int subscribeMax();

  void incrementStats(T rho, T uSqr);
  void incrementStats(Aggregatable& aggregatable);
  void gatherAverage(int whichAverage, T value);
  void gatherSum(int whichSum, T value);
  void gatherMin(int whichMin, T value);
  void gatherMax(int whichMax, T value);
  void incrementStats();
  T getAverageRho() const;
  T getAverageEnergy() const;
  T getMaxU() const;
  std::size_t getNumCells() const;

  T getAverage(int whichAverage) const;
  T getSum(int whichSum) const;
  T getMin(int whichMin) const;
  T getMax(int whichMax) const;

  std::vector<T>& getAverageVect();
  std::vector<T>& getSumVect();
  std::vector<T>& getMinVect();
  std::vector<T>& getMaxVect();

  void incrementTime();
  void resetTime(size_t value=0);
  std::size_t getTime() const;
  void print(int iterationStep, T physicalTime=-1) const;
  void initialize();

private:
  mutable OstreamManager clout;
  // variables for internal computations
  std::vector<T> tmpAv, tmpSum, tmpMin, tmpMax;
  std::size_t tmpNumCells;
  // variables containing the public result
  std::vector<T> averageVect, sumVect, minVect, maxVect;
  std::size_t numCells;
  std::size_t latticeTime;
  bool firstCall;

};

}  // namespace olb

#endif
