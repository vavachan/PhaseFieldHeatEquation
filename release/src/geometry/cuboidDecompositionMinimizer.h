/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
 *                2022 Adrian Kummerlaender
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

#ifndef CUBOID_GEOMETRY_MINIMIZER_H
#define CUBOID_GEOMETRY_MINIMIZER_H

#include "cuboid.h"

namespace olb {

template <typename T, unsigned D> class CuboidDecomposition;
template <typename T> class IndicatorF3D;

/// Splits largest cuboid by-volume until there are nC cuboids
template <typename T>
void continueMinimizeByVolume(CuboidDecomposition<T,3>& cGeometry, IndicatorF3D<T>& indicatorF, int nC)
{
  cGeometry.shrink(indicatorF);

  while (cGeometry.size() < nC) {
    // Search for the largest child cuboid
    T iCVolume[cGeometry.size()];
    int maxiC = 0;
    for (int iC = 0; iC < cGeometry.size(); iC++) {
      iCVolume[iC] = cGeometry.get(iC).getLatticeVolume();
      if ( iCVolume[iC] > iCVolume[maxiC] ) {
        maxiC = iC;
      }
    }

    // looking for largest extend, because halfing the cuboid by its largest extend will result in the smallest surface and therfore in the least comminication cells
    auto& largest = cGeometry.get(maxiC);
    if (largest.getNx() >= largest.getNy() && largest.getNx() >= largest.getNz()) {
      // clout << "Cut in x direction!" << std::endl;
      largest.divide({2,1,1}, cGeometry.cuboids());
    }
    else if (largest.getNy() >= largest.getNx() && largest.getNy() >= largest.getNz()) {
      // clout << "Cut in y direction!" << std::endl;
      largest.divide({1,2,1}, cGeometry.cuboids());
    }
    else {
      // clout << "Cut in z direction!" << std::endl;
      largest.divide({1,1,2}, cGeometry.cuboids());
    }
    cGeometry.remove(maxiC);
    // shrink the two new cuboids
    cGeometry.shrink(cGeometry.cuboids().size()-2, indicatorF);
    cGeometry.shrink(cGeometry.cuboids().size()-1, indicatorF);
  }
}

/// Splits into nC cuboids by-volume
template <typename T>
void minimizeByVolume(CuboidDecomposition<T,3>& cGeometry, IndicatorF3D<T>& indicatorF, int nC)
{
  // Search for the largest multiplier not dividable by two
  int initalNc = nC;
  while ( initalNc % 2 == 0 ) {
    initalNc /= 2;
  }

  // Split evenly in initalNc many cuboids and shrink all
  cGeometry.split(0, initalNc);
  cGeometry.shrink(indicatorF);

  continueMinimizeByVolume(cGeometry, indicatorF, nC);
}

template <typename T>
void minimizeByWeight(CuboidDecomposition<T,3>& cGeometry, IndicatorF3D<T>& indicatorF, int nC)
{
  cGeometry.setWeights(indicatorF);

  // conduct a prime factorisation for the number of cuboids nC
  std::vector<int> factors;
  int initalNc = nC;
  // iterate over the prime numbes from 0 to 100 (may have to be extended)
  for (int i : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71}) {
    while (initalNc % i == 0) {
      initalNc /= i;
      factors.push_back(i);
    }
  }

  // recursively split the current cuboids by each prime factor
  for (int i = factors.size() - 1; i >= 0; i--) {
    int currentNc = cGeometry.cuboids().size();
    for (int iC = 0; iC < currentNc; iC++) {
      // clout << "split cuboid number #" << iC << " in " << factors[i] << " parts" << std::endl;
      cGeometry.splitByWeight(iC, factors[i], indicatorF);
      cGeometry.shrink(indicatorF);
    }
    for (int iC = 0; iC < currentNc; iC++) {
      // clout << "delet cuboid number #" << iC << std::endl;
      cGeometry.remove(0);
    }
  }
}

}

#endif
