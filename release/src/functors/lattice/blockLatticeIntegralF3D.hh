/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Mathias J. Krause,
 *  Albert Mink, Adrian Kummeränder
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

#ifndef BLOCK_LATTICE_INTEGRAL_F_3D_HH
#define BLOCK_LATTICE_INTEGRAL_F_3D_HH

#include <vector>
#include "utilities/omath.h"

#include "blockLatticeIntegralF3D.h"
#include "blockGeometryFaces3D.h"
#include "blockCalcF3D.h" // for IdentityF
#include "core/olbDebug.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
BlockLatticePhysDrag3D<T,DESCRIPTOR>::BlockLatticePhysDrag3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF3D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3),
    _indicatorF(indicatorF),
    _facesF(indicatorF, converter.getConversionFactorLength()),
    _pBoundForceF(blockLattice, indicatorF, converter),
    _sumF(_pBoundForceF, indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDrag3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T faces[7] = { };
  T sum[4]   = { };
  _sumF(sum, input);
  _facesF(faces, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];
  output[2] = _factor * sum[2] / faces[2];

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::BlockLatticePhysCorrDrag3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF3D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3),
    _indicatorF(indicatorF),
    _facesF(indicatorF, converter.getConversionFactorLength()),
    _pBoundForceF(blockLattice, indicatorF, converter),
    _sumF(_pBoundForceF, indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physCorrDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T faces[7] = { };
  T sum[4]   = { };
  _facesF(faces, input);
  _sumF(sum, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];
  output[2] = _factor * sum[2] / faces[2];

  return true;
}


} // end namespace olb

#endif
