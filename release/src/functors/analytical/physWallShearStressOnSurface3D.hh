/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerländer, Christoph Gaul, Mathias J. Krause
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

#ifndef PHYS_WALL_SHEAR_STRESS_ON_SURFACE_3D_HH
#define PHYS_WALL_SHEAR_STRESS_ON_SURFACE_3D_HH

#include<cmath>

#include "functors/analytical/physWallShearStressOnSurface3D.h"
#include "core/vector.h"

namespace olb {
template<typename T, typename DESCRIPTOR>
PhysWallShearStressOnSurface3D<T,DESCRIPTOR>::PhysWallShearStressOnSurface3D(
  const UnitConverter<T,DESCRIPTOR>& converter,
  AnalyticalF3D<T,T>& densityF,
  AnalyticalF3D<T,T>& stressF,
  STLreader<T>& stlReader):
AnalyticalF<3,T,T>(3),
_converter(converter),
_densityF(densityF),
_stressF(stressF),
_stlReader(stlReader)
{
const T omega = 1. / _converter.getLatticeRelaxationTime();
const T dt = _converter.getConversionFactorTime();
_physFactor = -omega
* descriptors::invCs2<T,DESCRIPTOR>() / dt
* _converter.getPhysDensity() * _converter.getPhysViscosity();
}

template<typename T, typename DESCRIPTOR>
bool PhysWallShearStressOnSurface3D<T,DESCRIPTOR>::operator() (T output[], const T physR[]){
  Vector<T,3> origin(physR);
  auto normal = _stlReader.surfaceNormal(physR);
  normal = normalize(normal);

  T traction[3] { };
  T stress[6] { };
  T rho{};

  Vector<T,3> neighbor = origin - 0.5*_converter.getPhysDeltaX() * normal;
  _densityF(&rho, neighbor.data());
  _stressF(stress, neighbor.data());

  traction[0] = stress[0]/_physFactor*rho*normal[0] +
                stress[1]/_physFactor*rho*normal[1] +
                stress[2]/_physFactor*rho*normal[2];
  traction[1] = stress[1]/_physFactor*rho*normal[0] +
                stress[3]/_physFactor*rho*normal[1] +
                stress[4]/_physFactor*rho*normal[2];
  traction[2] = stress[2]/_physFactor*rho*normal[0] +
                stress[4]/_physFactor*rho*normal[1] +
                stress[5]/_physFactor*rho*normal[2];

  T traction_normal_SP;
  T tractionNormalComponent[3];
  // scalar product of traction and normal vector
  traction_normal_SP = traction[0] * normal[0] +
                       traction[1] * normal[1] +
                       traction[2] * normal[2];
  tractionNormalComponent[0] = traction_normal_SP * normal[0];
  tractionNormalComponent[1] = traction_normal_SP * normal[1];
  tractionNormalComponent[2] = traction_normal_SP * normal[2];

  output[0] = traction[0] - tractionNormalComponent[0];
  output[1] = traction[1] - tractionNormalComponent[1];
  output[2] = traction[2] - tractionNormalComponent[2];

  return true;
}

} // namespace olb
#endif
