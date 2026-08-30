/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_PHYS_VISCOSITY_2D_H
#define LATTICE_PHYS_VISCOSITY_2D_H

#include<vector>

#include "superBaseF2D.h"
#include "superCalcF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"

#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "indicator/blockIndicatorBaseF2D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"


/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor to get pointwise phys viscosity on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysViscosity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysViscosity2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                              const UnitConverter<T,DESCRIPTOR>& converter, bool logscale=false);
};

/// functor returns pointwise phys viscosity on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysViscosity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  const bool _logscale;
public:
  BlockLatticePhysViscosity2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                              const UnitConverter<T,DESCRIPTOR>& converter,
                              bool logscale);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
