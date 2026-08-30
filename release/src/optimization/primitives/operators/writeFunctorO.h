/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef WRITE_FUNCTOR_O_H
#define WRITE_FUNCTOR_O_H

#include <core/operator.h>
#include "../concept.h"
#include "../superLatticeO.h"

namespace olb {

namespace operators {

/// Operator which evaluate a functor and write its results into the passed field
template <typename FUNCTOR, typename TO>
struct WriteFunctorO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = FUNCTOR::parameters;

  int getPriority() const {
   return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    cell.template setField<TO>(FUNCTOR{}.compute(cell, parameters));
  }
};

}

// Returns a unique pointer to the SuperLatticeO with the writing operator
template <typename FUNCTOR, typename TO, typename T, typename DESCRIPTOR>
auto makeWriteFunctorO(SuperLattice<T,DESCRIPTOR>& lattice) {
  return makeSuperLatticeO<operators::WriteFunctorO<FUNCTOR,TO>>(lattice);
}

template <typename FUNCTOR, typename TO, typename T, typename DESCRIPTOR>
void writeFunctorTo(SuperLattice<T,DESCRIPTOR>& lattice) {
  auto superLatticeO = makeSuperLatticeO<operators::WriteFunctorO<FUNCTOR,TO>>(lattice);
  superLatticeO->execute();
}

template <typename FUNCTOR, typename TO, typename T, typename DESCRIPTOR>
void writeFunctorTo(SuperLattice<T,DESCRIPTOR>& lattice,
                           FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator) {
  auto superLatticeO = makeSuperLatticeO<operators::WriteFunctorO<FUNCTOR,TO>>(lattice);
  superLatticeO->restrictTo(std::forward<FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>>(indicator));
  superLatticeO->execute();
}

template <typename FUNCTOR, typename TO, typename T, typename DESCRIPTOR>
void writePhysFunctorTo(SuperLattice<T,DESCRIPTOR>& lattice,
                               FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                               T conversionFactor = 1.0) {
  auto superLatticeO = makeSuperLatticeO<operators::WriteFunctorO<FUNCTOR,TO>>(lattice);
  superLatticeO->template setParameter<descriptors::CONVERSION>(conversionFactor);
  superLatticeO->restrictTo(std::forward<FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>>(indicator));
  superLatticeO->execute();
}

}
#endif
