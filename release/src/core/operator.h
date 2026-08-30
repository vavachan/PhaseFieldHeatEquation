/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_OPERATOR_H
#define CORE_OPERATOR_H

#include "operatorScope.h"

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;
template <typename T, Platform PLATFORM> class ConcreteBlockMask;
template <typename T, typename DESCRIPTOR> struct Dynamics;

/// Base of any block operator
struct AbstractBlockO {
  virtual ~AbstractBlockO() = default;

  virtual std::type_index id() const = 0;
};

/// Base of block-wide operators such as post processors
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct BlockO : public AbstractBlockO {
  /// Set whether iCell is covered by the operator (optional)
  virtual void set(CellID iCell, bool state) = 0;
  /// Setup operator context
  virtual void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;
  /// Apply operator on block
  virtual void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;

  /// Get number of cells covered by operator (optional)
  virtual std::size_t weight() const = 0;
};

/// Block application of concrete OPERATOR called using SCOPE on PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename OPERATOR, OperatorScope SCOPE>
class ConcreteBlockO;

/// Base of collision operations performed by BlockDynamicsMap
template <typename T, typename DESCRIPTOR>
struct AbstractCollisionO : public AbstractBlockO {
  virtual Dynamics<T,DESCRIPTOR>* getDynamics() = 0;

  /// Returns number of assigned cells
  /**
   * Used to determine the dominant dynamics to choose e.g. which
   * collision operator to vectorize or to prefer in GPU kernels.
   **/
  virtual std::size_t weight() const = 0;
  /// Set whether iCell is covered by the present collision step
  /**
   * \param iCell   Cell index
   * \param state   (De)activate for this dynamics / collision
   * \param overlap Cell index in overlap (set dynamics but do not collide)
   **/
  virtual void set(CellID iCell, bool state, bool overlap) = 0;
};

/// Collision dispatch strategy
enum struct CollisionDispatchStrategy {
  /// Apply dominant dynamics using mask and fallback to virtual dispatch for others
  Dominant,
  /// Apply all dynamics individually (async for Platform::GPU_CUDA)
  Individual
};

/// Collision operation on concrete blocks of PLATFORM
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct BlockCollisionO : public AbstractCollisionO<T,DESCRIPTOR> {
  /// Setup collision on block
  virtual void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;
  /// Apply collision on subdomain of block using strategy
  /**
   * Subdomain is currently assumed to be the non-overlap area of the block
   **/
  virtual void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block,
                     ConcreteBlockMask<T,PLATFORM>&               subdomain,
                     CollisionDispatchStrategy                    strategy) = 0;
};

/// Collision operation of concrete DYNAMICS on concrete block lattices of PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
class ConcreteBlockCollisionO;

/// Base of block-wide coupling operators executed by SuperLatticeCoupling
template <typename COUPLEES>
struct AbstractCouplingO : public AbstractBlockO {
  /// Value type used for coupling parameters
  using value_t      = typename COUPLEES::values_t::template get<0>::value_t;
  /// Descriptor type used for coupling parameters
  using descriptor_t = typename COUPLEES::values_t::template get<0>::descriptor_t;

  using AbstractParameters = olb::AbstractParameters<value_t,descriptor_t>;
  using ParametersD = olb::ParametersD<value_t,descriptor_t>;

  using LatticeR = olb::LatticeR<descriptor_t::d>;
  template <typename FIELD>
  using FieldD = olb::FieldD<value_t,descriptor_t,FIELD>;

  /// Execute coupling operation
  virtual void execute() = 0;
  /// Return reference to parameters of coupling operator
  virtual AbstractParameters& getParameters() = 0;
  /// Set whether iCell is covered by the present coupling
  /**
   * By default the entire non-overlap block area is coupled.
   **/
  virtual void set(CellID iCell, bool state) = 0;

};

/// Coupling of COUPLEES using concrete OPERATOR with SCOPE on PLATFORM lattices
template<typename COUPLEES, Platform PLATFORM, typename OPERATOR, OperatorScope SCOPE>
class ConcreteBlockCouplingO;

/// Particle coupling of COUPLEES using concrete OPERATOR with SCOPE on PLATFORM lattices
template<typename COUPLEES, Platform PLATFORM, typename OPERATOR, OperatorScope SCOPE>
class ConcreteBlockPointCouplingO;

}

#endif
