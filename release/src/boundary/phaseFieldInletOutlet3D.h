/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

//This file contains all inlet and outlet boundary conditions of phase field models
//This boundary only contains free floating functions

#ifndef PHASE_FIELD_INLET_OUTLET_3D_H
#define PHASE_FIELD_INLET_OUTLET_3D_H

#include "setBoundary.h"
#include "phaseFieldInletOutlet.h"
#include "zouHeDynamics.h"
#include "dynamics/freeEnergyDynamics.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct FreeEnergyVelocity<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedVelocityBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedVelocityBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedVelocityBoundaryTuple
      >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct FreeEnergyPressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum3D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct FreeEnergyPressureConvective<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum3D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum3D>(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum3D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 3)
struct FreeEnergyOrderParameter<T,DESCRIPTOR> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationDynamics<T,DESCRIPTOR,FreeEnergyInletOutletDynamics>
    ::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletOrderParameter3D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 3)
struct FreeEnergyOrderParameterConvective<T,DESCRIPTOR> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationDynamics<T,DESCRIPTOR,FreeEnergyInletOutletDynamics>
    ::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletOrderParameter3D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

}

}

#endif
