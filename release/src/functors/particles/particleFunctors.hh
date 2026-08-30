/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause,
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

#ifndef PARTICLE_FUNCTORS_HH
#define PARTICLE_FUNCTORS_HH

namespace olb {


template <typename T, typename DESCRIPTOR>
ParticleCircumRadiusF<T,DESCRIPTOR>::ParticleCircumRadiusF(
  Container<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fields_t>>& container )
  : ContainerF<T,DESCRIPTOR,DynamicFieldGroupsD<T,typename DESCRIPTOR::fields_t>,T>( container, 1 )
{
  this->getName() = "ParticleCircumRadiusF";
}

template <typename T, typename DESCRIPTOR>
bool ParticleCircumRadiusF<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  using namespace olb::descriptors;
  using GROUP_EVAL = typename DESCRIPTOR::template derivedField<SURFACE>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<SINDICATOR>;
  auto sIndicatorPtr = this->getContainer().data().template get<GROUP_EVAL>()
                  .template getField<FIELD_EVAL>(input[0]);
  T circumRadius = sIndicatorPtr->getCircumRadius();
  output[0] = circumRadius;
  return true; //TODO: add Error handling using true/false
}


} //namespace olb

#endif
