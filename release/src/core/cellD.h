/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef CELL_D_H
#define CELL_D_H

#include <type_traits>

#include "fieldArrayD.hh"

namespace olb {

/// Minimal cell storing only population data
template<typename T, typename DESCRIPTOR>
class PopulationCellD {
private:
  FieldD<T,DESCRIPTOR,descriptors::POPULATION> _data;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  template <typename POPULATIONS>
  PopulationCellD(POPULATIONS&& pops) any_platform:
    _data{pops} { }

  const T& operator[](unsigned iPop) const any_platform {
    return _data[iPop];
  }

  T& operator[](unsigned iPop) any_platform {
    return _data[iPop];
  }
};

/// Single cell implementing the full field data interface
template<typename T, typename DESCRIPTOR>
class CellD {
private:
  MultiFieldArrayForDescriptorD<T,DESCRIPTOR,Platform::CPU_SISD> _fieldsD;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  CellD(): _fieldsD(1) { }

  const T& operator[](unsigned iPop) const {
    return _fieldsD.template getFieldComponent<descriptors::POPULATION>(0, iPop);
  }

  T& operator[](unsigned iPop) {
    return _fieldsD.template getFieldComponent<descriptors::POPULATION>(0, iPop);
  }

  template <typename FIELD>
  auto getField() const {
    return _fieldsD.template getField<FIELD>(0);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    _fieldsD.template setField<FIELD>(0, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _fieldsD.template getFieldPointer<FIELD>(0);
  }

  template <typename FIELD>
  auto getFieldPointer() const {
    return _fieldsD.template getFieldPointer<FIELD>(0);
  }

  template <typename FIELD>
  const auto& getFieldComponent(unsigned iD) const {
    return _fieldsD.template getFieldComponent<FIELD>(0, iD);
  }

  template <typename FIELD>
  auto& getFieldComponent(unsigned iD) {
    return _fieldsD.template getFieldComponent<FIELD>(0, iD);
  }

};

/// Single cell with full field data and dynamics interface
template<typename T, typename DESCRIPTOR, typename DYNAMICS>
class FullCellD : public CellD<T,DESCRIPTOR> {
public:
  FullCellD() : CellD<T,DESCRIPTOR>()
  { }

  T computeRho() const {
    return typename DYNAMICS::MomentaF().computeRho(*this);
  }
  void computeU(T* u) {
    typename DYNAMICS::MomentaF().computeU(*this, u);
  }
  void computeJ(T* j) {
    typename DYNAMICS::MomentaF().computeJ(*this, j);
  }
  void computeRhoU(T& rho, T* u) {
    typename DYNAMICS::MomentaF().computeRhoU(*this, rho, u);
  }
  void computeStress(T* pi) {
    T rho, u[DESCRIPTOR::d] { };
    typename DYNAMICS::MomentaF().computeRhoU(*this, rho, u);
    typename DYNAMICS::MomentaF().computeStress(*this, rho, u, pi);
  }
  void computeAllMomenta(T& rho, T* u, T* pi) {
    typename DYNAMICS::MomentaF().computeAllMomenta(*this, rho, u, pi);
  }

  void defineRho(T& rho) {
    typename DYNAMICS::MomentaF().defineRho(*this, rho);
  }
  void defineU(T* u) {
    typename DYNAMICS::MomentaF().defineU(*this, u);
  }
  void defineRhoU(T rho, T* u) {
    typename DYNAMICS::MomentaF().defineRhoU(*this, rho, u);
  }
  void defineAllMomenta(T rho, T* u, T* pi) {
    typename DYNAMICS::MomentaF().defineAllMomenta(*this, rho, u, pi);
  }

  void iniEquilibrium(T rho, T* u) {
    typename DYNAMICS::MomentaF().iniEquilibrium(*this, rho, u);
  }
  void iniRegularized(T rho, T* u, T* pi) {
    typename DYNAMICS::MomentaF().iniRegularized(*this, rho, u, pi);
  }
  void inverseShiftRhoU(T& rho, T* u) {
    typename DYNAMICS::MomentaF().inverseShiftRhoU(*this, rho, u);
  }

};

}

#endif
