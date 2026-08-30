/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Robin Trunk
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

#ifndef ADVECTION_DIFFUSION_FORCES_H
#define ADVECTION_DIFFUSION_FORCES_H

#include "core/unitConverter.h"

namespace olb {

/**
 * Refactored AdvectionDiffusion forces to fit
 * with new template style
 */
namespace ade_forces {

struct AdvDiffDragForce3D {

  struct INIT_ARG : public descriptors::FIELD_BASE<1> { };
  struct DRAG_COEFF : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<INIT_ARG,DRAG_COEFF>;

  template <typename V,typename PARAMETERS>
  V computeDragFromSt(V ST, V U, V dT, V dL, PARAMETERS& params) any_platform {
    V dragCoeff = (U*dT) / (ST*dL);
    params.template set<INIT_ARG>(8.);
    params.template set<DRAG_COEFF>(dragCoeff);
  }

  template <typename V, typename COUPLING, typename CONVERTER>
  static void computeParametersFromRhoAndRadius(V pRho, V pRadius, COUPLING& coupling, CONVERTER& converter) any_platform {
    V dragCoeff = (9.*converter.getPhysViscosity()*converter.getPhysDensity()*converter.getConversionFactorTime())
              / (2.*pRho*pRadius*pRadius);
    coupling.template setParameter<ade_forces::AdvDiffDragForce3D::INIT_ARG>(V(8.));
    coupling.template setParameter<ade_forces::AdvDiffDragForce3D::DRAG_COEFF>(dragCoeff);
  }

  template <typename V,typename PARAMETERS>
  V computeDragFromRhoAndRadius(V pRho, V pRadius, V dGamma, V dRho, V dT, PARAMETERS& params) any_platform {
    V dragCoeff = (9.*dGamma*dRho*dT) / (2.*pRho*pRadius*pRadius);
    params.template set<INIT_ARG>(8.);
    params.template set<DRAG_COEFF>(dragCoeff);
  }

  template <typename V, typename CELLS, typename PARAMETERS>
  static void applyForce(V force[], CELLS& cells, V velocity[], PARAMETERS& params) any_platform {
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    auto& cellNS = cells.template get<names::NavierStokes>();
    V dragCoeff = params.template get<DRAG_COEFF>();
    V velF[3] = {0.,0.,0.};
    cellNS.computeU(velF);
    for (int i=0; i < DESCRIPTOR::d; i++) {
      force[i] += dragCoeff*(velF[i]-velocity[i]);
    }
  }

};

}

template<typename T, typename DESCRIPTOR,
         typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvectionDiffusionForce3D {
public:
  AdvectionDiffusionForce3D()
  {
    initArg = 0;
  };
  virtual ~AdvectionDiffusionForce3D() {};
  virtual void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])=0;
  int getInitArg()
  {
    return initArg;
  }
private:
  int initArg;
};

template<typename T, typename DESCRIPTOR,
         typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffDragForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T St_);
  AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_);
  ~AdvDiffDragForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  int initArg;
  T dragCoeff;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffSNDragForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffSNDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_);
  ~AdvDiffSNDragForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  int initArg;
  T dragCoeff;
  T Re_pCoeff;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffBuoyancyForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffBuoyancyForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, Vector<T,3> g, T pRho_);
  ~AdvDiffBuoyancyForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  int initArg;
  T densDiff;
  Vector<T,3> gravity;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffRotatingForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffRotatingForce3D(SuperGeometry<T,3>& superGeometry_,
                         const UnitConverter<T,DESCRIPTOR>& converter_,
                         std::vector<T> axisPoint_,
                         std::vector<T> axisDirection_,
                         T w_, T* frac_,
                         bool centrifugeForceOn_ = true,
                         bool coriolisForceOn_ = true);
  AdvDiffRotatingForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_);
  virtual ~AdvDiffRotatingForce3D() {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]);

protected:
  SuperGeometry<T,3>& sg;
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T invMassLessForce;
  T w;
  T* frac;
  bool centrifugeForceOn;
  bool coriolisForceOn;

};

template<typename T, typename DESCRIPTOR,
         typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffMagneticWireForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffMagneticWireForce3D(SuperGeometry<T,3>& superGeometry_, UnitConverter<T,DESCRIPTOR> const& converter_, T pMass, AnalyticalF<3,T, T>& getMagForce);
  ~AdvDiffMagneticWireForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  SuperGeometry<T,3>& sg;
  int initArg;
  T _pMass;
  T _conversionVelocity;
  AnalyticalF<3,T, T>& _getMagForce;
};

}

#endif
