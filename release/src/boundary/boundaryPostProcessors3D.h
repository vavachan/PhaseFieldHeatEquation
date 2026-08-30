/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef BOUNDARY_POST_PROCESSORS_3D_H
#define BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"

#include "core/operator.h"

namespace olb {

/**
* This class computes the skordos BC
* on a plane wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class PlaneFdBoundaryProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

private:
  template <int deriveDirection, typename CELL, typename V=CELL::value_t>
  void interpolateGradients(CELL& cell, V velDeriv[DESCRIPTOR::d]) const any_platform;

};


template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  struct PREV_CELL : public descriptors::FIELD_BASE<
    util::populationsContributingToVelocity<DESCRIPTOR,direction,-orientation>().size()
  > { };

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void initialize(CELL& cell) any_platform;

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

};


/**
* This class computes the skordos BC
* on a convex edge wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template <typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

private:
  template <typename CELL, typename V=CELL::value_t>
  V getNeighborRho(CELL& cell, int step1, int step2) any_platform;

  template <int deriveDirection, int orientation, typename CELL, typename V=CELL::value_t>
  void interpolateGradients(CELL& cell, V velDeriv[DESCRIPTOR::d]) const any_platform;

};


template<typename T, typename DESCRIPTOR,
         int xNormal, int yNormal, int zNormal>
struct OuterVelocityCornerProcessor3D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes a slip BC in 3D
*/

template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
};


template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};

/**
* This class computes a partial slip BC in 3D
*/

template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessor3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
  T tuner;
};


template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessorGenerator3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
  T tuner;
};


/// PostProcessors for the chemical potential boundary condition in the free energy model.
/// The chemical potentials on the boundary are set equal to the chemical potential on the
/// fluid cell normal to the boundary. This is necessary because the coupling between the
/// lattices requires the calculation of the gradient of the chemical potential.
///
/// It would be preferable if these were implemented as a lattice coupling that ran
/// between the chemical potential and force lattice couplings. However there is no
/// access to the discrete normals in lattice couplings.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
struct FreeEnergyInletMomentum3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

    T rhoBoundary = cell.computeRho();
    T rhoBulk = cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).computeRho();

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundary - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
struct FreeEnergyInletOrderParameter3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
struct FreeEnergyOutletMomentum3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[3];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T, 3> normalVec({NORMAL_X,NORMAL_Y,NORMAL_Z});

    if (normalVec[2] == 0) {
      if (normalVec[1] == 0) {
        if (normalVec[0] < 0) {
          uPerp = -u[0];
        } else {
          uPerp = u[0];
        }
      } else if (normalVec[0] == 0) {
        if (normalVec[1] < 0) {
          uPerp = -u[1];
        } else {
          uPerp = u[1];
        }
      } else {
        uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1]);
      }
    } else if (normalVec[1] == 0) {
      if (normalVec[0] == 0) {
        if (normalVec[2] < 0) {
          uPerp = -u[2];
        } else {
          uPerp = u[2];
        }
      } else {
        uPerp = util::sqrt(u[0] * u[0] + u[2] * u[2]);
      }
    } else if (normalVec[0] == 0) {
      uPerp = util::sqrt(u[1] * u[1] + u[2] * u[2]);
    } else {
      uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundaryNew - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
struct FreeEnergyOutletOrderParameter3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[3];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T, 3> normalVec({NORMAL_X,NORMAL_Y,NORMAL_Z});

    if (normalVec[2] == 0) {
      if (normalVec[1] == 0) {
        if (normalVec[0] < 0) {
          uPerp = -u[0];
        } else {
          uPerp = u[0];
        }
      } else if (normalVec[0] == 0) {
        if (normalVec[1] < 0) {
          uPerp = -u[1];
        } else {
          uPerp = u[1];
        }
      } else {
        uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1]);
      }
    } else if (normalVec[1] == 0) {
      if (normalVec[0] == 0) {
        if (normalVec[2] < 0) {
          uPerp = -u[2];
        } else {
          uPerp = u[2];
        }
      } else {
        uPerp = util::sqrt(u[0] * u[0] + u[2] * u[2]);
      }
    } else if (normalVec[0] == 0) {
      uPerp = util::sqrt(u[1] * u[1] + u[2] * u[2]);
    } else {
      uPerp = util::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyChemPotBoundaryProcessor3DA {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyChemPotBoundaryProcessor3DB {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

};


/// PostProcessor for the wetting boundary condition in the free energy model. This is
/// required to set rho on the boundary (using the denisty of the neighbouring cell in
/// direction of inwards facing normal at the boundary), as the coupling between the
/// lattices requires the calculation of a density gradient.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyWallProcessor3D {
public:
  using parameters = meta::list<olb::descriptors::ADDEND>;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform;

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyWallMomentumProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X, -NORMAL_Y, -NORMAL_Z}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundary - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyWallOrderParameterProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X, -NORMAL_Y, -NORMAL_Z}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y,-NORMAL_Z}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

/// PostProcessor for the density / velocity outflow boundaries in the free energy model.
/// The density / order parameters are prescribed to the outflow nodes such that they obey
/// the local-velocity convective boundary condition given in Lou, Gou, Shi (2013).
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y, int NORMAL_Z>
class FreeEnergyConvectiveProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};


}

#endif
