/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_H

#include "utilities/omath.h"

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

namespace olb {

//======================================================================
// ======== Phase field coupling without bouancy 3D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class PhaseFieldCouplingPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PhaseFieldCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                    T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness,
                                    std::vector<BlockStructureD<3>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  using L = DESCRIPTOR;
  using PHI_CACHE = descriptors::FIELD_BASE<1,0,0>;

  int x0, x1, y0, y1, z0, z1;

  T _rho_L, _rho_H, _delta_rho;
  T _mu_L, _mu_H;
  T _surface_tension, _interface_thickness;
  T _beta, _kappa;

  BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::INTERPHASE_NORMAL>> *tPartner;

  std::vector<BlockStructureD<3>*> partners;
};

template<typename T, typename DESCRIPTOR>
class PhaseFieldCouplingGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  PhaseFieldCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<3>* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T _rho_L, _rho_H, _delta_rho;
  T _mu_L, _mu_H;
  T _surface_tension, _interface_thickness;
};


//======================================================================
// =============  SmagorinskyBoussinesqCoupling 3D ===================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_,
      std::vector<BlockStructureD<3>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>> *tPartner;
  T forcePrefactor[L::d];
  T tauTurbADPrefactor;
  T smagoPrefactor;

  std::vector<BlockStructureD<3>*> partners;
};

template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<3>* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  T smagoPrefactor;
};

//==================================================================================================
// ========Coupling 3D of Navier-Stokes on Advection-Diffusion====================//
//==================================================================================================

template<typename FORCE = void>
struct AdvectionDiffusionParticleCoupling3D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  struct LATTICE_TIME : public descriptors::FIELD_BASE<1> { };
  using parameters = meta::merge<meta::list<LATTICE_TIME>, typename FORCE::parameters>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELL::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELL::template value_t<names::NavierStokes>::descriptor_t;

    auto& cellNS = cells.template get<names::NavierStokes>();
    auto& cellAD = cells.template get<names::AdvectionDiffusion>();
    int iT = parameters.template get<LATTICE_TIME>();
    bool par = !(iT%2);

    auto vel = par ? cellAD.template getField<descriptors::VELOCITY2>() : cellAD.template getField<descriptors::VELOCITY>();
    V velocity[DESCRIPTOR::d] {};
    for (int i=0; i < DESCRIPTOR::d; i++) velocity[i] = vel[i];
    V velGrad[DESCRIPTOR::d] {};
    V forceValue[DESCRIPTOR::d] {};
    V velF[DESCRIPTOR::d] {};
    if constexpr (!std::is_void_v<FORCE>) {
      // calculating upwind Gradient
      // vel contains velocity information on ADlattice
      // velGrad contains upwind vel on ADlattice
      auto velXp = par ? cellAD.neighbor({1,0,0}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({1,0,0}).template getField<descriptors::VELOCITY>();
      auto velXn = par ? cellAD.neighbor({-1,0,0}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({-1,0,0}).template getField<descriptors::VELOCITY>();
      auto velYp = par ? cellAD.neighbor({0,1,0}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({0,1,0}).template getField<descriptors::VELOCITY>();
      auto velYn = par ? cellAD.neighbor({0,-1,0}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({0,-1,0}).template getField<descriptors::VELOCITY>();
      auto velZp = par ? cellAD.neighbor({0,0,1}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({0,0,1}).template getField<descriptors::VELOCITY>();
      auto velZn = par ? cellAD.neighbor({0,0,-1}).template getField<descriptors::VELOCITY2>() : cellAD.neighbor({0,0,-1}).template getField<descriptors::VELOCITY>();

      if (vel[0]<0.) {
        velGrad[0] = vel[0]*(velXp[0]-vel[0]);
        velGrad[1] = vel[0]*(velXp[1]-vel[1]);
        velGrad[2] = vel[0]*(velXp[2]-vel[2]);
      }
      else {
        velGrad[0] = vel[0]*(vel[0]-velXn[0]);
        velGrad[1] = vel[0]*(vel[1]-velXn[1]);
        velGrad[2] = vel[0]*(vel[2]-velXn[2]);
      }
      if (vel[1]<0.) {
        velGrad[0] += vel[1]*(velYp[0]-vel[0]);
        velGrad[1] += vel[1]*(velYp[1]-vel[1]);
        velGrad[2] += vel[1]*(velYp[2]-vel[2]);
      }
      else {
        velGrad[0] += vel[1]*(vel[0]-velYn[0]);
        velGrad[1] += vel[1]*(vel[1]-velYn[1]);
        velGrad[2] += vel[1]*(vel[2]-velYn[2]);
      }
      if (vel[2]<0.) {
        velGrad[0] += vel[2]*(velZp[0]-vel[0]);
        velGrad[1] += vel[2]*(velZp[1]-vel[1]);
        velGrad[2] += vel[2]*(velZp[2]-vel[2]);
      }
      else {
        velGrad[0] += vel[2]*(vel[0]-velZn[0]);
        velGrad[1] += vel[2]*(vel[1]-velZn[1]);
        velGrad[2] += vel[2]*(vel[2]-velZn[2]);
      }

      // writes force in forceValues, vel refers to ADlattice
      FORCE().applyForce(forceValue,cells,velocity,parameters);
      if (par) {
        cellAD.template setField<descriptors::VELOCITY2>(vel);
      }
      else {
        cellAD.template setField<descriptors::VELOCITY>(vel);
      }

      // compute new particle velocity
      Vector<V,DESCRIPTOR::d> newVel;
      for (int i=0; i < DESCRIPTOR::d; i++) {
        newVel[i] = vel[i] + forceValue[i] - velGrad[i];
      }

      if (par) {
        cellAD.template setField<descriptors::VELOCITY>(newVel);
      } else {
        cellAD.template setField<descriptors::VELOCITY2>(newVel);
      }
    }
    else {   // set particle velocity to fluid velocity
      //std::cout << "without force" << std::endl;
      cellNS.computeU(velF);
      Vector<V,DESCRIPTOR::d> newVel;
      for (int i = 0; i < DESCRIPTOR::d; i++) {
        newVel[i] = velF[i];
      }
      if (par) {
        cellAD.template setField<descriptors::VELOCITY>(newVel);
      } else {
        cellAD.template setField<descriptors::VELOCITY2>(newVel);
      }
    }
  }
};

//======================================================================
// ======== Porous Regularized NSDiffusion Coupling 3D =================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<3>* > partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;

  std::vector<BlockStructureD<3>*> partners;
};

template<typename T, typename DESCRIPTOR>
class PorousNavierStokesAdvectionDiffusionCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  PorousNavierStokesAdvectionDiffusionCouplingGenerator3D(
    int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<BlockStructureD<3>* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};

//=====================================================================================
//==============  MixedScaleBoussinesqCouplingPostProcessor3D ===============
//=====================================================================================
template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
      std::vector<BlockStructureD<3>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  typedef DESCRIPTOR L;
  using HEAT_FLUX_CACHE = descriptors::FIELD_BASE<1, 0, 0>;
  int x0, x1, y0, y1, z0, z1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *tPartner;
  Vector<T,L::d> forcePrefactor;
  T tauTurbADPrefactor;

  std::vector<BlockStructureD<3>*> partners;
};

template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                          T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<3>* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
};

//======================================================================
// ========VANS-ADE Particle Coupling 3D ====================//
//======================================================================
template<
  typename T,
  typename DESCRIPTOR,
  typename POROSITY,
  typename ADLattice,
  typename FIELD_A,
  typename FIELD_B
  >
class VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D :
  public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D(
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
    std::vector<BlockStructureD<3>* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR,ADLattice> > > forces_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_,  int z0_, int z1_) override;

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > _forces;
private:
  int x0, x1, y0, y1, z0, z1, iC;
  BlockLattice<T,ADLattice>* _partnerLattice;

  Cell<T,ADLattice> _cell;
  Cell<T,ADLattice> _cellXp;
  Cell<T,ADLattice> _cellXn;
  Cell<T,ADLattice> _cellYp;
  Cell<T,ADLattice> _cellYn;
  Cell<T,ADLattice> _cellZp;
  Cell<T,ADLattice> _cellZn;

  bool par = true;
};

template<
  typename T,
  typename DESCRIPTOR,
  typename POROSITY,
  typename ADLattice,
  typename FIELD_A,
  typename FIELD_B
  >
class VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D :
  public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D();
  PostProcessor3D<T,DESCRIPTOR>* generate(
    std::vector<BlockStructureD<3>* > partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
  void addForce(AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force);

protected:
  std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice> > > ADforces;
  };

}

#endif
