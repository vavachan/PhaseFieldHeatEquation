/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Mathias J. Krause, Geng Liu
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

/** \file
 * The description of optimization algorthims -- header file.
 */
/** \file
 * A collection of dynamics classes for dual LB methods
 * (e.g. dual MRT) with which a Cell object can be
 * instantiated -- header file.
 */

#ifndef DUAL_MRT_DYNAMICS_H
#define DUAL_MRT_DYNAMICS_H

//#include "optiStructure.h"
//#include "solver3D.h"



// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {


template<typename T> class Controller;
template<typename T, typename DESCRIPTOR> class DualController;




/// Implementation of the dual MRT collision step with external force
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class DualForcedMRTdynamics : public MRTdynamics<T,DESCRIPTOR,MOMENTA> {
public:
  /// Constructor
  DualForcedMRTdynamics(T omega_);

  /// Collision step
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

  virtual void defineRho(Cell<T,DESCRIPTOR>& cell, T rho);
private:

  T Mt_S_M[DESCRIPTOR::q][DESCRIPTOR::q];
  T invM_invMt[DESCRIPTOR::q][DESCRIPTOR::q];
  T Mt_S_invMt[DESCRIPTOR::q][DESCRIPTOR::q];
};

////////////////////// Class ForcedMRTdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
DualForcedMRTdynamics<T,DESCRIPTOR,MOMENTA>::DualForcedMRTdynamics (
  T omega_)
  : MRTdynamics<T,DESCRIPTOR,MOMENTA>(omega_)
{
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = DESCRIPTOR::S[iPop];
  }
  for (int iPop  = 0; iPop < DESCRIPTOR::shearIndexes; ++iPop) {
    rt[DESCRIPTOR::shearViscIndexes[iPop]] = omega_;
  }

  T tmp[DESCRIPTOR::q][DESCRIPTOR::q];

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      Mt_S_M[iPop][jPop] = T();
      Mt_S_invMt[iPop][jPop] = T();
      invM_invMt[iPop][jPop] = T();
      tmp[iPop][jPop] = T();
    }
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      for (int kPop=0; kPop < DESCRIPTOR::q; ++kPop) {
        if (jPop==kPop) {
          tmp[iPop][jPop] += DESCRIPTOR::M[kPop][iPop]*rt[kPop];
        }
      }
    }
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      for (int kPop=0; kPop < DESCRIPTOR::q; ++kPop) {
        Mt_S_M[iPop][jPop] += tmp[iPop][kPop]**DESCRIPTOR::M[kPop][jPop];
        Mt_S_invMt[iPop][jPop] += tmp[iPop][kPop]**DESCRIPTOR::invM[jPop][kPop];
        invM_invMt[iPop][jPop] += DESCRIPTOR::invM[iPop][kPop]**DESCRIPTOR::invM[jPop][kPop];
      }
    }
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      //std::cout << Mt_S_M[iPop][jPop];
      //  std::cout <<Mt_S_invMt[iPop][jPop];
      //    std::cout <<invM_invMt[iPop][jPop];
      //if (iPop == jPop) {
      // Mt_S_M[iPop][jPop] = omega_;
      // Mt_S_invMt[iPop][jPop] = omega_;
      // invM_invMt[iPop][jPop] = 1.;
      //}
    }
  }

  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void DualForcedMRTdynamics<T,DESCRIPTOR,MOMENTA>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop]=rho;
  }
}



template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> DualForcedMRTdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{

  cell.revert();


  T rho_phi_x = T();
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    rho_phi_x += cell[iPop];
  }

  // Who am I?
  T* ct = cell[DESCRIPTOR::q+3];
  int c[3];
  c[0] = (int)ct[0];
  c[1] = (int)ct[1];
  c[2] = (int)ct[2];



  // Preparation
  auto pop_f = cell.template getFieldPointer<opti::F>();
  auto dJdF = cell.template getFieldPointer<opti::DJDF>();

  T rho_f;
  T u_f[3];
  rho_f = T();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    u_f[iD] = T();
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    rho_f += pop_f[iPop];
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u_f[iD] += pop_f[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
  }
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    u_f[iD] /= rho_f;
  }

  // Computation of dual equillibrium
  T Meq_phi[DESCRIPTOR::q][DESCRIPTOR::q];
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T u_c = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u_c += u_f[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      T sum = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        sum += (u_f[iD]*u_f[iD] + 2.*descriptors::c<DESCRIPTOR>(jPop,iD)*(descriptors::c<DESCRIPTOR>(iPop,iD)-u_f[iD]))*descriptors::invCs2<T,DESCRIPTOR>()*0.5+(u_c*descriptors::c<DESCRIPTOR>(iPop,iD)*(2.*descriptors::c<DESCRIPTOR>(jPop,iD)-u_f[iD]))*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*0.5;
      }
      Meq_phi[iPop][jPop] = descriptors::t<T,DESCRIPTOR>(iPop)*(1.+sum);
    }
  }

  // Computation of dual force
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  T F1_phi[DESCRIPTOR::q][DESCRIPTOR::q];
  T F2_phi[DESCRIPTOR::q][DESCRIPTOR::q];
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T f_c = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      f_c += force[iD]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      T sum = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        sum += (f_c*descriptors::c<DESCRIPTOR>(iPop,iD)-force[iD]/(T)descriptors::invCs2<T,DESCRIPTOR>())*descriptors::c<DESCRIPTOR>(jPop,iD);
      }
      F1_phi[iPop][jPop] = descriptors::t<T,DESCRIPTOR>(iPop)*descriptors::invCs2<T,DESCRIPTOR>()*f_c;
      F2_phi[iPop][jPop] = descriptors::t<T,DESCRIPTOR>(iPop)*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*sum;
    }
  }


  // Computation of dual collision
  T dualMRTcollision[DESCRIPTOR::q];
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    dualMRTcollision[iPop] = T();
    for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      T I=T();
      if (iPop==jPop) {
        I=T(1);
      }
      dualMRTcollision[iPop] += (I - (I - Meq_phi[jPop][iPop])*Mt_S_invMt[iPop][iPop] + F1_phi[jPop][iPop] + F2_phi[jPop][iPop]*(1. - 0.5*Mt_S_invMt[iPop][iPop]))*cell[jPop];
    }
  }


  // Adding dirivitive of objective and writing back to cell
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = dualMRTcollision[iPop]  + dJdF[iPop];
  }

  // Incrementing statistic values for convergence
  T phi2 = T();
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    phi2 += cell[iPop]*cell[iPop];
  }

  statistics.incrementStats(1.0 + cell[0], phi2);
  cell.revert();
}







} // namespace opti

} // namespace olb

#endif
