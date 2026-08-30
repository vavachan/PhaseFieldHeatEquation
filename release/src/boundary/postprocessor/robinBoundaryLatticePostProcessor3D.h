/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
 *                2022 Nando Suntoyo, Adrian Kummerlaender, Shota Ito,
 *                2024 Marc Heinzelmann
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

#ifndef ROBIN_BOUNDARY_LATTICE_POST_PROCESSOR_3D_H
#define ROBIN_BOUNDARY_LATTICE_POST_PROCESSOR_3D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "core/util.h"
#include "descriptor/descriptor.h"
#include "utilities/omath.h"
#include <assert.h>

namespace olb {


//======================================================================
// ======== Robin Boundary for AD 3D ======//
//======================================================================


/**
 * First scheme adapted from
 * Xuhui Meng and Zhaoli Guo. “Boundary scheme for linear heterogeneous
 * surface reactions in the lattice Boltzmann method”.
 * In: Physical Review E 94.5 (2016), doi.org/10.1103/PhysRevE.94.053307
 * Eq.(20)

 * Second scheme from
 * Long Ju, Chunhua Zhang, and Zhaoli Guo. “Local reactive boundary scheme for
 * irregular geometries in lattice Boltzmann method”. In: International Journal of Heat
 * and Mass Transfer 150 (2020), p. 119314. doi: 10.1016/j.ijheatmasstransfer.2020.119314
 **/


// ======== MengCurvedCorr ======== //
template<typename T, typename DESCRIPTOR, int Normal1, int Normal2, int Normal3>
struct robinBoundaryLatticePostProcessor3D
{
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::OMEGA>;
  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    T omega = parameters.template get<descriptors::OMEGA>();
    auto v = cell.template getField<descriptors::VELOCITY>();
    auto a = cell.template getField<descriptors::G>(); //get robin boundary coefficients
    T a1 = a[0];
    T a2 = a[1];
    T a3 = a[2];

    const int direction = abs(Normal1*0 + Normal2*1 + Normal3*2); // ==0,1,2
    const int orientation = -(Normal1 + Normal2 + Normal3); // ==+-1
    Vector<int,3> n(-Normal1, -Normal2, -Normal3); //normal points into domain
    T NdotV = n[0]*v[0] + n[1]*v[1] + n[2]*v[2];
    T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
    T gamma = -descriptors::invCs2<T,DESCRIPTOR>() * omega; //a factor used in scheme

    Vector<T,1> unknownIndices;
    for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++){
      if(descriptors::c<DESCRIPTOR>(iPop, direction) == orientation){
        unknownIndices[0] = iPop;
      }
    }

    ///define variables, names match the characters in the formula in thesis
    T b1 = a1-a2*gamma*NdotV;
    T b2 = a2*gamma;

    ///calculation of the two sums used in the MengCurvedCorr scheme
    T sum1 = cell.computeRho();
    T sum2 = 0;
    for (unsigned iPop : unknownIndices) {
      sum1 -= (cell[iPop]+w + cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
      sum2 += cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w;
    }
    T beta = (a3-b1*sum1+2*b2*sum2)/(unknownIndices.size()*w*(b1+b2));

    for (unsigned iPop : unknownIndices) {
      cell[iPop] = w * beta - (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) - w; //set new populations
    }
  }
};

// ======== Ju2020 ======== //
template<typename T, typename DESCRIPTOR, int Normal1, int Normal2, int Normal3>
struct robinBoundaryLatticePostProcessor3Dother
{
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::OMEGA>;
  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    T omega = parameters.template get<descriptors::OMEGA>();
    auto v = cell.template getField<descriptors::VELOCITY>();
    auto a = cell.template getField<descriptors::G>(); //get robin boundary coefficients
    T a1 = a[0];
    T a2 = a[1];
    T a3 = a[2];

    const int direction = abs(Normal1*0 + Normal2*1 + Normal3*2); // ==0,1,2
    const int orientation = -(Normal1 + Normal2 + Normal3); // ==+-1
    Vector<int,3> n(-Normal1, -Normal2, -Normal3); //normal points into domain
    T NdotV = n[0]*v[0] + n[1]*v[1] + n[2]*v[2];
    T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
    T gamma = -descriptors::invCs2<T,DESCRIPTOR>() * omega; //a factor used in scheme

    Vector<T,1> unknownIndices;
    for(int iPop = 0; iPop<DESCRIPTOR::q; iPop++){
      if(descriptors::c<DESCRIPTOR>(iPop, direction) == orientation){
        unknownIndices[0] = iPop;
      }
    }

    ///define variables, names match the characters in the formula in thesis
    T chi = (a2!=0.) ? 1/ (gamma*a2) : (1/omega)/(1/omega-0.5);
    T k = a1; //reaction rate
    T A = unknownIndices.size()-1;
    T alpha = (-chi*k+NdotV+2*w*(1-A)) / (chi*k-NdotV+2*w*(1+A)); //bounced-back part

    ///necessary procedure to calculate B
    T B = 0;
    for (unsigned iPop : unknownIndices) {
      B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
    }

    ///final calculation of the new populations
    for (unsigned iPop : unknownIndices) {
      B -= 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //final value of B
      T beta = (2*w*chi*a3+2*w*B) / (chi*k-NdotV+2*w*(1+A)); //reverse reaction contribution
      cell[iPop] = alpha * (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) + beta - w; //set new populations
      B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //restore B
    }
  }
};

// ======== Edge ======== //
template<typename T, typename DESCRIPTOR, int Plane, int Normal1, int Normal2>
struct robinBoundaryExtendedPostProcessor3DEdges
{
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::OMEGA>;
  int getPriority() const {
    return 2;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    bool schemeSwitch = false; //switch between two schemes; true=MengCurvedCorr, false=Ju2020

    // get external values
    T omega = parameters.template get<descriptors::OMEGA>();
    auto v = cell.template getField<descriptors::VELOCITY>();
    auto a = cell.template getField<descriptors::G>();
    T a1 = a[0];
    T a2 = a[1];
    T a3 = a[2];

    constexpr auto unknownIndices = util::subIndexOutgoing3DonEdges<DESCRIPTOR,Plane,Normal1,Normal2>();
    assert(unknownIndices.size() == 2);
    auto c0=descriptors::c<DESCRIPTOR>(unknownIndices[0]);
    auto c1=descriptors::c<DESCRIPTOR>(unknownIndices[1]);
    Vector<int,3> n(c0[0]+c1[0], c0[1]+c1[1], c0[2]+c1[2]); //normal points into domain


    if(schemeSwitch){
      ///define variables, names match the characters in the formula
      T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
      T gamma = -descriptors::invCs2<T,DESCRIPTOR>() * omega;
      T b1 = a1-a2*gamma*(n[0]*v[0]+n[1]*v[1]+n[2]*v[2]);
      T b2 = a2*gamma;

      ///calculation of the two sums used in the MengCurvedCorr scheme
      T sum1 = cell.computeRho();
      T sum2 = 0;
      for (unsigned iPop : unknownIndices) {
        sum1 -= (cell[iPop]+w + cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
        sum2 += cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w;
      }
      T beta = (a3-b1*sum1+2*b2*sum2)/(unknownIndices.size()*w*(b1+b2));

      for (unsigned iPop : unknownIndices) {
        cell[iPop] = w * beta - (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) - w; //set new populations
      }
    } else{
      ///define variables, names match the characters in the formula
      T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
      T chi = (a2!=0.) ? -1/ (descriptors::invCs2<T,DESCRIPTOR>()*a2*omega) : (1/omega)/(1/omega-0.5);
      T k = a1; //reaction rate
      T NdotV = n[0]*v[0] + n[1]*v[1] + n[2]*v[2];
      T A = unknownIndices.size()-1;
      T alpha = (-chi*k+NdotV+2*w*(1-A)) / (chi*k-NdotV+2*w*(1+A)); //bounced-back part

      ///necessary procedure to calculate B
      T B = 0;
      for (unsigned iPop : unknownIndices) {
        B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
      }

      ///final calculation of the new populations
      for (unsigned iPop : unknownIndices) {
        B -= 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //final value of B
        T beta = (2*w*chi*a3+2*w*B) / (chi*k-NdotV+2*w*(1+A)); //reverse reaction contribution
        cell[iPop] = alpha * (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) + beta - w; //set new populations
        B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //restore B
      }
    }
  }
};


// ======== Corner ======== //
template<typename T, typename DESCRIPTOR, int Normal1, int Normal2, int Normal3>
struct robinBoundaryExtendedPostProcessor3DCorners
{
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::OMEGA>;
  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{
    bool schemeSwitch = false; //switch between two schemes; true=MengCurvedCorr, false=Ju2020

    // get external values
    T omega = parameters.template get<descriptors::OMEGA>();
    auto v = cell.template getField<descriptors::VELOCITY>(); //velocity field
    auto a = cell.template getField<descriptors::G>(); //field containing coefficients for the robin boundary condition
    T a1 = a[0]; T a2 = a[1]; T a3 = a[2];
    std::vector<int> n(3);
    n[0] = -Normal1; n[1] = -Normal2; n[2] = -Normal3; //change normal to point into domain

    constexpr auto unknownIndices = util::subIndexOutgoing3DonCorners<DESCRIPTOR,Normal1,Normal2,Normal3>(); //function returns the indices of missing populations
    assert(unknownIndices.size() == 3);

    if(schemeSwitch){
      //--- MengCurvedCorr ---//
      ///define variables, names match the characters in the formula
      T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
      T gamma = -descriptors::invCs2<T,DESCRIPTOR>() * omega;
      T b1 = a1-a2*gamma*(n[0]*v[0]+n[1]*v[1]+n[2]*v[2]);
      T b2 = a2*gamma;

      ///calculation of the two sums used in the MengCurvedCorr scheme
      T sum1 = cell.computeRho();
      T sum2 = 0;
      for (unsigned iPop : unknownIndices) {
        sum1 -= (cell[iPop]+w + cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
        sum2 += cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w;
      }
      T beta = (a3-b1*sum1+2*b2*sum2)/(unknownIndices.size()*w*(b1+b2));

      for (unsigned iPop : unknownIndices) {
        cell[iPop] = w * beta - (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) - w; //set new populations
      }
    } else{
      //--- Ju2020 ---//
      ///define variables, names match the characters in the formula
      T w = descriptors::t<T,DESCRIPTOR>(1); //lattice weight = 1/8
      T chi = (a2!=0.) ? -1/ (descriptors::invCs2<T,DESCRIPTOR>()*a2*omega) : (1/omega)/(1/omega-0.5);
      T k = a1; //reaction rate
      T NdotV = n[0]*v[0] + n[1]*v[1] + n[2]*v[2];
      T A = unknownIndices.size()-1;
      T alpha = (-chi*k+NdotV+2*w*(1-A)) / (chi*k-NdotV+2*w*(1+A)); //bounced-back part

      ///necessary procedure to calculate B
      T B = 0;
      for (unsigned iPop : unknownIndices) {
        B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w);
      }

      ///final calculation of the new populations
      for (unsigned iPop : unknownIndices) {
        B -= 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //final value of B
        T beta = (2*w*chi*a3+2*w*B) / (chi*k-NdotV+2*w*(1+A)); //reverse reaction part
        cell[iPop] = alpha * (cell[descriptors::opposite<DESCRIPTOR>(iPop)] + w) + beta - w; //set new populations
        B += 2*(cell[descriptors::opposite<DESCRIPTOR>(iPop)]+w); //restore B for next iteration
      }
    }
  }
};



}//end namespace

#endif
