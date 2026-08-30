/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerländer
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

#ifndef TURBULENT_WALL_MODEL_POST_PROCESSOR_H
#define TURBULENT_WALL_MODEL_POST_PROCESSOR_H

namespace olb {

  // interpolation of the velocity along the surface normal on the chosen distance from the boundary cell
  template <typename CELL, typename V = typename CELL::value_t, typename DESCRIPTOR = typename CELL::descriptor_t>
  void interpolateVelocity(CELL& cell, Vector<V,DESCRIPTOR::d> distance, V u_y2[DESCRIPTOR::d]) any_platform{
    Vector<V,DESCRIPTOR::d> floorV;
    for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
      floorV[iD] = util::floor(distance[iD]);
    }
    Vector<Vector<V,DESCRIPTOR::d>,(DESCRIPTOR::d==2)*4+(DESCRIPTOR::d==3)*8> surroundingPoints(floorV);
    surroundingPoints[1][0] += 1.;
    surroundingPoints[2][1] += 1.;
    surroundingPoints[3][0] += 1.; surroundingPoints[3][1] += 1.;
    if( DESCRIPTOR::d == 3) {
      surroundingPoints[4][2] += 1.;
      surroundingPoints[5][0] += 1.; surroundingPoints[5][2] += 1.;
      surroundingPoints[6][1] += 1.; surroundingPoints[6][2] += 1.;
      surroundingPoints[7][0] += 1.; surroundingPoints[7][1] += 1.; surroundingPoints[7][2] += 1.;
    }

    for (auto point : surroundingPoints) {
      const Vector<V,DESCRIPTOR::d> dist = distance - point;
      V weight = V(1);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        weight *=(V(1.) - util::abs(dist[iD]));
      }
      V uNP[DESCRIPTOR::d] = {0.};
      cell.neighbor(point).computeU(uNP);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        u_y2[iD] += uNP[iD]*weight;
      }
    }
  };

  template <typename CELL, typename V = typename CELL::value_t, typename DESCRIPTOR = typename CELL::descriptor_t>
  void interpolateStrainRate(CELL& cell, Vector<V,DESCRIPTOR::d> distance, V pi[util::TensorVal<DESCRIPTOR>::n]) any_platform{
    Vector<V,DESCRIPTOR::d> floorV;
    for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
      floorV[iD] = util::floor(distance[iD]);
    }
    Vector<Vector<V,DESCRIPTOR::d>,(DESCRIPTOR::d==2)*4+(DESCRIPTOR::d==3)*8> surroundingPoints(floorV);
    surroundingPoints[1][0] += 1.;
    surroundingPoints[2][1] += 1.;
    surroundingPoints[3][0] += 1.; surroundingPoints[3][1] += 1.;
    if( DESCRIPTOR::d == 3) {
      surroundingPoints[4][2] += 1.;
      surroundingPoints[5][0] += 1.; surroundingPoints[5][2] += 1.;
      surroundingPoints[6][1] += 1.; surroundingPoints[6][2] += 1.;
      surroundingPoints[7][0] += 1.; surroundingPoints[7][1] += 1.; surroundingPoints[7][2] += 1.;
    }

    for (auto point : surroundingPoints) {
      const Vector<V,DESCRIPTOR::d> dist = distance - point;
      V weight = V(1);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        weight *=(V(1.) - util::abs(dist[iD]));
      }
      V piNP[util::TensorVal<DESCRIPTOR>::n] {V(0)};
      cell.neighbor(point).computeStress(piNP);
      for ( int iN = 0; iN < util::TensorVal<DESCRIPTOR>::n; iN++ ) {
        pi[iN] += piNP[iN] * weight;
      }
    }
  };


  // Spalding wall function written in form of [f = (u+) + ... - (y+)] for the ongoing Newton approximation process
  // vector = {V u_tau, V u_2, V y_Norm, V nu, V u_plus, V y_plus, V spalding_u2}
  template <typename V, typename Par_f>
  V spaldingWallFunction(Par_f & params) any_platform{
    if(int (params[6]) == 2){
      params[4] = params[1]/params[0];
      params[5] = params[2]/params[3]*params[0];
    }
    V k = 0.4;
    V b = 5.5;
    V function = params[4] + util::exp(-k*b)*(util::exp(k*params[4]) - 1. - k*params[4] - 0.5*util::pow(k*params[4],2.) - 1./6.*util::pow(k*params[4],3.)) - params[5];
    return function;
  };

  // Spalding wall function derivatives [f'(u_2+),f'(u+)] for the ongoing Newton approximation process
  // vector = {V x, V u2, V y2, V nu, V k}
  template <typename V, typename Par_d>
  V spaldingDerivative(Par_d & param) any_platform{
    V derivative = 0.;
    V k = 0.4;
    V b = 5.5;
    if(int (param[4]) == 2){
      derivative = (-param[1]*(util::pow(param[0],-2.))) + util::exp(-k*b)*( (util::exp(k*param[1]/param[0]))*(-k*param[1]*(util::pow(param[0],-2.))) + k*param[1]/(util::pow(param[0],2.)) + (util::pow(k*param[1]/param[0],2.))/param[0] + 0.5*(util::pow(k*param[1]/param[0],3.))/param[0] ) - (param[2]/param[3]);
    }
    if(int (param[4]) == 1){
       derivative = 1. +  util::exp(-k*b)*(util::exp(k*param[0])*k - k - k*k*param[0] - 0.5*k*util::pow(k*param[0],2));
    }
    return derivative;
  };

  // Newton approximation process

  // Newton with derivatives, params=(f:function, d:derivative, args: params, optArg: variable, Crit: Convergence)
  template <typename V, typename Funcf,  typename Funcd, typename Par_f, typename Par_d>
  V newtonSpalding(Funcf f, Funcd d, Par_f & argsf, Par_d & argsd, int optArg, V convergeCrit, int iter) any_platform{
    V x = argsf[optArg];
    V x_old = argsf[optArg];
    for(int i = 0; i<iter; i++){
      argsf[optArg] = x;
      argsd[0] = x;
      x -= f(argsf)/d(argsd);
      argsf[optArg] = x;
      argsd[0] = x;
      if(util::abs(f(argsf)) < convergeCrit ){
        break;
      }
      if(util::isnan(x)){
         //printf("x is NAN!!,Applied x=x_old, iter i x_old=%d,%lf,%d\n",i,x_old,optArg);
         x = x_old;
         break;
      }
      x_old = x;
    }
    return x;
  };

  //Power Law profile
  //output = {u_tau, u_plus(x1)}
  template <typename V>
  Vector<V,2> powerLawWallFunction(V nu, V u2, V y2, V y1) any_platform{
    Vector<V,2> output(V(0), V(0));
    V u_tau = util::sqrt(V(0.0246384) * util::pow(nu, V(0.25)) * util::pow(u2, V(1.75)) / util::pow(y2, V(0.25)));
    V y_plus = y1 * u_tau / nu;

    if (y_plus > V(28.178)) {
      output[0] = u_tau;
      output[1] = V(5.5) + V(1./0.399) * util::log(y_plus);
    }
    else if (y_plus < V(28.178)  && y_plus > V(5.)) {
      output[0] = u_tau;
      output[1] = V(-3.05) + V(5.) * util::log(y_plus);
    }
    else {
      output[0] = util::sqrt(V(2.) * u2 * nu / y2);
      y_plus = y1 * output[0] / nu;
      output[1] = y_plus;
    }
    return output;
  };



//======================================================================
// ======== Wall function velocity and vsiscosity PostProcessor ======//
//======================================================================

template <bool bodyForce, bool interpolateSampleVelocity, bool useVanDriest, int wallFunctionProfile, bool movingWall>
class TurbulentWallModelPostProcessor {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::OMEGA, descriptors::SAMPLING_DISTANCE>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    // interpolate the velocity at the y2 point
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto y1 = cell.template getField<descriptors::Y1>();
    V y1Norm = util::norm<DESCRIPTOR::d>(y1);
    if( y1Norm > V(1.e-3) ){
      Vector<V,DESCRIPTOR::d> normal;
      Vector<V,DESCRIPTOR::d> uWall;
      //Vector<V,DESCRIPTOR::d> uWallTang;
      for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        normal[iD] = y1[iD] / y1Norm;
      }
      V samplingCellDistance = parameters.template get<descriptors::SAMPLING_DISTANCE>();
      Vector<V,DESCRIPTOR::d> y2 = samplingCellDistance*normal;
      V u_plusF {};
      auto y1y2 = y2 - y1;
      V u_y2[DESCRIPTOR::d] = {0.};
      if ( interpolateSampleVelocity ) {
        interpolateVelocity(cell, y1y2, u_y2);
      } else {
        Vector<V,DESCRIPTOR::d> roundedNormal;
        for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
          roundedNormal[iD] = int(util::round(y1y2[iD]));
        }
        cell.neighbor(roundedNormal).computeU(u_y2);
      }
      if (movingWall) {
        // I need to go backwards along the normal and take the solid velocity there
        uWall = cell.template getField<descriptors::VELOCITY>();
        for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
          u_y2[iD] -= uWall[iD];
        }
      }
      V u_y2Norm =  util::norm<DESCRIPTOR::d>(u_y2);
      if( u_y2Norm >= V(1.e-7)) {
        V u_1Actual[DESCRIPTOR::d] = {0.};
        cell.computeU(u_1Actual);
        if (movingWall) {
          for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
            u_1Actual[iD] -= uWall[iD];
          }
        }
        V u_y2normal = V(0);
        for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
          u_y2normal += u_y2[iD]*normal[iD];
        }
        auto u_t = u_y2 - u_y2normal*normal;
        auto tangent = u_t/util::norm<DESCRIPTOR::d>(u_t);

        // computing friction velocity u_tau and tangent velocity at y1 with Spalding function and Newton iteration
        const V omega = parameters.template get<descriptors::OMEGA>();
        V nu = (V(1)/omega - 0.5)/descriptors::invCs2<V,DESCRIPTOR>();
        V u_2 = V(0);
        for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
          u_2 += u_y2[iD]*tangent[iD];
        }
        V y2Norm = util::norm<DESCRIPTOR::d>(y2);
        V u_tau = cell.template getField<descriptors::U_TAU>();
        if(u_tau == V(0)){
          u_tau = util::sqrt(nu*util::abs(u_2)/y2Norm);
        }
        V u_plus = util::abs(u_2)/u_tau;
        V y_plus = y2Norm/nu*u_tau;
        V u_1t = V(0);
        for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
          u_1t += u_1Actual[iD]*tangent[iD];
        }
        if ( u_2 != 0. ) {
          if ( wallFunctionProfile == 0 ) {
            Vector<V,2> wf = powerLawWallFunction(nu, util::abs(u_2), y2Norm, y1Norm);
            u_tau = wf[0];
          } else {
            //// NewtonApprox. u_tau
            Vector<V,7> params_u2 = {u_tau, util::abs(u_2), y2Norm, nu, u_plus, y2Norm/nu*u_tau, 2.};  // Vector for spalding wall function
            Vector<V,5> params_d2 = {u_tau, util::abs(u_2), y2Norm, nu, 2.};  // Vector for spalding derivative
            u_tau = newtonSpalding( [](Vector<V,7> par1) { return spaldingWallFunction<V>(par1); }, [](Vector<V,5> par2) { return spaldingDerivative<V>(par2); }, params_u2, params_d2, 0, V{0.0001}, 10);
          }
          cell.template setField<descriptors::U_TAU>(u_tau);

          if (y2Norm/nu*u_tau > V(2.5)) {
            //// Newton Approx. u1_+
            u_plus = util::abs(u_1t)/u_tau;
            y_plus = y1Norm/nu*u_tau;
            if ( wallFunctionProfile == 0 ) {
              Vector<V,2> wf = powerLawWallFunction(nu, util::abs(u_2), y2Norm, y1Norm);
              u_plus = wf[1];
            } else {
              Vector<V,7> params_u1 = {u_tau, util::abs(u_1t), y1Norm, nu, u_plus, y_plus, 1.}; //Vector for spalding wall function
              Vector<V,5> params_d1 = {u_plus, 0., 0., 0., 1.};  // Vector for spalding derivative
              u_plus = newtonSpalding( [](Vector<V,7> par1) { return spaldingWallFunction<V>(par1); }, [](Vector<V,5> par2) { return spaldingDerivative<V>(par2); }, params_u1, params_d1, 4, V{0.0001},10);
            }
            u_plusF = u_plus;

            if(useVanDriest){
              // Van Driest damping
              V dUdY = (u_2 - u_plusF*u_tau)/util::norm<DESCRIPTOR::d>(y1y2);
              V nuTurb = util::pow((V(0.4)*y1Norm*(V(1)-util::exp(-y_plus/V(26)))),2) * util::abs(dUdY);
              cell.template setField<descriptors::VISCOSITY>(nu+nuTurb);
            }
          }
        }
        V u_1 = u_plusF*u_tau;
        // Exchange the computed tangential boundary velocity through the modelled one
        if(util::norm<DESCRIPTOR::d>(u_t) != 0. && u_1 > 0.){
          auto ub = u_1*tangent;
          if (movingWall) {
            ub = ub + uWall;
            auto ub_prev = cell.template getField<descriptors::WMVELOCITY>();
            ub = V{0.9}*ub_prev + V{0.1}*ub;
          }
          if (bodyForce) {
            auto force = cell.template getField<descriptors::FORCE>();
            V u_b[DESCRIPTOR::d] { };
            for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
              u_b[iD] = ub[iD]-V(0.5)*force[iD]*tangent[iD];
            }
            cell.template setField<descriptors::WMVELOCITY>(u_b);
          } else {
            cell.template setField<descriptors::WMVELOCITY>(ub);
          }
        }
      } else {
        cell.template setField<descriptors::WMVELOCITY>(0);
      }
    } else {
      cell.template setField<descriptors::WMVELOCITY>(0);
    }
  }
};

//======================================================================
// PostProcessor for fNeq from neighbor cell along the normal for wall modelling
//======================================================================

  template<int rhoMethod>
  class TurbulentWallModelFneqGuoPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
    using parameters = meta::list<descriptors::OMEGA>;

    int getPriority() const {
      return 2;
    }

    template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform {
      using DESCRIPTOR = typename CELL::descriptor_t;
      auto y1 = cell.template getField<descriptors::Y1>();
      V y1Norm = util::norm<DESCRIPTOR::d>(y1);
      if( y1Norm > V(1.e-3)){
        auto u = cell.template getField<descriptors::WMVELOCITY>();
        V uNorm =  util::norm<DESCRIPTOR::d>(u);
        if ( uNorm < V(1.e-7) ) {
          cell.template setField<descriptors::WMPOROSITY>(V(1));
          cell.template setField<collision::HYBRID>(V(1));
          if( rhoMethod == 1) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
          else if( rhoMethod == 2) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
        } else {
          cell.template setField<descriptors::WMPOROSITY>(V(0));
          cell.template setField<collision::HYBRID>(V(0));
          Vector<V,DESCRIPTOR::d> normal;
          for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
            normal[iD] = y1[iD] / y1Norm;
          }
          Vector<V,DESCRIPTOR::d> yNext = V(1.5)*normal;
          //Vector<V,DESCRIPTOR::d> yNextNext = V(3.)*normal;
          Vector<V,DESCRIPTOR::d> roundedNormal;
          for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
            roundedNormal[iD] = int(util::round(normal[iD]));
          }
          const V omega = parameters.template get<descriptors::OMEGA>();
          V kinVisc = (V(1)/omega - V(0.5))/descriptors::invCs2<V,DESCRIPTOR>();
          V invPreFactor = -V(1)/(omega * kinVisc * descriptors::invCs2<V,DESCRIPTOR>());
          V rho {};
          if( rhoMethod == 0) {
            rho = cell.computeRho();
          }
          else if( rhoMethod == 1) {
            rho = cell.neighbor(roundedNormal).computeRho();
            cell.template setField<descriptors::DENSITY>(rho);
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
          else if( rhoMethod == 2) {
            rho = V(1.);
            cell.template setField<descriptors::DENSITY>(rho);
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
          V factor = invPreFactor * V(2) * kinVisc * rho;
          V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
          V piNext[util::TensorVal<DESCRIPTOR>::n] {V(0)};
          //V piNextNext[util::TensorVal<DESCRIPTOR>::n] {V(0)};
          //cell.neighbor(roundedNormal).computeStress(pi);
          interpolateStrainRate(cell, yNext, piNext);
          //interpolateStrainRate(cell, yNextNext, piNextNext);
          for ( int i = 0; i < util::TensorVal<DESCRIPTOR>::n; i++ ) {
            pi[i] = piNext[i];// + (piNext[i] - piNextNext[i]);
            pi[i] /= factor;
          }
          cell.template setField<descriptors::TENSOR>(pi);
        }
      }
    }
  };

//======================================================================
// PostProcessor for fNeq from finite difference method (FDM) stress tensor for wall modelling
//======================================================================

  template<int rhoMethod>
  class TurbulentWallModelFneqFDMPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCell;

    int getPriority() const {
      return 2;
    }

    template <typename CELL, typename V = typename CELL::value_t>
    void apply(CELL& cell) any_platform{
      using DESCRIPTOR = typename CELL::descriptor_t;
      auto y1 = cell.template getField<descriptors::Y1>();
      V y1Norm = util::norm<DESCRIPTOR::d>(y1);
      if( y1Norm > V(1.e-3)){
        auto u = cell.template getField<descriptors::WMVELOCITY>();
        V uNorm =  util::norm<DESCRIPTOR::d>(u);
        if ( uNorm < V(1.e-7) ) {
          cell.template setField<descriptors::WMPOROSITY>(V(1));
          cell.template setField<collision::HYBRID>(V(1));
          if( rhoMethod == 1) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
          else if( rhoMethod == 2) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
        } else {
        cell.template setField<descriptors::WMPOROSITY>(V(0));
        cell.template setField<collision::HYBRID>(V(0));
          if constexpr ( DESCRIPTOR::d == 2 ) {
            using namespace olb::util::tensorIndices2D;
            V uxDx = 0, uyDx = 0, uxDy = 0, uyDy = 0;
            //----------- FDM in X-direction --------------------------------------------------//
            if (y1[0] > V(0)) {
              auto uXp = cell.neighbor({1,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = uXp[0] - u[0];
              uyDx = uXp[1] - u[1];
            }
            else if (y1[0] < V(0)) {
              auto uXm = cell.neighbor({-1,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = u[0] - uXm[0];
              uyDx = u[1] - uXm[1];
            }
            if (util::abs(y1[0]) < V(1.e-8)) {
              auto uXp = cell.neighbor({1,0}).template getField<descriptors::WMVELOCITY>();
              auto uXm = cell.neighbor({-1,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = V(0.5)*(uXp[0] - uXm[0]);
              uyDx = V(0.5)*(uXp[1] - uXm[1]);
            }
            //----------- FDM in Y-direction --------------------------------------------------//
            if (y1[1] > V(0)) {
              auto uYp = cell.neighbor({0,1}).template getField<descriptors::WMVELOCITY>();
              uxDy = uYp[0] - u[0];
              uyDy = uYp[1] - u[1];
            }
            else if (y1[1] < V(0)) {
              auto uYm = cell.neighbor({0,-1}).template getField<descriptors::WMVELOCITY>();
              uxDy = u[0] - uYm[0];
              uyDy = u[1] - uYm[1];
            }
            if (util::abs(y1[1]) < V(1.e-8)) {
              auto uYp = cell.neighbor({0,1}).template getField<descriptors::WMVELOCITY>();
              auto uYm = cell.neighbor({0,-1}).template getField<descriptors::WMVELOCITY>();
              uxDy = V(0.5)*(uYp[0] - uYm[0]);
              uyDy = V(0.5)*(uYp[1] - uYm[1]);
            }
            //----------- TENSOR tensor with FDM --------------------------------------------------//
            V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
            pi[xx] = V(0.5)*(V)2 * uxDx;
            pi[yy] = V(0.5)*(V)2 * uyDy;
            pi[xy] = V(0.5)*(uxDy + uyDx);
            cell.template setField<descriptors::TENSOR>(pi);
          }
          else {
            using namespace olb::util::tensorIndices3D;
            V uxDx = 0, uyDx = 0, uzDx = 0, uxDy = 0, uyDy = 0, uzDy = 0, uxDz = 0, uyDz = 0, uzDz = 0;
            //----------- FDM in X-direction --------------------------------------------------//
            if (y1[0] > V(0)) {
              auto uXp = cell.neighbor({1,0,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = uXp[0] - u[0];
              uyDx = uXp[1] - u[1];
              uzDx = uXp[2] - u[2];
            }
            else if (y1[0] < V(0)) {
              auto uXm = cell.neighbor({-1,0,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = u[0] - uXm[0];
              uyDx = u[1] - uXm[1];
              uzDx = u[2] - uXm[2];
            }
            if (util::abs(y1[0]) < V(1.e-8)) {
              auto uXp = cell.neighbor({1,0,0}).template getField<descriptors::WMVELOCITY>();
              auto uXm = cell.neighbor({-1,0,0}).template getField<descriptors::WMVELOCITY>();
              uxDx = V(0.5)*(uXp[0] - uXm[0]);
              uyDx = V(0.5)*(uXp[1] - uXm[1]);
              uzDx = V(0.5)*(uXp[2] - uXm[2]);
            }
            //----------- FDM in Y-direction --------------------------------------------------//
            if (y1[1] > V(0)) {
              auto uYp = cell.neighbor({0,1,0}).template getField<descriptors::WMVELOCITY>();
              uxDy = uYp[0] - u[0];
              uyDy = uYp[1] - u[1];
              uzDy = uYp[2] - u[2];
            }
            else if (y1[1] < V(0)) {
              auto uYm = cell.neighbor({0,-1,0}).template getField<descriptors::WMVELOCITY>();
              uxDy = u[0] - uYm[0];
              uyDy = u[1] - uYm[1];
              uzDy = u[2] - uYm[2];
            }
            if (util::abs(y1[1]) < V(1.e-8)) {
              auto uYp = cell.neighbor({0,1,0}).template getField<descriptors::WMVELOCITY>();
              auto uYm = cell.neighbor({0,-1,0}).template getField<descriptors::WMVELOCITY>();
              uxDy = V(0.5)*(uYp[0] - uYm[0]);
              uyDy = V(0.5)*(uYp[1] - uYm[1]);
              uzDy = V(0.5)*(uYp[2] - uYm[2]);
            }
            //----------- FDM in Z-direction --------------------------------------------------//
            if (y1[2] > V(0)) {
              auto uZp = cell.neighbor({0,0,1}).template getField<descriptors::WMVELOCITY>();
              uxDz = uZp[0] - u[0];
              uyDz = uZp[1] - u[1];
              uzDz = uZp[2] - u[2];
            }
            else if (y1[2] < V(0)) {
              auto uZm = cell.neighbor({0,0,-1}).template getField<descriptors::WMVELOCITY>();
              uxDz = u[0] - uZm[0];
              uyDz = u[1] - uZm[1];
              uzDz = u[2] - uZm[2];
            }
            if (util::abs(y1[2]) < V(1.e-8)) {
              auto uZp = cell.neighbor({0,0,1}).template getField<descriptors::WMVELOCITY>();
              auto uZm = cell.neighbor({0,0,-1}).template getField<descriptors::WMVELOCITY>();
              uxDz = V(0.5)*(uZp[0] - uZm[0]);
              uyDz = V(0.5)*(uZp[1] - uZm[1]);
              uzDz = V(0.5)*(uZp[2] - uZm[2]);
            }
            //----------- TENSOR tensor with FDM --------------------------------------------------//
            V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
            pi[xx] = V(0.5)*(V)2 * uxDx;
            pi[yy] = V(0.5)*(V)2 * uyDy;
            pi[zz] = V(0.5)*(V)2 * uzDz;
            pi[xy] = V(0.5)*(uxDy + uyDx);
            pi[xz] = V(0.5)*(uxDz + uzDx);
            pi[yz] = V(0.5)*(uyDz + uzDy);
            cell.template setField<descriptors::TENSOR>(pi);
          }

          if( rhoMethod == 1) {
            Vector<V,DESCRIPTOR::d> normal;
            for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
              normal[iD] = y1[iD] / y1Norm;
            }
            Vector<V,DESCRIPTOR::d> roundedNormal;
            for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
              roundedNormal[iD] = int(util::round(normal[iD]));
            }
            V extRho = cell.neighbor(roundedNormal).computeRho();
            cell.template setField<descriptors::DENSITY>(extRho);
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
          else if( rhoMethod == 2) {
            cell.template setField<descriptors::DENSITY>(V(1.));
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
        }
      }
    }
  };



//======================================================================
// PostProcessor for fNeq from finite difference method (FDM) stress tensor for wall modelling at porous domains
//======================================================================

  template<int rhoMethod>
  class TurbulentWallModelPorousFneqFDMPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCell;

    int getPriority() const {
      return 2;
    }

    template <typename CELL, typename V = typename CELL::value_t>
    void apply(CELL& cell) any_platform{
      using DESCRIPTOR = typename CELL::descriptor_t;
      auto y1 = cell.template getField<descriptors::Y1>();
      V y1Norm = util::norm<DESCRIPTOR::d>(y1);
      if( y1Norm > V(1.e-3)){
        auto u = cell.template getField<descriptors::WMVELOCITY>();
        V uNorm =  util::norm<DESCRIPTOR::d>(u);
        if ( uNorm < V(1.e-7) ) {
          cell.template setField<descriptors::WMPOROSITY>(V(1));
          cell.template setField<collision::HYBRID>(V(1));
          if( rhoMethod == 1) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
          else if( rhoMethod == 2) {
            cell.template setField<collision::HYBRID_RHO>(V(1));
          }
        } else {
        cell.template setField<descriptors::WMPOROSITY>(V(0));
        cell.template setField<collision::HYBRID>(V(0));
          if constexpr ( DESCRIPTOR::d == 2 ) {
            using namespace olb::util::tensorIndices2D;
            V uxDx = 0, uyDx = 0, uxDy = 0, uyDy = 0;
            //----------- FDM in X-direction --------------------------------------------------//
            if (y1[0] > V(0) || cell.neighbor({-1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXp[DESCRIPTOR::d] {};
              cell.neighbor({1,0}).computeU(uXp);
              uxDx = uXp[0] - u[0];
              uyDx = uXp[1] - u[1];
            }
            else if (y1[0] < V(0) || cell.neighbor({1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXm[DESCRIPTOR::d] {};
              cell.neighbor({-1,0}).computeU(uXm);
              uxDx = u[0] - uXm[0];
              uyDx = u[1] - uXm[1];
            }
            if (util::abs(y1[0]) < V(1.e-8) &&
                cell.neighbor({-1,0}).template getField<descriptors::POROSITY>() == V(0) &&
                cell.neighbor({1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXp[DESCRIPTOR::d], uXm[DESCRIPTOR::d] {};
              cell.neighbor({1,0}).computeU(uXp);
              cell.neighbor({-1,0}).computeU(uXm);
              uxDx = V(0.5)*(uXp[0] - uXm[0]);
              uyDx = V(0.5)*(uXp[1] - uXm[1]);
            }
            //----------- FDM in Y-direction --------------------------------------------------//
            if (y1[1] > V(0) || cell.neighbor({0,-1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYp[DESCRIPTOR::d] {};
              cell.neighbor({0,1}).computeU(uYp);
              uxDy = uYp[0] - u[0];
              uyDy = uYp[1] - u[1];
            }
            else if (y1[1] < V(0) || cell.neighbor({0,1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYm[DESCRIPTOR::d] {};
              cell.neighbor({0,-1}).computeU(uYm);
              uxDy = u[0] - uYm[0];
              uyDy = u[1] - uYm[1];
            }
            if (util::abs(y1[1]) < V(1.e-8) &&
                cell.neighbor({0,-1}).template getField<descriptors::POROSITY>() == V(0) &&
                cell.neighbor({0,1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYp[DESCRIPTOR::d], uYm[DESCRIPTOR::d] {};
              cell.neighbor({0,1}).computeU(uYp);
              cell.neighbor({0,-1}).computeU(uYm);
              uxDy = V(0.5)*(uYp[0] - uYm[0]);
              uyDy = V(0.5)*(uYp[1] - uYm[1]);
            }
            //----------- TENSOR tensor with FDM --------------------------------------------------//
            V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
            pi[xx] = V(0.5)*(V)2 * uxDx;
            pi[yy] = V(0.5)*(V)2 * uyDy;
            pi[xy] = V(0.5)*(uxDy + uyDx);
            cell.template setField<descriptors::TENSOR>(pi);
          }
          else {
            using namespace olb::util::tensorIndices3D;
            V uxDx = 0, uyDx = 0, uzDx = 0, uxDy = 0, uyDy = 0, uzDy = 0, uxDz = 0, uyDz = 0, uzDz = 0;
            //----------- FDM in X-direction --------------------------------------------------//
            if (y1[0] > V(0) || cell.neighbor({-1,0,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXp[DESCRIPTOR::d] {};
              cell.neighbor({1,0,0}).computeU(uXp);
              uxDx = uXp[0] - u[0];
              uyDx = uXp[1] - u[1];
              uzDx = uXp[2] - u[2];
            }
            else if (y1[0] < V(0) || cell.neighbor({1,0,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXm[DESCRIPTOR::d] {};
              cell.neighbor({-1,0,0}).computeU(uXm);
              uxDx = u[0] - uXm[0];
              uyDx = u[1] - uXm[1];
              uzDx = u[2] - uXm[2];
            }
            if (util::abs(y1[0]) < V(1.e-8) &&
                cell.neighbor({-1,0,0}).template getField<descriptors::POROSITY>() == V(0) &&
                cell.neighbor({1,0,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uXp[DESCRIPTOR::d], uXm[DESCRIPTOR::d] {};
              cell.neighbor({1,0,0}).computeU(uXp);
              cell.neighbor({-1,0,0}).computeU(uXm);
              uxDx = V(0.5)*(uXp[0] - uXm[0]);
              uyDx = V(0.5)*(uXp[1] - uXm[1]);
              uzDx = V(0.5)*(uXp[2] - uXm[2]);
            }
            //----------- FDM in Y-direction --------------------------------------------------//
            if (y1[1] > V(0) || cell.neighbor({0,-1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYp[DESCRIPTOR::d] {};
              cell.neighbor({0,1,0}).computeU(uYp);
              uxDy = uYp[0] - u[0];
              uyDy = uYp[1] - u[1];
              uzDy = uYp[2] - u[2];
            }
            else if (y1[1] < V(0) || cell.neighbor({0,1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYm[DESCRIPTOR::d] {};
              cell.neighbor({0,-1,0}).computeU(uYm);
              uxDy = u[0] - uYm[0];
              uyDy = u[1] - uYm[1];
              uzDy = u[2] - uYm[2];
            }
            if (util::abs(y1[1]) < V(1.e-8) &&
                cell.neighbor({0,-1,0}).template getField<descriptors::POROSITY>() == V(0) &&
                cell.neighbor({0,1,0}).template getField<descriptors::POROSITY>() == V(0)) {
              V uYp[DESCRIPTOR::d], uYm[DESCRIPTOR::d] {};
              cell.neighbor({0,1,0}).computeU(uYp);
              cell.neighbor({0,-1,0}).computeU(uYm);
              uxDy = V(0.5)*(uYp[0] - uYm[0]);
              uyDy = V(0.5)*(uYp[1] - uYm[1]);
              uzDy = V(0.5)*(uYp[2] - uYm[2]);
            }
            //----------- FDM in Z-direction --------------------------------------------------//
            if (y1[2] > V(0) || cell.neighbor({0,0,-1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uZp[DESCRIPTOR::d] {};
              cell.neighbor({0,0,1}).computeU(uZp);
              uxDz = uZp[0] - u[0];
              uyDz = uZp[1] - u[1];
              uzDz = uZp[2] - u[2];
            }
            else if (y1[2] < V(0) || cell.neighbor({0,0,1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uZm[DESCRIPTOR::d] {};
              cell.neighbor({0,0,-1}).computeU(uZm);
              uxDz = u[0] - uZm[0];
              uyDz = u[1] - uZm[1];
              uzDz = u[2] - uZm[2];
            }
            if (util::abs(y1[2]) < V(1.e-8) &&
                cell.neighbor({0,0,-1}).template getField<descriptors::POROSITY>() == V(0) &&
                cell.neighbor({0,0,1}).template getField<descriptors::POROSITY>() == V(0)) {
              V uZp[DESCRIPTOR::d], uZm[DESCRIPTOR::d] {};
              cell.neighbor({0,0,1}).computeU(uZp);
              cell.neighbor({0,0,-1}).computeU(uZm);
              uxDz = V(0.5)*(uZp[0] - uZm[0]);
              uyDz = V(0.5)*(uZp[1] - uZm[1]);
              uzDz = V(0.5)*(uZp[2] - uZm[2]);
            }
            //----------- TENSOR tensor with FDM --------------------------------------------------//
            V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
            pi[xx] = V(0.5)*(V)2 * uxDx;
            pi[yy] = V(0.5)*(V)2 * uyDy;
            pi[zz] = V(0.5)*(V)2 * uzDz;
            pi[xy] = V(0.5)*(uxDy + uyDx);
            pi[xz] = V(0.5)*(uxDz + uzDx);
            pi[yz] = V(0.5)*(uyDz + uzDy);
            cell.template setField<descriptors::TENSOR>(pi);
          }

          if( rhoMethod == 1) {
            Vector<V,DESCRIPTOR::d> normal;
            for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
              normal[iD] = y1[iD] / y1Norm;
            }
            Vector<V,DESCRIPTOR::d> roundedNormal;
            for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
              roundedNormal[iD] = int(util::round(normal[iD]));
            }
            V extRho = cell.neighbor(roundedNormal).computeRho();
            cell.template setField<descriptors::DENSITY>(extRho);
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
          else if( rhoMethod == 2) {
            cell.template setField<descriptors::DENSITY>(V(1.));
            cell.template setField<collision::HYBRID_RHO>(V(0));
          }
        }
      } else {
        cell.template setField<descriptors::WMPOROSITY>(V(1));
        cell.template setField<collision::HYBRID>(V(1));
        if( rhoMethod == 1) {
          cell.template setField<collision::HYBRID_RHO>(V(1));
        }
        else if( rhoMethod == 2) {
          cell.template setField<collision::HYBRID_RHO>(V(1));
        }
      }
    }
  };


//======================================================================
// PostProcessor for zero fNeq for wall modelling
//======================================================================

  template<int rhoMethod>
  class TurbulentWallModelFneqZeroPostProcessor {
  public:
    static constexpr OperatorScope scope = OperatorScope::PerCell;

    int getPriority() const {
      return 2;
    }

    template <typename CELL, typename V = typename CELL::value_t>
    void apply(CELL& cell) any_platform{
      using DESCRIPTOR = typename CELL::descriptor_t;
      auto u = cell.template getField<descriptors::WMVELOCITY>();
      V uNorm =  util::norm<DESCRIPTOR::d>(u);
      if ( uNorm < V(1.e-7) ) {
        cell.template setField<descriptors::WMPOROSITY>(V(1));
        cell.template setField<collision::HYBRID>(V(1));
        if( rhoMethod == 1) {
          cell.template setField<collision::HYBRID_RHO>(V(1));
        }
        else if( rhoMethod == 2) {
          cell.template setField<collision::HYBRID_RHO>(V(1));
        }
      } else {
        if( rhoMethod == 1) {
          auto y1 = cell.template getField<descriptors::Y1>();
          V y1Norm = util::norm<DESCRIPTOR::d>(y1);
          Vector<V,DESCRIPTOR::d> normal;
          for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
            normal[iD] = y1[iD] / y1Norm;
          }
          Vector<V,DESCRIPTOR::d> roundedNormal;
          for ( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
            roundedNormal[iD] = int(util::round(normal[iD]));
          }
          V extRho = cell.neighbor(roundedNormal).computeRho();
          cell.template setField<descriptors::DENSITY>(extRho);
          cell.template setField<collision::HYBRID_RHO>(V(0));
        }
        else if( rhoMethod == 2) {
          cell.template setField<descriptors::DENSITY>(V(1.));
          cell.template setField<collision::HYBRID_RHO>(V(0));
        }
        cell.template setField<descriptors::WMPOROSITY>(V(0));
        cell.template setField<collision::HYBRID>(V(0));
        V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};

        cell.template setField<descriptors::TENSOR>(pi);
      }
    }
  };
}
#endif
