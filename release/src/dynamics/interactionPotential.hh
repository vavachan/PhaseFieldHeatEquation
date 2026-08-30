/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Peter Weisbrod
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

#ifndef INTERACTION_POTENTIAL_HH
#define INTERACTION_POTENTIAL_HH


#include "dynamics/interactionPotential.h"


namespace olb {



template <typename T, typename S>
ShanChen93<T,S>::ShanChen93(T rhoZero) : AnalyticalF<1,T,S>(1), _rhoZero(rhoZero)
{
  this->getName() = "ShanChen93";
}

template <typename T, typename S>
bool ShanChen93<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=util::sqrt(_rhoZero)*(1-util::exp(-(rho[0]/_rhoZero)));
  return true;
}


template <typename T, typename S>
ShanChen94<T,S>::ShanChen94(T rhoZero, T psiZero) : AnalyticalF<1,T,S>(1), _rhoZero(rhoZero), _psiZero(psiZero)
{
  this->getName() = "ShanChen94";
}

template <typename T, typename S>
bool ShanChen94<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=_psiZero*util::exp(-_rhoZero/rho[0]);
  return true;
}


template <typename T, typename S>
PengRobinson<T,S>::PengRobinson(T G, T acentricFactor, T a, T b, T tr) : AnalyticalF<1,T,S>(1), _G(G), _acentricFactor(acentricFactor), _a(a), _b(b)
{
  _R = 1.;
  //a=0.45724*R*R*tc*tc/pc;
  //b=0.0778*R*tc/pc;
  _tc = 0.0778/0.45724*_a/_b/_R;
  //T pc = 0.0778*_R*tc/_b;
  //T rhoc = pc/0.307/_R/tc;
  _t = _tc*tr;
  //Zc=0.307 Tc=0.072922004 pc=0.059569985 rhoc=2.6609121
  _alpha = 1. + (0.37464+1.54226*_acentricFactor-0.26992*_acentricFactor*_acentricFactor)*(1.-util::sqrt(_t/_tc));
  _alpha = _alpha*_alpha;
  this->getName() = "PengRobinson";
}

template <typename T, typename S>
bool PengRobinson<T,S>::operator()(T psi[], const S rho[])
{
  T p = (rho[0]*_R*_t/(1.-_b*rho[0]))-(_a*_alpha*rho[0]*rho[0]/(1.+2.*_b*rho[0]-_b*_b*rho[0]*rho[0]));
  psi[0] = util::sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool PengRobinson<T,S>::operator()(T psi[], const S rho[], const S t[])
{
  _t = t[0];
  _alpha = 1. + (0.37464+1.54226*_acentricFactor-0.26992*_acentricFactor*_acentricFactor)*(1.-util::sqrt(_t/_tc));
  _alpha = _alpha*_alpha;
  T p = (rho[0]*_R*_t/(1.-_b*rho[0]))-(_a*_alpha*rho[0]*rho[0]/(1.+2.*_b*rho[0]-_b*_b*rho[0]*rho[0]));
  psi[0] = util::sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}


template <typename T, typename S>
CarnahanStarling<T,S>::CarnahanStarling(T G, T a, T b, T tr) : AnalyticalF<1,T,S>(1), _G(G), _a(a), _b(b)
{
  _R = 1.;
  //a=0.4963*tc*tc*R*R/pc;
  //b=0.18727*R*tc/pc;
  T tc = 0.18727/0.4963*_a/_b/_R;
  //T pc = 0.18727*_R*tc/_b;
  //T rhoc = pc/0.35930763/_R/tc;
  _t = tc*tr;
  //Zc=0.35930763 Tc=0.094333065 pc=0.0044164383 rhoc=0.13029921
  this->getName() = "CarnahanStarling";
}

template <typename T, typename S>
bool CarnahanStarling<T,S>::operator()(T psi[], const S rho[])
{
  T c = _b*rho[0]/4.;
  T p = rho[0]*_R*_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-_a*rho[0]*rho[0];
  psi[0] = util::sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool CarnahanStarling<T,S>::operator()(T psi[], const S rho[], const S t[])
{
  _t = t[0];
  T c = _b*rho[0]/4.;
  T p = rho[0]*_R*_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-_a*rho[0]*rho[0];
  psi[0] = util::sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}


template <typename T, typename S>
PsiEqualsRho<T,S>::PsiEqualsRho() : AnalyticalF<1,T,S>(1)
{
  this->getName() = "PsiEqualsRho";
}

template <typename T, typename S>
bool PsiEqualsRho<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=rho[0];
  return true;
}


template <typename T, typename S>
Krause<T,S>::Krause(T rhoZero, T psiZero) : AnalyticalF<1,T,S>(1), _rhoZero(rhoZero), _psiZero(psiZero)
{
  this->getName() = "Krause";
}

template <typename T, typename S>
bool Krause<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=_psiZero/1.77*1.414/_rhoZero*util::exp(-(_rhoZero-rho[0])*(_rhoZero-rho[0])/_rhoZero/_rhoZero);
  return true;
}


template <typename T, typename S>
WeisbrodKrause<T,S>::WeisbrodKrause(T rhoZero, T sigmu) : AnalyticalF<1,T,S>(1), _rhoZero(rhoZero), _sigmu(sigmu)
{
  _rhoZero=_rhoZero/1.005088;
  this->getName() = "WeisbrodKrause";
}

template <typename T, typename S>
bool WeisbrodKrause<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=util::sqrt(_rhoZero)*1.5179/_sigmu*util::exp(-(_sigmu-rho[0]/_rhoZero)*(_sigmu-rho[0]/_rhoZero)/_sigmu/_sigmu);
  return true;
}


template <typename T, typename S>
Normal<T,S>::Normal(T sigma, T mu) : AnalyticalF<1,T,S>(1), _sigma(sigma), _mu(mu)
{
  this->getName() = "Normal";
}

template <typename T, typename S>
bool Normal<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=1./2.507/_sigma*util::exp(-(rho[0]-_mu)*(rho[0]-_mu)/_sigma/_sigma/2.);
  return true;
}

MultiComponentPengRobinson::MultiComponentPengRobinson(double p_, double T_, std::vector<double> z_, std::vector<double> a_L, std::vector<double> b_L,
                                                       std::vector<double> M_L, std::vector<double> Tc_L, std::vector<double> pc_L, std::vector<double> omega_, std::vector<double> devi,
                                                       std::vector<double> alpha_, std::vector<double> gI_L, std::vector<double> gII_L) :
p(p_), T(T_), z(z_), num_comp(z.size()), a_c(a_L), b(b_L), M(M_L), T_c(Tc_L), p_c(pc_L), omega(omega_), m(devi), alpha(alpha_), g_I(gI_L), g_II(gII_L)
{ }
double MultiComponentPengRobinson::G_Excess(std::vector<double> x_i, std::vector<double> g_jk, std::vector<double> beta_jk){
  double G_E = 0;
  for(int i = 0; i < num_comp; i++){
    double sum1 = 0;
    double sum2 = 0;
    for(int j = 0; j < num_comp; j++){
      sum1 = sum1 + x_i[j]*g_jk[i+num_comp*j]*beta_jk[i+num_comp*j];
      sum2 = sum2 + x_i[j]*beta_jk[i+num_comp*j];
    }
    G_E = G_E + x_i[i] * sum1/sum2;
  }
  return G_E;
}
double MultiComponentPengRobinson::b_mix(std::vector<double> x_i){
  double b_m = 0;
  for(int i = 0; i < num_comp; i++){
    b_m = b_m + b[i]*x_i[i];
  }
  return b_m;
}
double MultiComponentPengRobinson::a_mix(std::vector<double> x_i, std::vector<double> a_i, double b_m, double G_E){
  double a_sum = 0;
  for(int i = 0; i < num_comp; i++){
    a_sum = a_sum + x_i[i]*a_i[i]/b[i];
  }
  /*double a_m = 0;
  for(int i = 0; i < num_comp; i++){
    a_i[i] = a_c[i] * util::pow((1 + m[i] * (1 - util::pow((_T/T_c[i]),0.5))),2);
    for(int j = 0; j < num_comp; j++){
      a_m += x_i[i]*x_i[j]*util::pow(a_i[i]*a_i[j],0.5);
    }
  }
  return a_m;*/
  return b_m * (a_sum - G_E / lambda);
}
double MultiComponentPengRobinson::p_PR(double v, double a_m, double b_m){
  return R*T/(v - b_m) - a_m/(util::pow(v,2) + 2*v*b_m - util::pow(b_m,2));
}
double MultiComponentPengRobinson::gamma_i(std::vector<double> x_i, double a_i, std::vector<double> g_jk, std::vector<double> beta_jk, double _T, int i){
  double sum1 = 0;
  double sum2 = 0;
  double sum3 = 0;
  double sum4 = 0;
  for(int j = 0; j < num_comp; j++){
    sum1 = sum1 + x_i[j]*g_jk[num_comp*j+i]*beta_jk[num_comp*j+i];
    sum2 = sum2 + x_i[j]*beta_jk[num_comp*j+i];
    double sub31 = 0;
    double sub41 = 0;
    double sub42 = 0;
    for(int k = 0; k < num_comp; k++){
      sub31 = sub31 + x_i[k]*beta_jk[j+num_comp*k];
      sub41 = sub41 + x_i[k]*g_jk[j+num_comp*k]*beta_jk[j+num_comp*k];
      sub42 = sub42 + x_i[k]*beta_jk[j+num_comp*k];
    }
    sum3 = sum3 + x_i[j]*g_jk[num_comp*i+j]*beta_jk[num_comp*i+j] / sub31;
    sum4 = sum4 + x_i[j]*beta_jk[num_comp*i+j]*sub41/util::pow(sub42,2);
  }
  return 1/(R*_T) * (a_i / b[i] - 1/lambda * (sum1/sum2 + sum3 - sum4));
}
double MultiComponentPengRobinson::lnfugacity_i(std::vector<double> x_i, double v, double a_m, double b_m, double gamma_i, double _T, int i){
  //TODO change inputs for non-isothermal cases!!!
  double term1 = b[i]/(v-b_m) + util::log(x_i[i]*R*_T/(v-b_m));
  double term2 = -gamma_i / util::pow(2,1.5) * util::log((v + (util::pow(2,0.5)+1)*b_m)/(v-(util::pow(2,0.5)-1)*b_m));
  double term3 = -b[i]*a_m/(util::pow(2,1.5)*b_m*R*_T)*
          ((util::pow(2,0.5)+1)/(v+(util::pow(2,0.5)+1)*b_m)+(util::pow(2,0.5)-1)/(v-(util::pow(2,0.5)-1)*b_m));
  /*std::vector<double> a_i(num_comp);
  a_m = 0;
  for(int i = 0; i < num_comp; i++){
    a_i[i] = a_c[i] * util::pow((1 + m[i] * (1 - util::pow((_T/T_c[i]),0.5))),2);
    for(int j = 0; j < num_comp; j++){
      a_m += x_i[i]*x_i[j]*util::pow(a_i[i]*a_i[j],0.5);
    }
  }
  gamma_i = -a_m/(b_m*R*_T)*(b[i]/b_m-2/a_m*(x_i[0]*util::pow(a_i[0]*a_i[i],0.5)+x_i[1]*util::pow(a_i[1]*a_i[i],0.5)));
  term1 = b[i]/(v-b_m) + util::log(x_i[i]*R*_T/(v-b_m));
  term2 = -gamma_i / util::pow(2,1.5) * util::log((v + (util::pow(2,0.5)+1)*b_m)/(v-(util::pow(2,0.5)-1)*b_m));
  term3 = -b[i]*a_m/(util::pow(2,1.5)*b_m*R*_T)*((util::pow(2,0.5)+1)/(v+(util::pow(2,0.5)+1)*b_m)+(util::pow(2,0.5)-1)/(v-(util::pow(2,0.5)-1)*b_m));*/
  return term1 + term2 + term3;
}
std::vector<double> MultiComponentPengRobinson::iterate_VLE(double residuum, double beta0){
  std::vector<double> a_i(num_comp), g_jk(num_comp*num_comp), beta_jk(num_comp*num_comp);
  for(int i = 0; i < num_comp; i++){
    a_i[i] = a_c[i] * util::pow((1 + m[i] * (1 - util::pow((T/T_c[i]),0.5))),2);
    for(int j = 0; j < num_comp; j++){
      g_jk[num_comp*i+j] = g_I[num_comp*i+j] + g_II[num_comp*i+j] * T;
      beta_jk[num_comp*i+j] = b[j] * util::exp(-alpha[num_comp*i+j]*g_jk[num_comp*i+j]/(R*T));
    }
  }
  std::vector<double> vx(2+2*num_comp), K(num_comp), x(num_comp), y(num_comp), G_E(2), b_m(2), a_m(2), gammas(2*num_comp);
  std::vector<double> deltaLogFugacity(num_comp), errorFugacity(num_comp);
  double totalError = 0;
  for(int i = 0; i < num_comp; i++){
    //adjust prefactor or initial values for partition coefficients, 0.9 for hydrocarbons
    K[i] = 1.0*(p_c[i]/(p))*util::exp(5.373*(1+omega[i])*(1-T_c[i]/T));
    errorFugacity[i] = 1;
  }
  for(int i = 0; i < 2; i++){
    vx[i] = init_v[i];
  }
  double beta = beta0;
  double f = 0;
  do{
    totalError = 0;
    f = 0;
    for(int i = 0; i < num_comp; i++){
      f += (K[i]-1)*z[i]/(1+(K[i]-1)*beta);
    }
    while(abs(f) > residuum){
      double dfdbeta = 0;
      for(int i = 0; i < num_comp; i++){
        dfdbeta += -1*(K[i]-1)*(K[i]-1)*z[i]/((1+(K[i]-1)*beta)*(1+(K[i]-1)*beta));
      }
      beta = beta - f / dfdbeta;
      f = 0;
      for(int i = 0; i < num_comp; i++){
        f += (K[i]-1)*z[i]/(1+(K[i]-1)*beta);
      }
    }
    double x_sum = 0;
    double y_sum = 0;
    for(int i = 0; i < num_comp; i++){
      x[i] = z[i]/(1+(K[i]-1)*beta);
      y[i] = K[i]*x[i];
      x_sum += x[i];
      y_sum += y[i];
    }
    if(x_sum != 1 || y_sum != 1){
      for(int i = 0; i < num_comp; i++){
        vx[i+2] = x[i]/x_sum;
        x[i] = vx[i+2];
        vx[i+2+num_comp] = y[i]/y_sum;
        y[i] = vx[i+2+num_comp];
      }
    }

    G_E[0] = G_Excess(x, g_jk, beta_jk);
    G_E[1] = G_Excess(y, g_jk, beta_jk);
    b_m[0] = b_mix(x);
    b_m[1] = b_mix(y);
    a_m[0] = a_mix(x, a_i, b_m[0], G_E[0]);
    a_m[1] = a_mix(y, a_i, b_m[1], G_E[1]);
    for(int i = 0; i < num_comp; i++){
      gammas[i] = gamma_i(x, a_i[i], g_jk, beta_jk, T, i);
      gammas[i+num_comp] = gamma_i(y, a_i[i], g_jk, beta_jk, T, i);
    }

    for(int i = 0; i < 2; i++){
      while(abs(p-p_PR(vx[i], a_m[i], b_m[i])) > residuum){
        double dpdv = -R*T/util::pow(vx[i]-b_m[i],2) +
                2*a_m[i]*(vx[i]+b_m[i])/util::pow(util::pow(vx[i],2) + 2*vx[i]*b_m[i] - util::pow(b_m[i],2),2);
        vx[i] = vx[i] - (p - p_PR(vx[i], a_m[i], b_m[i])) / (-dpdv);
      }
    }
    for(int i = 0; i < num_comp; i++){
      deltaLogFugacity[i] = lnfugacity_i(x, vx[0], a_m[0], b_m[0], gammas[i], T, i) -
              lnfugacity_i(y, vx[1], a_m[1], b_m[1], gammas[i+num_comp], T, i);
      errorFugacity[i] = util::exp(deltaLogFugacity[i]) - 1;
      K[i] = util::exp(deltaLogFugacity[i]) / (x[i]/y[i]);
      totalError += util::pow(errorFugacity[i],2);
    }
  }while(util::pow(totalError,0.5) > residuum);
  for(int i = 1; i<num_comp; i++){
    volatilities[i] = (M[i]*(vx[i+2]/vx[0] - vx[i+2+num_comp]/vx[1])) / (M[0]*(vx[2]/vx[0] - vx[2+num_comp]/vx[1]));
  }
  double vol_sum = 0;
  for(int i = 0; i<num_comp; i++){
    vol_sum += volatilities[i];
  }
  for(int i = 0; i<num_comp; i++){
    chis[i] = volatilities[i] / vol_sum;
  }
  return vx;
}
std::vector<double> MultiComponentPengRobinson::getChis(int n){
  std::vector<double> n_chis(n);
  for(int i = 0; i<n; i++){
    n_chis[i] = chis[i];
  }
  return n_chis;
}

} // end namespace olb

#endif
