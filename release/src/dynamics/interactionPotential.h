/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Czelusniak,
 *                     Maximilian Schecher, Adrian Kummerlaender
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

#ifndef INTERACTION_POTENTIAL_H
#define INTERACTION_POTENTIAL_H


/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

//Calculate VLEs with PR EoS with non-classical Huron-Vidal mixing rule for n-component mixtures
class MultiComponentPengRobinson {
public:
  MultiComponentPengRobinson(double p_, double T_, std::vector<double> z_, std::vector<double> a_L, std::vector<double> b_L,
                             std::vector<double> M_L, std::vector<double> Tc_L, std::vector<double> pc_L, std::vector<double> omega_, std::vector<double> devi,
                             std::vector<double> alpha_, std::vector<double> gI_L, std::vector<double> gII_L);
  double G_Excess(std::vector<double> x_i, std::vector<double> g_jk, std::vector<double> beta_jk);
  double b_mix(std::vector<double> x_i);
  double a_mix(std::vector<double> x_i, std::vector<double> a_i, double b, double G_E);
  double p_PR(double v, double a_m, double b_m);
  double gamma_i(std::vector<double> x_i, double a_i, std::vector<double> g_jk, std::vector<double> beta_jk, double _T, int i);
  double lnfugacity_i(std::vector<double> x_i, double v, double a_m, double b_m, double gamma_i, double _T, int i);
  std::vector<double> iterate_VLE(double residuum, double beta0);
  std::vector<double> getChis(int n);
private:
  double p, T;
  std::vector<double> z;
  const double G{-1.};
  const double R{1.};
  const double lambda{0.6232252401};
  int num_comp;
  std::vector<double> a_c;
  std::vector<double> b;
  std::vector<double> M;
  std::vector<double> T_c;
  std::vector<double> p_c;
  std::vector<double> omega;
  std::vector<double> m;
  std::vector<double> alpha;
  std::vector<double> g_I;
  std::vector<double> g_II;
  std::vector<double> init_v = {0.105,120.0};
  std::vector<double> volatilities = {1.,0.,0.,0.};
  std::vector<double> chis = {0.5,0.5,0.5,0.5};
};

namespace interaction {

// established -- only multicomponent flow
struct PsiEqualsRho {
  using parameters = meta::list<>;

  template <typename V, typename PARAMETERS>
  V compute(V rho, PARAMETERS& params) any_platform {
    return rho;
  };
};

struct ShanChen94 {
  using parameters = meta::list<>;

  template <typename V, typename PARAMETERS>
  V compute(V rho, PARAMETERS& params) any_platform {
    return 4 * util::exp(-200 / rho);
  };
};

struct CarnahanStarling {
  struct G : public descriptors::FIELD_BASE<1> {
    template <typename T, typename DESCRIPTOR,typename FIELD>
    static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
      return value != 0;
    }
  };
  struct A : public descriptors::FIELD_BASE<1> {
    template <typename T, typename DESCRIPTOR,typename FIELD>
    static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
      return value > 0;
    }
  };
  struct B : public descriptors::FIELD_BASE<1> {
    template <typename T, typename DESCRIPTOR,typename FIELD>
    static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
      return value > 0;
    }
   };
  struct T : public descriptors::FIELD_BASE<1> {
    template <typename T, typename DESCRIPTOR,typename FIELD>
    static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
      return value > 0;
    }
   };

  using parameters = meta::list<G,A,B,T>;

  template <typename V, typename PARAMETERS>
  V compute(V rho, PARAMETERS& params) any_platform {
    V c = params.template get<B>() * rho / V{4};
    V p = rho
        * (V(0.18727 / 0.4963) * params.template get<A>() / params.template get<B>())
        * params.template get<T>()
        * ((V{1} +c+c*c-c*c*c) / (V{1}-c) / (V{1}-c) / (V{1}-c))
        - params.template get<A>() * rho*rho;
    return util::sqrt(V{6}*(p-rho/V{3})/ params.template get<G>());
  };
};

// Polinomial EOS
struct Polinomial {
  struct G              : public descriptors::FIELD_BASE<1> { };
  struct KAPPAP         : public descriptors::FIELD_BASE<1> { };

  struct RHOV           : public descriptors::FIELD_BASE<1> { };
  struct RHOL           : public descriptors::FIELD_BASE<1> { };
  struct THICKNESS      : public descriptors::FIELD_BASE<1> { };
  struct SURFTENSION    : public descriptors::FIELD_BASE<1> { };

  struct RHOC           : public descriptors::FIELD_BASE<1> { };
  struct PC             : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<G,KAPPAP,RHOV,RHOL,THICKNESS,SURFTENSION,RHOC,PC>;

  template <typename V, typename COUPLING>
  static void computeParameters(COUPLING& coupling) {
    // Set G
    const V g = V(-1.);
    coupling.template setParameter<interaction::Polinomial::G>(g);
    // Compute and set rhoc (critical density)
    auto rhov = coupling.template getParameter<interaction::Polinomial::RHOV>()[0];
    auto rhol = coupling.template getParameter<interaction::Polinomial::RHOL>()[0];
    const V rhoc = 0.5 * ( rhov + rhol );
    coupling.template setParameter<interaction::Polinomial::RHOC>(rhoc);
    // Compute p_c
    const V rho1 = rhov + V(0.33) * ( rhol - rhov );
    const V rho2 = rhov + V(0.67) * ( rhol - rhov );
    // Coefficients
    const V A = 1. / ( rhol - rhov );
    const V B = 1. / 2. / sqrt(rhol) / ( sqrt(rhol) - sqrt(rhov) );
    const V C = 1. / 2. / sqrt(rhol) / ( sqrt(rhov) + sqrt(rhol) );
    // f1
    const V f1 = A * log( sqrt(rho1) - sqrt(rhov) )
               - B * log( sqrt(rhol) - sqrt(rho1) )
               - C * log( sqrt(rho1) + sqrt(rhol) );
    // f2
    const V f2 = A * log( sqrt(rho2) - sqrt(rhov) )
               - B * log( sqrt(rhol) - sqrt(rho2) )
               - C * log( sqrt(rho2) + sqrt(rhol) );
    // p_c
    auto thickness = coupling.template getParameter<interaction::Polinomial::THICKNESS>()[0];
    const V pc = rhoc * rhoc * rhoc / V(6.) * pow( ( f2 - f1 ) / thickness, 2 );
    coupling.template setParameter<interaction::Polinomial::PC>(pc);
    // Compute kappap
    const V a = V(1./4.);
    const V b = rhol / V(2.);
    const V c = sqrt(rhov) / V(3.);
    const V d = sqrt(rhov) * rhol;
    // g1
    const V g1 = a * rhov * rhov - b * rhov - c * pow( rhov, V(1.5) ) + d * sqrt( rhov );
    // g2
    const V g2 = a * rhol * rhol - b * rhol - c * pow( rhol, V(1.5) ) + d * sqrt( rhol );
    // kappaLi
    auto gamma = coupling.template getParameter<interaction::Polinomial::SURFTENSION>()[0];
    const V kappaP = V(9.) * sqrt( rhoc * rhoc * rhoc ) *gamma
                    / g / sqrt( V(6.) * pc ) / ( g2 - g1 ) ;
    coupling.template setParameter<interaction::Polinomial::KAPPAP>(kappaP);
  }

  template <typename V, typename PARAMETERS>
  V compute( V Rho, PARAMETERS& params ) any_platform {
    V p = computeP( Rho, params );
    auto g = params.template get<G>(); // critical density
    V psi = util::sqrt( V(2.)*( p - Rho/V(3.) )/g );
    return psi;
  };

  template <typename V, typename PARAMETERS>
  V computeP(V Rho, PARAMETERS& params ) any_platform {
    auto rho_c = params.template get<RHOC>(); // critical density
    auto p_c = params.template get<PC>(); // critical pressure
    auto rho_v = params.template get<RHOV>(); // vapor density
    auto rho_l = params.template get<RHOL>(); // liquid density
    const V sqrt_rho_v = sqrt(rho_v);
    const V sqrt_Rho = sqrt(Rho);
    V p = p_c / rho_c / rho_c / rho_c * (
          V(2.) * Rho * Rho * Rho + ( rho_v - V(2.) * rho_l ) * Rho * Rho
          - V(3.) * sqrt_rho_v * Rho* Rho * sqrt_Rho + V(2.) * sqrt_rho_v * rho_l * Rho * sqrt_Rho
          + sqrt_rho_v * rho_l * rho_l * sqrt_Rho );
    return p;
  };
};

template <unsigned N_COMPONENTS>
struct MCPRpseudoPotential {
  template <typename X_I, typename G_JK, typename BETA_JK, typename V=typename X_I::value_t>
  V G_Excess(const X_I& x_i, const G_JK& g_jk, const BETA_JK& beta_jk) any_platform {
    V G_E = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      V sum1 = 0;
      V sum2 = 0;
      for(unsigned j = 0; j < N_COMPONENTS; j++){
        sum1 = sum1 + x_i[j]*g_jk[i+N_COMPONENTS*j]*beta_jk[i+N_COMPONENTS*j];
        sum2 = sum2 + x_i[j]*beta_jk[i+N_COMPONENTS*j];
      }
      G_E = G_E + x_i[i] * sum1/sum2;
    }
    return G_E;
  };

  template <typename X_I, typename B, typename V=typename X_I::value_t>
  V b_mix(const X_I& x_i, const B& b) any_platform {
    V b_m = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      b_m += b[i]*x_i[i];
    }
    return b_m;
  };

  template <typename X_I, typename A_I, typename B, typename V=typename X_I::value_t>
  V a_mix(const X_I& x_i, const A_I& a_i, const B& b, V b_m, V G_E) any_platform {
    V a_sum = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      a_sum = a_sum + x_i[i]*a_i[i]/b[i];
    }
    /*double a_m = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      for(unsigned j = 0; j < N_COMPONENTS; j++){
        a_m += x_i[i]*x_i[j]*util::pow(a_i[i]*a_i[j],0.5);
      }
    }
    return a_m;*/
    return b_m * (a_sum - G_E / 0.6232252401);
  };

  template <typename RHO_FIELD, typename __T, typename K, typename A_C, typename B, typename T_C, typename M, typename ALPHA, typename G_I, typename G_II, typename BIG_M, typename V=typename RHO_FIELD::value_t>
  V compute(const RHO_FIELD& rhoField, const __T& _T, const K& k, const A_C& a_c, const B& b, const T_C& T_c, const M& m, const ALPHA& alpha, const G_I& g_I, const G_II& g_II, const BIG_M& _M) any_platform {
    Vector<V,N_COMPONENTS> molarRhoField{}, a_i{}, x_i{};
    Vector<V,N_COMPONENTS*N_COMPONENTS> g_jk{}, beta_jk{};
    V R = 1;
    V molarRho = 0;
    V Rho = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      a_i[i] = a_c[i] * util::pow((1 + m[i] * (1 - util::pow((_T/T_c[i]),0.5))),2);
      for(unsigned j = 0; j < N_COMPONENTS; j++){
        g_jk[N_COMPONENTS*i+j] = g_I[N_COMPONENTS*i+j] + g_II[N_COMPONENTS*i+j] * _T;
        beta_jk[N_COMPONENTS*i+j] = b[j] * util::exp(-alpha[N_COMPONENTS*i+j]*g_jk[N_COMPONENTS*i+j]/(R*_T));
      }
      molarRhoField[i] = rhoField[i]/_M[i];
      molarRho += molarRhoField[i];
      Rho += rhoField[i];
    }
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      x_i[i] = molarRhoField[i]/molarRho;
    }
    V g_E = G_Excess(x_i, g_jk, beta_jk);
    V b_m = b_mix(x_i, b);
    V a_m = a_mix(x_i, a_i, b, b_m, g_E);
    V p = molarRho*R*_T/(1-b_m*molarRho) - a_m*molarRho*molarRho/(1+2*b_m*molarRho-b_m*molarRho*b_m*molarRho);
    V psi = util::sqrt(-6.*(k*p - Rho/3.));
    return psi;
  };

  template <typename RHO_FIELD, typename A_C, typename B, typename T_C, typename M, typename ALPHA, typename G_I, typename G_II, typename BIG_M, typename V=typename RHO_FIELD::value_t>
  V computeP(const RHO_FIELD& rhoField, V _T, const A_C& a_c, const B& b, const T_C& T_c, const M& m, const ALPHA& alpha, const G_I& g_I, const G_II& g_II, const BIG_M& _M) any_platform {
    olb::Vector<V,N_COMPONENTS> molarRhoField{}, a_i{}, x_i{};
    olb::Vector<V,N_COMPONENTS*N_COMPONENTS> g_jk{}, beta_jk{};
    V R = 1;
    V molarRho = 0;
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      a_i[i] = a_c[i] * util::pow((1 + m[i] * (1 - util::pow((_T/T_c[i]),0.5))),2);
      for(unsigned j = 0; j < N_COMPONENTS; j++){
        g_jk[N_COMPONENTS*i+j] = g_I[N_COMPONENTS*i+j] + g_II[N_COMPONENTS*i+j] * _T;
        beta_jk[N_COMPONENTS*i+j] = b[j] * util::exp(-alpha[N_COMPONENTS*i+j]*g_jk[N_COMPONENTS*i+j]/(R*_T));
      }
      molarRhoField[i] = rhoField[i]/_M[i];
      molarRho += molarRhoField[i];
    }
    for(unsigned i = 0; i < N_COMPONENTS; i++){
      x_i[i] = molarRhoField[i]/molarRho;
    }
    V g_E = G_Excess(x_i, g_jk, beta_jk);
    V b_m = b_mix(x_i, b);
    V a_m = a_mix(x_i, a_i, b, b_m, g_E);
    V p = molarRho*R*_T/(1-b_m*molarRho) - a_m*molarRho*molarRho/(1+2*b_m*molarRho-b_m*molarRho*b_m*molarRho);
    return p;
  }

};

}

#ifndef USING_LEGACY_CODEGEN

// established -- original for both single- and multicomponent flow

template <typename T, typename S>
class ShanChen93 : public AnalyticalF<1,T,S> {
private:
  T _rhoZero;
public:
  ShanChen93(T rhoZero=1.);
  bool operator() (T psi[], const S rho[]);
};

// established -- only multicomponent flow

template <typename T, typename S>
class PsiEqualsRho : public AnalyticalF<1,T,S> {
private:
public:
  PsiEqualsRho();
  bool operator() (T psi[], const S rho[]) override;
};

// established -- only singlecomponent flow

template <typename T, typename S>
class ShanChen94 : public AnalyticalF<1,T,S> {
private:
  T _rhoZero;
  T _psiZero;
public:
  ShanChen94(T rhoZero=200., T psiZero=4.);
  bool operator() (T psi[], const S rho[]) override;
};

template <typename T, typename S>
class PengRobinson : public AnalyticalF<1,T,S> {
private:
  T _G;
  T _acentricFactor;
  T _a;
  T _b;
  T _R;
  T _alpha;
  T _t;
  T _tc;
public:
  PengRobinson(T G, T acentricFactor=0.334, T a=2./49., T b=2./21., T tr=.8);
  bool operator() (T psi[], const S rho[]);
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S t[]);
};

template <typename T, typename S>
class CarnahanStarling : public AnalyticalF<1,T,S> {
private:
  T _G;
  T _a;
  T _b;
  T _R;
  T _t;
public:
  CarnahanStarling(T G, T a=1., T b=4., T tr=.7);
  bool operator() (T psi[], const S rho[]) override;
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S t[]);
};

// under development -- for singlecomponent flow

// 0.5 -> psiZero=0.65
// 1 -> psiZero=1.9
// 1.5 -> psiZero=3.5
// 2. -> psiZero=5,45
template <typename T, typename S>
class Krause : public AnalyticalF<1,T,S> {
private:
  T _rhoZero;
  T _psiZero;
public:
  Krause(T rhoZero=1., T psiZero=1.9);
  bool operator() (T psi[], const S rho[]);
};

template <typename T, typename S> // density of liquid phase always equals rhoZero for G=-1
class WeisbrodKrause : public AnalyticalF<1,T,S> {
private:
  T _rhoZero;
  T _sigmu;
public:
  WeisbrodKrause(T rhoZero=1., T sigmu=1.);
  bool operator() (T psi[], const S rho[]);
};

template <typename T, typename S> // not very good
class Normal : public AnalyticalF<1,T,S> {
private:
  T _sigma;
  T _mu;
public:
  Normal(T sigma=1., T mu=1.);
  bool operator() (T psi[], const S rho[]);
};

#endif // not USING_LEGACY_CODEGEN

} // end namespace olb

#endif

#include "interactionPotential.cse.h"
