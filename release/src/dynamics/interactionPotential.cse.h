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

#ifndef DYNAMICS_INTERACTION_POTENTIAL_CSE_H
#define DYNAMICS_INTERACTION_POTENTIAL_CSE_H


#ifndef DISABLE_CSE

namespace olb {

namespace interaction {

template <>
struct MCPRpseudoPotential<3> {

template <typename _RHOFIELD, typename __T, typename _K, typename _A_C, typename _B, typename _T_C, typename _M, typename _ALPHA, typename _G_I, typename _G_II, typename _BIGM, typename V=typename _RHOFIELD::value_t>
auto compute(const _RHOFIELD& rhoField, const __T& _T, const _K& k, const _A_C& a_c, const _B& b, const _T_C& T_c, const _M& m, const _ALPHA& alpha, const _G_I& g_I, const _G_II& g_II, const _BIGM& bigM) any_platform
{
auto x0 = rhoField[0]/bigM[0];
auto x1 = x0*b[0];
auto x2 = rhoField[1]/bigM[1];
auto x3 = x2*b[1];
auto x4 = rhoField[2]/bigM[2];
auto x5 = x4*b[2];
auto x6 = x1 + x3 + x5;
auto x7 = _T*g_II[0] + g_I[0];
auto x8 = V{1}/_T;
auto x9 = x0*util::exp(-x7*x8*alpha[0]);
auto x10 = _T*g_II[3] + g_I[3];
auto x11 = x2*util::exp(-x10*x8*alpha[3]);
auto x12 = _T*g_II[6] + g_I[6];
auto x13 = x4*util::exp(-x12*x8*alpha[6]);
auto x14 = _T*g_II[1] + g_I[1];
auto x15 = x0*util::exp(-x14*x8*alpha[1]);
auto x16 = _T*g_II[4] + g_I[4];
auto x17 = x2*util::exp(-x16*x8*alpha[4]);
auto x18 = _T*g_II[7] + g_I[7];
auto x19 = x4*util::exp(-x18*x8*alpha[7]);
auto x20 = _T*g_II[2] + g_I[2];
auto x21 = x0*util::exp(-x20*x8*alpha[2]);
auto x22 = _T*g_II[5] + g_I[5];
auto x23 = x2*util::exp(-x22*x8*alpha[5]);
auto x24 = _T*g_II[8] + g_I[8];
auto x25 = x4*util::exp(-x24*x8*alpha[8]);
return V{2.44948974278318}*util::sqrt(k*(V{1}*_T*(x0 + x2 + x4)/(x6 + V{-1}) + x6*(x0*((-(util::sqrt(_T/T_c[0]) + V{-1})*m[0] + V{1})*(-(util::sqrt(_T/T_c[0]) + V{-1})*m[0] + V{1}))*a_c[0]/b[0] - V{1.60455694171447}*x0*(x10*x11 + x12*x13 + x7*x9)/(x11 + x13 + x9) + x2*((-(util::sqrt(_T/T_c[1]) + V{-1})*m[1] + V{1})*(-(util::sqrt(_T/T_c[1]) + V{-1})*m[1] + V{1}))*a_c[1]/b[1] - V{1.60455694171447}*x2*(x14*x15 + x16*x17 + x18*x19)/(x15 + x17 + x19) + x4*((-(util::sqrt(_T/T_c[2]) + V{-1})*m[2] + V{1})*(-(util::sqrt(_T/T_c[2]) + V{-1})*m[2] + V{1}))*a_c[2]/b[2] - V{1.60455694171447}*x4*(x20*x21 + x22*x23 + x24*x25)/(x21 + x23 + x25))*V{1} / (V{2}*x1 + V{2}*x3 + V{2}*x5 + V{1} - (x6*x6))) + V{0.333333333333333}*rhoField[0] + V{0.333333333333333}*rhoField[1] + V{0.333333333333333}*rhoField[2]);
}

template <typename _RHOFIELD, typename __T, typename _A_C, typename _B, typename _T_C, typename _M, typename _ALPHA, typename _G_I, typename _G_II, typename _BIGM, typename V=typename _RHOFIELD::value_t>
auto computeP(const _RHOFIELD& rhoField, const __T& _T, const _A_C& a_c, const _B& b, const _T_C& T_c, const _M& m, const _ALPHA& alpha, const _G_I& g_I, const _G_II& g_II, const _BIGM& bigM) any_platform
{
auto x0 = rhoField[0]/bigM[0];
auto x1 = x0*b[0];
auto x2 = rhoField[1]/bigM[1];
auto x3 = x2*b[1];
auto x4 = rhoField[2]/bigM[2];
auto x5 = x4*b[2];
auto x6 = x1 + x3 + x5;
auto x7 = _T*g_II[0] + g_I[0];
auto x8 = V{1}/_T;
auto x9 = x0*util::exp(-x7*x8*alpha[0]);
auto x10 = _T*g_II[3] + g_I[3];
auto x11 = x2*util::exp(-x10*x8*alpha[3]);
auto x12 = _T*g_II[6] + g_I[6];
auto x13 = x4*util::exp(-x12*x8*alpha[6]);
auto x14 = _T*g_II[1] + g_I[1];
auto x15 = x0*util::exp(-x14*x8*alpha[1]);
auto x16 = _T*g_II[4] + g_I[4];
auto x17 = x2*util::exp(-x16*x8*alpha[4]);
auto x18 = _T*g_II[7] + g_I[7];
auto x19 = x4*util::exp(-x18*x8*alpha[7]);
auto x20 = _T*g_II[2] + g_I[2];
auto x21 = x0*util::exp(-x20*x8*alpha[2]);
auto x22 = _T*g_II[5] + g_I[5];
auto x23 = x2*util::exp(-x22*x8*alpha[5]);
auto x24 = _T*g_II[8] + g_I[8];
auto x25 = x4*util::exp(-x24*x8*alpha[8]);
return -(V{1}*_T*(x0 + x2 + x4)/(x6 + V{-1}) + x6*(x0*((-(util::sqrt(_T/T_c[0]) + V{-1})*m[0] + V{1})*(-(util::sqrt(_T/T_c[0]) + V{-1})*m[0] + V{1}))*a_c[0]/b[0] - V{1.60455694171447}*x0*(x10*x11 + x12*x13 + x7*x9)/(x11 + x13 + x9) + x2*((-(util::sqrt(_T/T_c[1]) + V{-1})*m[1] + V{1})*(-(util::sqrt(_T/T_c[1]) + V{-1})*m[1] + V{1}))*a_c[1]/b[1] - V{1.60455694171447}*x2*(x14*x15 + x16*x17 + x18*x19)/(x15 + x17 + x19) + x4*((-(util::sqrt(_T/T_c[2]) + V{-1})*m[2] + V{1})*(-(util::sqrt(_T/T_c[2]) + V{-1})*m[2] + V{1}))*a_c[2]/b[2] - V{1.60455694171447}*x4*(x20*x21 + x22*x23 + x24*x25)/(x21 + x23 + x25))*V{1} / (V{2}*x1 + V{2}*x3 + V{2}*x5 + V{1} - (x6*x6)));
}

};

}

}

#endif

#endif