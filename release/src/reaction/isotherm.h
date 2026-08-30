/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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


#ifndef OLB_APPS_FLORIAN_ADSORPTION3D_PHOSPHATEREACTION_H_
#define OLB_APPS_FLORIAN_ADSORPTION3D_PHOSPHATEREACTION_H_

#include <vector>
#include <cmath>
#include <iostream>

namespace olb {

/**
 * @brief This is a new version of Isotherm made to be
 * suitable with the template format
 *
 * Returns loading in mg/g corresponding to concentration in mg/mL
 *
 * @param ISO_KONST_A first constant
 * @param ISO_KONST_B second constant or exponent
 */
namespace Isotherm {

// Linear Isotherm
struct LinearIsotherm {

  struct ISO_KONST_A : public descriptors::FIELD_BASE<1> { };
  struct ISO_KONST_B : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<ISO_KONST_A,ISO_KONST_B>;

  template <typename V, typename COUPLING>
  static void setParameters(V isoKonstA, V isoKonstB, COUPLING& coupling) any_platform {
    coupling.template setParameter<Isotherm::LinearIsotherm::ISO_KONST_A>(isoKonstA);
    coupling.template setParameter<Isotherm::LinearIsotherm::ISO_KONST_B>(isoKonstB);
    return;
  }

  template <typename V, typename COUPLING>
  static V getLoadingFromCoupling(V c, COUPLING& coupling) any_platform {
    auto isoKonstA = coupling.template getParameter<Isotherm::LinearIsotherm::ISO_KONST_A>();
    return isoKonstA[0] * c;
  }

  template <typename V, typename PARAMETERS>
  V getLoading(V c, PARAMETERS& params) any_platform {
    V isoKonstA = params.template get<ISO_KONST_A>();
    return isoKonstA * c;
  }

  template <typename V, typename COUPLING>
  static void print(std::ostream& clout, COUPLING& coupling) {
    auto isoKonstA = coupling.template getParameter<Isotherm::LinearIsotherm::ISO_KONST_A>();
    auto isoKonstB = coupling.template getParameter<Isotherm::LinearIsotherm::ISO_KONST_B>();
    clout << "----------------- Isotherm information -----------------" << std::endl;
    clout << "Isotherm exponent    n = " << isoKonstB << std::endl;
    clout << "Isotherm factor      K = " << isoKonstA << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;
  }

};

// LangmuirIsotherm
struct LangmuirIsotherm {

  struct ISO_KONST_A : public descriptors::FIELD_BASE<1> { };
  struct ISO_KONST_B : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<ISO_KONST_A,ISO_KONST_B>;

  template <typename V, typename COUPLING>
  static void setParameters(V isoKonstA, V isoKonstB, COUPLING& coupling) any_platform {
    coupling.template setParameter<Isotherm::LangmuirIsotherm::ISO_KONST_A>(isoKonstA);
    coupling.template setParameter<Isotherm::LangmuirIsotherm::ISO_KONST_B>(isoKonstB);
    return;
  }

  template <typename V, typename COUPLING>
  static V getLoadingFromCoupling(V c, COUPLING& coupling) any_platform {
    auto isoKonstA = coupling.template getParameter<Isotherm::LangmuirIsotherm::ISO_KONST_A>();
    auto isoKonstB = coupling.template getParameter<Isotherm::LangmuirIsotherm::ISO_KONST_B>();
    return V(isoKonstA[0]) * V(isoKonstB[0]) * c / ( 1 + V(isoKonstB[0]) * c );
  }

  template <typename V, typename PARAMETERS>
  V getLoading(V c, PARAMETERS& params) any_platform {
    V isoKonstA = params.template get<ISO_KONST_A>();
    V isoKonstB = params.template get<ISO_KONST_B>();
    return isoKonstA * isoKonstB * c / ( 1 + isoKonstB * c );
  }

  template <typename V, typename COUPLING>
  static void print(std::ostream& clout, COUPLING& coupling) {
    auto isoKonstA = coupling.template getParameter<Isotherm::LangmuirIsotherm::ISO_KONST_A>();
    auto isoKonstB = coupling.template getParameter<Isotherm::LangmuirIsotherm::ISO_KONST_B>();
    clout << "----------------- Isotherm information -----------------" << std::endl;
    clout << "Isotherm exponent    n = " << isoKonstB << std::endl;
    clout << "Isotherm factor      K = " << isoKonstA << std::endl;
    clout << "-------------------------------------------------------------" << std::endl;
  }

};

// FreundlichIsotherm
struct FreundlichIsotherm {

  struct ISO_KONST_A : public descriptors::FIELD_BASE<1> { };
  struct ISO_KONST_B : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<ISO_KONST_A,ISO_KONST_B>;

  template <typename V, typename PARAMETERS>
  V getLoading(V c, PARAMETERS& params) any_platform {
    V isoKonstA = params.template get<ISO_KONST_A>();
    V isoKonstB = params.template get<ISO_KONST_B>();
    return isoKonstA * std::pow(c, isoKonstB);
  }

};

}

}
#endif
