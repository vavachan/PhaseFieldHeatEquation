/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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


#ifndef NAMES_H
#define NAMES_H


#include "descriptor/tag.h"


namespace olb {

/// Define names as empty structs in order to enable calls like
/// lattice(NavierStokes()).
namespace names {

struct A { };
struct B { };
struct C { };
struct D { };

struct Component1 { };
struct Component2 { };
struct Component3 { };
struct Component4 { };

struct Lattice : public descriptors::DESCRIPTOR_TAG { };
struct Lattice2 : public descriptors::DESCRIPTOR_TAG { };
struct Primal : public descriptors::DESCRIPTOR_TAG { };
struct Dual : public descriptors::DESCRIPTOR_TAG { };
struct NavierStokes : public descriptors::DESCRIPTOR_TAG { };
struct AdvectionDiffusion : public descriptors::DESCRIPTOR_TAG { };
struct Temperature  : public descriptors::DESCRIPTOR_TAG { };
struct Points : public descriptors::DESCRIPTOR_TAG { };

struct TurbKineticEnergy : public descriptors::DESCRIPTOR_TAG { };
struct DissipationRate : public descriptors::DESCRIPTOR_TAG { };

template <unsigned DIM>
struct Concentration   : public descriptors::DESCRIPTOR_TAG { };

struct Concentration0  : public descriptors::DESCRIPTOR_TAG { };
struct Concentration1  : public descriptors::DESCRIPTOR_TAG { };
struct Concentration2  : public descriptors::DESCRIPTOR_TAG { };
struct VolumeRendering{};

struct Parameter { };
struct Opti                  : public Parameter { };
struct Output                : public Parameter { };
struct OutputOpti            : public Parameter { };
struct VisualizationVTK      : public Parameter { const std::string name {"VisualizationVTK"}; };
struct VisualizationGnuplot  : public Parameter { const std::string name {"VisualizationGnuplot"}; };
struct VisualizationImages   : public Parameter { const std::string name {"VisualizationImages"}; };
struct Simulation            : public Parameter { };
struct Stationarity          : public Parameter { };

struct Errors       : public Parameter { };
struct Results      : public Parameter { };
struct Reference    : public Parameter { };

struct OutputChannel { };
struct debug        : public OutputChannel { };
struct error        : public OutputChannel { };
struct file         : public OutputChannel { };
struct info         : public OutputChannel { };
struct performance  : public OutputChannel { };
struct results      : public OutputChannel { };

}


}




#endif
