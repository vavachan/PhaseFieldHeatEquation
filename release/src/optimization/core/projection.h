/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-13, 2016, 2023-24 Mathias J. Krause, Lukas Baron,
 *  Benjamin Förster, Julius Jeßberger
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


// This file contains projections which can be used to rescale control variables.

#ifndef PROJECTION_H
#define PROJECTION_H

#include <regex>

#include "core/unitConverter.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {

enum StartValueType {Control, ProjectedControl, Porosity, Permeability};

namespace projection {

template <typename T>
struct Base {

  virtual T project(T) const = 0;
  virtual T derivative(T) const = 0;
  virtual T inverse(T) const = 0;
};

template <typename T>
struct Identity : public Base<T> {

  T project(T x) const override { return x; }
  T derivative(T x) const override { return 1; }
  T inverse(T x) const override { return x; }
};

template <typename T>
struct Sigmoid : public Base<T> {

  T project(T x) const override { return util::exp(x) / ( util::exp(x) + T(1) ); }
  T derivative(T x) const override { return util::exp(x) / util::pow( util::exp(x) + T(1), 2 ); }
  T inverse(T x) const override { return -util::log(T(1)/x - T(1)); }
};

template <typename T>
struct SigmoidD : public Base<T> {

  T project(T x) const override { return util::exp(x) / util::pow( util::exp(x) + T(1), 2 ); }
  T derivative(T x) const override { return 0; }
  T inverse(T x) const override { return 0; }
};

/// Convert force to lattice units
template <typename T>
struct ForceFactor : public Base<T> {

  const T scale;

  template <typename DESCRIPTOR>
  ForceFactor(const UnitConverter<T,DESCRIPTOR>& converter)
   : scale {converter.getConversionFactorMass() / converter.getConversionFactorForce()}
  { }

  T project(T x) const override { return scale * x; }
  T derivative(T x) const override { return scale; }
  T inverse(T x) const override { return x / scale; }
};

template <typename T>
struct ForceFactorD : public Base<T> {

  const T scale;

  template <typename DESCRIPTOR>
  ForceFactorD(const UnitConverter<T,DESCRIPTOR>& converter)
   : scale {converter.getConversionFactorMass() / converter.getConversionFactorForce()}
  { }

  T project(T x) const override { return scale; }
  T derivative(T x) const override { return 0; }
  T inverse(T x) const override { return 0; }
};

template <typename T>
struct Rectifier : public Base<T> {

  T project(T x) const override { return T(1) - T(1)/util::max( T(1), x ); }
  T derivative(T x) const override { return (x >= T(1) ? T(1)/(x*x) : T(0)); }
  T inverse(T x) const override { return T(1)/(T(1)-x); }
};

template <typename T>
struct Softplus : public Base<T> {

  T project(T x) const override { return T(1) - T(1)/(util::log(T(1)+util::exp(x))+T(1)); }
  T derivative(T x) const override { return util::exp(x)
     / ((util::exp(x)+T(1))*(util::log(T(1)+util::exp(x))+T(1))*(util::log(T(1)+util::exp(x))+T(1))); }
  T inverse(T x) const override { return util::log(util::exp(x/(T(1)-x))-T(1)); }
};

template <typename T>
struct Baron : public Base<T> {

  const T pi      {T(4) * util::atan(T(1))};

  T project(T x) const override {
    const T y = T(0.5) + T(0.5) * util::sin((x-T(0.5))*pi);
    return y*y;
  }
  T derivative(T x) const override {
    return (T(1) + util::sin((x-T(0.5))*pi))
                     * util::cos((x-T(0.5))*pi)
                     * pi * T(0.5);
  }
  T inverse(T x) const override {
    return T(0.5) + util::asin(T(2)*x-T(1)) / pi;
  }
};

template <typename T>
struct Krause : public Base<T> {

  T project(T x) const override {
    return T(1) - T(1) / util::exp(x*x);
  }
  T derivative(T x) const override {
    return T(2)*x / util::exp(x*x);
  }
  T inverse(T x) const override {
    return util::sqrt(util::log(T(1)/(T(1)-x)));
  }
};


// ------------------------------------------ //
//             helper functions               //
// ------------------------------------------ //

template<typename T, typename DESCRIPTOR>
T gridTerm(const UnitConverter<T,DESCRIPTOR>& converter) {
  return converter.getPhysDeltaX()
    * converter.getPhysDeltaX()
    * converter.getLatticeViscosity()
    * converter.getLatticeRelaxationTime();
}

/// Get control value for given porosity
// Inverse projection function to get startValue from porosity
template<typename T>
T porosityToControl(T porosity, const Base<T>& projection)
{
  return projection.inverse(porosity);
}

/// Get porosity for given permeability
template<typename T, typename DESCRIPTOR>
T permeabilityToPorosity(T permeability,
  const UnitConverter<T,DESCRIPTOR>& converter)
{
  return (T(1) - gridTerm(converter) / permeability);
}

/// Get permeability for given porosity
template<typename T, typename DESCRIPTOR>
T porosityToPermeability(T porosity,
  const UnitConverter<T,DESCRIPTOR>& converter)
{
  return gridTerm(converter) / (T(1) - porosity);
}

/// Get inverse of permeability for given porosity
template<typename T, typename DESCRIPTOR>
T porosityToInvPermeability(T porosity,
  const UnitConverter<T,DESCRIPTOR>& converter)
{
  return (T(1) - porosity) / gridTerm(converter);
}

/// Get control for given permeability
template<typename T, typename DESCRIPTOR>
T permeabilityToControl(T permeability,
  const Base<T>& projection,
  const UnitConverter<T,DESCRIPTOR>& converter)
{
  const T porosity = permeabilityToPorosity(permeability, converter);
  return porosityToControl(porosity);
}

/// Transform porosity/ permeability/ other startValue into control
template<typename T, typename DESCRIPTOR>
T getInitialControl(T startValue,
  const Base<T>& projection,
  const UnitConverter<T,DESCRIPTOR>& converter,
  StartValueType type,
  bool verbose = true)
{
  T control {0};
  OstreamManager clout(std::cout, "getInitialControl");

  if (type == Porosity) {
    control = projection::porosityToControl(startValue, projection);
    if (verbose) {
      clout << "Transform porosity startValue into control startValue:\n";
      clout << "Porosity: " << startValue << std::endl;
      clout << "Control: " << control << std::endl;
    }
  }
  else if (type == Permeability) {
    const T tmpPorosity = projection::permeabilityToPorosity(startValue, converter);
    control = projection::porosityToControl(tmpPorosity, projection);
    if (verbose) {
      clout << "Transform permeability startValue into control startValue:\n";
      clout << "Permeability: " << startValue << std::endl;
      clout << "Porosity: " << tmpPorosity << std::endl;
      clout << "Control: " << control << std::endl;
    }
  }
  else if (type == ProjectedControl) {
    control = projection.inverse(startValue);
  } else if (type == Control) {
    control = startValue;
  } else {
    throw std::invalid_argument(
      "Error: unknown start value type.");
  }
  return control;
}

template<typename T, typename OptiCaseDual_Type>
T getInitialControl(T startValue, OptiCaseDual_Type& optiCase)
{
  return getInitialControl(startValue,
    *(optiCase._projection),
    *(optiCase._converter),
    optiCase._startValueType,
    optiCase._verbose);
}

// ------------------------------------------ //
//       grid-independent projections         //
// ------------------------------------------ //

/// Gridterm-dependent projection base class
/**
 * Projection = 1 - gridTerm / K(x)
 * with K(x) = subprojection(x) + gridTerm
 */
template <typename T, typename DESCRIPTOR>
struct GiBase : public Base<T> {

  const T _gridTerm;

  GiBase(const UnitConverter<T,DESCRIPTOR>& converter)
   : _gridTerm {gridTerm(converter)}
  { }

  virtual T subprojection(T x) const = 0;
  virtual T derivSubprojection(T x) const = 0;
  virtual T inverseSubprojection(T x) const = 0;

  T project(T x) const override {
    return subprojection(x) / (subprojection(x) + _gridTerm);
  }
  T derivative(T x) const override {
    const T y = subprojection(x) + _gridTerm;
    return _gridTerm*derivSubprojection(x) / (y*y);
  }
  T inverse(T x) const override {
    return inverseSubprojection(_gridTerm*x/(T(1)-x));
  }
};

template <typename T, typename DESCRIPTOR>
struct Foerster : public GiBase<T,DESCRIPTOR> {

  Foerster(const UnitConverter<T,DESCRIPTOR>& converter)
   : GiBase<T,DESCRIPTOR>{converter}
  { }

  T subprojection(T x) const override { return util::exp(x); }
  T derivSubprojection(T x) const override { return util::exp(x); }
  T inverseSubprojection(T x) const override { return util::log(x); }
};

/// FoersterProjection for arbitrary n
/**
 * subproj(a) = util::exp(a^(2n)) - 1
 */
template <typename T, typename DESCRIPTOR>
struct FoersterN : public GiBase<T,DESCRIPTOR> {

  const unsigned _n;

  FoersterN(const UnitConverter<T,DESCRIPTOR>& converter, unsigned n)
   : GiBase<T,DESCRIPTOR>(converter), _n(n)
  {
    OLB_PRECONDITION(_n >= 1);
  }

  T subprojection(T x) const override { return util::exp(util::pow(x, 2*_n)) - T(1); }
  T derivSubprojection(T x) const override {
    return T(2*_n) * util::pow(x, 2*_n-1) * subprojection(x);
  }
  T inverseSubprojection(T x) const override {
    return util::pow(util::log(x+T(1)), T(1)/T(2*_n));
  }
};

/// StasiusProjection for arbitrary n
/**
 * subproj(a) = a^(2n)
 */
template <typename T, typename DESCRIPTOR>
struct StasiusN : public GiBase<T,DESCRIPTOR> {

  const unsigned _n;

  StasiusN(const UnitConverter<T,DESCRIPTOR>& converter, unsigned n)
   : GiBase<T,DESCRIPTOR>(converter), _n(n)
  {
    OLB_PRECONDITION(_n >= 1);
  }

  T subprojection(T x) const override { return util::pow(x, 2*_n); }
  T derivSubprojection(T x) const override { return T(2*_n) * util::pow(x, 2*_n-1); }
  T inverseSubprojection(T x) const override {
    return util::pow(x, T(1)/T(2*_n));
  }
};



template <typename T, typename DESCRIPTOR>
std::shared_ptr<projection::Base<T>> construct(
  const UnitConverter<T,DESCRIPTOR>& converter, std::string& name)
{
  std::smatch match;
  if (name == "Identity") {
    return std::make_shared<Identity<T>>();
  } else if (name == "Sigmoid") {
    return std::make_shared<Sigmoid<T>>();
  } else if (name == "ForceFactor") {
    return std::make_shared<ForceFactor<T>>(converter);
  } else if (name == "Rectifier") {
    return std::make_shared<Rectifier<T>>();
  } else if (name == "Softplus") {
    return std::make_shared<Softplus<T>>();
  } else if (name == "Baron") {
    return std::make_shared<Baron<T>>();
  } else if (name == "Krause") {
    return std::make_shared<Krause<T>>();
  } else if (name == "Foerster") {
    return std::make_shared<Foerster<T,DESCRIPTOR>>(converter);
  } else if (std::regex_match(name, match, std::regex ("(Foerster)([0-9])"))) {
    const unsigned n = std::stoi(match[2]);
    return std::make_shared<FoersterN<T,DESCRIPTOR>>(converter, n);
  } else if (std::regex_match(name, match, std::regex ("(Stasius)([0-9])"))) {
    const unsigned n = std::stoi(match[2]);
    return std::make_shared<StasiusN<T,DESCRIPTOR>>(converter, n);
  } else {
    throw std::invalid_argument(
      "Error: unknown projection name selected.");
  }
}

}

template <typename PROJECTION, typename T>
std::vector<T> applyProjection(const std::vector<T>& vector) {
  std::vector<T> projected;
  for (auto& e : vector) {
    projected.push_back(PROJECTION().project(e));
  }
  return projected;
}

template <typename PROJECTION, typename T, typename ARG>
std::vector<T> applyProjection(const std::vector<T>& vector, ARG arg) {
  std::vector<T> projected;
  for (auto& e : vector) {
    projected.push_back(PROJECTION(arg).project(e));
  }
  return projected;
}

} // namespace opti

} // namespace olb
#endif // PROJECTION_H
