/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef BOUNDARIES_VORTEX_METHOD_H
#define BOUNDARIES_VORTEX_METHOD_H

#include "utilities/functorPtr.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

namespace olb {

struct VortexMethodPreProcessor;
struct VortexMethodPostProcessor;

struct U_PROFILE : public descriptors::FIELD_BASE<0,1> { };
struct VELOCITY_OLD : public descriptors::FIELD_BASE<0,1> { };

struct CONVERSION_FACTOR_LENGTH : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_FACTOR_VELOCITY : public descriptors::FIELD_BASE<1> { };

struct SEEDS_COUNT : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };
struct SEEDS : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };
struct SEEDS_VORTICITY : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };

struct AXIS_DIRECTION : public descriptors::FIELD_BASE<0,1> { };
struct SIGMA : public descriptors::FIELD_BASE<1> { };

namespace stage {

struct VortexMethod { };

}

template <typename T, typename DESCRIPTOR>
class VortexMethodTurbulentVelocityBoundary final {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _inletLatticeI;
  FunctorPtr<IndicatorF3D<T>> _inletPhysI;
  std::shared_ptr<AnalyticalF3D<T,T>> _velocityProfileF;
  std::shared_ptr<AnalyticalF3D<T,T>> _intensityProfileF;

  UnitConverter<T,DESCRIPTOR>& _converter;
  SuperLattice<T, DESCRIPTOR>& _sLattice;

  const int _nSeeds;
  T _nTime;
  T _inletArea;
  T _sigma;

  Vector<T,3> _axisDirection;

  std::vector<T> _AiC;
  std::vector<int> _NiC;

  SuperFieldArrayD<T,DESCRIPTOR,descriptors::LOCATION> _seeds;
  SuperFieldArrayD<T,DESCRIPTOR,descriptors::VORTICITY> _seedsVorticity;

  void generateSeeds();
  void updateSeedsVorticity(std::size_t iT);

public:
  VortexMethodTurbulentVelocityBoundary(
    FunctorPtr<SuperIndicatorF3D<T>>&& inletLatticeI,
    FunctorPtr<IndicatorF3D<T>>&& inletPhysI,
    UnitConverter<T,DESCRIPTOR>& converter,
    SuperLattice<T, DESCRIPTOR>& sLattice,
    int nSeeds,
    T nTime,
    T sigma,
    Vector<T,3> axisDirection);

  void setVelocityProfile(std::shared_ptr<AnalyticalF3D<T,T>> velocityProfileF);
  void setIntensityProfile(std::shared_ptr<AnalyticalF3D<T,T>> intensityProfileF);

  void apply(std::size_t iT);

};

template <typename T, typename DESCRIPTOR>
VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::VortexMethodTurbulentVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& inletLatticeI,
  FunctorPtr<IndicatorF3D<T>>&& inletPhysI,
  UnitConverter<T,DESCRIPTOR>& converter,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  int nSeeds,
  T nTime,
  T sigma,
  Vector<T,3> axisDirection)
  : _inletLatticeI(std::move(inletLatticeI)),
    _inletPhysI(std::move(inletPhysI)),
    _converter(converter),
    _sLattice(sLattice),
    _nSeeds(nSeeds),
    _nTime(nTime),
    _sigma(sigma),
    _axisDirection(axisDirection),
    _AiC(_sLattice.getLoadBalancer().size()),
    _NiC(_sLattice.getLoadBalancer().size()),
    _seeds(sLattice.getCuboidDecomposition(),
           sLattice.getLoadBalancer()),
    _seedsVorticity(sLattice.getCuboidDecomposition(),
                    sLattice.getLoadBalancer())
{
  OstreamManager clout(std::cout, "VortexMethod");

  sLattice.template addPostProcessor<stage::VortexMethod>(*_inletLatticeI,
                                                          meta::id<VortexMethodPreProcessor>{});
  SuperIndicatorLayer3D<T> extendedLatticeInletI(*_inletLatticeI);
  sLattice.template addPostProcessor<stage::VortexMethod>(extendedLatticeInletI,
                                                          meta::id<VortexMethodPostProcessor>{});

  if (_sigma < _converter.getPhysDeltaX()) {
    _sigma = _converter.getPhysDeltaX();
  }

  _inletArea = 0.;

  auto& cGeometry = sLattice.getCuboidDecomposition();
  Cuboid3D<T> indicatorCuboid(*inletPhysI, _converter.getPhysDeltaX());

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    _AiC[iC] = 0.;
    auto& cuboid = cGeometry.get(load.glob(iC));
    if (cuboid.intersects(indicatorCuboid)) {
      auto& block = _sLattice.getBlock(iC);
      for (int iX=0; iX < block.getNx(); ++iX) {
        for (int iY=0; iY < block.getNy(); ++iY) {
          for (int iZ=0; iZ < block.getNz(); ++iZ) {
            int latticeR[3] {iX, iY, iZ};
            auto physR = cuboid.getPhysR(latticeR);
            block.get(iX,iY,iZ).template setField<descriptors::LOCATION>(physR);
            bool inInlet{};
            _inletPhysI(&inInlet, physR.data());
            if (inInlet) {
              _AiC[iC] += util::pow(_converter.getPhysDeltaX(), 2.);
            }
          }
        }
      }
    }
    _inletArea += _AiC[iC];
  }

  #ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(_inletArea, MPI_SUM);
  #endif

  clout << "inletArea=" << _inletArea << std::endl;

  // initial seeding of random vortex points of number N
  for (int iC=0; iC < load.size(); ++iC) {
    _NiC[iC] = int(_nSeeds*_AiC[iC] / _inletArea);
    _seeds.getBlock(iC).resize(_NiC[iC]);
    _seedsVorticity.getBlock(iC).resize(_NiC[iC]);
  }

  for (int iC=0; iC < load.size(); ++iC) {
    if (_NiC[iC] > 0) {
      auto& block = _sLattice.getBlock(iC);
      block.template setParameter<SEEDS_COUNT>(_NiC[iC]);
      block.template setParameter<SEEDS>(_seeds.getBlock(iC));
      block.template setParameter<SEEDS_VORTICITY>(_seedsVorticity.getBlock(iC));
      block.template getField<VELOCITY_OLD>();
      block.template getField<U_PROFILE>();
    }
  }

  _sLattice.template setParameter<CONVERSION_FACTOR_LENGTH>(_converter.getConversionFactorLength());
  _sLattice.template setParameter<CONVERSION_FACTOR_VELOCITY>(_converter.getConversionFactorVelocity());
  _sLattice.template setParameter<SIGMA>(_sigma);
  _sLattice.template setParameter<AXIS_DIRECTION>(axisDirection);

  _sLattice.setProcessingContext(ProcessingContext::Simulation);
}


// random points seed each time step
template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::generateSeeds()
{
  std::random_device seed;
  std::mt19937 generator(seed());
  std::uniform_real_distribution<T> distribution(0, 1);
  auto randomness = [&distribution,&generator]() -> T {
    return distribution(generator);
  };

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    auto& seeds = _seeds.getBlock(iC);
    if (_NiC[iC] > 0) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int j = 0; j < _NiC[iC]; j++) {
        auto seed = _inletPhysI->getSample(randomness);
        seeds.set(j, seed);
      }
      seeds.setProcessingContext(ProcessingContext::Simulation);
    }
  }
}

// vorticity calculation
template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::updateSeedsVorticity(std::size_t iT)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(0, 1);
  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    if (_NiC[iC] > 0) {
      auto& seeds = _seeds.getBlock(iC);
      auto& seedsVorticity = _seedsVorticity.getBlock(iC);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for(int i = 0; i<_NiC[iC]; i++) {
        Vector<T,3> inlet = seeds.get(i);
        Vector<T,3> velocity;
        _velocityProfileF->operator()(velocity.data(), inlet.data());
        Vector<T,1> intensity;
        _intensityProfileF->operator()(intensity.data(), inlet.data());
        T normVel = olb::norm(velocity);
        T kinE = 3./2. * util::pow(normVel * intensity[0], 2. );
        T vorticity = 4. * util::sqrt( M_PI * _AiC[iC] * kinE / 3. / _NiC[iC] / (2.*util::log(3.) - 3.*util::log(2.)) );

        // random vorticity sign change
        T sign = 1.;
        if(iT % _converter.getLatticeTime(_nTime) == 0) {
          if(distrib(gen) == 0) {
            sign *= -1;
          }
        }
        vorticity *= sign;

        seedsVorticity.set(i,vorticity);
      }
      seedsVorticity.setProcessingContext(ProcessingContext::Simulation);
    }
  }
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setVelocityProfile(std::shared_ptr<AnalyticalF3D<T,T>> velocityProfileF)
{
  _velocityProfileF = velocityProfileF;
  _sLattice.template defineField<U_PROFILE>(*_inletLatticeI, *_velocityProfileF);
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setIntensityProfile(std::shared_ptr<AnalyticalF3D<T,T>> intensityProfileF)
{
  _intensityProfileF = intensityProfileF;
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::apply(std::size_t iT)
{
  generateSeeds();
  updateSeedsVorticity(iT);

  {
    auto& load = _sLattice.getLoadBalancer();
    for (int iC=0; iC < load.size(); ++iC) {
      if (_NiC[iC] > 0) {
        auto& block = _sLattice.getBlock(iC);
        block.template setParameter<SEEDS_COUNT>(_NiC[iC]);
        block.template setParameter<SEEDS>(_seeds.getBlock(iC));
        block.template setParameter<SEEDS_VORTICITY>(_seedsVorticity.getBlock(iC));
      }
    }
    _sLattice.template setParameter<SIGMA>(_sigma);
    _sLattice.template setParameter<AXIS_DIRECTION>(_axisDirection);
  }
  _sLattice.executePostProcessors(stage::VortexMethod{});
}


struct VortexMethodPreProcessor {
  using parameters = meta::list<
    SEEDS_COUNT,
    SEEDS,
    SEEDS_VORTICITY,
    CONVERSION_FACTOR_LENGTH,
    CONVERSION_FACTOR_VELOCITY,
    AXIS_DIRECTION,
    SIGMA
  >;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  Vector<V,3> computeVelocityGradient(CELL& cell, PARAMETERS& parameters) any_platform
  {
    const V dx = parameters.template get<CONVERSION_FACTOR_LENGTH>();
    Vector<V,3> actVel = cell.template getField<VELOCITY_OLD>();
    V aXPlus = util::max(actVel[0], 0); V aXMinus = util::min(actVel[0], 0);
    V aYPlus = util::max(actVel[1], 0); V aYMinus = util::min(actVel[1], 0);
    V aZPlus = util::max(actVel[2], 0); V aZMinus = util::min(actVel[2], 0);
    Vector<V,3> uXPlus  = cell.neighbor({ 1,  0,  0}).template getField<VELOCITY_OLD>();
    Vector<V,3> uXMinus = cell.neighbor({-1,  0,  0}).template getField<VELOCITY_OLD>();
    Vector<V,3> uYPlus  = cell.neighbor({ 0,  1,  0}).template getField<VELOCITY_OLD>();
    Vector<V,3> uYMinus = cell.neighbor({ 0, -1,  0}).template getField<VELOCITY_OLD>();
    Vector<V,3> uZPlus  = cell.neighbor({ 0,  0,  1}).template getField<VELOCITY_OLD>();
    Vector<V,3> uZMinus = cell.neighbor({ 0,  0, -1}).template getField<VELOCITY_OLD>();
    Vector<V,3> grad;
    const V actVelNorm = olb::norm(actVel);
    grad[0] = -aXPlus*(actVelNorm - olb::norm(uXMinus))/dx - aXMinus*(olb::norm(uXPlus) - actVelNorm)/dx;
    grad[1] = -aYPlus*(actVelNorm - olb::norm(uYMinus))/dx - aYMinus*(olb::norm(uYPlus) - actVelNorm)/dx;
    grad[2] = -aZPlus*(actVelNorm - olb::norm(uZMinus))/dx - aZMinus*(olb::norm(uZPlus) - actVelNorm)/dx;
    return grad;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    const auto seeds = parameters.template get<SEEDS>();
    const auto seedsVorticity = parameters.template get<SEEDS_VORTICITY>();
    const auto axisDirection = parameters.template get<AXIS_DIRECTION>();
    const auto sigma = parameters.template get<SIGMA>();

    Vector<V,DESCRIPTOR::d> x = cell.template getField<descriptors::LOCATION>();
    Vector<V,DESCRIPTOR::d> output = cell.template getField<U_PROFILE>();
    const V conversionFactorVelocity = parameters.template get<CONVERSION_FACTOR_VELOCITY>();
    output /= conversionFactorVelocity;

    // calculation of velocities from the placed vortexes
    Vector<V,3> uVortex {0, 0, 0};
    for (std::size_t j=0; j < parameters.template get<SEEDS_COUNT>(); ++j) {
      Vector<V,DESCRIPTOR::d> diffX = {0,0,0};
      for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
        diffX[iD] = seeds[iD][j] - x[iD];
      }
      Vector<V,3> crossP = crossProduct3D(diffX, axisDirection);
      V norm2 = diffX[0]*diffX[0] + diffX[1]*diffX[1] + diffX[2]*diffX[2];
      V expTerm = util::exp(V{-0.5}*norm2/(sigma*sigma));
      uVortex += (V{0.5}/M_PI * seedsVorticity[j] * crossP / norm2 * (V{1} - expTerm)*expTerm) / conversionFactorVelocity;
    }

    // calculation of the fluctuation streamwise with Langevin equation
    auto grad = computeVelocityGradient(cell, parameters);
    V normGrad = olb::norm(grad);
    if (normGrad == V{0}) {
      normGrad = 1;
    }
    V streamFluct = -(uVortex*grad);
    streamFluct /= normGrad;

    output += uVortex + streamFluct * axisDirection;

    cell.defineU(output.data());
  }

};

struct VortexMethodPostProcessor {
  using parameters = meta::list<CONVERSION_FACTOR_VELOCITY>;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return std::numeric_limits<int>::max();
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    Vector<V,DESCRIPTOR::d> u{};
    cell.computeU(u.data());
    u *= parameters.template get<CONVERSION_FACTOR_VELOCITY>();
    cell.template setField<VELOCITY_OLD>(u);
  }

};

}

#endif
