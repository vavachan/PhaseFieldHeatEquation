/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

#ifndef DUAL_FUNCTORS_H
#define DUAL_FUNCTORS_H


#include "functors/analytical/analyticalBaseF.h"
#include "functors/lattice/reductionF3D.h"
#include "functors/lattice/superBaseF3D.h"

namespace olb {

namespace opti {

/// functor to get the pointwise dual dissipation density on local lattices, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class BlockLatticeDphysDissipationDf : public BlockLatticeF<T,DESCRIPTOR> {
protected:
  const int _overlap;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeDphysDissipationDf(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                   int overlap,
                                   const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator()(T output[], const int input[]);
};

/// functor to get pointwise dual dissipation density on local lattices, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeDphysDissipationDf : public SuperLatticePhysF<T,DESCRIPTOR> {
public:
  SuperLatticeDphysDissipationDf(SuperLattice<T,DESCRIPTOR>& sLattice,
                                   const UnitConverter<T,DESCRIPTOR>& converter);
};


/// functor to get pointwise dual velocity density on local lattices, if globIC is not on
/// the local processor, the returned vector is empty
// this is made a BGK-type momenta computation.
template <typename T, typename DESCRIPTOR>
class BlockLatticeDphysVelocityDf3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const int _overlap;
  const UnitConverter<T,DESCRIPTOR>& _converter;
private:
  const int _nDim;
  const int _extractDim;
public:
  BlockLatticeDphysVelocityDf3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                int overlap,
                                const UnitConverter<T,DESCRIPTOR>& converter,
                                int nDim, int extractDim);
  bool operator()(T output[], const int input[]);
};

/// functor to get pointwise dual velocity density on local lattices, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeDphysVelocityDf3D : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeDphysVelocityDf3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticeDphysVelocityDf3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter,
                                SuperExtractComponentF3D<T,T>& sExtract);
};


/// functor to compute 0.5*L2Norm(f-f_wanted)^2 on a lattice
template <typename T, typename DESCRIPTOR>
class DifferenceObjective3D : public SuperIdentity3D<T,T> {
public:
  DifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&     sLattice,
    FunctorPtr<SuperF3D<T,T>>&&        f,
    FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  DifferenceObjective3D(
    FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>>&& f,
    FunctorPtr<AnalyticalF3D<T,T>>&&            wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&&          indicatorF);

  DifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&   sLattice,
    FunctorPtr<SuperF3D<T,T>>&&      f,
    FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
    SuperGeometry<T,3>& geometry,
    int                 material);
};


/// functor to compute 0.5*(f-f_wanted)^2 on a lattice
template <typename T, typename DESCRIPTOR>
class BlockDdifferenceObjectiveDf3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>&          _f;
  BlockF3D<T>&          _dFdF;
  BlockF3D<T>&          _wantedF;
  BlockIndicatorF3D<T>& _indicatorF;

  const T _weight;
public:
  BlockDdifferenceObjectiveDf3D(
    BlockLattice<T,DESCRIPTOR>& blockLattice,
    BlockF3D<T>&          f,
    BlockF3D<T>&          dFdF,
    BlockF3D<T>&          wantedF,
    BlockIndicatorF3D<T>& indicatorF,
    T weight);

  bool operator()(T output[], const int input[]);
};

/// functor to compute 0.5*(f-f_wanted)^2 on a lattice
template <typename T, typename DESCRIPTOR>
class DdifferenceObjectiveDf3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>> _f;
  FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>> _dFdF;
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR>  _wantedF;
  FunctorPtr<SuperIndicatorF3D<T>>           _indicatorF;
public:
  DdifferenceObjectiveDf3D(
    FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>>&& f,
    FunctorPtr<SuperLatticePhysF3D<T,DESCRIPTOR>>&& dFdF,
    FunctorPtr<AnalyticalF3D<T,T>>&&             wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&&           indicatorF);

  bool operator()(T output[], const int input[]);
};


/// functor to compute 0.5*L2Norm(f-f_wanted)^2/L2Norm(f_wanted)^2 on a lattice
template<typename T, typename DESCRIPTOR>
class RelativeDifferenceObjective3D : public SuperIdentity3D<T,T> {
public:
  RelativeDifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&     sLattice,
    FunctorPtr<SuperF3D<T,T>>&&        f,
    FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  RelativeDifferenceObjective3D(
    FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>>&& f,
    FunctorPtr<AnalyticalF3D<T,T>>&&            wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&&          indicatorF);

  RelativeDifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&   sLattice,
    FunctorPtr<SuperF3D<T,T>>&&      f,
    FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
    SuperGeometry<T,3>& geometry,
    std::vector<int>    materials);
  RelativeDifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&   sLattice,
    FunctorPtr<SuperF3D<T,T>>&&      f,
    FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
    SuperGeometry<T,3>& geometry,
    int                 material);

  // with two lattice functors
  RelativeDifferenceObjective3D(
    SuperLattice<T, DESCRIPTOR>&     sLattice,
    FunctorPtr<SuperF3D<T,T>>&&        f,
    FunctorPtr<SuperF3D<T,T>>&&   wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
};


/// functor to compute 0.5(f-f_wanted)^2/f_wanted^2 on a lattice
template <typename T, typename DESCRIPTOR>
class BlockDrelativeDifferenceObjectiveDf3D : public BlockF3D<T> {
private:
  BlockF3D<T>&          _f;
  BlockF3D<T>&          _dFdF;
  BlockF3D<T>&          _wantedF;
  BlockIndicatorF3D<T>& _indicatorF;

  const T   _globalValue;
  const T   _weight;
public:
  BlockDrelativeDifferenceObjectiveDf3D(
    BlockStructureD<3>& blockStructure,
    BlockF3D<T>&          f,
    BlockF3D<T>&          dFdF,
    BlockF3D<T>&          wantedF,
    BlockIndicatorF3D<T>& indicatorF,
    T globalValue,
    T weight);

  bool operator()(T output[], const int input[]);
};

/// functor to compute 0.5(f-f_wanted)^2/f_wanted^2 on a lattice
template <typename T, typename DESCRIPTOR>
class DrelativeDifferenceObjectiveDf3D : public SuperF3D<T,T> {
private:
  FunctorPtr<SuperF3D<T,T>>                 _f;
  FunctorPtr<SuperF3D<T,T>>                 _dFdF;
  FunctorPtr<SuperF3D<T,T>>                 _wantedF;
  FunctorPtr<SuperIndicatorF3D<T>>          _indicatorF;
public:
  DrelativeDifferenceObjectiveDf3D(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   FunctorPtr<SuperF3D<T,T>>&&        f,
                                   FunctorPtr<SuperF3D<T,T>>&&        dFdF,
                                   FunctorPtr<SuperF3D<T,T>>&&   wantedF,
                                   FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  DrelativeDifferenceObjectiveDf3D(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   FunctorPtr<SuperF3D<T,T>>&&        f,
                                   FunctorPtr<SuperF3D<T,T>>&&        dFdF,
                                   FunctorPtr<AnalyticalF3D<T,T>>&&   wantedF,
                                   FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  DrelativeDifferenceObjectiveDf3D(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   FunctorPtr<SuperF3D<T,T>>&&      f,
                                   FunctorPtr<SuperF3D<T,T>>&&      dFdF,
                                   FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
                                   SuperGeometry<T,3>& geometry,
                                   std::vector<int>    materials);
  DrelativeDifferenceObjectiveDf3D(SuperLattice<T, DESCRIPTOR>& sLattice,
                                   FunctorPtr<SuperF3D<T,T>>&&      f,
                                   FunctorPtr<SuperF3D<T,T>>&&      dFdF,
                                   FunctorPtr<AnalyticalF3D<T,T>>&& wantedF,
                                   SuperGeometry<T,3>& geometry,
                                   int                 material);

  bool operator()(T output[], const int input[]);
};

/// functor to compute 0.5*(f[extractDim]-f_wanted[0])^2/f_wanted^2 on a lattice
template <typename T, typename DESCRIPTOR>
class BlockDrelativeDifferenceObjectiveComponentDf3D : public BlockF3D<T> {
private:
  BlockF3D<T>&          _f;
  BlockF3D<T>&          _dFdF;
  BlockF3D<T>&          _wantedF;
  BlockIndicatorF3D<T>& _indicatorF;

  const int _extractDim;
  const T   _globalValue;
  const T   _weight;
public:
  BlockDrelativeDifferenceObjectiveComponentDf3D(
    BlockStructureD<3>& blockStructure,
    BlockExtractComponentF3D<T>& f,
    BlockF3D<T>&                 dFdF,
    BlockF3D<T>&                 wantedF,
    BlockIndicatorF3D<T>&        indicatorF,
    T globalValue,
    T weight);

  bool operator()(T output[], const int input[]);
};

/// functor to compute 0.5*(f[extractDim]-f_wanted[0])^2/f_wanted^2 on a lattice
template <typename T, typename DESCRIPTOR>
class DrelativeDifferenceObjectiveComponentDf3D : public SuperF3D<T,T> {
private:
  SuperExtractComponentF3D<T,T>             _f;
  FunctorPtr<SuperF3D<T,T>>                 _dFdF;
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> _wantedF;
  FunctorPtr<SuperIndicatorF3D<T>>          _indicatorF;
public:
  DrelativeDifferenceObjectiveComponentDf3D(
    SuperLattice<T, DESCRIPTOR>& sLattice,
    FunctorPtr<SuperF3D<T,T>>&&         f, int extractDim,
    FunctorPtr<SuperF3D<T,T>>&&         dFdF,
    FunctorPtr<AnalyticalF3D<T,T>>&&    wantedF,
    FunctorPtr<SuperIndicatorF3D<T>>&&  indicatorF);

  bool operator()(T output[], const int input[]);
};

} // namespace opti

} // namespace olb

#endif
