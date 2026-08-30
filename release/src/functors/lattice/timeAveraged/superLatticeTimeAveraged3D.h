/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Mathias J. Krause, Benedict Hasenauer
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

#ifndef SUPER_LATTICE_TIME_AVERAGED_F3_D_H
#define SUPER_LATTICE_TIME_AVERAGED_F3_D_H

#include <vector>

#include <functors/lattice/euklidNorm3D.h>

namespace olb {


// Averages a functor value about a timespan and gives back the averaged value(TA)
// and in 3*Dim the root mean square value(RMS) of the functorvalue in Dim in the operator
// TA  = SUM(functorvalue(iT)/SUM(iT))
// RMS = SQRT( SUM( (functorvalue(iT) - TA)^2 / SUM(iT) ) )
template <typename T>
class SuperLatticeTimeAveragedF3D : public SuperF3D<T,T> {
protected:
  int _ensembles;
  SuperF3D<T,T>& _sFunctor;
  SuperData<3,T, T> _sData;
  SuperData<3,T, T> _sDataP2;

  SuperData<3,T,T> _sDataHelp;     // for Kahan summation
  SuperData<3,T,T> _sDataP2help;   // for Kahan summation

public:
  SuperLatticeTimeAveragedF3D(SuperF3D<T,T>& sFunctor);

  bool operator() (T output[], const int input[]);

  void addEnsemble();

  int getEnsembles();
  int getBlockFSize() const;

};

template <typename T>
class SuperLatticeTimeAveragedMagnitudesF3D : public SuperF3D<T,T> {
protected:
  int _ensembles;
  SuperF3D<T,T>& _sFunctor;
  SuperData<3,T, T> _sData;
  SuperData<3,T, T> _sDataP2;

  SuperData<3,T,T> _sDataHelp;     // for Kahan summation
  SuperData<3,T,T> _sDataP2help;   // for Kahan summation

public:
  SuperLatticeTimeAveragedMagnitudesF3D(SuperF3D<T,T>& sFunctor);

  bool operator() (T output[], const int input[]);

  void addEnsemble();

  int getEnsembles();
  int getBlockFSize() const;

};

// The functor calculates the crosscorrelation(CC) of two functorvalues averaged above the Time
// CC = SUM((functorM[iT] - functorMAverage)*(functorN[iT] - functorNAverage))
// the dimesion of the functor is the product of the given functor dimensions
// the output if the functor M and N have two dimesnions is {m0*n0,m0*n1,m1*n0,m1*n0}
template <typename T>
class SuperLatticeTimeAveragedCrossCorrelationF3D final:  public SuperF3D<T,T> {
private:
  int _ensembles;
  SuperF3D<T,T>& _sFunctorM;
  SuperF3D<T,T>& _sFunctorN;
  SuperData<3,T, T> _sDataM;
  SuperData<3,T, T> _sDataN;
  SuperData<3,T, T> _sDataMN;

  // for Kahan summation
  SuperData<3,T,T> _sDataMhelp;
  SuperData<3,T,T> _sDataNhelp;
  SuperData<3,T,T> _sDataMNhelp;

public:
  SuperLatticeTimeAveragedCrossCorrelationF3D(SuperF3D<T,T>& sFunctorM, SuperF3D<T,T>& sFunctorN);

  bool operator() (T output[], const int input[]);

  void addEnsemble();

};

}

#endif
