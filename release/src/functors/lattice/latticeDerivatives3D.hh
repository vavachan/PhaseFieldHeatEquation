/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen, Mathias J. Krause
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

#ifndef LATTICE_DERIVATIVES_3D_HH
#define LATTICE_DERIVATIVES_3D_HH

#include "utilities/finiteDifference.h"
#include "latticeDerivatives3D.h"

namespace olb {


////////////////////////BlockFiniteDifference3D//////////////////////////////////
template <typename T>
BlockFiniteDifference3D<T>::BlockFiniteDifference3D
(BlockGeometry<T,3>& blockGeometry, BlockF3D<T>& blockFunctor, std::list<int>& matNumber)
  : BlockF3D<T>(blockFunctor.getBlockStructure(), 3*blockFunctor.getTargetDim()), _blockGeometry(blockGeometry), _blockFunctor(blockFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  _targetDim = _blockFunctor.getTargetDim();
  _n[0] = this-> _blockGeometry.getNx()-1;
  _n[1] = this-> _blockGeometry.getNy()-1;
  _n[2] = this-> _blockGeometry.getNz()-1;

}

template <typename T>
bool BlockFiniteDifference3D<T>::operator() (T output[], const int input[])
{
//  // derivation tensor
  std::vector<std::vector<T>> fdGrad;

  fdGrad.resize(_targetDim);
  for (int i = 0; i < _targetDim; i++) {
    fdGrad[i].resize(3);
  }

  for (int i = 0; i < 3; i++) {
    int fInput_p[3];
    fInput_p[0] = input[0];
    fInput_p[1] = input[1];
    fInput_p[2] = input[2];
    fInput_p[i]+=1;

    int fInput_2p[3];
    fInput_2p[0] = input[0];
    fInput_2p[1] = input[1];
    fInput_2p[2] = input[2];
    fInput_2p[i]+=2;

    int fInput_3p[3];
    fInput_3p[0] = input[0];
    fInput_3p[1] = input[1];
    fInput_3p[2] = input[2];
    fInput_3p[i]+=3;

    int fInput_4p[3];
    fInput_4p[0] = input[0];
    fInput_4p[1] = input[1];
    fInput_4p[2] = input[2];
    fInput_4p[i]+=4;

    int fInput_n[3];
    fInput_n[0] = input[0];
    fInput_n[1] = input[1];
    fInput_n[2] = input[2];
    fInput_n[i]-=1;

    int fInput_2n[3];
    fInput_2n[0] = input[0];
    fInput_2n[1] = input[1];
    fInput_2n[2] = input[2];
    fInput_2n[i]-=2;

    int fInput_3n[3];
    fInput_3n[0] = input[0];
    fInput_3n[1] = input[1];
    fInput_3n[2] = input[2];
    fInput_3n[i]-=3;

    int fInput_4n[3];
    fInput_4n[0] = input[0];
    fInput_4n[1] = input[1];
    fInput_4n[2] = input[2];
    fInput_4n[i]-=4;

    T fOutput[_targetDim];
    _blockFunctor(fOutput,input);

    if (input[i] < 3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2p[0], fInput_2p[1], fInput_2p[2]})) == _matNumber.end()) {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
        }
      }
      else {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
        }
      }
    }
    else if (input[i] > _n[i]-3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2n[0], fInput_2n[1], fInput_2n[2]})) == _matNumber.end()) {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
        }
      }
      else {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
        }
      }
    }
    else {
      if ( std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_n[0], fInput_n[1], fInput_n[2]})) == _matNumber.end()  &&
          std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_p[0], fInput_p[1], fInput_p[2]})) == _matNumber.end() ) {
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=0.;
        }
        // boundary treatment with Second-order asymmetric gradient
      }
      else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_n[0], fInput_n[1], fInput_n[2]})) == _matNumber.end()) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2p[0], fInput_2p[1], fInput_2p[2]})) == _matNumber.end()) {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
          }
        }
        else {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          T fOutput_2p[_targetDim];
          _blockFunctor(fOutput_2p,fInput_2p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
          }
        }
      }
      else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_p[0], fInput_p[1], fInput_p[2]})) == _matNumber.end() ) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get({fInput_2n[0], fInput_2n[1], fInput_2n[2]})) == _matNumber.end()) {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
          }
        }
        else {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          T fOutput_2n[_targetDim];
          _blockFunctor(fOutput_2n,fInput_2n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
          }
        }
      }
      else {
        //inner domain 8th order central difference
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);

        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);

        T fOutput_3n[_targetDim];
        _blockFunctor(fOutput_3n,fInput_3n);

        T fOutput_4n[_targetDim];
        _blockFunctor(fOutput_4n,fInput_4n);

        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);

        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);

        T fOutput_3p[_targetDim];
        _blockFunctor(fOutput_3p,fInput_3p);

        T fOutput_4p[_targetDim];
        _blockFunctor(fOutput_4p,fInput_4p);
        for (int j=0; j < _targetDim; j++) {
          //fdGrad[j][i]=fd::centralGradient(fOutput_p[j], fOutput_n[j]);
          fdGrad[j][i]=((T)672*(fOutput_p[j]-fOutput_n[j])+(T)168*(fOutput_2n[j]-fOutput_2p[j])
                        +(T)32*(fOutput_3p[j]-fOutput_3n[j])+(T)3*(fOutput_4n[j]-fOutput_4p[j])) / 840.;
        }
      }
    }
    for (int i=0; i < 3; i++) {
      for (int j=0; j < _targetDim; j++) {
        output[j*3+i] = fdGrad[j][i];
      }
    }
  }
  return true;
}

////////////////////////SuperFiniteDifference3D//////////////////////////////////
template <typename T>
SuperFiniteDifference3D<T>::SuperFiniteDifference3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sGeometry(sGeometry),_sFunctor(sFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockFiniteDifference3D<T> ( _sGeometry.getBlockGeometry(iC), _sFunctor.getBlockF(iC), _matNumber ));
  }
}

////////////////////////BlockPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockPhysFiniteDifference3D<T,DESCRIPTOR>::BlockPhysFiniteDifference3D
(BlockF3D<T>& blockFinDiff, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockF3D<T>(blockFinDiff.getBlockStructure(), 3*blockFinDiff.getTargetDim()), _blockFinDiff(blockFinDiff), _converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  _targetDim = _blockFinDiff.getTargetDim();

}

template <typename T, typename DESCRIPTOR>
bool BlockPhysFiniteDifference3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockFinDiff(output,input);
  for (int i = 0; i < _targetDim; i++) {
    output[i] /= _converter.getConversionFactorLength();
  }
  return true;
}

////////////////////////SuperPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperPhysFiniteDifference3D<T,DESCRIPTOR>::SuperPhysFiniteDifference3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sFinDiff(sGeometry,sFunctor,matNumber),_converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockPhysFiniteDifference3D<T,DESCRIPTOR> (_sFinDiff.getBlockF(iC), _converter ));
  }
}


////////////////////////BlockLaplacian3D//////////////////////////////////
template <typename T>
BlockLaplacian3D<T>::BlockLaplacian3D
(BlockGeometry<T,3>& blockGeometry, BlockF3D<T>& blockFunctor,
  bool forthOrder)
  : BlockF3D<T>(blockFunctor.getBlockStructure(), blockFunctor.getTargetDim()),
    _blockGeometry(blockGeometry), _blockFunctor(blockFunctor),
    _forthOrder(forthOrder)
{
  this->getName() = "Laplacian(" + blockFunctor.getName() + ")";
  _n[0] = this-> _blockGeometry.getNx()-1;
  _n[1] = this-> _blockGeometry.getNy()-1;
  _n[2] = this-> _blockGeometry.getNz()-1;
}

template <typename T>
bool BlockLaplacian3D<T>::operator() (T output[], const int input[])
{
  for (int i=0; i<this->getTargetDim(); ++i) {
    output[i] = 0;
  }
  int inputMod[3];
  for (int j=0; j<3; ++j) {
    inputMod[j] = input[j];
  }
  T u_0[this->getTargetDim()];
  T u_m[this->getTargetDim()];
  T u_p[this->getTargetDim()];

  // only used for forth order dq
  T u_2m[this->getTargetDim()];
  T u_2p[this->getTargetDim()];

  // compute discrete second derivatives for each direction and add them
  for (int j=0; j<3; ++j) {
    if (input[j] < 1) {
      // use forward difference quotient at the boundary
      _blockFunctor.operator()(u_m, inputMod);
      inputMod[j] += 1;
      _blockFunctor.operator()(u_0, inputMod);
      inputMod[j] += 1;
      _blockFunctor.operator()(u_p, inputMod);
    } else if (input[j] > _n[j]-1) {
      // use backward difference quotient at the boundary
      _blockFunctor.operator()(u_p, inputMod);
      inputMod[j] -= 1;
      _blockFunctor.operator()(u_0, inputMod);
      inputMod[j] -= 1;
      _blockFunctor.operator()(u_m, inputMod);
    } else {
      // use central difference quotient
      _blockFunctor.operator()(u_0, inputMod);
      inputMod[j] -= 1;
      _blockFunctor.operator()(u_m, inputMod);
      inputMod[j] += 2;
      _blockFunctor.operator()(u_p, inputMod);
    }
    inputMod[j] = input[j];  // reset

    if ((_forthOrder) && (input[j] > 1) && (input[j] < _n[j]-1)) {
      // additional evaluations
      inputMod[j] -= 2;
      _blockFunctor.operator()(u_2m, inputMod);
      inputMod[j] += 4;
      _blockFunctor.operator()(u_2p, inputMod);
      inputMod[j] = input[j];  // reset

      // forth order dq
      for (int i=0; i<this->getTargetDim(); ++i) {
        output[i] += fd::centralSecondDeriv(u_2m[i], u_m[i], u_0[i], u_p[i], u_2p[i]);
      }
    } else {
      // second order dq
      for (int i=0; i<this->getTargetDim(); ++i) {
        output[i] += fd::centralSecondDeriv(u_m[i], u_0[i], u_p[i]);
      }
    }
  }
  return true;
/*
  // todo: enable treating all spatial directions at once?
  _blockFunctor.operator()(u_0, inputMod);

  inputMod[0] -= 1;
  _blockFunctor.operator()(u_xm1, inputMod);
  inputMod[0] += 2;
  _blockFunctor.operator()(u_xp1, inputMod);
  inputMod[0] -= 1;

  inputMod[1] -= 1;
  _blockFunctor.operator()(u_ym1, inputMod);
  inputMod[1] += 2;
  _blockFunctor.operator()(u_yp1, inputMod);
  inputMod[1] -= 1;

  inputMod[2] -= 1;
  _blockFunctor.operator()(u_zm1, inputMod);
  inputMod[2] += 2;
  _blockFunctor.operator()(u_zp1, inputMod);
  inputMod[2] -= 1;

  for (int i=0; i<_targetDim; ++i) {
    output[i] = fd::laplacian3D(u_xm1[i], u_ym1[i], u_zm1[i], u_0[i], u_xp1[i], u_yp1[i], u_zp1[i]);
  }
*/
}


////////////////////////SuperLaplacian3D//////////////////////////////////
template <typename T>
SuperLaplacian3D<T>::SuperLaplacian3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor,
  bool forthOrder)
  : SuperF3D<T>(sFunctor.getSuperStructure(),sFunctor.getTargetDim()),
    _sGeometry(sGeometry),_sFunctor(sFunctor)
{
  this->getName() = "Laplacian(" + sFunctor.getName() + ")";
  const int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLaplacian3D<T> (
      _sGeometry.getBlockGeometry(iC), _sFunctor.getBlockF(iC), forthOrder ));
  }
}


////////////////////////BlockPhysLaplacian3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockPhysLaplacian3D<T,DESCRIPTOR>::BlockPhysLaplacian3D
(BlockF3D<T>& blockLaplacian, const UnitConverter<T,DESCRIPTOR>& converter,
  bool forthOrder)
  : BlockF3D<T>(blockLaplacian.getBlockStructure(), blockLaplacian.getTargetDim()),
    _blockLaplacian(blockLaplacian), _converter(converter),
    _factor(T{1} / converter.getConversionFactorLength())
{
  this->getName() = "Phys" + _blockLaplacian.getName();
}

template <typename T, typename DESCRIPTOR>
bool BlockPhysLaplacian3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockLaplacian(output,input);
  for (int i = 0; i < this->getTargetDim(); i++) {
    output[i] *= _factor * _factor;
  }
  return true;
}

////////////////////////SuperPhysLaplacian3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperPhysLaplacian3D<T,DESCRIPTOR>::SuperPhysLaplacian3D
(SuperGeometry<T,3>& sGeometry, SuperF3D<T>& sFunctor,
  const UnitConverter<T,DESCRIPTOR>& converter,
  bool forthOrder)
  : SuperF3D<T>(sFunctor.getSuperStructure(),sFunctor.getTargetDim()),
    _laplacian(sGeometry, sFunctor, forthOrder), _converter(converter)
{
  this->getName() = "Phys" + _laplacian.getName();
  const int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);

  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockPhysLaplacian3D<T,DESCRIPTOR> (
      _laplacian.getBlockF(iC), _converter, forthOrder ));
  }
}


} // end namespace olb
#endif
