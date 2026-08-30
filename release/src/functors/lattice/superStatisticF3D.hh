/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Jakob Mangold, Mathias J. Krause, 2024 Julius Jeßberger
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

#ifndef SUPER_STATISTIC_F3D_HH
#define SUPER_STATISTIC_F3D_HH

#include "superStatisticF3D.h"
#include "indicator/superIndicatorF3D.h"

namespace olb {


template <typename T, typename W>
SuperVarianceF3D<T,W>::SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                        FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                                        T expectedValue)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()+1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF)),
    _expectedValue(expectedValue),
    _expectedValueF(this->getSuperStructure(), _expectedValue)
{
  this->getName() = "Variance("+_f->getName()+")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockVarianceF3D<T,W>(_f->getBlockF(iC),
                                  _indicatorF->getBlockIndicatorF(iC),
                                  _expectedValue)
      );
    }
  }
}

template <typename T, typename W>
SuperVarianceF3D<T,W>::SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                        SuperGeometry<T,3>& superGeometry,
                                        const int material,
                                        T expectedValue)
  : SuperVarianceF3D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material),
      expectedValue)
{ }

template <typename T, typename W>
bool SuperVarianceF3D<T,W>::operator() (W output[], const int input[])
{
  return SuperAverage3D<T,W>((*_f-_expectedValueF)*(*_f-_expectedValueF),*_indicatorF)(output, input);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template <typename T, typename W>
SuperStdDeviationF3D<T,W>::SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
    T expectedValue)
  : SuperVarianceF3D<T,W>(std::forward<decltype(f)>(f), std::forward<decltype(indicatorF)>(indicatorF), expectedValue)
{
  this->getName() = "StdDeviation("+this->_f->getName()+")";

  LoadBalancer<T>& load = this->_f->getSuperStructure().getLoadBalancer();

  if ( this->_f->getBlockFSize()          == load.size() &&
       this->_indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockStdDeviationF3D<T,W>(this->_f->getBlockF(iC),
                                      this->_indicatorF->getBlockIndicatorF(iC),
                                      this->_expectedValue)
      );
    }
  }
}

template <typename T, typename W>
SuperStdDeviationF3D<T,W>::SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
    SuperGeometry<T,3>& superGeometry,
    const int material,
    T expectedValue)
  : SuperStdDeviationF3D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material),
      expectedValue)
{ }

template <typename T, typename W>
bool SuperStdDeviationF3D<T,W>::operator() (W output[], const int input[])
{
  const bool res = SuperVarianceF3D<T,W>::operator()(output, input);
  for (int i=0; i<this->getTargetDim(); ++i) {
    output[i] = util::sqrt(output[i]);
  }
  return res;
}



}

#endif

