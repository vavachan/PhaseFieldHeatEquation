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

#ifndef BLOCK_STATISTIC_F3D_HH
#define BLOCK_STATISTIC_F3D_HH

#include "blockStatisticF3D.h"

namespace olb {


template <typename T, typename W>
BlockVarianceF3D<T,W>::BlockVarianceF3D(BlockF3D<W>&          f,
                                        BlockIndicatorF3D<T>& indicatorF,
                                        T expectedValue)
  : BlockF3D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _expectedValue(BlockConst3D<T>(this->getBlockStructure(), expectedValue))
{
  this->getName() = "BlockVariance("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockVarianceF3D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  return BlockAverage3D<T,W>((_f-_expectedValue)*(_f-_expectedValue),_indicatorF)(output, input);
}


template <typename T, typename W>
BlockStdDeviationF3D<T,W>::BlockStdDeviationF3D(BlockF3D<W>&          f,
    BlockIndicatorF3D<T>& indicatorF,
    T expectedValue)
  : BlockF3D<W>(f.getBlockStructure(), f.getTargetDim()),
    _variance(BlockVarianceF3D<T,W>(f,indicatorF,expectedValue))
{
  this->getName() = "BlockStdDeviation("+f.getName()+")";
}

template <typename T, typename W>
bool BlockStdDeviationF3D<T,W>::operator() (W output[], const int input[])
{
  const bool res = _variance(output, input);
  for (int i=0; i<_variance.getTargetDim(); ++i) {
    output[i] = util::sqrt(output[i]);
  }
  return res;
}

}

#endif
