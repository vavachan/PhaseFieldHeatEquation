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

#ifndef BLOCK_STATISTIC_F3D_H
#define BLOCK_STATISTIC_F3D_H

#include "geometry/cuboid.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "blockAverage3D.h"
#include "blockBaseF3D.h"


namespace olb {


/// BlockVarianceF3D returns the variance in each component of f on a indicated subset
template <typename T, typename W = T>
class BlockVarianceF3D final : public BlockF3D<W> {
private:
  BlockF3D<W>&          _f;
  BlockIndicatorF3D<T>& _indicatorF;
  BlockConst3D<T> _expectedValue;

public:
  BlockVarianceF3D(BlockF3D<W>&          f,
                   BlockIndicatorF3D<T>& indicatorF,
                   T expectedValue);
  bool operator() (W output[], const int input[]) override;
};

/// BlockStdDeviationF3D returns the standard deviation in each component of f on a indicated subset
template <typename T, typename W = T>
class BlockStdDeviationF3D final : public BlockF3D<W> {
private:
  BlockVarianceF3D<T,W> _variance;
public:
  BlockStdDeviationF3D(BlockF3D<W>&          f,
                       BlockIndicatorF3D<T>& indicatorF,
                       T expectedValue);
  bool operator() (W output[], const int input[]) override;
};

}

#endif
