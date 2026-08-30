/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef INDICATOR_FROM_BLOCKDATA_F_3D_HH
#define INDICATOR_FROM_BLOCKDATA_F_3D_HH

namespace olb {

template <typename S>
IndicatorBlockData3D<S>::IndicatorBlockData3D(
  BlockData<3,S,S>& blockData,
  Vector<S,3> extend, Vector<S,3> origin,
  S deltaR, bool invert)
  : _blockData(blockData), _deltaR(deltaR), _invert(invert)
{
  this->_myMin = origin;
  this->_myMax = origin + extend;
}


template <typename S>
bool IndicatorBlockData3D<S>::operator()(bool output[], const S input[])
{

  // Translation
  S xDist = input[0] - this->_myMin[0];
  S yDist = input[1] - this->_myMin[1];
  S zDist = input[2] - this->_myMin[2];

  int x = ((this->_myMin[0] + xDist)/_deltaR)+0.5;
  int y = ((this->_myMin[1] + yDist)/_deltaR)+0.5;
  int z = ((this->_myMin[2] + zDist)/_deltaR)+0.5;

  if (x >= 0 && x < _blockData.getNx() && y >= 0 && y < _blockData.getNy() && z >= 0 && z < _blockData.getNz()) {
    LatticeR<3> input(x,y,z);
    if (this->_blockData.get(input) > std::numeric_limits<S>::epsilon()) {
      if (!_invert){
        output[0] = S(this->_blockData.get(input));
        return true;
      } else {
        output[0] = S(0);
        return false;
      }
    }
  }
  if (!_invert){
    output[0] = S(0);
    return false;
  } else {
    output[0] = 1.-S(this->_blockData.get(input));
    return true;
  }
}


} // namespace olb


#endif
