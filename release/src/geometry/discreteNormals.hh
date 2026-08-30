/**  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Mathias Krause, Simon Zimny, Adrian Kummerlaender
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

#ifndef GEOMETRY_DISCRETE_NORMAL_HH
#define GEOMETRY_DISCRETE_NORMAL_HH

namespace olb {

/// Returns type (e.g. edge / corner) and discrete normal in 2D
template<concepts::BaseType T>
std::pair<DiscreteNormalType,Vector<int,2>> computeBoundaryTypeAndNormal(
  BlockIndicatorF2D<T>& fluidI,
  BlockIndicatorF2D<T>& outsideI,
  Vector<int,2> latticeR)
{
  const int iX = latticeR[0];
  const int iY = latticeR[1];

  std::vector<int> discreteNormal(3, 0);

  if (!fluidI({iX, iY})
      && !outsideI({iX, iY})) {

    ///boundary0N and boundary 0P
    if (!fluidI({iX, iY + 1})
        && !outsideI({iX, iY + 1})
        && !fluidI({iX, iY - 1})
        && !outsideI({iX, iY - 1})) {

      if (fluidI({iX + 1, iY})) {
        discreteNormal[0] = 0;
        discreteNormal[1] = -1;
        discreteNormal[2] = 0;
      }

      if (fluidI({iX - 1, iY})) {
        discreteNormal[0] = 0;
        discreteNormal[1] = 1;
        discreteNormal[2] = 0;
      }
    }

    /// boundary1N and boundary1P
    if (!fluidI({iX + 1, iY})
        && !outsideI({iX + 1, iY})
        && !fluidI({iX - 1, iY})
        && !outsideI({iX - 1, iY})) {

      if (fluidI({iX, iY + 1})) {
        discreteNormal[0] = 0;
        discreteNormal[1] = 0;
        discreteNormal[2] = -1;
      }

      if (fluidI({iX, iY - 1})) {
        discreteNormal[0] = 0;
        discreteNormal[1] = 0;
        discreteNormal[2] = 1;
      }
    }

    /// externalCornerNN and externalCornerNP
    if (   !fluidI({iX + 1, iY})
        && !outsideI({iX + 1, iY})) {

      if (!fluidI({iX, iY + 1})
          && !outsideI({iX, iY + 1})
          && fluidI({iX + 1, iY + 1})) {
        discreteNormal[0] = 1;
        discreteNormal[1] = -1;
        discreteNormal[2] = -1;
      }

      if (!fluidI({iX, iY - 1})
          && !outsideI({iX, iY - 1})
          && fluidI({iX + 1, iY - 1})) {
        discreteNormal[0] = 1;
        discreteNormal[1] = -1;
        discreteNormal[2] = 1;
      }
    }

    /// externalCornerPN and externalCornerPP
    if (!fluidI({iX - 1, iY})
        && !outsideI({iX - 1, iY})) {

      if (!fluidI({iX, iY + 1})
          && !outsideI({iX, iY + 1})
          && fluidI({iX - 1, iY + 1})) {
        discreteNormal[0] = 1;
        discreteNormal[1] = 1;
        discreteNormal[2] = -1;
      }

      if (!fluidI({iX, iY - 1})
          && !outsideI({iX, iY - 1})
          && fluidI({iX - 1, iY - 1})) {
        discreteNormal[0] = 1;
        discreteNormal[1] = 1;
        discreteNormal[2] = 1;
      }
    }

    /// internalCornerNN and internalCornerNP
    if (!fluidI({iX - 1, iY})
        && !outsideI({iX - 1, iY})) {

      if (!fluidI({iX, iY - 1})
          && !outsideI({iX, iY - 1})
          && (outsideI({iX - 1, iY - 1}) || !fluidI({iX - 1, iY - 1}))
          && fluidI({iX + 1, iY + 1})) {
        discreteNormal[0] = 2;
        discreteNormal[1] = -1;
        discreteNormal[2] = -1;
      }

      if (!fluidI({iX, iY + 1})
          && !outsideI({iX, iY + 1})
          && (outsideI({iX - 1, iY + 1}) || !fluidI({iX - 1, iY + 1}))
          && fluidI({iX + 1, iY - 1})) {
        discreteNormal[0] = 2;
        discreteNormal[1] = -1;
        discreteNormal[2] = 1;
      }
    }

    /// internalCornerPN and internalCornerPP
    if (!fluidI({iX + 1, iY})
        && !outsideI({iX + 1, iY})) {

      if (!fluidI({iX, iY - 1})
          && !outsideI({iX, iY - 1})
          && (outsideI({iX + 1, iY - 1}) || !fluidI({iX + 1, iY - 1}))
          && fluidI({iX - 1, iY + 1})) {
        discreteNormal[0] = 2;
        discreteNormal[1] = 1;
        discreteNormal[2] = -1;
      }

      if (!fluidI({iX, iY + 1})
          && !outsideI({iX, iY + 1})
          && (outsideI({iX + 1, iY + 1}) || !fluidI({iX + 1, iY + 1}))
          && fluidI({iX - 1, iY - 1})) {
        discreteNormal[0] = 2;
        discreteNormal[1] = 1;
        discreteNormal[2] = 1;
      }
    }
  }

  return {static_cast<DiscreteNormalType>(discreteNormal[0]), discreteNormal.data()+1};
}

/// Returns type (e.g. edge / corner) and discrete normal in 3D
template<concepts::BaseType T>
std::pair<DiscreteNormalType,Vector<int,3>> computeBoundaryTypeAndNormal(
  BlockIndicatorF3D<T>& fluidI,
  BlockIndicatorF3D<T>& outsideI,
  Vector<int,3> latticeR)
{
  const int iX = latticeR[0];
  const int iY = latticeR[1];
  const int iZ = latticeR[2];

  std::vector<int> discreteNormal(4,0);
  std::vector<int> discreteNormal2(4,0);
  std::vector<int> nullVector(4,0);

  if (   !fluidI({iX, iY, iZ})
      && !outsideI({iX, iY, iZ})) {

    //boundary0N and boundary 0P
    if ( !fluidI({iX, iY, iZ + 1})
      && !outsideI({iX, iY, iZ + 1})
      && !fluidI({iX, iY, iZ - 1})
      && !outsideI({iX, iY, iZ - 1})
      && !fluidI({iX, iY + 1, iZ})
      && !outsideI({iX, iY + 1, iZ})
      && !fluidI({iX, iY - 1, iZ})
      && !outsideI({iX, iY - 1, iZ})) {

        if (fluidI({iX + 1, iY, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 0;
        }
      }

        if (fluidI({iX - 1, iY, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 0;
        }
      }
    }

    // boundary1N and boundary1P
    if (   !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})) {

      if (fluidI({iX, iY + 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }

      if (fluidI({iX, iY - 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // boundary2N and boundary2P
    if (!fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})) {

      if (fluidI({iX, iY, iZ + 1})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }

      if (fluidI({iX, iY, iZ - 1})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 0;
          discreteNormal[1] = 0;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 0;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerNNN and externalCornerNPN
    if (   !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX + 1, iY, iZ + 1})
        && !outsideI({iX + 1, iY, iZ + 1})) {

      if (   !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX + 1, iY + 1, iZ})
          && !outsideI({iX + 1, iY + 1, iZ})
          && !fluidI({iX, iY + 1, iZ + 1})
          && !outsideI({iX, iY + 1, iZ + 1})
          && fluidI({iX + 1, iY + 1, iZ + 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (   !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX + 1, iY - 1, iZ})
          && !outsideI({iX + 1, iY - 1, iZ})
          && !fluidI({iX, iY - 1, iZ + 1})
          && !outsideI({iX, iY - 1, iZ + 1})
          && fluidI({iX + 1, iY - 1, iZ + 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalCornerNPP and externalCornerNNP
    if (   !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX + 1, iY, iZ - 1})
        && !outsideI({iX + 1, iY, iZ - 1})) {

      if (   !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX + 1, iY - 1, iZ})
          && !outsideI({iX + 1, iY - 1, iZ})
          && !fluidI({iX, iY - 1, iZ - 1})
          && !outsideI({iX, iY - 1, iZ - 1})
          && fluidI({iX + 1, iY - 1, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (   !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX + 1, iY + 1, iZ})
          && !outsideI({iX + 1, iY + 1, iZ})
          && !fluidI({iX, iY + 1, iZ - 1})
          && !outsideI({iX, iY + 1, iZ - 1})
          && fluidI({iX + 1, iY + 1, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerPPP and externalCornerPNP
    if (   !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX - 1, iY, iZ - 1})
        && !outsideI({iX - 1, iY, iZ - 1})) {

      if (   !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX, iY - 1, iZ - 1})
          && !outsideI({iX, iY - 1, iZ - 1})
          && !fluidI({iX - 1, iY - 1, iZ})
          && !outsideI({iX - 1, iY - 1, iZ})
          && fluidI({iX - 1, iY - 1, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (   !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX, iY + 1, iZ - 1})
          && !outsideI({iX, iY + 1, iZ - 1})
          && !fluidI({iX - 1, iY + 1, iZ})
          && !outsideI({iX - 1, iY + 1, iZ})
          && fluidI({iX - 1, iY + 1, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalCornerPNN and externalCornerPPN
    if (   !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX - 1, iY, iZ + 1})
        && !outsideI({iX - 1, iY, iZ + 1})) {

      if (   !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX, iY + 1, iZ + 1})
          && !outsideI({iX, iY + 1, iZ + 1})
          && !fluidI({iX - 1, iY + 1, iZ})
          && !outsideI({iX - 1, iY + 1, iZ})
          && fluidI({iX - 1, iY + 1, iZ + 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (   !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX, iY - 1, iZ + 1})
          && !outsideI({iX, iY - 1, iZ + 1})
          && !fluidI({iX - 1, iY - 1, iZ})
          && !outsideI({iX - 1, iY - 1, iZ})
          && fluidI({iX - 1, iY - 1, iZ + 1})) {

        if (discreteNormal == nullVector) {

          discreteNormal[0] = 1;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }

        else {
          discreteNormal2[0] = 1;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // internalCornerPPP and internalCornerPNP
    if (   fluidI({iX - 1, iY, iZ})
        && fluidI({iX, iY, iZ - 1})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})) {

      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // internalCornerPNN and InternalCornerPPN
    if (   fluidI({iX - 1, iY, iZ})
        && fluidI({iX, iY, iZ + 1})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})) {

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }

      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // internalCornerNPP and internalCornerNNP
    if (   fluidI({iX + 1, iY, iZ})
        && fluidI({iX, iY, iZ - 1})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})) {

      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // internalCornerNPN and internalCornerNNN
    if (   fluidI({iX + 1, iY, iZ})
        && fluidI({iX, iY, iZ + 1})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})) {

      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }

        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 2;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 2;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge0PN and externalEdge0NN
    if (   !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX + 1, iY, iZ + 1})
        && !fluidI({iX - 1, iY, iZ + 1})) {

      if (   fluidI({iX, iY - 1, iZ + 1})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = -1;
        }
      }

      if (   fluidI({iX, iY + 1, iZ + 1})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge0NP and externalEdge0PP
    if (   !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX + 1, iY, iZ - 1})
        && !fluidI({iX - 1, iY, iZ - 1})) {

      if (   fluidI({iX, iY + 1, iZ - 1})
          && !fluidI({iX, iY + 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = -1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 1;
        }
      }

      if (   fluidI({iX, iY - 1, iZ - 1})
          && !fluidI({iX, iY - 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 0;
          discreteNormal[2] = 1;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 0;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalEdge1NN and externalEdge1NP
    if (   !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})) {

      if (   fluidI({iX + 1, iY, iZ + 1})
          && !fluidI({iX + 1, iY, iZ})
          && !outsideI({iX + 1, iY, iZ})
          && !fluidI({iX - 1, iY, iZ})
          && !fluidI({iX, iY, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }

      if (   fluidI({iX - 1, iY, iZ + 1})
          && !fluidI({iX - 1, iY, iZ})
          && !outsideI({iX - 1, iY, iZ})
          && !fluidI({iX + 1, iY, iZ})
          && !fluidI({iX, iY, iZ - 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = -1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = -1;
        }
      }
    }

    // externalEdge1PN and externalEdge1PP
    if (   !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})) {

      if (   fluidI({iX + 1, iY, iZ - 1})
          && !fluidI({iX + 1, iY, iZ})
          && !outsideI({iX + 1, iY, iZ})
          && !fluidI({iX - 1, iY, iZ})
          && !fluidI({iX, iY, iZ + 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }

      if (   fluidI({iX - 1, iY, iZ - 1})
          && !fluidI({iX - 1, iY, iZ})
          && !outsideI({iX - 1, iY, iZ})
          && !fluidI({iX + 1, iY, iZ})
          && !fluidI({iX, iY, iZ + 1})) {

        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 0;
          discreteNormal[3] = 1;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 0;
          discreteNormal2[3] = 1;
        }
      }
    }

    // externalEdge2NN and externalEdge2PN
    if (   !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})) {

      if (   !fluidI({iX + 1, iY, iZ})
          && fluidI({iX + 1, iY + 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }

      if (   !fluidI({iX - 1, iY, iZ})
          && fluidI({iX - 1, iY + 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = -1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = -1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // externalEdge2PP and externalEdge2NP
    if (   !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})) {

      if (   !fluidI({iX - 1, iY, iZ})
          && fluidI({iX - 1, iY - 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = 1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = 1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }

      if (   !fluidI({iX + 1, iY, iZ})
          && fluidI({iX + 1, iY - 1, iZ})) {
        if (discreteNormal == nullVector) {
          discreteNormal[0] = 3;
          discreteNormal[1] = -1;
          discreteNormal[2] = 1;
          discreteNormal[3] = 0;
        }
        else {
          discreteNormal2[0] = 3;
          discreteNormal2[1] = -1;
          discreteNormal2[2] = 1;
          discreteNormal2[3] = 0;
        }
      }
    }

    // internalEdge0NN and internalEdge0PN
    if (   !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && fluidI({iX, iY, iZ + 1})) {

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX - 1, iY - 1, iZ})
          && !outsideI({iX - 1, iY - 1, iZ})
          && !fluidI({iX + 1, iY - 1, iZ})
          && !outsideI({iX + 1, iY - 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = -1;
        discreteNormal[3] = -1;
      }
      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX - 1, iY + 1, iZ})
          && !outsideI({iX - 1, iY + 1, iZ})
          && !fluidI({iX + 1, iY + 1, iZ})
          && !outsideI({iX + 1, iY + 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = 1;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge0NP and internalEdge0PP
    if (   !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && fluidI({iX, iY, iZ - 1})) {

      if (   fluidI({iX, iY + 1, iZ})
          && !fluidI({iX - 1, iY - 1, iZ})
          && !outsideI({iX - 1, iY - 1, iZ})
          && !fluidI({iX, iY - 1, iZ})
          && !outsideI({iX, iY - 1, iZ})
          && !fluidI({iX + 1, iY - 1, iZ})
          && !outsideI({iX + 1, iY - 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = -1;
        discreteNormal[3] = 1;
      }

      if (   fluidI({iX, iY - 1, iZ})
          && !fluidI({iX - 1, iY + 1, iZ})
          && !outsideI({iX - 1, iY + 1, iZ})
          && !fluidI({iX, iY + 1, iZ})
          && !outsideI({iX, iY + 1, iZ})
          && !fluidI({iX + 1, iY + 1, iZ})
          && !outsideI({iX + 1, iY + 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 0;
        discreteNormal[2] = 1;
        discreteNormal[3] = 1;
      }
    }

    // internalEdge1PP and internalEdge 1NP
    if (   fluidI({iX - 1, iY, iZ})
        && !fluidI({iX + 1, iY, iZ})
        && !outsideI({iX + 1, iY, iZ})
        && !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})) {

      if (   fluidI({iX, iY, iZ - 1})
          && !fluidI({iX, iY, iZ + 1})
          && !outsideI({iX, iY, iZ + 1})
          && !fluidI({iX, iY + 1, iZ + 1})
          && !outsideI({iX, iY + 1, iZ + 1})
          && !fluidI({iX, iY - 1, iZ + 1})
          && !outsideI({iX, iY - 1, iZ + 1})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 0;
        discreteNormal[3] = 1;
      }

      if (   fluidI({iX, iY, iZ + 1})
          && !fluidI({iX, iY, iZ - 1})
          && !outsideI({iX, iY, iZ - 1})
          && !fluidI({iX, iY + 1, iZ - 1})
          && !outsideI({iX, iY + 1, iZ - 1})
          && !fluidI({iX, iY - 1, iZ - 1})
          && !outsideI({iX, iY - 1, iZ - 1})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 0;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge1PN and internalEdge1NN
    if (   fluidI({iX + 1, iY, iZ})
        && !fluidI({iX - 1, iY, iZ})
        && !outsideI({iX - 1, iY, iZ})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})
        && !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})) {

      if (   fluidI({iX, iY, iZ - 1})
          && !fluidI({iX, iY, iZ + 1})
          && !outsideI({iX, iY, iZ + 1})
          && !fluidI({iX, iY + 1, iZ + 1})
          && !outsideI({iX, iY + 1, iZ + 1})
          && !fluidI({iX, iY - 1, iZ + 1})
          && !outsideI({iX, iY - 1, iZ + 1})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 0;
        discreteNormal[3] = 1;
      }

      if (   fluidI({iX, iY, iZ + 1})
          && !fluidI({iX, iY, iZ - 1})
          && !outsideI({iX, iY, iZ - 1})
          && !fluidI({iX, iY + 1, iZ - 1})
          && !outsideI({iX, iY + 1, iZ - 1})
          && !fluidI({iX, iY - 1, iZ - 1})
          && !outsideI({iX, iY - 1, iZ - 1})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 0;
        discreteNormal[3] = -1;
      }
    }

    // internalEdge2PP and internalEdge2NP
    if (   fluidI({iX, iY - 1, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX, iY + 1, iZ})
        && !outsideI({iX, iY + 1, iZ})) {

      if (fluidI({iX - 1, iY, iZ})
          && fluidI({iX - 1, iY - 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = 1;
        discreteNormal[3] = 0;
      }

      if (   fluidI({iX + 1, iY, iZ})
          && fluidI({iX + 1, iY - 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = 1;
        discreteNormal[3] = 0;
      }
    }

    // internalEdge2PN and internalEdge2NN
    if (   fluidI({iX, iY + 1, iZ})
        && !fluidI({iX, iY, iZ - 1})
        && !outsideI({iX, iY, iZ - 1})
        && !fluidI({iX, iY, iZ + 1})
        && !outsideI({iX, iY, iZ + 1})
        && !fluidI({iX, iY - 1, iZ})
        && !outsideI({iX, iY - 1, iZ})) {

      if (   fluidI({iX - 1, iY, iZ})
          && fluidI({iX - 1, iY + 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = 1;
        discreteNormal[2] = -1;
        discreteNormal[3] = 0;
      }

      if (   fluidI({iX + 1, iY, iZ})
          && fluidI({iX + 1, iY + 1, iZ})) {

        discreteNormal[0] = 4;
        discreteNormal[1] = -1;
        discreteNormal[2] = -1;
        discreteNormal[3] = 0;
      }
    }
  }
  return {static_cast<DiscreteNormalType>(discreteNormal[0]), discreteNormal.data()+1};
}

}

#endif
