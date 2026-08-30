/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt
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

#ifndef INDIC_MOD_H
#define INDIC_MOD_H

#include "utilities/aliases.h"
#include "utilities/functorPtr.h"

namespace olb {

template <typename S, unsigned D>
class IndicInverse : public IndicatorF<S, D> {
protected:
  FunctorPtr<IndicatorF<S, D>> _f;

public:
  IndicInverse(FunctorPtr<IndicatorF<S, D>> f, PhysR<S, D> min,
               PhysR<S, D> max);
  IndicInverse(FunctorPtr<IndicatorF<S, D>> f);
  S signedDistance(const Vector<S, D>& input) override;
  S signedDistanceExact(const Vector<S, D>& input) override;
  Vector<S,D> surfaceNormal(const Vector<S,D>& pos, const S meshSize) override;
  Vector<S,D> surfaceNormalExact(const Vector<S,D>& pos, const S meshSize) override;
};

template <typename S, unsigned D>
class IndicScale : public IndicatorF<S, D> {
protected:
  FunctorPtr<IndicatorF<S, D>> _f;
  S                            _scalingFactor;

public:
  IndicScale(FunctorPtr<IndicatorF<S, D>> f, S scalingFactor = S {1});
  void setScalingFactor(const S scalingFactor);
  S getScalingFactor();
  Vector<S,D> getEstimatedCenter();
  S    signedDistance(const Vector<S, D>& input) override;
  S    signedDistanceExact(const Vector<S, D>& input) override;
};

template <typename S, unsigned D>
class IndicElongation : public IndicatorF<S, D> {
protected:
  FunctorPtr<IndicatorF<S, D>> _f;
  Vector<S,D>                  _elongation;

public:
  IndicElongation(FunctorPtr<IndicatorF<S, D>> f, const Vector<S,D>& elongation = Vector<S,D>(0.));
  void setElongation(const Vector<S,D>& elongation);
  Vector<S,D> getElongation();
  Vector<S,D> getEstimatedCenter();
  S    signedDistance(const Vector<S, D>& input) override;
  S    signedDistanceExact(const Vector<S, D>& input) override;
};

} // namespace olb
#endif
