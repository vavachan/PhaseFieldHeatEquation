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

#ifndef INDIC_MOD_HH
#define INDIC_MOD_HH

#include "indicMod.h"

namespace olb {

template <typename S, unsigned D>
IndicInverse<S, D>::IndicInverse(FunctorPtr<IndicatorF<S, D>> f,
                                 PhysR<S, D> min, PhysR<S, D> max)
    : _f(std::move(f))
{
  this->_myMin = min;
  this->_myMax = max;
}


/// Alternative constructor enabling quick inversion
/**
 * \warning Do to the neglection of specifying min and max, no volume but rather a plane is created!
 * \warning Might lead to errors, when using min/max checks assuming a volume.
 **/
template <typename S, unsigned D>
IndicInverse<S, D>::IndicInverse(FunctorPtr<IndicatorF<S, D>> f)
  : _f(std::move(f))
{
  this->_myMin = f->getMin();
  this->_myMax = f->getMax();
}

template <typename S, unsigned D>
S IndicInverse<S, D>::signedDistance(const Vector<S, D>& input)
{
  return -this->_f->signedDistance(input);
}

template <typename S, unsigned D>
S IndicInverse<S, D>::signedDistanceExact(const Vector<S, D>& input)
{
  return -this->_f->signedDistanceExact(input);
}

template <typename S, unsigned D>
Vector<S,D> IndicInverse<S, D>::surfaceNormal(const Vector<S,D>& pos, const S meshSize)
{
  return this->_f->surfaceNormal(pos, meshSize) * -1;
}

template <typename S, unsigned D>
Vector<S,D> IndicInverse<S, D>::surfaceNormalExact(const Vector<S,D>& pos, const S meshSize)
{
  return this->_f->surfaceNormalExact(pos, meshSize) * -1;
}

template <typename S, unsigned D>
IndicScale<S, D>::IndicScale(FunctorPtr<IndicatorF<S, D>> f, S scalingFactor)
    : _f(std::move(f))
{
  setScalingFactor(scalingFactor);
}

template <typename S, unsigned D>
void IndicScale<S, D>::setScalingFactor(S scalingFactor)
{
  _scalingFactor = scalingFactor;
  this->_myMin = _f->getMin() - scalingFactor*_f->getMin();
  this->_myMax = _f->getMax() + scalingFactor*_f->getMax();
}

template <typename S, unsigned D>
S IndicScale<S, D>::getScalingFactor()
{
  return _scalingFactor;
}

template <typename S, unsigned D>
Vector<S,D> IndicScale<S, D>::getEstimatedCenter()
{
  return 0.5 * (_f->getMin() + _f->getMax());
}


template <typename S, unsigned D>
S IndicScale<S, D>::signedDistance(const Vector<S, D>& input)
{
  std::function<S(const Vector<S, D>&)> sdf =
      [this](const Vector<S, D>& input) {
        return this->_f->signedDistance(input);
      };

  return sdf::scale(sdf, input, _scalingFactor);
}

template <typename S, unsigned D>
S IndicScale<S, D>::signedDistanceExact(const Vector<S, D>& input)
{
  std::function<S(const Vector<S, D>&)> sdf =
      [this](const Vector<S, D>& input) {
        return this->_f->signedDistanceExact(input);
      };

  return sdf::scale(sdf, input, _scalingFactor);
}

template <typename S, unsigned D>
IndicElongation<S, D>::IndicElongation(FunctorPtr<IndicatorF<S, D>> f, const Vector<S,D>& elongation)
    : _f(std::move(f))
{
  setElongation(elongation);
}

template <typename S, unsigned D>
void IndicElongation<S, D>::setElongation(const Vector<S,D>& elongation)
{
  _elongation = elongation;
  this->_myMin = _f->getMin() - elongation;
  this->_myMax = _f->getMax() + elongation;
}

template <typename S, unsigned D>
Vector<S,D> IndicElongation<S, D>::getElongation()
{
  return _elongation;
}

template <typename S, unsigned D>
Vector<S,D> IndicElongation<S, D>::getEstimatedCenter()
{
  return 0.5 * (_f->getMin() + _f->getMax());
}

template <typename S, unsigned D>
S IndicElongation<S, D>::signedDistance(const Vector<S, D>& input)
{
  std::function<S(const Vector<S, D>&)> sdf =
      [this](const Vector<S, D>& input) {
        return this->_f->signedDistance(input);
      };

  return sdf::elongation(sdf, input, _elongation, getEstimatedCenter());
}

template <typename S, unsigned D>
S IndicElongation<S, D>::signedDistanceExact(const Vector<S, D>& input)
{
  std::function<S(const Vector<S, D>&)> sdf =
      [this](const Vector<S, D>& input) {
        return this->_f->signedDistanceExact(input);
      };

  return sdf::elongation(sdf, input, _elongation, getEstimatedCenter());
}

} // namespace olb
#endif
