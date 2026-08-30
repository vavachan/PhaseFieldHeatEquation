/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2024 Dennis Teutscher
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
#ifndef ATMBOUNDARY_H
#define ATMBOUNDARY_H


namespace olb {
/// @brief Analytical functor for setting a wind profile at the inlet. The function used is based on the equation used by OpenFoam.
/// @tparam T
template <typename T>
class AnalyticalWindProfileF3D : public AnalyticalF3D<T,T> {

protected:
  T _u0;
  T _z0;
  T _z_ref;
  T _kappa;
  T _d =0;
  int _windDirection;
  int _normalToEarth;
  T _angleRad;


public:
  AnalyticalWindProfileF3D(T u0, T z0, T z_ref,T kappa, T d, int windDirection = 0, int normalToEarth =2, T angleRad=0) : AnalyticalF3D<T,T>( 3 )
  ,_u0(u0), _z0(z0), _z_ref(z_ref), _kappa(kappa), _d(d),_windDirection(windDirection), _normalToEarth(normalToEarth), _angleRad(angleRad) {
    this->getName() = "AnalyticalWindProfileF3D";
  };

  bool operator()( T output[3], const T input[3] ) override {
    output[0] = util::cos(_angleRad)*(_u0*_kappa)/util::log((_z_ref+_z0)/_z0) / _kappa *util::log(((util::abs(input[_normalToEarth])-_d)+_z0)/_z0);
    output[1] = util::sin(_angleRad)*(_u0*_kappa)/util::log((_z_ref+_z0)/_z0) / _kappa *util::log(((util::abs(input[_normalToEarth])-_d)+_z0)/_z0);
    return true;
  };

};
/// @brief Analytical functor for setting a wind profile based on the power law at the inlet.
/// @tparam T
template <typename T>
class AnalyticalWindProfilePowerLawF3D : public AnalyticalF3D<T,T> {

protected:
  T _u_r;            /// Reference wind speed
  T _z_r;            /// Reference height
  T _alpha;          /// Power-law exponent
  T _alphaT;
  int _windDirection;
  int _normalToEarth;
  T _angleRad;
  T _I_u_r;          /// Reference turbulence intensity
  T _randomFactorStdDev;  /// Standard deviation for random wind fluctuation

  /// Random number generation
  std::default_random_engine _generator;
  std::normal_distribution<T> _distribution;

public:
  AnalyticalWindProfilePowerLawF3D(T u_r, T z_r, T alpha, T I_u_r,T alphaT , int windDirection = 0, int normalToEarth = 2, T angleRad = 0)
  : AnalyticalF3D<T,T>( 3 ), _u_r(u_r), _z_r(z_r), _alpha(alpha), _I_u_r(I_u_r), _alphaT(alphaT),
    _windDirection(windDirection), _normalToEarth(normalToEarth), _angleRad(angleRad),
    _distribution(-1, 1) {
    this->getName() = "AnalyticalWindProfileF3D";
  }

  bool operator()( T output[3], const T input[3] ) override {
    OstreamManager clout(std::cout, "Test");
      T z = util::abs(input[_normalToEarth]);
      T windSpeed = _u_r * util::pow(z / _z_r, _alpha);
      T turbulenceIntensity = _I_u_r * util::pow(z / _z_r, -_alphaT);
      output[0] = util::cos(_angleRad) * windSpeed + windSpeed*(_distribution(_generator)*turbulenceIntensity);
      output[1] = util::sin(_angleRad) * windSpeed + windSpeed*(_distribution(_generator)*turbulenceIntensity);
      output[2] = windSpeed * (_distribution(_generator)*turbulenceIntensity);
      return true;
  }
};



template <typename T>
class AnalyticalTurbulentWindProfileF3D : public AnalyticalF3D<T, T> {
protected:
    T _u0;
    T _z0;
    T _z_ref;
    T _kappa;
    T _d = 0;
    T _intensity_0;
    int _windDirection;
    int _normalToEarth;
    T _angleRad;
    std::default_random_engine generator;
    std::normal_distribution<T> distribution;

public:
    AnalyticalTurbulentWindProfileF3D(T u0, T z0, T z_ref, T kappa, T d, T intensity_0 = 0.1, int windDirection = 0, int normalToEarth = 2, T angleRad = 0)
        : AnalyticalF3D<T, T>(3),
          _u0(u0), _z0(z0), _z_ref(z_ref), _kappa(kappa), _d(d), _intensity_0(intensity_0),
          _windDirection(windDirection), _normalToEarth(normalToEarth), _angleRad(angleRad),
          distribution(0.0, 1.0) {
        this->getName() = "AnalyticalTurbulentWindProfileF3D";
    }

    bool operator()(T output[3], const T input[3]) override {
        T height = util::abs(input[_normalToEarth]) - _d + _z0;
        T log_height_ratio = util::log(height / _z0);
        T mean_wind_speed = (_u0 * _kappa) / util::log((_z_ref + _z0) / _z0) / _kappa * log_height_ratio;
        T turbulence_intensity = _intensity_0 * log_height_ratio;

        // Adding random turbulence perturbation
        T random_perturbation_x = turbulence_intensity * distribution(generator);
        T random_perturbation_y = turbulence_intensity * distribution(generator);

        output[0] = cos(_angleRad) * mean_wind_speed * (1 + random_perturbation_x);
        output[1] = sin(_angleRad) * mean_wind_speed * (1 + random_perturbation_y);
        output[2] = mean_wind_speed*random_perturbation_y;

        return true;
    }
};

}
#endif
