/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
 *                2025 Adrian Kummerlaender
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

#ifndef CUBOID_3D_HH
#define CUBOID_3D_HH

#include "cuboid.h"

namespace olb {

template <typename T, unsigned D>
Cuboid<T,D>& Cuboid<T,D>::operator=(const Cuboid<T,D>& rhs) {
  _origin = rhs._origin;
  _extent = rhs._extent;
  _delta = rhs._delta;
  _weight = rhs._weight;
  return *this;
}

template <typename T, unsigned D>
std::optional<LatticeR<D>> Cuboid<T,D>::getCloseLatticeR(Vector<T,D> physR, T eps) const {
  auto globR = _origin;
  auto physLatticeR = (physR - globR) / _delta;
  if (util::abs(util::round(physLatticeR) - physLatticeR) <= eps) {
    return LatticeR<D>(util::round(physLatticeR));
  } else {
    return std::nullopt;
  }
}

template <typename T, unsigned D>
T Cuboid<T,D>::getPhysVolume() const {
  T physVolume = util::pow(_delta,D);
  for (unsigned iD=0; iD < D; ++iD) {
    physVolume *= _extent[iD];
  }
  return physVolume;
}

template <typename T, unsigned D>
std::size_t Cuboid<T,D>::getLatticeVolume() const {
  std::size_t latticeVolume = _extent[0];
  for (unsigned iD=1; iD < D; ++iD) {
    latticeVolume *= _extent[iD];
  }
  return latticeVolume;
}

template <typename T, unsigned D>
T Cuboid<T,D>::getPhysPerimeter() const requires (D == 2) {
  return 2*_extent[1]*_delta + 2*_extent[0]*_delta;
}

template <typename T, unsigned D>
T Cuboid<T,D>::getPhysPerimeter() const requires (D == 3) {
  const T area = util::pow(_delta,2);
  return 2*area*(_extent[0]*_extent[1] + _extent[1]*_extent[2] + _extent[2]*_extent[0]);
}

template <typename T, unsigned D>
std::size_t Cuboid<T,D>::getLatticePerimeter() const requires (D == 2) {
  return 2*_extent[1] + 2*_extent[0] - 4;
}

template <typename T, unsigned D>
std::size_t Cuboid<T,D>::getLatticePerimeter() const requires (D == 3) {
  return 2*((_extent[0]-1)*(_extent[1]-1) + (_extent[1]-1)*(_extent[2]-1) + (_extent[2]-1)*(_extent[0]-1));
}

template <typename T, unsigned D>
void Cuboid<T,D>::refine(int factor) {
  if (factor < 1) {
    throw std::invalid_argument("refinement factor must be >= 1");
  } else if (factor == 2) {
    _delta /= factor;
    _extent *= factor;
    _extent -= 1;
    _weight *= util::pow(factor,D);
  } else if (factor != 1) {
    throw std::invalid_argument("TBD refinement factor must be == 2");
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::write(std::ostream& cout) const
{
  cout << "--------Cuboid Details----------" << std::endl;
  cout << " Corner: " << "\t" << _origin << std::endl;
  cout << " Delta: " << "\t" << "\t" << getDeltaR() << std::endl;
  cout << " Perimeter: " << "\t" << "\t" << getPhysPerimeter() << std::endl;
  cout << " Volume: " << "\t" << "\t" << getPhysVolume() << std::endl;
  cout << " Extent: " << "\t" << _extent << std::endl;
  cout << " Nodes at Perimeter: " << "\t" << getLatticePerimeter() << std::endl;
  cout << " Nodes in Volume: " << "\t" << getLatticeVolume() << std::endl;
  cout << " Nodes in Indicator: " << "\t" << getWeight() << std::endl;
  cout << " Other Corner: "  << "\t"<< _origin + (_extent-0.5)*_delta << std::endl;
  cout << "--------------------------------" << std::endl;
}

template <typename T, unsigned D>
void Cuboid<T,D>::print() const
{
  OstreamManager clout(std::cout, "cuboid");
  write(clout);
}

template <typename T, unsigned D>
bool Cuboid<T,D>::operator==(const Cuboid<T,D>& rhs) const
{
  return    util::nearZero(_origin - rhs._origin)
         && _extent == rhs._extent
         && _weight == rhs._weight;
}

template <typename T, unsigned D>
bool Cuboid<T,D>::isInside(Vector<T,D> pos, int overlap) const
{
  return _origin <= pos + overlap*_delta + _delta / 2
      && _origin + (_extent+overlap) * _delta > pos + _delta / 2;
}

template <typename T, unsigned D>
bool Cuboid<T,D>::intersects(Vector<T,D> min, Vector<T,D> max, int overlap) const
{
  auto loc0 = maxv(_origin-overlap*_delta, min);
  auto loc1 = minv(_origin+(_extent+overlap-1)*_delta, max);
  return loc1 >= loc0;

}

template <typename T, unsigned D>
bool Cuboid<T,D>::intersects(const Cuboid<T,D>& cuboid) const
{
  return intersects(cuboid.getOrigin(),
                    cuboid.getOrigin() + cuboid.getDeltaR()*(cuboid.getExtent()-1),
                    0);
}

template <typename T, unsigned D>
void Cuboid<T,D>::divide(Vector<int,D> division, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 2)
{
  T globPosX_child, globPosY_child;
  int xN_child = 0;
  int yN_child = 0;

  globPosX_child = _origin[0];
  globPosY_child = _origin[1];

  for (int iX=0; iX < division[0]; iX++) {
    for (int iY=0; iY < division[1]; iY++) {
      xN_child       = (_extent[0]+division[0]-iX-1)/division[0];
      yN_child       = (_extent[1]+division[1]-iY-1)/division[1];
      Cuboid2D<T> child({globPosX_child, globPosY_child}, _delta, {xN_child, yN_child});
      childrenC.push_back(child);
      globPosY_child += yN_child*_delta;
    }
    globPosY_child = _origin[1];
    globPosX_child += xN_child*_delta;
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::divide(Vector<int,D> division, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 3)
{
  int xN_child = 0;
  int yN_child = 0;
  int zN_child = 0;

  T globPosX_child = _origin[0];
  T globPosY_child = _origin[1];
  T globPosZ_child = _origin[2];

  for (int iX=0; iX < division[0]; iX++) {
    for (int iY=0; iY < division[1]; iY++) {
      for (int iZ=0; iZ < division[2]; iZ++) {
        xN_child       = (_extent[0]+division[0]-iX-1)/division[0];
        yN_child       = (_extent[1]+division[1]-iY-1)/division[1];
        zN_child       = (_extent[2]+division[2]-iZ-1)/division[2];
        Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                          _delta, {xN_child, yN_child, zN_child});
        childrenC.push_back(child);
        globPosZ_child += zN_child*_delta;
      }
      globPosZ_child = _origin[2];
      globPosY_child += yN_child*_delta;
    }
    globPosY_child = _origin[1];
    globPosX_child += xN_child*_delta;
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::divideP(int p, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 2)
{
  int nX = 0;
  int nY = 0;
  T ratio;
  T bestRatio = (T)_extent[0]/(T)_extent[1];
  T difRatio = util::fabs(bestRatio - 1) + 1;
  for (int i=1; i<= p; i++) {
    int j = p / i;
    if (i*j<=p) {
      if ( util::fabs(bestRatio - (T)i/(T)j) <= difRatio) {
        difRatio = util::fabs(bestRatio - (T)i/(T)j);
        nX = i;
        nY = j;
      }
    }
  }

  ratio = T(nX)/(T)nY;
  int rest = p - nX*nY;

  if (rest==0) {
    divide({nX,nY},childrenC);
    return;
  }

  if (ratio < bestRatio && (nY-rest) >= 0) {
    int n_QNoInsertions = nX*(nY-rest);
    T bestVolume_QNoInsertions = (T)_extent[0]*_extent[1] * n_QNoInsertions/(T)p;
    int yN_QNoInsertions = (int)(bestVolume_QNoInsertions / (T)_extent[0]);
    int xN_QNoInsertions = _extent[0];
    int yN_QInsertions = _extent[1]-yN_QNoInsertions;
    int xN_QInsertions = _extent[0];
    Cuboid2D<T> firstChildQ({_origin[0], _origin[1]}, _delta, {xN_QNoInsertions, yN_QNoInsertions});
    Cuboid2D<T> secondChildQ({_origin[0], _origin[1]+yN_QNoInsertions*_delta}, _delta, {xN_QInsertions, yN_QInsertions});
    firstChildQ.divide({nX, nY-rest}, childrenC);
    secondChildQ.divide({nX+1,rest}, childrenC);
  }
  else {
    int n_QNoInsertions = nY*(nX-rest);
    T bestVolume_QNoInsertions = (T)_extent[0]*_extent[1] * n_QNoInsertions/(T)p;
    int xN_QNoInsertions = (int)(bestVolume_QNoInsertions / (T)_extent[1] + 0.9999);
    int yN_QNoInsertions = _extent[1];
    int xN_QInsertions = _extent[0]-xN_QNoInsertions;
    int yN_QInsertions = _extent[1];
    Cuboid2D<T> firstChildQ({_origin[0], _origin[1]}, _delta, {xN_QNoInsertions, yN_QNoInsertions});
    Cuboid2D<T> secondChildQ({_origin[0]+xN_QNoInsertions*_delta, _origin[1]}, _delta, {xN_QInsertions, yN_QInsertions});
    firstChildQ.divide({nX-rest, nY}, childrenC);
    secondChildQ.divide({rest,nY+1}, childrenC);
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::divideP(int p, std::vector<Cuboid<T,D>>& childrenC) const requires (D == 3)
{
  int iXX = 1;
  int iYY = 1;
  int iZZ = p;
  int nX = _extent[0]/iXX;
  int bestIx = iXX;
  int nY = _extent[1]/iYY;
  int bestIy = iYY;
  int nZ = _extent[2]/iZZ;
  int bestIz = iZZ;
  T bestRatio = ((T)(_extent[0]/iXX)/(T)(_extent[1]/iYY)-1)*((T)(_extent[0]/iXX)/(T)(_extent[1]/iYY)-1)
                + ((T)(_extent[1]/iYY)/(T)(_extent[2]/iZZ)-1)*((T)(_extent[1]/iYY)/(T)(_extent[2]/iZZ)-1)
                + ((T)(_extent[2]/iZZ)/(T)(_extent[0]/iXX)-1)*((T)(_extent[2]/iZZ)/(T)(_extent[0]/iXX)-1);

  for (int iX=1; iX<=p; iX++) {
    for (int iY=1; iY*iX<=p; iY++) {
      for (int iZ=p/(iX*iY); iZ*iY*iX<=p; iZ++) {
        if ((iX+1)*iY*iZ>p && iX*(iY+1)*iZ>p ) {
          T ratio = ((T)(_extent[0]/iX)/(T)(_extent[1]/iY)-1)*((T)(_extent[0]/iX)/(T)(_extent[1]/iY)-1)
                    + ((T)(_extent[1]/iY)/(T)(_extent[2]/iZ)-1)*((T)(_extent[1]/iY)/(T)(_extent[2]/iZ)-1)
                    + ((T)(_extent[2]/iZ)/(T)(_extent[0]/iX)-1)*((T)(_extent[2]/iZ)/(T)(_extent[0]/iX)-1);
          if (ratio<bestRatio) {
            bestRatio = ratio;
            bestIx = iX;
            bestIy = iY;
            bestIz = iZ;
            nX = _extent[0]/iX;
            nY = _extent[1]/iY;
            nZ = _extent[2]/iZ;
          }
        }
      }
    }
  }

  int rest = p - bestIx*bestIy*bestIz;

  // split in one cuboid
  if (rest==0) {
    divide({bestIx, bestIy, bestIz}, childrenC);
    return;
  }
  else {

    // add in z than in y direction
    if (nZ>nY && nZ>nX) {

      int restY = rest%bestIy;
      // split in two cuboid
      if (restY==0) {
        int restX = rest/bestIy;
        CuboidDecomposition<T,2> helpG({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restX);

        int yN_child = 0;
        T globPosY_child = _origin[1];

        for (int iY=0; iY<bestIy; iY++) {
          yN_child         = (_extent[1]+bestIy-iY-1)/bestIy;
          for (int iC=0; iC<helpG.size(); iC++) {
            int xN_child     = helpG.get(iC).getNx();
            int zN_child     = helpG.get(iC).getNy();
            T globPosX_child = helpG.get(iC).getOrigin()[0];
            T globPosZ_child = helpG.get(iC).getOrigin()[1];

            Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                              _delta, {xN_child, yN_child, zN_child});
            childrenC.push_back(child);

          }
          globPosY_child += yN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restX = rest/bestIy+1;
      int yN_child = 0;
      T globPosY_child = _origin[1];
      int splited_nY = (int) (_extent[1] * (T)((bestIx*bestIz+restX)*restY)/(T)p);
      CuboidDecomposition<T,2> helpG0({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restX);

      for (int iY=0; iY<restY; iY++) {
        yN_child         = (splited_nY+restY-iY-1)/restY;
        for (int iC=0; iC<helpG0.size(); iC++) {
          int xN_child     = helpG0.get(iC).getNx();
          int zN_child     = helpG0.get(iC).getNy();
          T globPosX_child = helpG0.get(iC).getOrigin()[0];
          T globPosZ_child = helpG0.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }

      splited_nY = _extent[1] - splited_nY;
      restX = rest/bestIy;
      CuboidDecomposition<T,2> helpG1({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restX);
      yN_child = 0;

      for (int iY=0; iY<bestIy-restY; iY++) {
        yN_child         = (splited_nY+bestIy-restY-iY-1)/(bestIy-restY);
        for (int iC=0; iC<helpG1.size(); iC++) {
          int xN_child     = helpG1.get(iC).getNx();
          int zN_child     = helpG1.get(iC).getNy();
          T globPosX_child = helpG1.get(iC).getOrigin()[0];
          T globPosZ_child = helpG1.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }
      return;
    }

    // add in x than in y direction
    else if (nX>nY && nX>nZ) {
      int restY = rest%bestIy;
      // split in two cuboid
      if (restY==0) {
        int restZ = rest/bestIy;
        CuboidDecomposition<T,2> helpG({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restZ);

        int yN_child = 0;
        T globPosY_child = _origin[1];

        for (int iY=0; iY<bestIy; iY++) {
          yN_child         = (_extent[1]+bestIy-iY-1)/bestIy;
          for (int iC=0; iC<helpG.size(); iC++) {
            int xN_child     = helpG.get(iC).getNx();
            int zN_child     = helpG.get(iC).getNy();
            T globPosX_child = helpG.get(iC).getOrigin()[0];
            T globPosZ_child = helpG.get(iC).getOrigin()[1];

            Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                              _delta, {xN_child, yN_child, zN_child});
            childrenC.push_back(child);

          }
          globPosY_child += yN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest/bestIy+1;

      int yN_child = 0;
      T globPosY_child = _origin[1];
      int splited_nY = (int) (_extent[1] * (T)((bestIx*bestIz+restZ)*restY)/(T)p);
      CuboidDecomposition<T,2> helpG0({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restZ);

      for (int iY=0; iY<restY; iY++) {
        yN_child         = (splited_nY+restY-iY-1)/restY;
        for (int iC=0; iC<helpG0.size(); iC++) {
          int xN_child     = helpG0.get(iC).getNx();
          int zN_child     = helpG0.get(iC).getNy();
          T globPosX_child = helpG0.get(iC).getOrigin()[0];
          T globPosZ_child = helpG0.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }

      splited_nY = _extent[1] - splited_nY;
      restZ = rest/bestIy;

      CuboidDecomposition<T,2> helpG1({_origin[0], _origin[2]}, _delta, {_extent[0], _extent[2]}, bestIx*bestIz+restZ);
      yN_child = 0;

      for (int iY=0; iY<bestIy-restY; iY++) {
        yN_child         = (splited_nY+bestIy-restY-iY-1)/(bestIy-restY);
        for (int iC=0; iC<helpG1.size(); iC++) {
          int xN_child     = helpG1.get(iC).getNx();
          int zN_child     = helpG1.get(iC).getNy();
          T globPosX_child = helpG1.get(iC).getOrigin()[0];
          T globPosZ_child = helpG1.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosY_child += yN_child*_delta;
      }
      return;
    }

    // add in y than in x direction
    else {
      int restX = rest%bestIx;
      // split in two cuboid
      if (restX==0) {
        int restZ = rest/bestIx;
        CuboidDecomposition<T,2> helpG({_origin[2], _origin[1]}, _delta, {_extent[2], _extent[1]}, bestIz*bestIy+restZ);


        int xN_child = 0;
        T globPosX_child = _origin[0];

        for (int iX=0; iX<bestIx; iX++) {
          xN_child         = (_extent[0]+bestIx-iX-1)/bestIx;
          for (int iC=0; iC<helpG.size(); iC++) {
            int zN_child     = helpG.get(iC).getNx();
            int yN_child     = helpG.get(iC).getNy();
            T globPosZ_child = helpG.get(iC).getOrigin()[0];
            T globPosY_child = helpG.get(iC).getOrigin()[1];

            Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                              _delta, {xN_child, yN_child, zN_child});
            childrenC.push_back(child);
          }
          globPosX_child += xN_child*_delta;
        }
        return;
      }

      // split in four cuboid

      int restZ = rest/bestIx+1;
      int xN_child = 0;
      T globPosX_child = _origin[0];
      int splited_nX = (int) (_extent[0] * (T)((bestIz*bestIy+restZ)*restX)/(T)p);
      CuboidDecomposition<T,2> helpG0({_origin[2], _origin[1]}, _delta, {_extent[2], _extent[1]}, bestIz*bestIy+restZ);

      for (int iX=0; iX<restX; iX++) {
        xN_child         = (splited_nX+restX-iX-1)/restX;
        for (int iC=0; iC<helpG0.size(); iC++) {
          int zN_child     = helpG0.get(iC).getNx();
          int yN_child     = helpG0.get(iC).getNy();
          T globPosZ_child = helpG0.get(iC).getOrigin()[0];
          T globPosY_child = helpG0.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosX_child += xN_child*_delta;
      }

      splited_nX = _extent[0] - splited_nX;
      restZ = rest/bestIx;
      CuboidDecomposition<T,2> helpG1({_origin[2], _origin[1]}, _delta, {_extent[2], _extent[1]}, bestIz*bestIy+restZ);
      xN_child = 0;

      for (int iX=0; iX<bestIx-restX; iX++) {
        xN_child         = (splited_nX+bestIx-restX-iX-1)/(bestIx-restX);
        for (int iC=0; iC<helpG1.size(); iC++) {
          int zN_child     = helpG1.get(iC).getNx();
          int yN_child     = helpG1.get(iC).getNy();
          T globPosZ_child = helpG1.get(iC).getOrigin()[0];
          T globPosY_child = helpG1.get(iC).getOrigin()[1];

          Cuboid3D<T> child({globPosX_child, globPosY_child, globPosZ_child},
                            _delta, {xN_child, yN_child, zN_child});
          childrenC.push_back(child);
        }
        globPosX_child += xN_child*_delta;
      }
      return;
    }
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::resize(Vector<int,D> offset, Vector<int,D> extent) {
  _origin += offset*_delta;
  _extent = extent;
}

template <typename T, unsigned D>
void Cuboid<T,D>::divideFractional(int iD, std::vector<T> fractions, std::vector<Cuboid<T,D>>& childrenC) const
{
  auto delta = Vector<T,D>([&](int i) -> T {
    return i == iD ? _delta : 0;
  });
  auto base = Vector<int,D>([&](int i) -> T {
    return i == iD ? 1 : 0;
  });

  std::vector<int> fractionWidths;
  int totalWidth = 0;
  for (T f : fractions) {
    fractionWidths.emplace_back(f * (getExtent()*base));
    totalWidth += fractionWidths.back();
  }
  fractionWidths.back() += getExtent()*base - totalWidth;

  auto origin = getOrigin();
  auto extent = Vector<int,D>([&](int i) -> T {
    return i == iD ? 0 : getExtent()[i];
  });

  for (int width : fractionWidths) {
    Cuboid<T,D> child(origin, _delta, extent + (width)*base);
    child.print();
    origin += width*delta;
    childrenC.push_back(child);
  }
}

template <typename T, unsigned D>
void Cuboid<T,D>::writeAsXML(std::ostream& ss) const {
  ss << " extent=\"";
  for (int i = 0; i<3; i++) {
    ss << getExtent()[i] << " ";
  }

  ss << "\" origin=\"";
  for (int i = 0; i<3; i++) {
    ss << getOrigin()[i] << " ";
  }

  ss << "\" deltaR=\"" << getDeltaR();
  ss << "\" weight=\"" << getWeight();
}

}

#endif
