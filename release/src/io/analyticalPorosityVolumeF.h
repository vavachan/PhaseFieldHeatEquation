/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerländer, Dennis Teutscher
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

#ifndef ANALYTICAL_POROSITY_VOLUME_F_H
#define ANALYTICAL_POROSITY_VOLUME_F_H

#ifdef FEATURE_VDB
#include <openvdb/openvdb.h>
#include <openvdb/tree/ValueAccessor.h>
#include <openvdb/Grid.h>
#endif

namespace olb {

template <typename T>
class AnalyticalPorosityVolumeF final : public AnalyticalF<3,T,T> {
private:
#ifdef FEATURE_VDB
  openvdb::GridPtrVecPtr _grids;
  openvdb::FloatGrid::Ptr _grid;
  openvdb::math::Transform::Ptr _rotation;
#endif

  Vector<int,3> _shape;
  Vector<T, 3> _origin;
  T _volumeDeltaX;

public:
  /// Reads first grid of VDB file with spacing volumeDeltaX and applies rotation [rad]
  AnalyticalPorosityVolumeF(std::string fileName, T volumeDeltaX, T rotation=0)
      : AnalyticalF<3,T,T>(1), _volumeDeltaX{volumeDeltaX}
  {
    #ifdef FEATURE_VDB
    openvdb::io::File _vdbFile(fileName);
    // Read the OpenVDB file
    _vdbFile.open();
    _grids = _vdbFile.getGrids();

    _grid = openvdb::gridPtrCast<openvdb::FloatGrid>(*_grids->begin());
    auto bbox = _grid->evalActiveVoxelBoundingBox();

    if (rotation != 0) {
      openvdb::Vec3d minIdx(bbox.min().x(), bbox.min().y(), bbox.min().z());
      openvdb::Vec3d extents(bbox.extents().x(), bbox.extents().y(), bbox.extents().z());
      openvdb::Vec3d centerIdx = minIdx + extents * 0.5;

      _rotation = openvdb::math::Transform::createLinearTransform(
        openvdb::math::Mat4d::identity());
      _rotation->preTranslate( centerIdx);
      _rotation->preRotate(rotation, openvdb::math::Z_AXIS);
      _rotation->preTranslate(-centerIdx);
      _grid->setTransform(_rotation);
    } else {
      _rotation = openvdb::math::Transform::createLinearTransform(
        openvdb::math::Mat4d::identity());
      _grid->setTransform(_rotation);
    }

    bbox = _grid->evalActiveVoxelBoundingBox();
    _shape[0] = bbox.max().x() - bbox.min().x();
    _shape[1] = bbox.max().y() - bbox.min().y();
    _shape[2] = bbox.max().z() - bbox.min().z();
    _origin[0] = bbox.min().x();
    _origin[1] = bbox.min().y();
    _origin[2] = bbox.min().z();

    _vdbFile.close();
    #else
    throw std::runtime_error("VDB support not enabled, set FEATURE := VDB");
    #endif
  }

  Vector<T,3> getPhysShape() const {
    return _shape * _volumeDeltaX;
  }
  Vector<T,3> getOrigin() const{
    return _origin * _volumeDeltaX;
  }

  bool operator()(T output[], const T physR[]) override{
    #ifdef FEATURE_VDB
    openvdb::Vec3d worldPos(
      util::floor(physR[0] / _volumeDeltaX),
      util::floor(physR[1] / _volumeDeltaX),
      util::floor(physR[2] / _volumeDeltaX)
    );

    openvdb::Vec3d indexPos = _grid->transform().worldToIndex(worldPos);
    openvdb::Coord location(
      static_cast<int>(std::round(indexPos.x())),
      static_cast<int>(std::round(indexPos.y())),
      static_cast<int>(std::round(indexPos.z()))
    );

    openvdb::FloatGrid::Accessor accessor = _grid->getAccessor();
    if(accessor.isValueOn(location)){
      output[0] = _grid->tree().getValue(location);
    } else {
      output[0] = 1;
    }
    return true;
    #endif
    return false;
  }

};

}

#endif
