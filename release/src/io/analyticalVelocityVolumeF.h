/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerländer, Dennis Teutscher, Shota Ito
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

#ifndef ANALYTICAL_VELOCITY_VOLUME_F_H
#define ANALYTICAL_VELOCITY_VOLUME_F_H

#if defined(FEATURE_VTK)
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkCellLocator.h>
#endif

namespace olb {
template <typename T>
class AnalyticalVelocityVolumeF final : public AnalyticalF<3,T,T> {
private:
#if defined(FEATURE_VTK)
  vtkSmartPointer<vtkXMLImageDataReader> _reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  vtkSmartPointer<vtkImageData> _data;
  vtkSmartPointer<vtkCellLocator> _cellLocator = vtkSmartPointer<vtkCellLocator>::New();
#endif

  std::string _fileName;
  Vector<double,3> _origin{};
  Vector<int,3> _shape{};
  T _spacing = 0.0;

  bool CellData = false;
  bool PointData = false;

public:
  AnalyticalVelocityVolumeF(std::string fileName)
      : AnalyticalF<3,T,T>(3)
  {
    this->getName() = "VelocityVolumeImporter";

    OstreamManager clout(std::cout, "VelocityVolumeImporter");
    _fileName = fileName;

    #if not defined(FEATURE_VTK)
    if (isVTIFile(_fileName))
    {
      std::cerr << "To use the VTK format, add VTK to the FEATURES list in config.mk";
      exit(1);
    }
    #endif
    #if defined(FEATURE_VTK)
    if (!isVTIFile(_fileName)) {
      std::cerr << "Please convert your file to .vti to use the functor.";
      exit(1);
    } else {
      clout << "Read .vti file..." << std::endl;
      _reader->SetFileName(fileName.c_str());
      _reader->Update();
      auto* _data = _reader->GetOutput();

      _spacing = static_cast<T>(_data->GetSpacing()[0]);
      clout << "Read spacing : " << _spacing << std::endl;
      _data->GetDimensions(_shape.data());
      clout << "Read shape : " << _shape << std::endl;
      _data->GetOrigin(_origin.data());
      clout << "Read origin : " << _origin << std::endl;

      vtkCellData* _cellData = _data->GetCellData();
      vtkPointData* _pointData = _data->GetPointData();
      if (_cellData->GetNumberOfArrays() > 0) {
        clout << "Found CellData. " << std::endl;
        CellData = true;
        _cellLocator->SetDataSet(_data);
        _cellLocator->BuildLocator();
      } else if (_pointData->GetNumberOfArrays() > 0) {
        clout << "Found PointData." << std::endl;
        clout << "Warning: If functor reads only 0, consider converting PointData to CellData." << std::endl;
        PointData = true;
      } else {
        std::cerr << "Neither cell data nor point data found.";
        exit(1);
      }
    }
    #endif
  }

  Vector<T,3> getPhysShape() const {
    return _shape * _spacing;
  }

  T getSpacing() const {
    return _spacing;
  }

  bool operator()(T output[], const T physR[]) override{
    #if defined(FEATURE_VTK)
    if (!isInside(physR)) {
      //std::cout << "Point outside." << std::endl;
      output[0] = 0;
      output[1] = 0;
      output[2] = 0;
      return true;
    }

    auto* _data = _reader->GetOutput();

    // cell data (preferred)
    if (CellData) {
      vtkCellData* _cellData = _data->GetCellData();
      // Assuming cellData holds the velocity data in the first array
      vtkDataArray* _vectors = _cellData->GetVectors(_cellData->GetArrayName(0));
      if (_vectors == nullptr) {
        std::cerr << "[AnalyticalVelocityVolumeF]: ";
        std::cerr << "Cell data: No vector data was found." << std::endl;
        exit(1);
      }

      double x[3] = {physR[0], physR[1], physR[2]};
      vtkIdType cellId = _cellLocator->FindCell(x);
      if (cellId >= 0) {
        double velocity[3] = {};
        _vectors->GetTuple(cellId, velocity);
        for (int dim = 0; dim < 3; ++dim) {
          output[dim] = static_cast<T>(velocity[dim]);
        }
      } else {
        //std::cout << "Cell not found." << std::endl;
        output[0] = 0;
        output[1] = 0;
        output[2] = 0;
      }
    }

    // point data (can cause problems due to round-off)
    else if (PointData) {
      vtkPointData* _pointData = _data->GetPointData();
      // Assuming pointData holds the velocity data in the first array
      vtkDataArray* _vectors = _pointData->GetVectors(_pointData->GetArrayName(0));
      if (_vectors == nullptr) {
        std::cerr << "[AnalyticalVelocityVolumeF]: ";
        std::cerr << "Point data: No vector data was found." << std::endl;
        return false;
      }

      Vector<int,3> x = getVtiLatticeR(physR);
      vtkIdType pointId = _data->FindPoint(x[0], x[1], x[2]);
      if (pointId >= 0) {
        double velocity[3] = {};
        _vectors->GetTuple(pointId, velocity);
        for (int dim = 0; dim < 3; ++dim) {
          output[dim] = static_cast<T>(velocity[dim]);
        }
      } else {
        //std::cout << "Point not found." << std::endl;
        output[0] = 0;
        output[1] = 0;
        output[2] = 0;
      }
    }
    #endif
    return true;
  }

private:
  bool isVTIFile(const std::string& fileName) {
    return (fileName.size() >= 4 && fileName.substr(fileName.size() - 4) == ".vti");
  }

  Vector<int,3> getVtiLatticeR(const T physR[]) {
    int iX;
    int iY;
    int iZ;
    iX = util::floor(physR[0] / _spacing);
    iY = util::floor(physR[1] / _spacing);
    iZ = util::floor(physR[2] / _spacing);
    return {iX, iY, iZ};
  }

  bool isInside(const T physR[]) {
    Vector<int,3> vtiLatticeR = getVtiLatticeR(physR);
    int iX = vtiLatticeR[0];
    int iY = vtiLatticeR[1];
    int iZ = vtiLatticeR[2];
    return (iX >= 0 && iY >= 0 && iZ >= 0 && iX < _shape[0] && iY < _shape[1] && iZ < _shape[2]);
  }
};
}

#endif
