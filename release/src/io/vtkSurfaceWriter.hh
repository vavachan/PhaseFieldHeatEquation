/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef IO_VTK_SURFACE_WRITER_HH
#define IO_VTK_SURFACE_WRITER_HH

#ifdef FEATURE_VTK

#include "vtkSurfaceWriter.h"

#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

namespace olb {

template<typename T>
void vtkSurfaceWriter<T>::init()
{
  _grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  _points = vtkSmartPointer<vtkPoints>::New();
  _cells = vtkSmartPointer<vtkCellArray>::New();
  _localTriangles.clear();

  std::vector<Cuboid3D<T>*> localCuboids;
  for (int globC=0; globC < _cuboidDecomposition.cuboids().size(); ++globC) {
    if (_loadBalancer.isLocal(globC)) {
      localCuboids.emplace_back(&_cuboidDecomposition.cuboids()[globC]);
    }
  }

  std::size_t pointId = 0;
  for (STLtriangle<T>& triangle : _surfaceI.getMesh().getTriangles()) {
    // Check whether triangle is fully contained in local cuboids (incl. overlap)
    bool isTriangleContained = true;
    for (int iD=0; iD < 3; ++iD) {
      bool isPointContained = false;
      for (auto* cuboid : localCuboids) {
        if (cuboid->isInside(triangle.point[iD], 1)) {
          isPointContained = true;
          break;
        }
      }
      isTriangleContained &= isPointContained;
    }

    if (isTriangleContained) {
      _localTriangles.emplace_back(&triangle);
      for (const auto& point : triangle.point) {
        _points->InsertNextPoint(point[0], point[1], point[2]);
      }

      vtkSmartPointer<vtkTriangle> cell = vtkSmartPointer<vtkTriangle>::New();
      cell->GetPointIds()->SetId(0, pointId++);
      cell->GetPointIds()->SetId(1, pointId++);
      cell->GetPointIds()->SetId(2, pointId++);
      _cells->InsertNextCell(cell);
    }
  }

  _grid->SetPoints(_points);
  _grid->SetCells(VTK_TRIANGLE, _cells);
}

template<typename T>
void vtkSurfaceWriter<T>::write(int iT)
{
  _grid->GetPointData()->Reset();
  std::vector<vtkSmartPointer<vtkFloatArray>> data(_f.size());
  for (std::size_t iF=0; iF < _f.size(); ++iF) {
    data[iF] = vtkSmartPointer<vtkFloatArray>::New();
    data[iF]->SetName(_f[iF]->getName().c_str());
    data[iF]->SetNumberOfComponents(_f[iF]->getTargetDim());
    _grid->GetPointData()->AddArray(data[iF]);
  }

  std::string filePath = singleton::directories().getVtkOutDir()
                       + createFileName(_fileName, iT)
                       + ".pvtu";

  if (singleton::mpi().isMainProcessor()) {
    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> multiWriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
    multiWriter->SetInputData(_grid);
    multiWriter->SetNumberOfPieces(singleton::mpi().getSize());
    multiWriter->SetFileName(filePath.c_str());
    multiWriter->SetUseSubdirectory(true);
    multiWriter->SetStartPiece(0);
    multiWriter->SetEndPiece(singleton::mpi().getSize()-1);
    multiWriter->SetWriteSummaryFile(1);
    multiWriter->Write();
  }

  for (const STLtriangle<T>* triangle : _localTriangles) {
    for (const auto& point : triangle->point) {
      for (std::size_t iF=0; iF < _f.size(); ++iF) {
        std::vector<T> result(_f[iF]->getTargetDim(), 0);
        _f[iF](result.data(), point.data());
        data[iF]->InsertNextTuple(result.data());
      }
    }
  }

  vtkSmartPointer<vtkXMLPUnstructuredGridWriter> multiWriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  multiWriter->SetWriteSummaryFile(0);
  multiWriter->SetInputData(_grid);
  multiWriter->SetNumberOfPieces(singleton::mpi().getSize());
  multiWriter->SetFileName(filePath.c_str());
  multiWriter->SetUseSubdirectory(true);
  multiWriter->SetStartPiece(singleton::mpi().getRank());
  multiWriter->SetEndPiece(singleton::mpi().getRank());
  multiWriter->Write();

  _f.clear();
}

}

#endif

#endif
