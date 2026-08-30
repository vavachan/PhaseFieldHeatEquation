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

#ifndef IO_VTK_SURFACE_WRITER_H
#define IO_VTK_SURFACE_WRITER_H

#ifdef FEATURE_VTK

#include "stlReader.h"
#include "fileName.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

namespace olb {

template<typename T>
class vtkSurfaceWriter {
private:
  STLreader<T>& _surfaceI;
  CuboidDecomposition3D<T>& _cuboidDecomposition;
  LoadBalancer<T>& _loadBalancer;

  const std::string _fileName;
  std::vector<FunctorPtr<AnalyticalF3D<T,T>>> _f;

  vtkSmartPointer<vtkUnstructuredGrid> _grid;
  vtkSmartPointer<vtkPoints> _points;
  vtkSmartPointer<vtkCellArray> _cells;
  std::vector<const STLtriangle<T>*> _localTriangles;

  void init();

public:
  vtkSurfaceWriter(STLreader<T>& surfaceI,
                   CuboidDecomposition3D<T>& cuboidDecomposition,
                   LoadBalancer<T>& loadBalancer,
                   const std::string& fileName):
    _surfaceI(surfaceI),
    _cuboidDecomposition(cuboidDecomposition),
    _loadBalancer(loadBalancer),
    _fileName(fileName)
  {
    init();
  }

  void addFunctor(FunctorPtr<AnalyticalF3D<T,T>>&& f) {
    _f.emplace_back(std::move(f));
  }

  void write(int iT);

};

}

#endif

#endif
