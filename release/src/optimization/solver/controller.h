/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
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

/** \file
 * The description of a Controller -- header file.
 * A controller manages an array of control variables.
 */


#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "io/ostreamManager.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

/// All optimization code is contained in this namespace.
namespace opti {


template<typename S>
class Controller {

protected:
  int _dimControl;
  std::vector<S> _control;

public:
  Controller (int dimControl)
   : _dimControl (dimControl),
     _control (dimControl, 0.0)
  { }

  Controller (std::vector<S>& initialValues)
   : _dimControl(initialValues.size()), _control(initialValues)
  { }

  S getControl(int i) const
  {
    return _control[i];
  };
  const std::vector<S>& getControl() const
  {
    return _control;
  };

  const S* getControlPointer(int i) const
  {
    return &_control[i];
  }
  S* getControlPointer(int i)
  {
    return &_control[i];
  }

  int getDimControl() const
  {
    return _dimControl;
  };

  void setControl(const std::vector<S>& control, int dimControl)
  {
    _control = control;
    _dimControl = dimControl;
  };

  void print()
  {
    static OstreamManager clout(std::cout,"Controller");

    clout << " control={ ";
    for (int i=0; i<_dimControl; i++) {
      clout << _control[i] << " ";
    }
    clout << "}" << std::endl;

  };

  void writeToFile(const char* fname=nullptr)
  {
    static OstreamManager clout(std::cout,"Controller");
    std::ofstream file;
    if (fname==nullptr) {
      static int iter(0);
      std::stringstream number;
      number << std::setw(3) << std::setfill('0') << iter++;
      file.open(("Controller" + number.str() + ".txt").c_str());
    }
    else {
      file.open(fname);
    }

    file << _dimControl << std::endl;
    for (int i=0; i<_dimControl; i++) {
      file << _control[i] << std::endl;
    }
    file.close();
  };

  void readFromFile(const char* fname)
  {
    static OstreamManager clout(std::cout,"Controller");
    std::ifstream file;
    //    if (file.open(fname)) {
    file.open(fname);
    int dim;
    file >> dim;
    if (dim!=_dimControl) {
      clout << "Error: dimensions do not match! dim_controller=" << _dimControl << "; dim_file=" << dim << std::endl;
      assert(false);
    }
    else {
      for (int i=0; i<_dimControl; i++) {
        double tmpVal;
        file >> tmpVal;
        _control[i]  = tmpVal;
      }
    }
    file.close();
    //    }
  }
};


} // namespace opti

} // namespace olb

#endif
