/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef GPU_CUDA_DYNAMICS_HH
#define GPU_CUDA_DYNAMICS_HH

#include "dynamics.h"
#include "dynamics/lbm.h"

namespace olb {

namespace gpu {

namespace cuda {

template <typename T, typename DESCRIPTOR> class DeviceContext;
template <typename T, typename DESCRIPTOR> class DataOnlyCell;

/// Virtual interface for device-side dynamically-dispatched dynamics access
template <typename T, typename DESCRIPTOR>
struct Dynamics {
  virtual CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) __device__ = 0;

  virtual T    computeRho (DataOnlyCell<T,DESCRIPTOR>& cell              ) __device__ = 0;
  virtual void computeU   (DataOnlyCell<T,DESCRIPTOR>& cell,         T* u) __device__ = 0;
  virtual void computeJ   (DataOnlyCell<T,DESCRIPTOR>& cell,         T* j) __device__ = 0;
  virtual void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) __device__ = 0;

  virtual void computeStress    (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) __device__ = 0;
  virtual void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) __device__ = 0;

  virtual void defineRho       (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho             ) __device__ = 0;
  virtual void defineU         (DataOnlyCell<T,DESCRIPTOR>& cell,         T* u       ) __device__ = 0;
  virtual void defineRhoU      (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u       ) __device__ = 0;
  virtual void defineAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) __device__ = 0;

  virtual void computeEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* fEq) __device__ = 0;

  virtual T getOmegaOrFallback(T fallback) __device__ = 0;

  void iniEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T rho, T* u) __device__ {
    T fEq[DESCRIPTOR::q] { };
    computeEquilibrium(cell, rho, u, fEq);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop];
    }
  }

  void iniRegularized(DataOnlyCell<T,DESCRIPTOR>& cell,
                      T rho, T* u, T* pi) __device__ {
    iniEquilibrium(cell, rho, u);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
    }
  }

  virtual void inverseShiftRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u)  __device__ = 0;

};

/// Implementation of gpu::cuda::Dynamics for concrete DYNAMICS
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public Dynamics<T,DESCRIPTOR> {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* _parameters;

public:
  ConcreteDynamics(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* parameters) __device__:
    _parameters{parameters} {
  }

  CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) override __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    return DYNAMICS().collide(cell, *_parameters);
  }

  T computeRho(DataOnlyCell<T,DESCRIPTOR>& cell) override __device__ {
    return DYNAMICS::MomentaF().computeRho(cell);
  }
  void computeU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) override __device__ {
    DYNAMICS::MomentaF().computeU(cell, u);
  }
  void computeJ(DataOnlyCell<T,DESCRIPTOR>& cell, T* j) override __device__ {
    DYNAMICS::MomentaF().computeJ(cell, j);
  }
  void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override __device__ {
    DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  }
  void computeStress(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override __device__ {
    DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  }
  void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override __device__ {
    DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
  }

  T getOmegaOrFallback(T fallback) override __device__ {
    if constexpr (DYNAMICS::parameters::template contains<descriptors::OMEGA>()) {
      return _parameters->template get<descriptors::OMEGA>();
    } else {
      return fallback;
    }
  }

  void computeEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* fEq) override __device__ {
    DYNAMICS::EquilibriumF().compute(cell, rho, u, fEq);
  }

  void defineRho(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho) override __device__ {
    typename DYNAMICS::MomentaF().defineRho(cell, rho);
  }

  void defineU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) override __device__ {
    typename DYNAMICS::MomentaF().defineU(cell, u);
  }

  void defineRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override __device__ {
    typename DYNAMICS::MomentaF().defineRhoU(cell, rho, u);
  }

  void defineAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override __device__ {
    typename DYNAMICS::MomentaF().defineAllMomenta(cell, rho, u, pi);
  }

  void inverseShiftRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override __device__ {
    typename DYNAMICS::MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};

/// Last node in a MaskedDynamics chain in kernel::call_operators
struct DynamicDispatchCollision {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
      collisionO->collide(lattice, iCell);
      return true;
    }
    return false;
  }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
      statistic = collisionO->collide(lattice, iCell);
      return true;
    }
    return false;
  }

};

}

}


}

#endif
