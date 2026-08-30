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

#ifndef GPU_CUDA_CONTEXT_H
#define GPU_CUDA_CONTEXT_H

#include "device.h"

namespace olb {

/// Representation of (Dynamics,Operator)Parameters<DYNAMICS> for CUDA block lattice
template <typename T, typename DESCRIPTOR, typename PARAMETERS>
class ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,PARAMETERS> final
  : public AbstractedConcreteParameters<T,DESCRIPTOR>
  , public Serializable {
private:
  using ParametersD = typename olb::ParametersD<T,DESCRIPTOR>::template include<PARAMETERS>;
  gpu::cuda::device::unique_ptr<ParametersD> _deviceParameters;

public:
  ParametersD parameters;

  ConcreteParametersD(std::size_t); // TODO: Implement more generic non-cellwise field allocation in Data

  AbstractParameters<T,DESCRIPTOR>& asAbstract() override {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) override;

  ParametersD* deviceData() {
    return _deviceParameters.get();
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

}

#endif
