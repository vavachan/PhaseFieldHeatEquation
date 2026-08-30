/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
 *
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

#ifndef SUPER_D_H
#define SUPER_D_H

#include "blockD.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class SuperD {
private:
  std::vector<std::unique_ptr<BlockD<T,DESCRIPTOR>>> _block;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  SuperD(const std::vector<std::pair<Vector<std::size_t,DESCRIPTOR::d>,Platform>>& blocks) {
    for (auto [extent, platform] : blocks) {
      _block.emplace_back(constructUsingConcretePlatform<ConcretizableBlockD<T,DESCRIPTOR>>(
        platform, extent, 0));
    }
  }

  SuperD(const LoadBalancer<T>& load) {
    for (int iB = 0; iB < load.size(); ++iB) {
      _block.emplace_back(constructUsingConcretePlatform<ConcretizableBlockD<T,DESCRIPTOR>>(
        load.platform(iB), 0, 0));
    }
  }

  std::size_t size() const {
    return _block.size();
  }

  BlockD<T,DESCRIPTOR>& getBlock(unsigned iB) {
    return *_block.at(iB);
  }

  template <typename OPERATOR>
  void apply() {
    OPERATOR{}.apply(*this);
  }

  void setProcessingContext(ProcessingContext context) {
    for (auto& block : _block) {
      block->setProcessingContext(context);
    }
  }

};

}

#endif
