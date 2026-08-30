/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerlaender
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

#ifndef SUPER_FIELD_ARRAY_D_H
#define SUPER_FIELD_ARRAY_D_H

namespace olb {

/// Maintains per-block arrays of FIELD
/**
 * Useful for storing additional data in the correct per-block platforms
 * for usage in operators (by passing pointer as parameter).
 **/
template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperFieldArrayD final : public SuperStructure<T,DESCRIPTOR::d> {
private:
  std::vector<std::unique_ptr<ColumnVectorBase>> _block;

public:
  SuperFieldArrayD(CuboidDecomposition<T,DESCRIPTOR::d>& cGeometry,
                   LoadBalancer<T>& loadBalancer):
    SuperStructure<T,DESCRIPTOR::d>(cGeometry, loadBalancer, 0)
  {
    auto& load = this->getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      _block.emplace_back(
        constructUsingConcretePlatform<ConcretizableFieldArrayD<T,DESCRIPTOR,FIELD>>(
          load.platform(iC), 0));
    }
  }

  AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& getBlock(int iC) {
    auto& load = this->getLoadBalancer();
    return *callUsingConcretePlatform<ConcretizableFieldArrayD<T,DESCRIPTOR,FIELD>>(
      load.platform(iC),
      _block[iC].get(),
      [&](auto* field) -> auto* {
        return &(field->asAbstract());
      });
  }

  template <Platform PLATFORM>
  FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& getBlock(int iC) {
    if (auto* ptr = dynamic_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>*>(_block[iC].get());
        ptr != nullptr) {
      return *ptr;
    } else {
      throw std::invalid_argument("FieldArrayD is not of PLATFORM");
    }
  }

};

}

#endif
