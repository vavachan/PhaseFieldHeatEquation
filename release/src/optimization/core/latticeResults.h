/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef LATTICE_RESULTS_H
#define LATTICE_RESULTS_H

namespace olb {

/// Wrapper of SuperVTMwriter allowing to specify which functors are
/// dumped as .vtm, without having the actual SuperLattice instance.
/// This is useful when handling multiple simulations as in adjoint flow control.
template <typename T, typename DESCRIPTOR>
class LatticeResults {
private:
  std::string _name;
  std::vector<std::function<
    std::unique_ptr<SuperF<DESCRIPTOR::d,T,T>>(const UnitConverter<T,DESCRIPTOR>&, SuperLattice<T,DESCRIPTOR>&)>
  > _creators;

  void prepare() {
    SuperVTMwriter<T,DESCRIPTOR::d> writer(_name);
    writer.createMasterFile();
  }

public:
  LatticeResults(std::string name):
    _name(name)
  {
    this->prepare();
  };

  /// Used to specify list of functors which should be written
  template <typename... FUNCTOR>
  void add() {
    (_creators.push_back([](const UnitConverter<T,DESCRIPTOR>& converter, SuperLattice<T,DESCRIPTOR>& lattice)
      -> std::unique_ptr<SuperF<DESCRIPTOR::d,T,T>> {
      if constexpr (std::is_base_of_v<SuperLatticePhysF<T,DESCRIPTOR>, FUNCTOR>) {
        return std::make_unique<FUNCTOR>(lattice, converter);
      } else {
        return std::make_unique<FUNCTOR>(lattice);
      }
    }), ...);
  }

  template <typename... FIELD>
  void addFields() {
    (_creators.push_back([](const UnitConverter<T,DESCRIPTOR>& converter, SuperLattice<T,DESCRIPTOR>& lattice)
      -> std::unique_ptr<SuperF<DESCRIPTOR::d,T,T>> {
      return std::make_unique<SuperLatticeField<T,DESCRIPTOR,FIELD>>(lattice);
    }), ...);
  }

  /// Evaluate functors added for IO
  void write(LatticeData<T,DESCRIPTOR>& data, std::size_t arg) {
    SuperVTMwriter<T,DESCRIPTOR::d> writer(_name);
    data.getSuperLattice().setProcessingContext(ProcessingContext::Evaluation);
    std::vector<std::unique_ptr<SuperF<DESCRIPTOR::d,T,T>>> functors;
    for (auto& f : _creators) {
      functors.push_back(f(data.getUnitConverter(), data.getSuperLattice()));
      writer.addFunctor(*functors.back());
    }
    writer.write(arg);
  }

  /// Evaluate functors added for IO
  void write(const UnitConverter<T,DESCRIPTOR>& converter, SuperLattice<T,DESCRIPTOR>& lattice, std::size_t arg) {
    SuperVTMwriter<T,DESCRIPTOR::d> writer(_name);
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    std::vector<std::unique_ptr<SuperF<DESCRIPTOR::d,T,T>>> functors;
    for (auto& f : _creators) {
      functors.push_back(f(converter, lattice));
      writer.addFunctor(*functors.back());
    }
    writer.write(arg);
  }
};

}

#endif
