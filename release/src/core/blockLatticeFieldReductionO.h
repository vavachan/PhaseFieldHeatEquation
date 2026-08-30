/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Yuji (Sam) Shimojima
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

#ifndef BLOCK_LATTICE_FIELD_REDUCTION_O_H
#define BLOCK_LATTICE_FIELD_REDUCTION_O_H

namespace olb {

namespace field::reduction {

//this filed is tag for extracting the 'core' grid area(no padding area)  for the reduction, especially in GPU.
struct TAG_CORE : public descriptors::TYPED_FIELD_BASE<int, 1> {};

}

namespace reduction {
struct TAGS_U             : public descriptors::TYPED_FIELD_BASE<int, 1> {};

struct ConditionTrue {
  template <typename CELL>
  bool operator()(CELL& cell) any_platform {
    return true;
  }
};

template <typename TAG_FIELD>
  struct checkBulkTag {
    template <typename CELL>
    bool operator()(CELL& cell) any_platform
    {
      return cell.template getField<TAG_FIELD>() == (int) 1;
    }
  };

struct SumO {
  template <typename FIELDD>
  FIELDD operator()(FIELDD lhs, FIELDD rhs) any_platform
  {
    return lhs + rhs;
  }

  template <typename T>
  T reset(T lhs) any_platform
  {
    T reset {};
    return reset;
  }

#ifdef PARALLEL_MODE_MPI
  MPI_Op forMPI() { return MPI_SUM; }
#endif // PARALLEL_MODE_MPI
};

struct MaxO {
  template <typename FIELDD>
  FIELDD operator()(FIELDD lhs, FIELDD rhs) any_platform
  {
    //      return maxv(lhs, rhs);
    return lhs > rhs ? lhs : rhs;
  }

  template <typename T>
  T reset(T lhs) any_platform
  {
    return std::numeric_limits<T>::lowest();
  }

#ifdef PARALLEL_MODE_MPI
  MPI_Op forMPI() { return MPI_MAX; }
#endif // PARALLEL_MODE_MPI
};

struct MinO {
  template <typename FIELDD>
  FIELDD operator()(FIELDD lhs, FIELDD rhs) any_platform
  {
    //     return minv(lhs, rhs);
    return lhs < rhs ? lhs : rhs;
  }

  template <typename T>
  T reset(T lhs) any_platform
  {
    return std::numeric_limits<T>::max();
  }

#ifdef PARALLEL_MODE_MPI
  MPI_Op forMPI() { return MPI_MIN; }
#endif // PARALLEL_MODE_MPI
};

struct ResetO {
  template <typename FIELDD>
  FIELDD operator()(FIELDD lhs) any_platform
  {
    FIELDD reset{};
    return reset;
  }
};

}

template <typename FIELD, typename REDUCTION_OP, typename CONDITION>
struct BlockLatticeFieldReductionO {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  using parameters = meta::list<fields::array_of<FIELD>>;

  int getPriority() const { return 2; }

  template <typename BLOCK>
  struct type {
    void setup(BLOCK& blockLattice) {}
    void apply(BLOCK& blockLattice) { throw std::runtime_error("BlockLatticeFieldReductionO not implemented"); }
  };

  template <typename BLOCK>
  void setup(BLOCK& blockLattice)
  {
    type<BLOCK> {}.setup(blockLattice);
  }

  template <typename BLOCK>
  void apply(BLOCK& blockLattice)
  {
    type<BLOCK> {}.apply(blockLattice);
  }
};

template <typename FIELD, typename REDUCTION_OP, typename CONDITION>
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
requires(isPlatformCPU(PLATFORM)) struct BlockLatticeFieldReductionO<FIELD, REDUCTION_OP, CONDITION>::type<
    ConcreteBlockLattice<T, DESCRIPTOR, PLATFORM>> {

  void setup(ConcreteBlockLattice<T, DESCRIPTOR, PLATFORM>& blockLattice)
  {
    blockLattice.template getData<OperatorParameters<BlockLatticeFieldReductionO>>();
  }

  void apply(ConcreteBlockLattice<T, DESCRIPTOR, PLATFORM>& blockLattice)
  {

    OstreamManager clout(std::cout, "Apply function in BlockLatticeFieldReductionO on CPU");

    auto& parameters = blockLattice.template getData<OperatorParameters<BlockLatticeFieldReductionO>>().parameters;

    const auto&                                    blockField = blockLattice.template getField<FIELD>();
    FieldD<T, DESCRIPTOR, fields::array_of<FIELD>> elementField =
        parameters.template get<fields::array_of<FIELD>>(); //not auto. because of considering 1 dimensional field

    for (unsigned iD = 0; iD < blockField.d; iD++) {
      if (elementField[iD] == nullptr) {
        return;
      }
      else {
        elementField[iD][0] = REDUCTION_OP {}.reset(elementField[iD][0]);
      }
    }

    blockLattice.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
      const std::size_t              iCell = blockLattice.getCellId(latticeR);
      const ConstCell<T, DESCRIPTOR> cell  = blockLattice.get(iCell);
      if (cell.template getField<field::reduction::TAG_CORE>() == (int)1) {
        if (CONDITION {}(cell)) {
          FieldD<T, DESCRIPTOR, FIELD> blockCellField = blockField.getRow(iCell);
          for (unsigned iD = 0; iD < blockField.d; iD++) {
            elementField[iD][0] = REDUCTION_OP {}(elementField[iD][0], blockCellField[iD]);
          }
        }
      }
    });
  }
};

} // namespace olb

#endif //BLOCK_LATTICE_FIELD_REDUCTION_O_H
