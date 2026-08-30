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

#ifndef BLOCK_D_H
#define BLOCK_D_H

namespace olb {

template<typename T, typename DESCRIPTOR> class BlockD;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockD;

/// Curried ConcreteBlockD eemplate for use in callUsingConcretePlatform
template<typename T, typename DESCRIPTOR>
struct ConcretizableBlockD {

using value_t = T;
using descriptor_t = DESCRIPTOR;
using base_t = BlockD<T,DESCRIPTOR>;

template <Platform PLATFORM>
using type = ConcreteBlockD<T,DESCRIPTOR,PLATFORM>;

};

template<typename T, typename DESCRIPTOR>
class BlockLatticeRowView {
private:
  BlockD<T,DESCRIPTOR>& _lattice;
  std::size_t _iRow;

public:
  BlockLatticeRowView(BlockD<T,DESCRIPTOR>& lattice,
                      std::size_t iRow):
    _lattice(lattice),
    _iRow(iRow) { }

  template<typename FIELD>
  auto getField() const {
    return _lattice.template getField<FIELD>().get(_iRow);
  }

  template<typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& data) const {
    _lattice.template getField<FIELD>().set(_iRow, data);
  }

};

template <typename OPERATOR>
struct CellOperatorO;

template <typename OPERATOR>
requires (OPERATOR::scope == OperatorScope::PerCell)
struct CellOperatorO<OPERATOR> {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  template <typename BLOCK>
  void apply(BLOCK& block) {
    for (std::size_t iCell=0; iCell < block.getNcells(); ++iCell) {
      cpu::PlainCell<typename BLOCK::value_t,
                     typename BLOCK::descriptor_t,
                     BLOCK::platform>
        cell(block, iCell);
      OPERATOR{}.apply(cell);
    }
  }
};

template <typename OPERATOR>
requires (OPERATOR::scope == OperatorScope::PerCellWithParameters)
struct CellOperatorO<OPERATOR> {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  template <typename BLOCK>
  void apply(BLOCK& block) {
    auto params = block.template getData<OperatorParameters<OPERATOR>>()
                       .parameters
                       .template copyAs<typename BLOCK::value_t>();
    for (std::size_t iCell=0; iCell < block.getNcells(); ++iCell) {
      cpu::PlainCell<typename BLOCK::value_t,
                     typename BLOCK::descriptor_t,
                     BLOCK::platform>
        cell(block, iCell);
      OPERATOR{}.apply(cell, params);
    }
  }
};

template<typename T, typename DESCRIPTOR>
class BlockD : public BlockStructure<DESCRIPTOR>
             , public Serializable {
protected:
  const Platform _platform;

  /// Calls f on concrete implementation of the BlockD interface
  template<typename F>
  auto dispatch(F&& f) const {
    return callUsingConcretePlatform<ConcretizableBlockD<T,DESCRIPTOR>>(
      _platform,
      this,
      f);
  }
  template<typename F>
  auto dispatch(F&& f) {
    return callUsingConcretePlatform<ConcretizableBlockD<T,DESCRIPTOR>>(
      _platform,
      this,
      f);
  }

public:
  BlockD(Vector<int,DESCRIPTOR::d> size, int padding, Platform platform);
  virtual ~BlockD();

  /// Set processing context
  /**
   * This is currently used to trigger data transfers between host
   * and GPU data for Platform::GPU_CUDA.
   **/
  virtual void setProcessingContext(ProcessingContext) = 0;

  virtual void resize(Vector<int,DESCRIPTOR::d> size) = 0;

  /// Return platform used to process lattice
  Platform getPlatform() const {
    return _platform;
  }

  template <Platform PLATFORM>
  ConcreteBlockD<T,DESCRIPTOR,PLATFORM>& asConcrete() {
    if (auto* ptr = dynamic_cast<ConcreteBlockD<T,DESCRIPTOR,PLATFORM>*>(this)) [[likely]] {
      return *ptr;
    } else {
      throw std::runtime_error("Invalid PLATFORM");
    }
  }

  /// Return whether FIELD_TYPE is available / has been allocated
  template<typename FIELD_TYPE>
  bool hasData() {
    return dispatch([&](const auto* lattice) -> bool {
      return lattice->template hasData<FIELD_TYPE>();
    });
  }

  /// Return abstract interface for concrete FIELD_TYPE data
  template<typename FIELD_TYPE>
  const auto& getData(FIELD_TYPE field = FIELD_TYPE{}) const
  {
    return *dispatch([&](const auto* lattice) -> const auto* {
      return &(lattice->template getData<FIELD_TYPE>().asAbstract());
    });
  }
  /// Return abstract interface for concrete FIELD_TYPE data
  template<typename FIELD_TYPE>
  auto& getData(FIELD_TYPE field = FIELD_TYPE{})
  {
    return *dispatch([&](auto* lattice) -> auto* {
      return &(lattice->template getData<FIELD_TYPE>().asAbstract());
    });
  }

  /// Return abstract interface for FIELD array
  template<typename FIELD>
  const auto& getField(FIELD field = FIELD{}) const
  {
    return getData(Array<FIELD>{});
  }
  /// Return value of FIELD at loc
  template<typename FIELD>
  auto getField(LatticeR<DESCRIPTOR::d> loc) const
  {
    return getField<FIELD>().get(this->getCellId(loc));
  }

  /// Return abstract interface for FIELD array
  template<typename FIELD>
  auto& getField(FIELD field = FIELD{})
  {
    return getData(Array<FIELD>{});
  }
  /// Set FIELD to fieldD at loc
  template<typename FIELD>
  void setField(LatticeR<DESCRIPTOR::d> loc, const FieldD<T,DESCRIPTOR,FIELD>& fieldD)
  {
    getField<FIELD>().set(this->getCellId(loc), fieldD);
  }

  /// Set value of parameter FIELD
  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD> value) {
    dispatch([&](auto* lattice) {
      lattice->template setParameter<FIELD>(value);
    });
  }
  template <typename PARAMETER, typename FIELD>
  void setParameter(AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& fieldArray) {
    dispatch([&](auto* lattice) {
      lattice->template setParameter<PARAMETER>(fieldArray);
    });
  }
  template <typename PARAMETER, Platform PLATFORM, typename FIELD>
  void setParameter(FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& fieldArray) {
    dispatch([&](auto* lattice) {
      lattice->template setParameter<PARAMETER>(fieldArray);
    });
  }

  BlockLatticeRowView<T,DESCRIPTOR> get(std::size_t iCell) {
    return {*this, iCell};
  }

  template <typename OPERATOR>
  void apply() {
    dispatch([&](auto* lattice) {
      lattice->template apply<OPERATOR>();
    });
  }

  virtual bool hasCommunicatable(std::type_index) const = 0;
  virtual Communicatable& getCommunicatable(std::type_index) = 0;

};

/// Implementation of BlockD on a concrete PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockD final : public BlockD<T,DESCRIPTOR> {
private:
  /// Field data
  ConcreteData<T,DESCRIPTOR,PLATFORM> _data;
  /// Index of DESCRIPTOR-declared field arrays
  utilities::FixedTypeIndexedMap<typename DESCRIPTOR::fields_t, ColumnVectorBase*> _descriptorFields;
  /// Pointers to Communicatable-casted FieldArrayD instances for overlap communication
  std::map<std::type_index, std::unique_ptr<Communicatable>> _communicatables;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  static constexpr Platform platform = PLATFORM;

  ConcreteBlockD(Vector<int,DESCRIPTOR::d> size, int padding=0);

  void setProcessingContext(ProcessingContext context) override {
    _data.setProcessingContext(context);
  }

  void resize(Vector<int,DESCRIPTOR::d> size) override;

  template<typename FIELD_TYPE>
  bool hasData() const;
  template<typename FIELD_TYPE>
  const auto& getData() const;
  template<typename FIELD_TYPE>
  auto& getData();

  template<typename FIELD>
  const auto& getField(FIELD field = FIELD()) const;
  template<typename FIELD>
  auto& getField(FIELD field = FIELD());

  bool hasCommunicatable(std::type_index field) const override {
    return _communicatables.find(field) != _communicatables.end();
  }
  Communicatable& getCommunicatable(std::type_index field) override {
    return *_communicatables.at(field).get();
  }

  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD> value)
  {
    _data.template forEachCastable<AbstractedConcreteParameters<T,DESCRIPTOR>>([&](auto* parameters) {
      auto& params = parameters->asAbstract();
      if (params.template provides<FIELD>()) {
        params.template set<FIELD>(value);
        parameters->setProcessingContext(ProcessingContext::Simulation);
      }
    });
  }

  template <typename PARAMETER, Platform _PLATFORM, typename FIELD>
  void setParameter(FieldArrayD<T,DESCRIPTOR,_PLATFORM,FIELD>& fieldArray)
  {
    if constexpr (PLATFORM == Platform::GPU_CUDA) {
      static_assert(PLATFORM == _PLATFORM, "FieldArrayD must be available on PLATFORM");
    }
    FieldD<T,DESCRIPTOR,PARAMETER> fieldArrayPointers;
    for (unsigned iD=0; iD < fieldArray.d; ++iD) {
      if constexpr (PLATFORM == Platform::GPU_CUDA) {
        fieldArrayPointers[iD] = fieldArray[iD].deviceData();
      } else {
        fieldArrayPointers[iD] = fieldArray[iD].data();
      }
    }
    setParameter<PARAMETER>(std::move(fieldArrayPointers));
  }

  template <typename PARAMETER, typename FIELD>
  void setParameter(AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& abstractFieldArray)
  {
    auto& fieldArray = dynamic_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(abstractFieldArray);
    setParameter<PARAMETER>(fieldArray);
  }

  /// Return reference to Data's FieldTypeRegistry
  /**
   * Only of interest for implementing specific device support
   **/
  auto& getDataRegistry() {
    return _data.getRegistry();
  }

  template <typename OPERATOR>
  void apply() {
    static_assert(OPERATOR::scope == OperatorScope::PerBlock,
                  "OPERATOR is not PerBlock scoped");
    OPERATOR{}.apply(*this);
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  /// Reinit population structure after deserialization
  void postLoad() override;

};


/// Constructs BlockD accessible on all locally used platforms
/**
 * For data that is shared between blocks, e.g. reference porosity fields
 **/
template <typename T, typename DESCRIPTOR, typename... ARGS>
std::unique_ptr<BlockD<T,DESCRIPTOR>> makeSharedBlockD(
  LoadBalancer<T>& loadBalancer, ARGS&&... args)
{
  #if defined(PLATFORM_GPU_CUDA)
  if (loadBalancer.isLocal(Platform::GPU_CUDA)) {
    return std::unique_ptr<BlockD<T,DESCRIPTOR>>(
      new ConcreteBlockD<T,DESCRIPTOR,Platform::GPU_CUDA>(std::forward<ARGS&&>(args)...));
  } else {
    return std::unique_ptr<BlockD<T,DESCRIPTOR>>(
      new ConcreteBlockD<T,DESCRIPTOR,Platform::CPU_SISD>(std::forward<ARGS&&>(args)...));
  }
  #else
  return std::unique_ptr<BlockD<T,DESCRIPTOR>>(
    new ConcreteBlockD<T,DESCRIPTOR,Platform::CPU_SISD>(std::forward<ARGS&&>(args)...));
  #endif
}

}

#endif
