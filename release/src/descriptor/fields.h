/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
 *                2024 Dennis Teutscher
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

#ifndef DESCRIPTOR_FIELDS_H
#define DESCRIPTOR_FIELDS_H

#include <type_traits>
#include <stdexcept>
#include <optional>
#include "core/meta.h"
#include "core/vector.h"
#include "utilities/aliases.h"

#include "core/platform/column.h"
#include "core/matrixView.h"

namespace olb {

/// Forward declaration of Expr type used for code generation
class Expr;

namespace descriptors {

// *INDENT-OFF*

/// \defgroup descriptor
//@{

/// Base of a field whose size is defined by [C_0,C_1,C_2]^T * [1,D,Q]
template <unsigned C_0, unsigned C_1=0, unsigned C_2=0>
struct FIELD_BASE {
  FIELD_BASE() = default;

  /// Return value type of field
  /**
   * Most fields are stored using the same value type as the T type
   * parameter of their associated lattice. However this template
   * offers the possibility of declaring a different value type per
   * field. See TYPED_FIELD_BASE.
   **/
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractColumn<T>;

  /// Get size of field for descriptor
  template <typename DESCRIPTOR>
  static constexpr std::size_t size() {
    return C_0
         + C_1 * DESCRIPTOR::d
         + C_2 * DESCRIPTOR::q;
  }

  // Initial value used for allocation of field data in FieldArrayD (optional)
  /**
   * Return value must be a correctly sized and typed olb::Vector.
   **/
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<C_0,C_1,C_2>>()>{};
  }

  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return true;
  }

  static constexpr bool isSerializable() {
    return true;
  }
};

struct SCALAR_FIELD     : public FIELD_BASE<1> { };
struct SPATIAL_FIELD    : public FIELD_BASE<0,1> { };
struct POPULATION_FIELD : public FIELD_BASE<0,0,1> { };

struct FIELD_BASE_CUSTOM_SIZE {
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractColumn<T>;

  static constexpr bool isSerializable() {
    return true;
  }
};

struct NEIGHBOR_FIELD : public FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
           if constexpr (DESCRIPTOR::d == 1) {
      return 2;
    } else if constexpr (DESCRIPTOR::d == 2) {
      return 8;
    } else if constexpr (DESCRIPTOR::d == 3) {
      return 26;
    }
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<NEIGHBOR_FIELD>()>{};
  }
};


/// Base of a descriptor field of scalar TYPE
template <typename TYPE, unsigned... Cs>
struct TYPED_FIELD_BASE : public FIELD_BASE<Cs...> {
  template <typename T>
  // TODO Find better solution for this hacky workaround to instantiate the entire lattice with Expr
  using value_type = std::conditional_t<std::is_same_v<T,Expr>,
                                        std::conditional_t<std::is_pointer_v<TYPE> || std::is_enum_v<TYPE>,
                                                           TYPE,
                                                           Expr>,
                                        TYPE>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<Cs...>>()>{};
  }

  static constexpr bool isSerializable() {
    return !std::is_pointer_v<TYPE>;
  }
};

/// Base of a field of scalar TYPE<BASE_TYPE> with dimensions same as FIELD_BASE
template <template<typename> typename TYPE, unsigned... Cs>
struct TEMPLATE_FIELD_BASE : public FIELD_BASE<Cs...> {
  template <typename T>
  using value_type = TYPE<T>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<Cs...>>()>{};
  }

  static constexpr bool isSerializable() {
    return false;
  }
};

/// Base of a field of pointers to the individual columns of FIELD
template <typename FIELD>
struct POINTER_FIELD_BASE {
  template <typename T>
  using value_type = std::add_pointer_t<typename FIELD::template value_type<T>>;

  template <typename T>
  using column_type = AbstractColumn<value_type<T>>;

  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return FIELD::template size<DESCRIPTOR>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<POINTER_FIELD_BASE>()>{};
  }

  static constexpr bool isSerializable() {
    return false;
  }
};

/// Base of a implicitly propagatable descriptor field
struct PROPAGATABLE_FIELD_BASE : public POPULATION_FIELD {
  template <typename T>
  using value_type = T;

  template <typename T>
  using column_type = AbstractCyclicColumn<T>;
};

/// Evaluates to true iff FIELD is propagatable (e.g. POPULATION)
template <typename FIELD>
using is_propagatable_field = typename std::is_base_of<PROPAGATABLE_FIELD_BASE, FIELD>::type;

/// Base of a tensor-valued descriptor field
struct TENSOR : public FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return (DESCRIPTOR::d * (DESCRIPTOR::d+1)) / 2; // see `TensorVal` in `core/util.h`
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};

/// Base of a matrix-valued descriptor field
template <unsigned ROWS, unsigned COLS>
struct MATRIX : public FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return ROWS*COLS;
  }

  static constexpr unsigned rows() {
    return ROWS;
  }

  static constexpr unsigned cols() {
    return COLS;
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<ROWS, COLS>()>{};
  }
};

/// Base of a matrix-valued descriptor field based on two FIELDS
template <typename FIELD, typename FIELD2>
struct FIELD_MATRIX : public FIELD_BASE_CUSTOM_SIZE  {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return DESCRIPTOR::template size<FIELD>() * DESCRIPTOR::template size<FIELD2>();
  }

  template <typename DESCRIPTOR>
  static constexpr unsigned rows() {
    return DESCRIPTOR::template size<FIELD>();
  }

  template <typename DESCRIPTOR>
  static constexpr unsigned cols() {
    return DESCRIPTOR::template size<FIELD2>();
  }

  template <typename T, typename DESCRIPTOR>
  static auto getMatrixView(FieldD<T,DESCRIPTOR,FIELD_MATRIX<FIELD,FIELD2>>& field) {
    return MatrixD<T,DESCRIPTOR,FIELD_MATRIX<FIELD,FIELD2>>(field);
  }
  template <typename T, typename DESCRIPTOR>
  static auto getTransposedMatrixView(FieldD<T,DESCRIPTOR,FIELD_MATRIX<FIELD,FIELD2>>& field) {
    return TransposedMatrixD<T,DESCRIPTOR,FIELD_MATRIX<FIELD,FIELD2>>(field);
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};

/// Base of degrees-of-rotational-freedom-dimensionalized fields
struct ROTATION_FIELD_BASE : public FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    if constexpr (DESCRIPTOR::d == 2) {
      return 1;
    } else {
      return DESCRIPTOR::d;
    }
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }

};

/// Base of a descriptor field of pointer type
template <typename TYPE>
struct OBJECT_POINTER_FIELD_BASE : public TYPED_FIELD_BASE<std::add_pointer_t<TYPE>,1> {
  static_assert(std::is_object_v<TYPE>, "TYPE must be object");
};

/// \defgroup descriptor_fields Set of common descriptor fields
/// \ingroup descriptor
//@{

struct CELL_ID      : public TYPED_FIELD_BASE<std::size_t,1> { };
struct MATERIAL     : public TYPED_FIELD_BASE<int,        1> { };
struct LATTICE_TIME : public TYPED_FIELD_BASE<std::size_t,1> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0;
  }
};

struct POPULATION : public PROPAGATABLE_FIELD_BASE { };

struct STATISTIC_GENERATED : public TYPED_FIELD_BASE<int,1> { };
struct STATISTIC           : public FIELD_BASE<2> { };

// Field types need to be distinct (i.e. not aliases)
// (Field size parametrized: Cs + Ds*D + Qs*Q)  Cs Ds Qs
struct VELOCITY             : public FIELD_BASE<0,  1, 0> { };
struct WMVELOCITY           : public FIELD_BASE<0,  1, 0> { };
struct WALL_VELOCITY        : public FIELD_BASE<0,  1, 0> { };
struct EXTERNAL_VELOCITY    : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY2            : public FIELD_BASE<0,  1, 0> { };
struct DAMPING              : public FIELD_BASE<1,  0, 0> { };
struct UX                   : public FIELD_BASE<1,  0, 0> { };
struct UY                   : public FIELD_BASE<1,  0, 0> { };
struct UZ                   : public FIELD_BASE<1,  0, 0> { };
struct TEMPGRADIENT         : public FIELD_BASE<0,  1, 0> { };
struct AVERAGE_VELOCITY     : public FIELD_BASE<0,  1, 0> { };
struct MAX_VELOCITY         : public FIELD_BASE<1,  0, 0> { };
struct AVERAGE_DENSITY      : public FIELD_BASE<1,  0, 0> { };
struct AVERAGE_TKE          : public FIELD_BASE<0,  1, 0> { };
struct SOURCE               : public FIELD_BASE<1,  0, 0> { };
struct U_TAU                : public FIELD_BASE<1,  0, 0> { };
struct PRESSCORR            : public FIELD_BASE<1,  0, 0> { };
struct FORCE                : public FIELD_BASE<0,  1, 0> { };
struct EXTERNAL_FORCE       : public FIELD_BASE<0,  1, 0> { };
struct NABLARHO             : public FIELD_BASE<0,  1, 0> { };
struct TAU_EFF              : public FIELD_BASE<1,  0, 0> { };
struct RHO                  : public FIELD_BASE<1,  0, 0> { };
struct VISCOSITY            : public FIELD_BASE<1,  0, 0> { };
struct GAMMA                : public FIELD_BASE<1,  0, 0> { };
struct CUTOFF_KIN_ENERGY    : public FIELD_BASE<1,  0, 0> { };
struct CUTOFF_HEAT_FLUX     : public FIELD_BASE<1,  0, 0> { };
struct CHEM_POTENTIAL       : public FIELD_BASE<1,  0, 0> { };
struct ADDEND               : public FIELD_BASE<1,  0, 0> { };
struct V6                   : public FIELD_BASE<6,  0, 0> { };
struct V12                  : public FIELD_BASE<12, 0, 0> { };
struct OMEGA                : public FIELD_BASE<1,  0, 0> {
  template <typename T, typename DESCRIPTOR,typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};
struct PREV_PHI             : public FIELD_BASE<1,  0, 0> { };
struct TAU                  : public FIELD_BASE<1,  0, 0> { };
struct INTERFACE_WIDTH      : public FIELD_BASE<1,  0, 0> { };
struct MAGIC                : public FIELD_BASE<1,  0, 0> { };
struct G                    : public FIELD_BASE<0,  1, 0> { };
struct EPSILON              : public FIELD_BASE<1,  0, 0> { };
struct BODY_FORCE           : public FIELD_BASE<0,  1, 0> { };
struct K                    : public FIELD_BASE<1,  0, 0> { };
struct NU                   : public FIELD_BASE<1,  0, 0> { };
struct VELOCITY_NUMERATOR   : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY_DENOMINATOR : public FIELD_BASE<1,  0, 0> { };
struct ZETA                 : public FIELD_BASE<0,  0, 1> { };
struct LOCAL_DRAG           : public FIELD_BASE<0,  1, 0> { };
struct VELOCITY_SOLID       : public FIELD_BASE<0,  1, 0> { };
struct COORDINATE           : public FIELD_BASE<0,  1, 0> { };
struct STRAINRATE           : public FIELD_BASE<6,  0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<STRAINRATE>()>(T(0));
  }
};
struct DENSITY           : public FIELD_BASE<1, 0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<DENSITY>()>(T(1.));
  }
};
struct NEIGHBOR             : public FIELD_BASE<1,  0, 0>{
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,1>(0.);
  }
};
struct CONVERSION           : public FIELD_BASE<1,  0, 0> { };
struct NORMALIZE            : public FIELD_BASE<1,  0, 0> { };
struct DX                   : public FIELD_BASE<1,  0, 0> { };
struct AV_SHEAR             : public FIELD_BASE<1,  0, 0> { };
struct SHEAR_RATE_MAGNITUDE : public FIELD_BASE<1,  0, 0> { };
struct TAU_W                : public FIELD_BASE<1,  0, 0> { };
struct SCALAR               : public FIELD_BASE<1,  0, 0> { };
struct SMAGO_CONST          : public FIELD_BASE<1,  0, 0> { };
struct EFFECTIVE_OMEGA      : public FIELD_BASE<1,  0, 0> { };
struct VELO_GRAD            : public FIELD_BASE<0,  3, 0> { };
struct FIL_RHO              : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_X      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Y      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Z      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_DISS        : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_TKE         : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_SIGMA_ADM      : public FIELD_BASE<1,  0, 0> { };
struct LOCAL_NU_EDDY        : public FIELD_BASE<1,  0, 0> { };
struct FILTERED_VEL_GRAD    : public FIELD_BASE<0,  3, 0> { };
struct ERROR_COVARIANCE     : public FIELD_BASE<1,  0, 0> { };
struct VARIANCE             : public FIELD_BASE<1,  0, 0> { };
struct TAU_SGS              : public FIELD_BASE<1,  0, 0> { };
struct FILTERED_POPULATION  : public FIELD_BASE<0,  0, 1> { };
struct INDICATE             : public FIELD_BASE<1,  0, 0> { };
struct BIOGAS_INSTANT       : public FIELD_BASE<1,  0, 0> { };
struct BIOGAS_CUMULATIVE    : public FIELD_BASE<1,  0, 0> { };
struct METHANE_INSTANT      : public FIELD_BASE<1,  0, 0> { };
struct METHANE_CUMULATIVE   : public FIELD_BASE<1,  0, 0> { };
struct CO2_INSTANT          : public FIELD_BASE<1,  0, 0> { };
struct CO2_CUMULATIVE       : public FIELD_BASE<1,  0, 0> { };
struct TEMPERATURE          : public FIELD_BASE<1,  0, 0> { };
struct SIGMA                : public FIELD_BASE<1,  0, 0> { };
struct PSI                  : public FIELD_BASE<1,  0, 0> { };
struct THETA                : public FIELD_BASE<1,  0, 0> { };
struct NORMGRADPSI          : public FIELD_BASE<1,  0, 0> { };
struct PSI0                 : public FIELD_BASE<1,  0, 0> { };
struct INTERPHASE_NORMAL    : public FIELD_BASE<0,  1, 0> { };
struct MASS                 : public FIELD_BASE<1,  0, 0> { };
struct CELL_TYPE            : public FIELD_BASE<1,  0, 0> { };
struct BOUNDARY             : public FIELD_BASE<1,  0, 0> { };
struct CONV_POPS            : public FIELD_BASE<0,  0, 1> { };
struct SOURCE_OLD           : public FIELD_BASE<1,  0, 0> { };
struct TOP                  : public FIELD_BASE<1,  0, 0> { };
struct BOTTOM               : public FIELD_BASE<1,  0, 0> { };
struct OLD_PHIU             : public FIELD_BASE<0,  1, 0> { };
struct PHIWETTING           : public FIELD_BASE<1, 0, 0> { };
struct CONTACT_DETECTION    : public TYPED_FIELD_BASE<size_t, 1,  0, 0> { };
struct POROSITY             : public FIELD_BASE<1,  0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<POROSITY>()>{1};
  }

  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0 && value[0]<=1;
  }
};
struct WMPOROSITY             : public FIELD_BASE<1,  0, 0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<POROSITY>()>{1};
  }

  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value >= 0 && value[0]<=1;
  }
};
struct SAMPLING_DISTANCE : public descriptors::FIELD_BASE<1,0,0> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<SAMPLING_DISTANCE>()>(3);
  }
};
struct POROSITY2            : public FIELD_BASE<1,  0, 0> { };
struct EUL2LAGR             : public FIELD_BASE<1,  0, 0> { };
struct LOCATION             : public FIELD_BASE<0,  1, 0> { };
struct Y1                   : public FIELD_BASE<0,  1, 0> { };
struct Y2                   : public FIELD_BASE<0,  1, 0> { };
struct VORTICITY            : public FIELD_BASE<1,  0, 0> { };

//TODO: This expression has been removed on master lately. As no obvious equivalent could be found immediately,
//      it is added back in to enable some functionality on feature/unifiedParticleFramework.
//      In case a proper equivalent exists, this expression can be removed for good!
template <
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD equals the head of our field list, terminate recursion
  std::enable_if_t<std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
  >
constexpr unsigned getIndexInFieldList()
{
  return 0;
}

template <
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD doesn't equal the head of our field list
  std::enable_if_t<!std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
  >
constexpr unsigned getIndexInFieldList()
{
  // Break compilation when WANTED_FIELD is not provided by list of fields
  static_assert(sizeof...(FIELDS) > 0, "Field not found.");

  return 1 + getIndexInFieldList<WANTED_FIELD,FIELDS...>();
}

//@}

//@}

}

namespace fields {

struct PHYS_R : public descriptors::FIELD_BASE<0,1> { };
struct BLOCK_LOWER : public descriptors::FIELD_BASE<0,1> { };
struct BLOCK_UPPER : public descriptors::FIELD_BASE<0,1> { };

namespace membrane {

struct VELOCITY : public descriptors::FIELD_BASE<0,1> { };
struct OLD_VELOCITY : public descriptors::FIELD_BASE<0,1> { };
struct FORCE : public descriptors::FIELD_BASE<0,1> { };

struct STENCIL_WIDTH : public descriptors::FIELD_BASE<1> { };

}

namespace moments {

struct U : public descriptors::FIELD_BASE<0,1> { };

}

}

namespace cse {

struct NEIGHBOR_PREFIX : public descriptors::FIELD_BASE<1> { };

class SymbolGenerator;
struct SCALAR_SYMBOL_GENERATOR : public descriptors::OBJECT_POINTER_FIELD_BASE<cse::SymbolGenerator> { };
struct VECTOR_SYMBOL_GENERATOR : public descriptors::OBJECT_POINTER_FIELD_BASE<cse::SymbolGenerator> { };
struct TENSOR_SYMBOL_GENERATOR : public descriptors::OBJECT_POINTER_FIELD_BASE<cse::SymbolGenerator> { };
struct POPULATION_SYMBOL_GENERATOR : public descriptors::OBJECT_POINTER_FIELD_BASE<cse::SymbolGenerator> { };
struct FUNCTION_CALLS : public descriptors::OBJECT_POINTER_FIELD_BASE<std::vector<Expr>> { };

}

namespace opti {

/// Objective functional results, always scalar
struct J                    : public descriptors::FIELD_BASE<1,  0, 0> { };
/// Stores populations of the primal problems for adjoint simulations
struct F                    : public descriptors::FIELD_BASE<0,  0, 1> { };
/// Derivative of Objective functional regarding populations
struct DJDF                 : public descriptors::FIELD_BASE<0,  0, 1> { };
/// Derivative of control projection regarding control variable
template <typename CONTROLS>
struct DPROJECTIONDALPHA : public descriptors::FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return DESCRIPTOR::template size<CONTROLS>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};
/// Total derivative of the objective functional regarding controls
template <typename CONTROLS>
struct SENSITIVITY : public descriptors::FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return DESCRIPTOR::template size<CONTROLS>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};
/// Objective functional derivative regarding optimization controls
template <typename CONTROLS>
struct DJDALPHA : public descriptors::FIELD_BASE_CUSTOM_SIZE {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return DESCRIPTOR::template size<CONTROLS>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, size<DESCRIPTOR>()>{};
  }
};
/// Collision operator derivative regarding optimization controls
template <typename CONTROLS>
struct DCDALPHA : public descriptors::FIELD_MATRIX<descriptors::POPULATION,CONTROLS> {
  template <typename DESCRIPTOR>
  static constexpr unsigned size() {
    return descriptors::FIELD_MATRIX<descriptors::POPULATION,CONTROLS>::template size<DESCRIPTOR>();
  }

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    using V = typename descriptors::FIELD_MATRIX<descriptors::POPULATION,CONTROLS>::template value_type<T>;
    constexpr unsigned DIM = descriptors::FIELD_MATRIX<descriptors::POPULATION,CONTROLS>::template size<DESCRIPTOR>();
    return Vector<V,DIM>{};
  }
};

}

}

// *INDENT-ON*

#endif
