/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_DESCRIPTORS_H
#define PARTICLE_DESCRIPTORS_H

#include "utilities/dimensionConverter.h"
#include "particles/descriptor/particleDescriptorUtilities.h"

namespace olb {

namespace descriptors {

//INFO: Initial values can be defined by providing getInitialValue()

//Parent fields ensuring generalized access
struct ANG_VELOCITY {};
struct ANG_ACC_STRD {};
struct ANGLE {};
struct SINDICATOR {};
struct TORQUE {};
struct MOFI {};           //Momentum of inertia
struct ROT_MATRIX {};

//Common
struct POSITION          : public FIELD_BASE<0,  1, 0> { };
struct INVALID           : public FIELD_BASE<1,  0, 0> {  //Storage invalidation, filtered in particleManager
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<INVALID>()>{false};
  }
};
struct YOUNG_MODULUS     : public FIELD_BASE<1,  0, 0> { };
struct SHEAR_MODULUS     : public FIELD_BASE<1,  0, 0> { };
struct POISSON_RATIO     : public FIELD_BASE<1,  0, 0> { };
struct ADHESION          : public FIELD_BASE<2,  0, 0> { };
struct DYNAMICS_ID       : public TYPED_FIELD_BASE<unsigned short, 1,  0, 0> { };
struct COMPUTE_MOTION    : public TYPED_FIELD_BASE<bool, 1, 0, 0> { };
struct COMPUTE_CONTACT   : public TYPED_FIELD_BASE<bool, 1, 0, 0> { };
//Resolved
struct ACCELERATION_STRD : public FIELD_BASE<0,  1, 0> { };
struct SURFACE_ID        : public FIELD_BASE<1,  0, 0> { };
struct COR_OFFSET        : public FIELD_BASE<0,  1, 0> { }; //Centre of rotation offset (considering untransformed coordinate system)
struct IC                : public TYPED_FIELD_BASE<int, 1,  0, 0> { };
struct DETACHING         : public TYPED_FIELD_BASE<bool, 1, 0, 0> { };
struct ELONGATION        : public FIELD_BASE<0,  1, 0> { };
//Subgrid
struct FORCE_STRD        : public FIELD_BASE<0,  1, 0> { };
struct ACTIVE            : public TYPED_FIELD_BASE<bool, 1, 0, 0> { //Disables movement, usually used in dynamics
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,DESCRIPTOR::template size<ACTIVE>()>{true};
  }
};

struct ID                : public TYPED_FIELD_BASE<size_t, 1, 0, 0> { };
struct MASS_ADDED        : public FIELD_BASE<1,  0, 0> { };
struct RADIUS            : public FIELD_BASE<1,  0, 0> { };
struct SPECIES           : public FIELD_BASE<1,  0, 0> { };
struct FLUIDVEL          : public FIELD_BASE<0,  1, 0> { };

// Contact model
struct ENLARGEMENT_FOR_CONTACT  : public FIELD_BASE<1,  0, 0> { };
struct IS_IN_CONTACT            : public FIELD_BASE<1,  0, 0> { };

//Dimension sensitive fields
template<unsigned D>
struct ANG_VELOCITY_XD   : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANG_VELOCITY { };
template<unsigned D>
struct ANG_ACC_STRD_XD   : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANG_ACC_STRD { };
template<unsigned D>
struct ANGLE_XD          : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public ANGLE { };
template<unsigned D>
struct TORQUE_XD         : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public TORQUE { };
template<unsigned D>
struct MOFI_XD           : public FIELD_BASE<utilities::dimensions::convert<D>::rotation,  0, 0>, public MOFI { };
template<unsigned D>
struct ROT_MATRIX_XD     : public FIELD_BASE<utilities::dimensions::convert<D>::matrix,  0, 0>, public ROT_MATRIX { };


namespace Haider_Levenspiel
{
//Haider Levenspiel constatns
struct A                 : public FIELD_BASE<1,0,0> {};
struct B                 : public FIELD_BASE<1,0,0> {};
struct C                 : public FIELD_BASE<1,0,0> {};
struct D                 : public FIELD_BASE<1,0,0> {};
}
//TRAN CONG
struct CIRCULARITY          : public FIELD_BASE<1,0,0> {};
struct PROJECTED_SURFACE    : public FIELD_BASE<1,0,0> {};
struct DA_DN                : public FIELD_BASE<1,0,0> {};

//EULER LAGRANGE EULER ROTATION SPHEROID
struct QUATERNION : public FIELD_BASE <4,0,0> {};
struct MOMENT_OF_INERTIA : public FIELD_BASE <3,0,0> {};
struct ORIENTATION : public FIELD_BASE <3,0,0> {};
struct ORIENTATION_ANGLE : public FIELD_BASE <1,0,0> {};
struct TRANSLATIONAL_DYADIC : public FIELD_BASE <9,0,0> {};
struct BETA : public FIELD_BASE <1,0,0> {};
///dimensionless constant for Hydrodynamic Torque computation
struct ABG0 : public FIELD_BASE <3,0,0> {}; ///dimensionless constant for Hydrodynamic Torque computation
///used for visualisation purposes in Paraview, scales the spheroid in main axis directions ELER
struct SCALING : public FIELD_BASE <3,0,0> {};
///Paraview visualisation, different notation of the Quaternions than in ELER
struct PARAQUATERNION: public FIELD_BASE <4,0,0> {};
struct LABEL : public FIELD_BASE <1,0,0> {};



//Indicator
template<unsigned D>
struct SINDICATOR_XD : public SINDICATOR, public FIELD_BASE<1> {
  template <typename T>
  using value_type = typename utilities::dimensions::convert<D>::template surfaceType<T>*;

  template <typename T>
  using column_type = AbstractColumn<typename utilities::dimensions::convert<D>::template surfaceType<T>*>;

  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>, DESCRIPTOR::template size<FIELD_BASE<1>>()>{};
  }
};

//Counter
template<typename NAME>
struct COUNTER          : public TYPED_FIELD_BASE<size_t,1,0,0> {};

//Parent descriptors ensuring generalized access
struct GENERAL {};
struct MOBILITY {};
struct SURFACE {};
struct FORCING {};
struct PHYSPROPERTIES {};
struct DYNBEHAVIOUR {};
/// Mechanical properties
struct MECHPROPERTIES {};
struct NUMERICPROPERTIES {};
/// Communication
struct PARALLELIZATION {};

//Field wrapping descriptors acting as GROUP
// -like the field list itself, this list can be extended by custom GROUP descriptors
// -in order to allow generalized access via getField<GROUP,FIELD> those should always inherit from the parent descriptors
template<unsigned D>
struct GENERAL_TMP                : public PARTICLE_DESCRIPTOR<D,POSITION>, public GENERAL {};

template<unsigned D>
struct GENERAL_EXTENDABLE         : public PARTICLE_DESCRIPTOR<D,POSITION,INVALID>, public GENERAL {};

template<unsigned D>
struct MOBILITY_VERLET            : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD,ANG_VELOCITY_XD<D>,ANG_ACC_STRD_XD<D>>, public MOBILITY {};

template<unsigned D>
struct MOBILITY_VERLET_NO_ANGLE   : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD>, public MOBILITY {};

template<unsigned D>
struct MOBILITY_EULER_NO_ANGLE    : public PARTICLE_DESCRIPTOR<D,VELOCITY>, public MOBILITY {};

struct DYNBEHAVIOUR_BASIC         : public PARTICLE_DESCRIPTOR<1,ACTIVE>, public DYNBEHAVIOUR {};

struct DYNBEHAVIOUR_MULTI_DYN     : public PARTICLE_DESCRIPTOR<1,DYNAMICS_ID,COMPUTE_MOTION,COMPUTE_CONTACT>, public DYNBEHAVIOUR {};

struct DYNBEHAVIOUR_DETACHABLE    : public PARTICLE_DESCRIPTOR<1,DETACHING,ACTIVE,COUNTER<ACTIVE>,
                                      COMPUTE_CONTACT>, public DYNBEHAVIOUR {};

template<unsigned D>
struct SURFACE_RESOLVED           : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,ROT_MATRIX_XD<D>,SINDICATOR_XD<D>>, public SURFACE {};

template<unsigned D>
struct SURFACE_RESOLVED_COR       : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,ROT_MATRIX_XD<D>,SINDICATOR_XD<D>,COR_OFFSET>, public SURFACE {};

template<unsigned D>
struct SURFACE_RESOLVED_CIRCULAR  : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,SINDICATOR_XD<D>>, public SURFACE {};

template<unsigned D>
struct SURFACE_RESOLVED_PARALLEL  : public PARTICLE_DESCRIPTOR<D,ANGLE_XD<D>,ROT_MATRIX_XD<D>,SINDICATOR_XD<D>,SURFACE_ID>, public SURFACE {};

template<unsigned D>
struct NUMERICPROPERTIES_RESOLVED_CONTACT
                                  : public PARTICLE_DESCRIPTOR<D,ENLARGEMENT_FOR_CONTACT>, public NUMERICPROPERTIES {};

template<unsigned D>
struct NUMERICPROPERTIES_RESOLVED_LUBRICATION_FORCE_CORRECTION
                                  : public PARTICLE_DESCRIPTOR<D,IS_IN_CONTACT>, public NUMERICPROPERTIES {};

template<unsigned D>
struct FORCING_RESOLVED           : public PARTICLE_DESCRIPTOR<D,FORCE,TORQUE_XD<D>>, public FORCING {};

template<unsigned D>
struct FORCING_ADHESIVE           : public PARTICLE_DESCRIPTOR<D,FORCE,TORQUE_XD<D>,ADHESION>, public FORCING {};

template<unsigned D>
struct FORCING_SUBGRID            : public PARTICLE_DESCRIPTOR<D,FORCE,FORCE_STRD>, public FORCING {};

template<unsigned D>
struct PHYSPROPERTIES_RESOLVED    : public PARTICLE_DESCRIPTOR<D,MASS,MOFI_XD<D>>, public PHYSPROPERTIES {};

template<unsigned D>
struct PHYSPROPERTIES_RESOLVED_PERMEABLE : public PARTICLE_DESCRIPTOR<D,MASS,MOFI_XD<D>,POROSITY>, public PHYSPROPERTIES {};

template<unsigned D>
struct PHYSPROPERTIES_SUBGRID     : public PARTICLE_DESCRIPTOR<D,MASS,MASS_ADDED,MOFI_XD<D>,RADIUS>, public PHYSPROPERTIES {};

template<unsigned D>
struct PHYSPROPERTIES_SUBGRID_REACTIVE  : public PARTICLE_DESCRIPTOR<D,MASS,MASS_ADDED,MOFI_XD<D>,RADIUS,SPECIES>, public PHYSPROPERTIES {};

template<unsigned D>
struct MECHPROPERTIES_COLLISION   : public PARTICLE_DESCRIPTOR<D,MATERIAL,YOUNG_MODULUS,SHEAR_MODULUS,POISSON_RATIO>, public MECHPROPERTIES {};

struct PARALLELIZATION_SUBGRID    : public PARTICLE_DESCRIPTOR<1,ID>, public PARALLELIZATION {};

struct PARALLELIZATION_RESOLVED   : public PARTICLE_DESCRIPTOR<1,ID,IC>, public PARALLELIZATION {};

template<unsigned D>
struct MOBILITY_VERLET_NO_ANGLE_FLUIDVEL   : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD, FLUIDVEL>, public MOBILITY {};

template<unsigned D>
struct MOBILITY_ANALYTIC_NO_ANGLE_FLUIDVEL   : public PARTICLE_DESCRIPTOR<D,VELOCITY,ACCELERATION_STRD, FLUIDVEL>, public MOBILITY {};

template<unsigned D>
struct PHYSPROPERTIES_SUBGRID_DENSITY : public PARTICLE_DESCRIPTOR<D,MASS,MASS_ADDED,MOFI_XD<D>,RADIUS, DENSITY>, public PHYSPROPERTIES {};


struct HAIDER_LEVENSPIEL          : public PARTICLE_DESCRIPTOR<1,Haider_Levenspiel::A,Haider_Levenspiel::B,Haider_Levenspiel::C,Haider_Levenspiel::D,SCALING,LABEL>, public NUMERICPROPERTIES {};

struct TRAN_CONG : public PARTICLE_DESCRIPTOR<1,CIRCULARITY, PROJECTED_SURFACE, DA_DN,SCALING,LABEL>, public NUMERICPROPERTIES {};

struct EULER_ROTATION : public PARTICLE_DESCRIPTOR<3, MOMENT_OF_INERTIA, ROT_MATRIX_XD<3>, ANG_VELOCITY_XD<3>, ORIENTATION, QUATERNION,TRANSLATIONAL_DYADIC,BETA,ABG0,SCALING,PARAQUATERNION,ORIENTATION_ANGLE, LABEL>, public NUMERICPROPERTIES {};


} //namespace descriptors

} //namespace olb


#endif
