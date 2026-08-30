/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023-24 Julius Jeßberger
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

/** \file Serialization of data on super geometry level.
 * Please do not use this in performance-relevant situations since it breaks the
 * parallelization concept (!)
 */

 #ifndef SERIALIZATION_H
 #define SERIALIZATION_H

#include "optimization/solver/controller.h"
#include "geometry/superGeometry.h"
#include "geometry/cuboidDecomposition.h"
#include "utilities/aliases.h"
#include "utilities/vectorHelpers.h"
#include "optimization/core/optiCaseDual.h"

namespace olb {

namespace opti {

template<typename S>
class Controller;



/////////////////////////// SERIALIZATION CLASSES /////////////////////////////


/** This class serializes the cells inside the geometry.
 * Intended for creating one long vector out of field data on super level.
 * Methods for indexing cells as well as field components are offered.
 */
template<typename S, unsigned dim>
class GeometrySerializer
{
public:
  /// Bundle for lattice coordinates + field component
  struct LatticeAndFieldR {
    LatticeR<dim+1> latticeR;
    std::size_t     fieldComponent;
  };

  /// Compute serialized cell index from lattice coordinates
  virtual std::size_t getSerializedCellIndex(const int latticeR[]) const = 0;

  /// Compute serialized cell index from lattice coordinates
  std::size_t getSerializedCellIndex(LatticeR<dim+1> latticeR) const {
    return getSerializedCellIndex(latticeR.data());
  }

  /// Get lattice coordinates from serialized cell index
  virtual LatticeR<dim+1> getLatticeR(std::size_t index) const = 0;

  /// Get index of field component from lattice coordinates and component index
  std::size_t getSerializedComponentIndex(const int latticeR[],
    unsigned iD, unsigned fieldDim) const {
    return fieldDim * getSerializedCellIndex(latticeR) + iD;
  }

  /// Get index of field component from lattice coordinates and component index
  std::size_t getSerializedComponentIndex(LatticeR<dim+1> latticeR,
    unsigned iD, unsigned fieldDim) const {
    return fieldDim * getSerializedCellIndex(latticeR) + iD;
  }

  /// Get index of field component from lattice coordinates and component index
  std::size_t getSerializedComponentIndex(LatticeAndFieldR coords,
    unsigned fieldDim) const {
    return getSerializedComponentIndex(coords.latticeR, coords.fieldComponent, fieldDim);
  }

  /// Get lattice coordinates and field component from serialized field index
  LatticeAndFieldR getLatticeAndFieldR(std::size_t index, unsigned fieldDim) const {
    // index = fieldDim * cellIndex + fieldComponent
    const auto fieldCoord = std::div((long int) index, (long int) fieldDim);
    const auto geomCoord = getLatticeR(fieldCoord.quot);
    return {geomCoord, (std::size_t) fieldCoord.rem};
  }

  virtual unsigned getNoCells() const = 0;
};

/** This class serializes the cells inside the geometry.
 *
 * Formula for serialized cell index (3d): cuboid_offset + NX*NY*z + NX*y + x
 * Formula for serialized field component index: fieldDim * cellIndex + fieldComponent
 *
 * This class is not intended to be used in performance-relevant context.
 * Index computations could be accellerated if desired.
 */
template<typename S, unsigned dim>
class SimpleGeometrySerializer : public GeometrySerializer<S,dim>
{
protected:
  CuboidDecomposition<S,dim>&     _cGeometry;
  const unsigned                  _noCuboids;

private:
  std::vector<unsigned>           _cuboidSizes;   // number of cells per cuboid
  std::vector<unsigned>           _offsets;       // accumulated cuboid sizes
  unsigned                        _noCells;

public:
  SimpleGeometrySerializer(CuboidDecomposition<S,dim>& cGeometry)
   : _cGeometry(cGeometry),
     _noCuboids(_cGeometry.size())
  {
    _cuboidSizes.reserve(_noCuboids);
    _offsets.reserve(_noCuboids);
    for (unsigned i = 0; i < _noCuboids; ++i) {
      _cuboidSizes.push_back(_cGeometry.get(i).getLatticeVolume());
    }
    _offsets.push_back(0);
    for (unsigned i = 1; i < _noCuboids; ++i) {
      _offsets.push_back(_offsets[i-1] + _cuboidSizes[i-1]);
    }

    const auto mc = _cGeometry.getMotherCuboid();
    if constexpr (dim == 2) {
      _noCells = (mc.getNx() + 1) * (mc.getNy() + 1);
    } else {
      _noCells = (mc.getNx() + 1) * (mc.getNy() + 1) * (mc.getNz() + 1);
    }
  }

  SimpleGeometrySerializer(SuperGeometry<S,dim>& sGeometry)
   : SimpleGeometrySerializer(sGeometry.getCuboidDecomposition())
  { }

  /// Compute serialized cell index from lattice coordinates
  std::size_t getSerializedCellIndex(const int latticeR[]) const override {
    const Cuboid<S,dim>& c = _cGeometry.get(latticeR[0]);
    const int nX = c.getNx();
    const int nY = c.getNy();
    std::size_t res;
    if constexpr (dim == 2) {
      // cuboid_offset + NX*y + x
      res = _offsets[latticeR[0]] + nX*latticeR[2] + latticeR[1];
    } else {
      // cuboid_offset + NX*NY*z + NX*y + x
      res = _offsets[latticeR[0]] + nX*nY*latticeR[3] + nX*latticeR[2] + latticeR[1];
    }
    return res;
  }

  using GeometrySerializer<S,dim>::getSerializedCellIndex;

  /// Get lattice coordinates from serialized cell index
  LatticeR<dim+1> getLatticeR(std::size_t index) const override {
    const auto cuboidIt = std::upper_bound(_offsets.begin(), _offsets.end(), index) - 1;
    LatticeR<dim+1> res;
    res[0] = std::distance(_offsets.begin(), cuboidIt);

    const Cuboid<S,dim>& c = _cGeometry.get(res[0]);
    const std::size_t nX = c.getNx();
    if constexpr (dim == 2) {
      // index = cuboid_offset + NX*y + x
      const auto divByNx = std::ldiv((long int) index - *cuboidIt, (long int) nX);
      res[2] = divByNx.quot;
      res[1] = divByNx.rem;
    } else {
      // index = cuboid_offset + NX*NY*z + NX*y + x
      const std::size_t nY = c.getNy();
      const auto divByNxNy = std::ldiv((long int) index - *cuboidIt, (long int) nX*nY);
      res[3] = divByNxNy.quot;
      const auto divByNx = std::ldiv(divByNxNy.rem, (long int) nX);
      res[2] = divByNx.quot;
      res[1] = divByNx.rem;
    }
    return res;
  }

  unsigned getNoCells() const override {
    return _noCells;
  }
};


/** This class serializes the cells which are marked by indicator.
 * It produces a smaller control vector, but access is more expensive compared
 * to SimpleGeometrySerializer.
 *
 * @tparam S Floating point type
 * @tparam dim Geometrical dimension
 */
template<typename S, unsigned dim>
class SparseGeometrySerializer : public GeometrySerializer<S,dim>
{
private:
  using L = LatticeR<dim+1>;
  std::vector<L>                             _coords;
  mutable FunctorPtr<SuperIndicatorF<S,dim>> _indicator;  // indicates which cells are serialized
      // the _indicator should never be changed, otherwise this class loses its
      // correctness. The "mutable" keyword is only necessary because FunctorPtr
      // does not support const arithmetic and should be removed if possible.

public:
  SparseGeometrySerializer(SuperGeometry<S,dim>& superGeometry,
    FunctorPtr<SuperIndicatorF<S,dim>>&& indicator)
   : _indicator(std::move(indicator))
  {
    const auto cGeometry = superGeometry.getCuboidDecomposition();
    L latticeR;
    for (unsigned iC = 0; iC < cGeometry.size(); ++iC) {
      latticeR[0] = iC;
      const int nX = cGeometry.get(iC).getNx();
      const int nY = cGeometry.get(iC).getNy();
      for (int iX=0; iX<nX; iX++) {
        latticeR[1] = iX;
        for (int iY=0; iY<nY; iY++) {
          latticeR[2] = iY;
          if constexpr (dim == 2) {
            if (_indicator(latticeR.data())) {
              _coords.push_back(latticeR);
            }
          } else {
            const int nZ = cGeometry.get(iC).getNz();
            for (int iZ=0; iZ<nZ; iZ++) {
              latticeR[3] = iZ;

              if (_indicator(latticeR.data())) {
                _coords.push_back(latticeR);
              }
            }
          }
        }
      }
    }
  }

  SparseGeometrySerializer(SuperGeometry<S,dim>& superGeometry,
    FunctorPtr<IndicatorF<S,dim>>&& indicator)
   : SparseGeometrySerializer(superGeometry,
     new SuperIndicatorFfromIndicatorF<S,dim>(std::move(indicator), superGeometry))
  { }

  /// Compute serialized cell index from lattice coordinates
  std::size_t getSerializedCellIndex(const int latticeR[]) const override {
#ifdef OLB_DEBUG
    // check if cell lies in indicated domain
    bool isInside;
    _indicator(&isInside, latticeR);
    if (! isInside) {
      OstreamManager clout (std::cout, "SparseGeometrySerializer");
      clout << "Warning: the passed cell does not lie in the indicated domain "
        << "and can therefore not be accessed." << std::endl;
      std::exit(1);
    }
#endif
    const L point (latticeR);
    const auto upper = std::upper_bound(
      _coords.begin(), _coords.end(), point, [](const L& a, const L& b){
      return lex_smaller(a, b);
    });
    return std::distance(_coords.begin(), upper - 1);
  }

  using GeometrySerializer<S,dim>::getSerializedCellIndex;

  /// Get lattice coordinates from serialized cell index
  LatticeR<dim+1> getLatticeR(std::size_t index) const override {
#ifdef OLB_DEBUG
    // check if index lies in accessible range
    if (index >= _coords.size()) {
      OstreamManager clout (std::cout, "SparseGeometrySerializer");
      clout << "Warning: the passed index is higher than the number of cells "
        << "which are marked by indicator." << std::endl;
      std::exit(1);
    }
#endif
    return _coords[index];
  }

  unsigned getNoCells() const override {
    return _coords.size();
  }
};





/////////////////////// FUNCTORS BASED ON SERIAL DATA /////////////////////////

template <typename T, typename DESCRIPTOR>
class BlockLatticeSerialDataF : public BlockLatticeF<T,DESCRIPTOR> {
protected:
  static constexpr unsigned d = DESCRIPTOR::d;
  Controller<T>&                    _controller;
  const GeometrySerializer<T,d>&    _serializer;
  const int                         _iC;
  std::function<T(T)>               _projection;
public:
  BlockLatticeSerialDataF(BlockLattice<T,DESCRIPTOR>& blockLattice,
    Controller<T>& controller, const GeometrySerializer<T,d>& serializer,
    int iC, unsigned dataDim, std::function<T(T)> projection)
  : BlockLatticeF<T,DESCRIPTOR>(blockLattice, dataDim),
    _controller(controller), _serializer(serializer), _iC(iC),
    _projection(projection)
  {
    this->getName() = "BlockLatticeSerialDataF";
  }

  bool operator() (T output[], const int input[])
  {
    if (this->getBlock().isInsideCore(input)) {
      int latticeR[4] = {_iC, input[0], input[1], 0};
      if constexpr (d == 3) {
        latticeR[3] = input[2];
      }
      const auto index = _serializer.getSerializedComponentIndex(latticeR, 0, this->getTargetDim());
      for (int i = 0; i < this->getTargetDim(); ++i) {
        output[i] = _projection(_controller.getControl(index + i));
      }
      return true;
    }
    return false;
  }
};


/// @brief A data field whose values are managed by a controller
/// @tparam T data type
/// @tparam DESCRIPTOR lattice design
/// @param dataDim dimension of field
// Controller stores field values, which are (component-wise, point by point)
// transformed by the projection function and then returned.
// The controller is global, it is expected to have data for the entire domain.
// Hence, every process must have access to the ''same'' controller.
template <typename T, typename DESCRIPTOR>
class SuperLatticeSerialDataF : public SuperLatticeF<T,DESCRIPTOR> {
protected:
  std::shared_ptr<const GeometrySerializer<T,DESCRIPTOR::d>> _serializer;
public:
  SuperLatticeSerialDataF(SuperLattice<T,DESCRIPTOR>& superLattice,
    Controller<T>& controller, unsigned dataDim,
    std::shared_ptr<const GeometrySerializer<T,DESCRIPTOR::d>> serializer,
    std::function<T(T)> projection = [](T x){ return x; })
  : SuperLatticeF<T,DESCRIPTOR>(superLattice, dataDim),
    _serializer(serializer)
  {
    this->getName() = "SuperLatticeSerialDataF";

    const int maxC = superLattice.getLoadBalancer().size();
    this->_blockF.reserve(maxC);

    for (int iC = 0; iC < maxC; ++iC) {
      this->_blockF.emplace_back(
        new BlockLatticeSerialDataF<T,DESCRIPTOR>(
          superLattice.getBlock(iC),
          controller,
          *_serializer,
          superLattice.getLoadBalancer().glob(iC),
          dataDim,
          projection)
      );
    }
  }

  SuperLatticeSerialDataF(SuperLattice<T,DESCRIPTOR>& superLattice,
    Controller<T>& controller, unsigned dataDim,
    std::function<T(T)> projection = [](T x){ return x; })
  : SuperLatticeSerialDataF(superLattice, controller, dataDim,
      std::make_shared<const SimpleGeometrySerializer<T,DESCRIPTOR::d>>(superLattice),
      projection)
  { }
};





/////////////////////////// GLOBAL FUNCTOR ACCESS /////////////////////////////


/// Helper that gives global access to material numbers
// warning: this function breaks the parallelization concept. It is expensive
// and should be used with care
template<typename S, unsigned dim>
int getMaterialGlobally (SuperGeometry<S,dim>& sGeometry, LatticeR<dim+1> latticeR)
{
  int material = 0;
  if ( sGeometry.getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank() ) {
    material = sGeometry.get(latticeR);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(&material, 1, sGeometry.getLoadBalancer().rank(latticeR[0]));
#endif
  return material;
}

/// Helper that gives global access to the values of an indicator
// warning: this function breaks the parallelization concept. It is expensive
// and should be used with care
template <unsigned D, typename T>
bool evaluateSuperIndicatorFglobally (SuperIndicatorF<T,D>& f, const int input[])
{
  bool result = false;
  const auto loadBalancer = f.getSuperGeometry().getLoadBalancer();
  if ( loadBalancer.rank(input[0]) == singleton::mpi().getRank() ) {
    result = f(input);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(&result, 1, loadBalancer.rank(input[0]));
#endif
  return result;
}

/// Helper that gives global access to the values of a functor
// warning: this function breaks the parallelization concept. It is expensive
// and should be used with care
template <unsigned D, typename T, typename U=T>
bool evaluateSuperFglobally (SuperF<D,T,U>& f, U* output, const int input[])
{
  const auto loadBalancer = f.getSuperStructure().getLoadBalancer();
  if ( loadBalancer.rank(input[0]) == singleton::mpi().getRank() ) {
    f(output, input);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(output, f.getTargetDim(), loadBalancer.rank(input[0]));
#endif
  return true;
}





/////////////////////////// SERIALIZE FROM FIELDS /////////////////////////////

/// @brief Take values of a field and put them into a long vector
/// Idea: FIELD[cartesianCoordinates] = result[serialIndex]
/// -> we want to get the vector ''result''.
/// @tparam T floating point type
/// @tparam DESCRIPTOR lattice descriptor, must provide the field FIELD
/// @tparam C container type, where result is written
/// @tparam FIELD name of the field
/// @param sLattice super lattice
/// @param serializer implements the transformation from cartesian to serial data
/// @param indicator describes, at which positions FIELD is evaluated (at other positions, the default value of C remains)
/// @param controlDim length of result vector
/// @return vector with the field data
template <typename FIELD, typename T, typename DESCRIPTOR, typename C=std::vector<T>>
C serialDataFromField(SuperLattice<T,DESCRIPTOR>& sLattice,
  const GeometrySerializer<T,DESCRIPTOR::d>& serializer,
  SuperIndicatorF<T,DESCRIPTOR::d>& indicator,
  unsigned controlDim
  )
{
  constexpr unsigned dim = DESCRIPTOR::d;
  constexpr unsigned fieldDim = DESCRIPTOR::template size<FIELD>();
  auto result = util::ContainerCreator<C>::create(controlDim);
  const auto& cGeometry = sLattice.getCuboidDecomposition();
  const auto& loadBalancer = sLattice.getLoadBalancer();

  LatticeR<dim+1> latticeR;
  for (int iC=0; iC<cGeometry.size(); iC++) {
    latticeR[0] = iC;
    const int nX = cGeometry.get(iC).getNx();
    const int nY = cGeometry.get(iC).getNy();
    const int nZ = cGeometry.get(iC).getNz();
    for (int iX=0; iX<nX; iX++) {
      latticeR[1] = iX;
      for (int iY=0; iY<nY; iY++) {
        latticeR[2] = iY;
        for (int iZ=0; iZ<nZ; iZ++) {
          latticeR[3] = iZ;

          if (evaluateSuperIndicatorFglobally<dim,T>(indicator, latticeR.data())) {
            Vector<T,fieldDim> dataHelp;
            if (loadBalancer.rank(iC) == singleton::mpi().getRank()) {
              for (unsigned iDim=0; iDim<fieldDim; iDim++) {
                const auto cell = sLattice.get(latticeR);
                dataHelp[iDim]
                  = cell.template getFieldComponent<FIELD>(iDim);
              }
            }
#ifdef PARALLEL_MODE_MPI
            singleton::mpi().bCast(&dataHelp[0], fieldDim, loadBalancer.rank(iC));
#endif
            for (unsigned iDim=0; iDim<fieldDim; iDim++) {
              const auto index = serializer.getSerializedComponentIndex(latticeR, iDim, fieldDim);
              result[index] = dataHelp[iDim];
            }
          }
        }
      }
    }
  }
  return result;
}

// forward declaration
template<
  typename S,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C>
class OptiCaseDual;


/// @brief Take values of FIELD and put them into a long vector
/// Idea: FIELD[cartesianCoordinates] = result[serialIndex]
/// -> we want to get the vector ''result''.
/// @tparam T floating point type
/// @tparam C container type, where result is written
/// @tparam FIELD name of the field
/// @param optiCase dual opti case, provides the information on serialization etc.
/// @param solver the solver, whose FIELD data is taken
/// @return vector with the field data
template<
  typename T,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C=std::vector<T>>
C serialDataFromField(OptiCaseDual<T,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>& optiCase,
  std::shared_ptr<SOLVER<T,SolverMode::Reference>> solver)
{
  using descriptor = typename SOLVER<T,SolverMode::Reference>::AdjointLbSolver::DESCRIPTOR;
  return serialDataFromField<
    CONTROLLED_FIELD, T, descriptor, C>(
    solver->lattice(),
    *(optiCase._serializer),
    *(optiCase._controlIndicator),
    optiCase._dimCtrl
  );
}


/// @brief Get control values of some simulation in the context of adjoint optimization
/// Take values of FIELD, apply inverse projection and put them into a long vector
/// Idea: FIELD[cartesianCoordinates] = projection(result[serialIndex])
/// -> we want to get the vector ''result''.
/// @tparam T floating point type
/// @tparam C container type, where result is written
/// @tparam FIELD name of the field
/// @param optiCase dual opti case, provides the information on serialization etc.
/// @param solver the solver, whose FIELD data is taken
/// @return vector with the field data
template<
  typename T,
  template<typename,SolverMode> typename SOLVER,
  concepts::Field CONTROLLED_FIELD,
  template<typename...> typename PRIMAL_DYNAMICS,
  typename C=std::vector<T>>
C getControl(OptiCaseDual<T,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>& optiCase,
  std::shared_ptr<SOLVER<T,SolverMode::Reference>> solver)
{
  C result = serialDataFromField<T,SOLVER,CONTROLLED_FIELD,PRIMAL_DYNAMICS,C>(optiCase, solver);
  std::transform(result.begin(), result.end(), result.begin(), [&](auto r){
    return optiCase._projection->inverse(r);
  });
  return result;
}


}
}
#endif
