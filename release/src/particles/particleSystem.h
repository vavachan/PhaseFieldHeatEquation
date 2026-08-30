/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen
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


#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H


namespace olb {

namespace particles {


template<typename T, typename PARTICLETYPE>
class ParticleSystem {
private:
  /// Main data storage
  using DATA = DynamicFieldGroupsD<T, typename PARTICLETYPE::fields_t>;
  Container<T,PARTICLETYPE,DATA> _fieldGroupsContainer;
  /// Danamics Vector
  std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>> _dynamicsVector;
  //Additional storage for multiple complex data types linked to particle (extendable by user...)
  using associatedTypes = meta::list<
    std::vector<std::unique_ptr<SmoothIndicatorF2D<T,T,true>>>,
    std::vector<std::unique_ptr<SmoothIndicatorF3D<T,T,true>>>
  >;
  typename associatedTypes::template decompose_into<std::tuple> _associatedData; //Storage for complex types
  //Serial size of DATA
  const std::size_t _serialSize;

  /// Iterator (Should be reworked for C++20)
  struct Iterator {
    //Iterator tags
    using iterator_category = std::random_access_iterator_tag;
    using pType             = Particle<T,PARTICLETYPE>;
    //Constructor
    Iterator(ParticleSystem<T,PARTICLETYPE>* pSystemPtr, std::size_t iParticle)
      : _pSystemPtr(pSystemPtr), _iParticle(iParticle) {}
    //Reference operator
    pType operator*() { return _pSystemPtr->get(_iParticle); } //is called for "auto p : particleSystem"
    // Prefix increment
    Iterator& operator++() { _iParticle++; return *this; }
    // Postfix increment
    Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
    //Comparison operators
    friend bool operator== (const Iterator& a, const Iterator& b) { return a._iParticle == b._iParticle; };
    friend bool operator!= (const Iterator& a, const Iterator& b) { return a._iParticle != b._iParticle; };
  private:
    ParticleSystem<T,PARTICLETYPE>* _pSystemPtr;
    std::size_t _iParticle;
  };

public:
  static_assert(PARTICLETYPE::d == 2 || PARTICLETYPE::d == 3, "Only D=2 and D=3 are supported");

  /// Default constructor
  ParticleSystem();

  /// Constructor with initial particle size
  ParticleSystem( std::size_t count );

  /// Constructor with initial particle size and associated Data
  ParticleSystem( std::size_t count,
                  typename associatedTypes::template decompose_into<std::tuple> associatedData );

  //Iterator instantiation
  Iterator begin()  { return Iterator(this,0); }
  Iterator end()    { return Iterator(this,size()); }

  /// Size of ParticleSystem
  template<typename PCONDITION=conditions::all_particles>
  constexpr std::size_t size();

  /// Size of ParticleSystem (enables passing of globiC, if used in SuperParticleSystem)
  template<typename PCONDITION=conditions::all_particles>
  constexpr std::size_t size(int globiC);

  /// Get specific dynamics (optional boundsCheck via .at())
  template<bool boundsCheck=false>
  dynamics::ParticleDynamics<T,PARTICLETYPE>* getDynamics(unsigned iDyn=0);

  /// Add dynamics (from shared_pointer of dynamics)
  void addDynamics(std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>& dynamicsSPtr );

  /// Define dynamics (factory method)
  template <typename DYNAMICS, typename ...Args>
  void defineDynamics(Args&& ...args);

  /// Extend particle system by one particle
  void extend();

  /// Swap particles by index
  void swapParticles(std::size_t iP, std::size_t jP);

  /// Upate/process particle quantities on a subset (mostly local equivalent to collide)
  void process( std::size_t p0, std::size_t p1, T timeStepSize, unsigned iDyn=0 );

  /// Upate/process particle quantities on all particles
  void process(T timeStepSize, unsigned iDyn=0);

  /// Update particles (equivalent to stream())
  void update();

  /// Expose container
  auto& get();

  /// Expose dynamics vector
  std::vector<std::shared_ptr<dynamics::ParticleDynamics<T,PARTICLETYPE>>>& getDynamicsVector();

  /// Create and return particle (analogous to cell TODO: make same baseclass) (optional boundsCheck)
  template<bool boundsCheck=false>
  Particle<T,PARTICLETYPE> get(std::size_t iParticle);

  /// Create and return particle with operator[]
  Particle<T,PARTICLETYPE> operator[](std::size_t iParticle);

  /// Get whole Field by GROUP (or base of GROUP) and FIELD
  template <typename GROUP, typename FIELD>
  auto& getFieldD();

  /// Get FieldPointer by GROUP (or base of GROUP) and FIELD for specific iParticle
  template <typename GROUP, typename FIELD>
  auto getFieldPointer( std::size_t iParticle );

  /// Get associated data by specifying data type
  template<typename TYPE>
  auto& getAssociatedData();

  /// Get serial size of dynamic fields group
  std::size_t getSerialSize() const;

  /// Print relevant infos
  void print();

  /// Check for erroes
  void checkForErrors();
};


} //namespace particles

} //namespace olb


#endif
