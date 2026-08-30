/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito, Timm Krüger
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

#ifndef SUPER_MEMBRANE_3D_H
#define SUPER_MEMBRANE_3D_H

#include "io/ostreamManager.h"

namespace olb {

namespace membrane {

template<typename T, typename DESCRIPTOR>
class MembraneParticleSystem3D {
private:
  std::vector<std::shared_ptr<MembraneParticle3D<T>>> _particle;
  // serialized data container for ibm communication
  std::shared_ptr<SuperD<T,DESCRIPTOR>> _data;
  std::vector<std::shared_ptr<CMesh<T>>> _meshes;
  std::vector<unsigned> _meshNumParticles;
  std::vector<unsigned> _meshForParticle;
  unsigned _numTotalNodes{0};
  unsigned _numTotalFaces{0};

  mutable OstreamManager clout;

  // Mesh initialization data.
  std::vector<std::string> _meshFilename;
  std::vector<T> _meshRadius; // Radius of particles in mesh
  std::vector<T> _meshKV; // Volume modulus of particles in mesh
  std::vector<T> _meshKSph; // Sphere modulus of particles in mesh
  std::vector<T> _meshKA; // Surface modulus of particles in mesh
  std::vector<T> _meshKAlpha; // Dilation modulus of particles in mesh
  std::vector<T> _meshKS; // Shear modulus of particles in mesh
  std::vector<T> _meshKB; // Bending modulus of particles in mesh

  // Particle initialisation data.
  struct initParticleStruct {
    int nodeOffset; // Node offset
    int meshType; // Mesh index
    T radius_pos; // Radius of particle
    T angle; // Rotation angle
    p_vec<> centre; // Centre position
    p_vec<> axis; // Rotation axis
  };
  std::vector<initParticleStruct> _initParticles;

  // Type of normal vector correction
  int _normalMethod{0};

  bool _checkInitialParticleOverlap = true;

  void readMeshParameters(std::string filename);
  void readPositionParameters(std::string filename);
  /// Checks if there are overlapping particles and if yes abort.
  /// Can be disabled via flags.
  void checkParticleParticleOverlap();
  /// Check if the particle moving speed exceeds 0.5
  void checkParticles();
  /// The particles are positioned according to data given in parameter file.
  void positionParticlesManually();
  /// The particle volume fraction is computed.
  void computeVolumeFraction();
  /// Preparation for the IBM-coupling
  void downloadData();
  void uploadData();
public:
  MembraneParticleSystem3D();
  void initializeMeshesFromXML(std::string filename);
  void initializeParticlesFromXML(std::string filename);
  void initializeData(LoadBalancer<T>& loadBalancer);
  /// Trigger internal forces for all membranes
  void computeForcesInternal();
  /// All membrane force contributions are combined for IBM spreading
  void combineForces();
  /// The particle configurations are updated.
  void updateParticles();
  /// One discrete time step solving FEM
  void updateForcesFromDisplacements();
  /// Update membrane momenta, only required for IO
  void updateMomenta();
  /// Access to membrane particle via index
  std::shared_ptr<MembraneParticle3D<T>> get(int c_i);
  /// Access to particle data via index required for ibm
  SuperD<T,DESCRIPTOR>& getData();
  unsigned getNumberMeshes();
  unsigned getNumberParticles();
  unsigned getTotalNumberNodes();
  unsigned getTotalNumberFaces();
};

}

}

#endif
