/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Timm Krüger, Shota Ito
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

#ifndef MEMBRANE_3D_H
#define MEMBRANE_3D_H

#include <vector>
#include "membraneMesh3D.h"
#include "utilities/pvecten.h"

namespace olb {

namespace membrane {

template<typename T>
struct MembraneParticle3D {
  /// Constructor.
  explicit MembraneParticle3D(std::shared_ptr<CMesh<T>> amesh) :
    mesh(amesh),
    pos(amesh->numNodes),
    posOld(amesh->numNodes),
    vel(amesh->numNodes),
    area(amesh->numFaces),
    weight(amesh->numNodes),
    normal(amesh->numFaces),
    forceStrain(amesh->numNodes),
    forceBending(amesh->numNodes),
    forceVolume(amesh->numNodes),
    forceSurface(amesh->numNodes),
    forceInteraction(amesh->numNodes),
    force(amesh->numNodes) {

    // Set initial individual properties.
    numNodes = mesh->numNodes;
    numFaces = mesh->numFaces;
    nodeOffset = 0;
    modShear = 0;
    modAlpha = 0;
    modBending = 0;
    modVolume = 0;
    modSurface = 0;
    modSphere = 0;

    // Set initial node positions and velocities.
    for(int n_i = 0; n_i < mesh->numNodes; ++n_i) {
      pos[n_i] = mesh->pos[n_i];
      posOld[n_i] = mesh->pos[n_i];
      vel[n_i].reset();
    }

    // Set initial centre and radius.
    setCentre(0, 0, 0);
    radius = 1;
    setRadius(radius, true);
  }

  /// Destructor
  ~MembraneParticle3D() {}

  /// Variables.
  int numNodes; // Number of nodes
  int numFaces; // Number of faces
  int nodeOffset; // Offset for local node to global node index transition
  std::shared_ptr<CMesh<T>> mesh;
  int meshType; // Mesh belonging to particle
  T radiusRef; // Reference radius
  T radius; // Current radius
  T volumeRef; // Reference volume
  T volume; // Current volume
  T surfaceRef; // Reference surface
  T surface; // Current surface
  T modShear; // Shear modulus
  T modAlpha; // Dilation modulus
  T modBending; // Bending modulus
  T modVolume; // Volume modulus
  T modSurface; // Surface modulus
  T modSphere; // Sphere modulus
  T ergStrain; // Strain and area energy
  T ergBending; // Bending energy
  T ergVolume; // Volume energy
  T ergSurface; // Surface energy
  std::vector<p_vec<> > pos; // Node positions
  std::vector<p_vec<> > posOld; // Old node positions
  std::vector<p_vec<> > vel; // Node velocities
  std::vector<T> area; // Face area
  std::vector<T> weight; // Node weight (equal to its Voronoi area)
  std::vector<p_vec<> > normal; // Unit face normal vector
  std::vector<p_vec<> > forceStrain; // Strain force
  std::vector<p_vec<> > forceBending; // Bending force
  std::vector<p_vec<> > forceVolume; // Volume force
  std::vector<p_vec<> > forceSurface; // Surface force
  std::vector<p_vec<> > forceInteraction; // Interaction force
  std::vector<p_vec<> > force; // Total force
  p_vec<> centre; // Current centre position
  p_vec<> centreOld; // Old centre position
  p_vec<> angularMomentum; // Angular momentum
  p_vec<> angularVelocity; // Angular velocity
  p_vec<> velocity; // Velocity
  p_vec<> Force; // Total force
  p_vec<> torque; // Torque

  /// Methods.
  T computeAngleEdges(int, int, int); // Compute angle defined by three nodes
  T computeAngleNormals(int, int); // Compute angle between two face normals
  p_vec<> computeNodeNormal(int); // Compute node normal vector
  p_vec<> computeFaceCentre(int); // Compute face centroid position
  p_ten2 computeInertiaTensorVolume(); // Compute inertia tensor (volume)
  p_ten2 computeInertiaTensorSurface(); // Compute inertia tensor (surface)
  p_vec<> computeAngularVelocity(); // Compute angular velocity
  p_vec<> computeOrientationVector(); // Compute orientation vector
  void computePosMinMax(p_vec<>, p_vec<>&, p_vec<>&); // Compute minimum and maximum node position with respect to a given plane
  T computeVolumeSlice(p_vec<>, p_vec<>); // Compute volume separated by slice
  void updateNodes(); // Update node positions
  void updateAreas(); // Update areas, normals and surface
  void updateVolume(); // Update volume
  void updateWeights(); // Update node weights
  void updateCentre(); // Update centre position
  void updateMomenta(); // Update momenta
  void setCentre(const p_vec<>&); // Set centre position (pass vector)
  void setCentre(T, T, T); // Set centre position
  void rotate(T, T, T, T); // Rotate particle
  void resize(T); // Resize particle without changing the equilibrium shape
  void setRadius(T, bool); // Set particle radius
  void forcesCompute(); // Call force computation methods
  void forcesPurge(); // Reset all force arrays
  void forcesVolume(); // Compute volume forces
  void forcesSurface(); // Compute surface forces
  void forcesSphere(); // Compute spherical shape recovery forces
  void forcesStrain(); // Compute strain forces
  void forcesBending(); // Compute bending forces
};

}

}

#endif
