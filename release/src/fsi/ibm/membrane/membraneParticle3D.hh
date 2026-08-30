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

#ifndef MEMBRANE_3D_HH
#define MEMBRANE_3D_HH

#include <iomanip>
#include <fstream>
#include "membraneParticle3D.h"

namespace olb {

namespace membrane {

using namespace std;

/*********************************************************************************
************* COMPUTE EDGE ANGLE *************************************************
*********************************************************************************/

/// The angle between two edges within the same face is computed. The edges are
/// specified by three nodes where the angle is taken at the second node. To
/// increase numerical stability, extremely small angles are neglected.
template<typename T>
inline T MembraneParticle3D<T>::computeAngleEdges(int n_i, int n_j, int n_k) {
  const p_vec<> edge1 = pos[n_k] - pos[n_j];
  const p_vec<> edge2 = pos[n_i] - pos[n_j];
  const T cosine = edge1 * edge2 / (edge1.length() * edge2.length());
  T angle;
  if(cosine >= 0.999999) {
    angle = 0.;
  }
  else if(cosine <= -0.999999) {
    angle = M_PI;
  }
  else {
    angle = acos(cosine);
  }

  return angle;
}

/*********************************************************************************
************* COMPUTE NORMAL ANGLE ***********************************************
*********************************************************************************/

/// The angle between two face normals is computed. To increase numerical
/// stability, extremely small angles are neglected.
template<typename T>
T MembraneParticle3D<T>::computeAngleNormals(int f_i, int f_j) {
  const T cosine = normal[f_i] * normal[f_j];
  T angle;
  if(cosine >= 0.999999) {
    angle = 0;
  }
  else if(cosine <= -0.999999) {
    angle = M_PI;
  }
  else {
    angle = acos(cosine);
  }

  return angle;
}

/*********************************************************************************
************* COMPUTE NODE NORMAL VECTOR *****************************************
*********************************************************************************/

/// Compute node normal vector as average of neighbouring face normals.
/// NOTE: More accurate implementations are possible.
template<typename T>
p_vec<> MembraneParticle3D<T>::computeNodeNormal(int n_i) {
  p_vec<> normalVec(0, 0, 0);
  for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
   const int f_i = mesh->nodeToFace[n][n_i];
   if(f_i >= 0) {
     normalVec += normal[f_i];
   }
 }
 normalVec.norm_to(1);

 return normalVec;
}

/*********************************************************************************
************* COMPUTE FACE CENTRE ************************************************
*********************************************************************************/

/// The centroid position of a face is computed.
template<typename T>
p_vec<> MembraneParticle3D<T>::computeFaceCentre(int f_i) {
  const int n_j = mesh->faceToNode[0][f_i];
  const int n_k = mesh->faceToNode[1][f_i];
  const int n_l = mesh->faceToNode[2][f_i];

  return ((pos[n_j] + pos[n_k] + pos[n_l]) / 3);
}

/*********************************************************************************
************* COMPUTE INERTIA TENSOR (VOLUME) ************************************
*********************************************************************************/

/// The inertia tensor of the membrane is computed based on its volume. It is
/// assumed that the density is constant throughout the volume. This approach is
/// usually employed for the computation of the equivalent ellipsoidal axes of the
/// membrane.
template<typename T>
p_ten2 MembraneParticle3D<T>::computeInertiaTensorVolume() {
  p_ten2 inertiaTensor;
  inertiaTensor.reset();
  /// Compute independent components.
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const p_vec<> pv = computeFaceCentre(f_i) - centre;
    const T sp = pv * normal[f_i];
    inertiaTensor.xx += (area[f_i] / 5 * ((pv * pv * sp) - util::pow(pv.x,2.) * sp));
    inertiaTensor.yy += (area[f_i] / 5 * ((pv * pv * sp) - util::pow(pv.y,2.) * sp));
    inertiaTensor.zz += (area[f_i] / 5 * ((pv * pv * sp) - util::pow(pv.z,2.) * sp));
    inertiaTensor.xy += (area[f_i] / 5 * (-pv.x * pv.y * sp));
    inertiaTensor.xz += (area[f_i] / 5 * (-pv.x * pv.z * sp));
    inertiaTensor.yz += (area[f_i] / 5 * (-pv.y * pv.z * sp));
  }
  /// Symmetrise inertia tensor.
  inertiaTensor.yx = inertiaTensor.xy;
  inertiaTensor.zx = inertiaTensor.xz;
  inertiaTensor.zy = inertiaTensor.yz;

  return inertiaTensor;
}

/*********************************************************************************
************* COMPUTE INERTIA TENSOR (SURFACE) ***********************************
*********************************************************************************/

/// The inertia tensor of the membrane is computed based on its surface. It is
/// assumed that the mass is concentrated on the membrane. This approach is
/// important for the rotational behaviour.
template<typename T>
p_ten2 MembraneParticle3D<T>::computeInertiaTensorSurface() {
  p_ten2 inertiaTensor;
  inertiaTensor.reset();
  /// Compute independent components.
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    const p_vec<> pv = pos[n_i] - centre;
    inertiaTensor.xx += (weight[n_i] * (util::pow(pv.y,2.) + util::pow(pv.z,2.)));
    inertiaTensor.yy += (weight[n_i] * (util::pow(pv.x,2.) + util::pow(pv.z,2.)));
    inertiaTensor.zz += (weight[n_i] * (util::pow(pv.x,2.) + util::pow(pv.y,2.)));
    inertiaTensor.xy -= (weight[n_i] * pv.x * pv.y);
    inertiaTensor.xz -= (weight[n_i] * pv.x * pv.z);
    inertiaTensor.yz -= (weight[n_i] * pv.y * pv.z);
  }
  /// Symmetrise inertia tensor.
  inertiaTensor.yx = inertiaTensor.xy;
  inertiaTensor.zx = inertiaTensor.xz;
  inertiaTensor.zy = inertiaTensor.yz;
  /// Normalise tensor to account for total surface area (= total mass).
  inertiaTensor *= surface;

  return inertiaTensor;
}

/*********************************************************************************
************* COMPUTE ANGULAR VELOCITY *******************************************
*********************************************************************************/

/// The angular velocity of the membrane is computed.
/// NOTE: The important assumption is that the angular velocity is defined by
/// L = I * omega where L is the angular momentum and I is the inertia tensor of
/// the thin membrane shell. This approach is only correct for a rigid object.
/// For deformable objects, another method is required.
template<typename T>
p_vec<> MembraneParticle3D<T>::computeAngularVelocity() {
  p_ten2 inertiaTensor = computeInertiaTensorSurface();
  const p_ten2 inertiaTensorInv = inertiaTensor.invert();

  return (inertiaTensorInv * angularMomentum);
}

/*********************************************************************************
************* COMPUTE ORIENTATION VECTOR OF MEMBRANE *****************************
*********************************************************************************/

/// The orientation vector of the membrane is computed: It is an approximate
/// measure for distinguishing the top and bottom of a particle. This additional
/// information is not accessible from the inertia tensor. The orientation vector
/// is defined to point always in the direction of the top of the particle. The
/// top of the particle is defined at the start of the simulation by assigning a
/// positive index to the top nodes and a negative index to the bottom nodes.
template<typename T>
p_vec<> MembraneParticle3D<T>::computeOrientationVector() {
  p_vec<> centrePos(0, 0, 0);
  p_vec<> centreNeg(0, 0, 0);
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    if(mesh->nodeOrientation[n_i] > 0) {
      centrePos += pos[n_i];
    }
    else if(mesh->nodeOrientation[n_i] < 0) {
      centreNeg += pos[n_i];
    }
  }
  centrePos /= mesh->numPosNodeOrientation;
  centreNeg /= mesh->numNegNodeOrientation;
  const p_vec<> orientation = centrePos - centreNeg;

  return orientation / orientation.length();
}

/*********************************************************************************
************* COMPUTE MINIMUM AND MAXIMUM POSITIONS ******************************
*********************************************************************************/

/// For an arbitrarily given plane (defined by its normal vector), the positions
/// of the two nodes with extremal distances to the plane are computed. If the
/// plane cuts the membrane, the two nodes have the maximum distance in either
/// direction of the plane. If the plane does not cut the membrane, the closest
/// and the farthest positions are computed. This information can be used to
/// compute the extension of the membrane along any axis.
template<typename T>
void MembraneParticle3D<T>::computePosMinMax(p_vec<> planeNormal, p_vec<>& posMin, p_vec<>& posMax) {
  T distMax = 0;
  T distMin = 0;
  planeNormal.norm_to(1);
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    const T dist = (pos[n_i] - centre) * planeNormal;
    if(dist < distMin) {
      distMin = dist;
      posMin = pos[n_i];
    }
    if(dist > distMax) {
      distMax = dist;
      posMax = pos[n_i];
    }
  }

  return;
}

/*********************************************************************************
************* COMPUTE VOLUME SEPARATED BY SLICE **********************************
*********************************************************************************/

/// Any plane separates the membrane in two parts (where one part can have zero
/// size). The plane is defined by an origin and its normal vector. The volume of
/// the part of the membrane in normal direction of the plane is computed. This
/// function is useful for finding the volume of a membrane between two parallel
/// planes. This allows to compute the local volume fraction in a planar geometry.
template<typename T>
T MembraneParticle3D<T>::computeVolumeSlice(p_vec<> sliceOrigin, p_vec<> sliceNormal) {
  T volumeSlice = 0.;
  sliceNormal.norm_to(1.);
  /// Distinguish cases.
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_1 = mesh->faceToNode[0][f_i];
    const int n_2 = mesh->faceToNode[1][f_i];
    const int n_3 = mesh->faceToNode[2][f_i];
    const T dist1 = (pos[n_1] - sliceOrigin) * sliceNormal;
    const T dist2 = (pos[n_2] - sliceOrigin) * sliceNormal;
    const T dist3 = (pos[n_3] - sliceOrigin) * sliceNormal;
    /// If all three nodes of the face are in normal direction of the plane,
    /// the face volume is completely taken.
    if(dist1 > 0 && dist2 > 0 && dist3 > 0) {
      volumeSlice += ((computeFaceCentre(f_i) - sliceOrigin) * normal[f_i] * area[f_i] / 3);
    }
    /// If only one node is in normal direction of the plane,
    /// the part of the face volume in normal direction is taken.
    else if(dist1 > 0 && dist2 <= 0 && dist3 <= 0) {
      const p_vec<> intsec12 = pos[n_1] + (pos[n_2] - pos[n_1]) * dist1 / (dist1 - dist2);
      const p_vec<> intsec13 = pos[n_1] + (pos[n_3] - pos[n_1]) * dist1 / (dist1 - dist3);
      const T faceArea = ((intsec12 - pos[n_1]) % (intsec13 - pos[n_1])).length() / 2;
      volumeSlice += (pos[n_1] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
    else if(dist1 <= 0 && dist2 > 0 && dist3 <= 0) {
      const p_vec<> intsec21 = pos[n_2] + (pos[n_1] - pos[n_2]) * dist2 / (dist1 - dist2);
      const p_vec<> intsec23 = pos[n_2] + (pos[n_3] - pos[n_2]) * dist2 / (dist2 - dist3);
      const T faceArea = ((intsec21 - pos[n_2]) % (intsec23 - pos[n_2])).length() / 2;
      volumeSlice += (pos[n_2] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
    else if(dist1 <= 0 && dist2 <= 0 && dist3 > 0) {
      const p_vec<> intsec31 = pos[n_3] + (pos[n_1] - pos[n_3]) * dist3 / (dist1 - dist3);
      const p_vec<> intsec32 = pos[n_3] + (pos[n_2] - pos[n_3]) * dist3 / (dist2 - dist3);
      const T faceArea = ((intsec31 - pos[n_3]) % (intsec32 - pos[n_3])).length() / 2;
      volumeSlice += (pos[n_3] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
    /// If two nodes are in normal direction of the plane, the part of the face volume
    /// in normal direction is taken as complement to the part in negative normal direction.
    else if(dist1 > 0 && dist2 > 0 && dist3 <= 0) {
      const p_vec<> intsec31 = pos[n_3] + (pos[n_1] - pos[n_3]) * dist3 / (dist1 - dist3);
      const p_vec<> intsec32 = pos[n_3] + (pos[n_2] - pos[n_3]) * dist3 / (dist2 - dist3);
      const T faceArea = area[f_i] - ((intsec31 - pos[n_3]) % (intsec32 - pos[n_3])).length() / 2;
      volumeSlice += (pos[n_3] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
    else if(dist1 > 0 && dist2 <= 0 && dist3 > 0) {
      const p_vec<> intsec21 = pos[n_2] + (pos[n_1] - pos[n_2]) * dist2 / (dist1 - dist2);
      const p_vec<> intsec23 = pos[n_2] + (pos[n_3] - pos[n_2]) * dist2 / (dist2 - dist3);
      const T faceArea = area[f_i] - ((intsec21 - pos[n_2]) % (intsec23 - pos[n_2])).length() / 2;
      volumeSlice += (pos[n_2] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
    else if(dist1 <= 0 && dist2 > 0 && dist3 > 0) {
      const p_vec<> intsec12 = pos[n_1] + (pos[n_2] - pos[n_1]) * dist1 / (dist1 - dist2);
      const p_vec<> intsec13 = pos[n_1] + (pos[n_3] - pos[n_1]) * dist1 / (dist1 - dist3);
      const T faceArea = area[f_i] - ((intsec12 - pos[n_1]) % (intsec13 - pos[n_1])).length() / 2;
      volumeSlice += (pos[n_1] - sliceOrigin) * normal[f_i] * faceArea / 3;
    }
  }

  return volumeSlice;
}

/*********************************************************************************
************* UPDATE NODE POSITIONS **********************************************
*********************************************************************************/

/// Node positions are updated. The old positions are stored for data analysis
/// (not required for the simulation itself). The periodicity of the simulation
/// is not taken into account here.
/// NOTE: The current implementation is forward Euler, but other schemes could
/// be used.
template<typename T>
void MembraneParticle3D<T>::updateNodes() {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    posOld[n_i] = pos[n_i];
    pos[n_i] += vel[n_i];
  }

  return;
}

/*********************************************************************************
************* UPDATE AREAS, NORMALS AND SURFACE **********************************
*********************************************************************************/

/// Face areas, normal vectors and membrane surface are updated.
template<typename T>
void MembraneParticle3D<T>::updateAreas() {
  surface = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_j = mesh->faceToNode[0][f_i];
    const int n_k = mesh->faceToNode[1][f_i];
    const int n_l = mesh->faceToNode[2][f_i];
    const p_vec<> norm = ((pos[n_j] - pos[n_k]) % (pos[n_l] - pos[n_k]));
    const T normLen = norm.length();
    area[f_i] = normLen / 2;
    surface += area[f_i];
    normal[f_i] = norm / normLen;
  }

  return;
}

/*********************************************************************************
************* UPDATE MEMBRANE VOLUME *********************************************
*********************************************************************************/

/// The membrane volume is updated.
/// NOTE: The volume-averaged centre position could be updated here as well.
template<typename T>
void MembraneParticle3D<T>::updateVolume() {
  volume = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const p_vec<> pos = computeFaceCentre(f_i);
//     const p_ten2 posTensor(SQ2(pos.x), pos.x * pos.y, pos.x * pos.z, pos.y * pos.x, SQ2(pos.y), pos.y * pos.z, pos.z * pos.x, pos.z * pos.y, SQ2(pos.z));
    volume += ((normal[f_i] * computeFaceCentre(f_i)) * area[f_i]) / 3;
  }

  return;
}

/*********************************************************************************
************* UPDATE NODE WEIGHTS ************************************************
*********************************************************************************/

/// The node weights are updated. They are proportional to the Voronoi area of
/// the node and normalised: the sum of all weights is unity. To find the weight
/// of a node, all neighbouring face areas have to be considered.
template<typename T>
void MembraneParticle3D<T>::updateWeights() {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    weight[n_i] = 0;
  }
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_j = mesh->faceToNode[0][f_i];
    const int n_k = mesh->faceToNode[1][f_i];
    const int n_l = mesh->faceToNode[2][f_i];
    weight[n_j] += area[f_i] / (3 * surface);
    weight[n_k] += area[f_i] / (3 * surface);
    weight[n_l] += area[f_i] / (3 * surface);
  }

  return;
}

/*********************************************************************************
************* UPDATE CENTRE POSITION *********************************************
*********************************************************************************/

/// The particle centre is updated. Based on the advection of the individual
/// nodes, the centre position also changes.
template<typename T>
void MembraneParticle3D<T>::updateCentre() {
  centreOld = centre;
  centre.reset();
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    centre += (pos[n_i] * weight[n_i]);
  }
}

/*********************************************************************************
************* UPDATE MOMENTA *****************************************************
*********************************************************************************/

/// The membrane momenta are updated:
/// - angular momentum
/// - linear and angular velocity
/// - total membrane force
/// - total mambrane torque
/// The surface takes the role of the mass.
/// NOTE: Torque requires old positions (since they correspond to the old forces).
template<typename T>
void MembraneParticle3D<T>::updateMomenta() {
  velocity.reset();
  angularMomentum.reset();
  Force.reset();
  torque.reset();
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    velocity += (vel[n_i] * weight[n_i]);
    angularMomentum += ((pos[n_i] - centre) % vel[n_i] * weight[n_i]);
    Force += force[n_i];
    torque += (posOld[n_i] - centreOld) % force[n_i];
  }
  angularMomentum *= surface;
  angularVelocity = computeAngularVelocity();

  return;
}

/*********************************************************************************
************* SET CENTRE POSITION ************************************************
*********************************************************************************/

/// The position of the particle centre is set to a given value. The node
/// positions are updated accordingly. Periodicity is not checked here since this
/// function is only intended to be used at the beginning of the simulation for
/// setup purposes.
/// There are two versions of this function: p_vec and Ts as parameters.
template<typename T>
void MembraneParticle3D<T>::setCentre(const p_vec<>& position) {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    pos[n_i] += (position - centre);
  }
  centre = position;

  return;
}

template<typename T>
void MembraneParticle3D<T>::setCentre(T X, T Y, T Z) {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    pos[n_i].x += (X - centre.x);
    pos[n_i].y += (Y - centre.y);
    pos[n_i].z += (Z - centre.z);
  }
  centre.set_to(X, Y, Z);

  return;
}

/*********************************************************************************
************* ROTATE MEMBRANE ****************************************************
*********************************************************************************/

/// The orientation of the particle is rotated about the axis (X, Y, Z) by angle
/// in degrees. A rotation matrix is computed from the rotation axis and the angle.
/// The particle centre and other properties are not changed.
template<typename T>
void MembraneParticle3D<T>::rotate(T angle, T X, T Y, T Z) {
  p_vec<> axis(X, Y, Z);
  axis.norm_to(1);
  angle *= (M_PI / 180);
  p_ten2 rotMatrix;
  rotMatrix.set_to(
    cos(angle) + axis.x * axis.x * (1 - cos(angle)),
    axis.x * axis.y * (1 - cos(angle)) - axis.z * sin(angle),
    axis.x * axis.z * (1 - cos(angle)) + axis.y * sin(angle),
    axis.y * axis.x * (1 - cos(angle)) + axis.z * sin(angle),
    cos(angle) + axis.y * axis.y * (1 - cos(angle)),
    axis.y * axis.z * (1 - cos(angle)) - axis.x * sin(angle),
    axis.z * axis.x * (1 - cos(angle)) - axis.y * sin(angle),
    axis.z * axis.y * (1 - cos(angle)) + axis.x * sin(angle),
    cos(angle) + axis.z * axis.z * (1 - cos(angle))
  );
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    pos[n_i] = rotMatrix * (pos[n_i] - centre) + centre;
  }

  return;
}

/*********************************************************************************
************* RESIZE SHAPE *******************************************************
*********************************************************************************/

/// The particle is radially resized by a factor (unity corresponds to no change).
/// The centre and the equilibrium properties of the membrane are not modified.
template<typename T>
void MembraneParticle3D<T>::resize(T a) {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    pos[n_i] = (pos[n_i] - centre) * a + centre;
  }

  return;
}

/*********************************************************************************
************* SET MEMBRANE RADIUS ************************************************
*********************************************************************************/

/// The particle radius is set to a new value. Consequently, the equilibrium
/// surface and volume are updated. The user has to decide whether the reference
/// radius has to be changed as well or not. This is particularly important for
/// the particle growth at the beginning of a simulation.
template<typename T>
void MembraneParticle3D<T>::setRadius(T radiusNew, bool switchRadius) {
  const T radiusTemp = radius;
  radius = radiusNew;
  if(switchRadius) {
    radiusRef = radiusNew;
  }
  surfaceRef = mesh->surface * util::pow(radiusNew,2.);
  volumeRef = mesh->volume * util::pow(radiusNew,3.);
  resize(radiusNew / radiusTemp);
  updateAreas();
  updateWeights();
  updateVolume();

  return;
}

/*********************************************************************************
************* COMPUTE FORCES *****************************************************
*********************************************************************************/

/// The particle forces are computed. This method calls all other force
/// computation methods. The forces must always be purged at the beginning in
/// order to reset the force arrays.
template<typename T>
void MembraneParticle3D<T>::forcesCompute() {
  forcesPurge();
  forcesVolume();
  forcesSurface();
  forcesSphere();
  forcesStrain();
  forcesBending();

  return;
}

/*********************************************************************************
************* RESET FORCES *******************************************************
*********************************************************************************/

/// All forces and energies are set to zero. This function must be called at each
/// time step before new forces and energies are computed.
template<typename T>
void MembraneParticle3D<T>::forcesPurge() {
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    force[n_i].reset();
    forceStrain[n_i].reset();
    forceBending[n_i].reset();
    forceVolume[n_i].reset();
    forceSurface[n_i].reset();
    forceInteraction[n_i].reset();
  }
  ergStrain = 0;
  ergBending = 0;
  ergVolume = 0;
  ergSurface = 0;

  return;
}

/*********************************************************************************
************* VOLUME FORCE *******************************************************
*********************************************************************************/

/// The volume force is computed. It is proportional to the relative volume
/// deviation. If the volume modulus is non-positive, the computation is skipped.
template<typename T>
void MembraneParticle3D<T>::forcesVolume() {
  if(modVolume < 1.e-15) {
    return;
  }
  const T strength = modVolume * (volume / volumeRef - 1) / 6;
  for(int f_j = 0; f_j < numFaces; ++f_j) {
    const int n_i = mesh->faceToNode[0][f_j];
    const int n_j = mesh->faceToNode[1][f_j];
    const int n_k = mesh->faceToNode[2][f_j];
    forceVolume[n_i] += (((pos[n_j] - centre) % (pos[n_k] - centre)) * strength);
    forceVolume[n_j] += (((pos[n_k] - centre) % (pos[n_i] - centre)) * strength);
    forceVolume[n_k] += (((pos[n_i] - centre) % (pos[n_j] - centre)) * strength);
  }
  ergVolume = modVolume * util::pow(volume - volumeRef,2.) / (2 * volumeRef);

  return;
}

/*********************************************************************************
************* SURFACE FORCE ******************************************************
*********************************************************************************/

/// The surface force is computed. It is proportional to the relative surface
/// deviation. If the surface modulus is non-positive, the computation is skipped.
template<typename T>
void MembraneParticle3D<T>::forcesSurface() {
  if(modSurface < 1.e-15) {
    return;
  }
  const T strength = modSurface * (surface / surfaceRef - 1) / 2;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_i = mesh->faceToNode[0][f_i];
    const int n_j = mesh->faceToNode[1][f_i];
    const int n_k = mesh->faceToNode[2][f_i];
    forceSurface[n_i] += ((normal[f_i] % (pos[n_k] - pos[n_j])) * strength);
    forceSurface[n_j] += ((normal[f_i] % (pos[n_i] - pos[n_k])) * strength);
    forceSurface[n_k] += ((normal[f_i] % (pos[n_j] - pos[n_i])) * strength);
  }
  ergSurface = modSurface * util::pow(surface - surfaceRef,2.) / (2 * surfaceRef);

  return;
}

/*********************************************************************************
************* SPHERICAL SHAPE RECOVERY FORCE *************************************
*********************************************************************************/

/// The spherical shape recovery force is computed. If the sphere modulus is
/// non-positive, the computation is skipped.
/// NOTE: This is an experimental function meant to be used for rigid spheres
/// and will likely change in the future or be replaced.
template<typename T>
void MembraneParticle3D<T>::forcesSphere() {
  if(modSphere < 1.e-15) {
    return;
  }
  /// Compute elastic and friction force contributions.
  p_vec<> forceTot(0., 0., 0.);
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    const p_vec<> posRel = pos[n_i] - centre;
    const T posLen = posRel.length();
    const T strength = modSphere * (posLen - radius) * surface;
    const p_vec<> direction = posRel / posLen;
    const p_vec<> forceElast = -direction * strength;
    const p_vec<> velRad = direction * ((vel[n_i] - velocity) * direction);
    const p_vec<> forceFric = -velRad * modSphere * surface * 3.;
    forceVolume[n_i] += (forceElast + forceFric) * weight[n_i];
    forceTot += forceVolume[n_i];
  }
  /// Correct for momentum conservation.
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    forceVolume[n_i] -= (forceTot) * weight[n_i];
  }

  return;
}

/*********************************************************************************
************* STRAIN FORCE *******************************************************
*********************************************************************************/

/// The strain force is computed via a simplified FEM approach. If the strain and
/// dilation moduli are non-positive, the computation is skipped.
/// NOTE: Selection of different strain models could be implemented.
template<typename T>
void MembraneParticle3D<T>::forcesStrain() {
  if(modShear < 1.e-15 && modAlpha < 1.e-15) {
    return;
  }
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_0 = mesh->faceToNode[0][f_i];
    const int n_1 = mesh->faceToNode[1][f_i];
    const int n_2 = mesh->faceToNode[2][f_i];
    /// Reference geometry.
    /// the angles phi and phiEq shall be at node n_1.
    /// The aligned edges l1 and l1Eq shall be between nodes n_1 and n_2.
    /// The edges l2 and l2Eq shall be between nodes n_1 and n_0.
    const T l1Eq = mesh->l1Eq[f_i] * radius;
    const T l2Eq = mesh->l2Eq[f_i] * radius;
    const T cosEq = mesh->cosEq[f_i];
    const T sinInvEq = mesh->sinInvEq[f_i];
    /// Deformed geometry.
    const T l1 = (pos[n_2] - pos[n_1]).length();
    const T l2 = (pos[n_0] - pos[n_1]).length();
    const T phi = computeAngleEdges(n_0, n_1, n_2);
    const T cosine = cos(phi);
    const T sine = sin(phi);
    /// Shape functions.
    /// The area in the denominator has been cancelled for numerical efficiency.
    /// The shape functions are stored in memory instead of recomputing them.
    const T b0 = mesh->b0[f_i] * radius;
    const T a1 = mesh->a1[f_i] * radius;
    const T b1 = mesh->b1[f_i] * radius;
    /// Deformation gradient tensor.
    const T Dxx = l1 / l1Eq;
    const T Dyy = l2 / l2Eq * sine * sinInvEq;
    const T Dxy = (l2 / l2Eq * cosine - l1 / l1Eq * cosEq) * sinInvEq;
    const T Gxx = util::pow(Dxx,2.);
    const T Gyy = util::pow(Dxy,2.) + util::pow(Dyy,2.);
    const T Gxy = Dxx * Dxy;
    /// Principal stresses.
    const T I1 = Gxx + Gyy - 2;
    const T I2 = Gxx * Gyy - util::pow(Gxy,2.) - 1;
    /// Neo-Hookean model.
    // const T w = modShear / 6. * (I1 - 1. + 1. / (I2 + 1.));
    // const T dw_dI1 = modShear / 6.;
    // const T dw_dI2 = -modShear / 6. / SQ2(I2 + 1.);
    /// Skalak model.
    const T w = modShear / 12 * (util::pow(I1,2.) + 2 * I1 - 2 * I2) + modAlpha / 12 * util::pow(I2,2.);
    const T dw_dI1 = modShear / 6 * (I1 + 1);
    const T dw_dI2 = -modShear / 6 + modAlpha / 6 * I2;
    /// Derivatives of strain invariants.
    const T dI1_dGxx = 1;
    const T dI1_dGyy = 1;
    const T dI2_dGxx = Gyy;
    const T dI2_dGyy = Gxx;
    const T dI2_dGxy = -2 * Gxy;
    /// Derivatives of squared deformation tensor.
    const T dGxx_du1x = 2 * a1 * Dxx;
    const T dGxy_du0x = b0 * Dxx;
    const T dGxy_du1x = a1 * Dxy + b1 * Dxx;
    const T dGxy_du1y = a1 * Dyy;
    const T dGyy_du0x = 2 * b0 * Dxy;
    const T dGyy_du0y = 2 * b0 * Dyy;
    const T dGyy_du1x = 2 * b1 * Dxy;
    const T dGyy_du1y = 2 * b1 * Dyy;
    /// Compute force components in common plane.
    const T force0x = dw_dI1 * (dI1_dGyy * dGyy_du0x) + dw_dI2 * (dI2_dGyy * dGyy_du0x + dI2_dGxy * dGxy_du0x);
    const T force0y = dw_dI1 * (dI1_dGyy * dGyy_du0y) + dw_dI2 * (dI2_dGyy * dGyy_du0y);
    const T force1x = dw_dI1 * (dI1_dGxx * dGxx_du1x + dI1_dGyy * dGyy_du1x) + dw_dI2 * (dI2_dGxx * dGxx_du1x + dI2_dGyy * dGyy_du1x + dI2_dGxy * dGxy_du1x);
    const T force1y = dw_dI1 * (dI1_dGyy * dGyy_du1y) + dw_dI2 * (dI2_dGyy * dGyy_du1y + dI2_dGxy * dGxy_du1y);
    /// Coordinate system.
    p_vec<> ex = pos[n_2] - pos[n_1];
    ex.norm_to(1.);
    p_vec<> ez = (pos[n_2] - pos[n_1]) % (pos[n_0] - pos[n_1]);
    ez.norm_to(1.);
    p_vec<> ey = ez % ex;
    /// Compute final forces.
    const p_vec<> force_0 = ex * force0x + ey * force0y;
    const p_vec<> force_1 = ex * force1x + ey * force1y;
    const p_vec<> force_2 = (force_0 + force_1) * -1.;
    /// Update forces and energy.
    /// Face areas are already contained in the shape functions.
    forceStrain[n_0] -= (force_0);
    forceStrain[n_1] -= (force_1);
    forceStrain[n_2] -= (force_2);
    ergStrain += (w * mesh->area[f_i] * util::pow(radius,2.));
  }

  return;
}

/*********************************************************************************
************* BENDING FORCE (ENERGY) *********************************************
*********************************************************************************/

/// The bending force is computed. If the bending modulus is non-positive, the
/// computation is skipped.
/// NOTE: Implement different bending models.

template<typename T>
void MembraneParticle3D<T>::forcesBending() {
  if(modBending < 1e-15) {
    return;
  }
  const T lambda = sqrt(3) * modBending;
  for(int f_i = 0; f_i < numFaces - 1; ++f_i) {
    for(int n = 0; n < 3; ++n) {
      const int f_j = mesh->faceNeighbours[n][f_i];
      /// Run through each pair of faces only once.
      /// Each pair of faces is considered only once, hence the force is twice as large.
      /// This additional factor 2 is cancalled with the factor 2 in the denominator of the derivative.
      if(f_j > f_i) {
        const T scalarProduct = normal[f_i] * normal[f_j];
        p_vec<> n_ij = normal[f_i] - normal[f_j] * scalarProduct;
        p_vec<> n_ji = normal[f_j] - normal[f_i] * scalarProduct;
        n_ij.norm_to(1);
        n_ji.norm_to(1);
        const T iar_i = 1 / area[f_i];
        const T iar_j = 1 / area[f_j];
        /// Identify common and single nodes
        /// n_2 is node in face f_i, but not in f_j.
        /// n_4 is node in face f_j, but not in f_i.
        /// n_1 and n_3 are the two common nodes in this specific order.
        const int n_2 = mesh->nodeDiff(f_i, f_j);
        const int n_4 = mesh->nodeDiff(f_j, f_i);
        int n_1, n_3;
        mesh->nodesSame(f_i, f_j, n_1, n_3);
        /// Compute current angle between faces (including sign).
        /// convex is a boolean for convex/concave angles; convex = 1 means that the pair of faces is convex.
        const T angleCurrent = computeAngleNormals(f_i, f_j);
        int convex;
        if((normal[f_i] % normal[f_j]) * (pos[n_1] - pos[n_3]) >= 0) {
          convex = 1;
        }
        else {
          convex = -1;
        }
        const T strength = lambda * (angleCurrent - convex * mesh->angleEq[n][f_i]);
        //  const T strength = lambda * (angleCurrent);
        /// Update force and energy.
        forceBending[n_1] += (((pos[n_2] - pos[n_3]) % n_ji * iar_i + (pos[n_3] - pos[n_4]) % n_ij * iar_j) * strength);
        forceBending[n_2] += (((pos[n_3] - pos[n_1]) % n_ji * iar_i) * strength);
        forceBending[n_3] += (((pos[n_4] - pos[n_1]) % n_ij * iar_j + (pos[n_1] - pos[n_2]) % n_ji * iar_i) * strength);
        forceBending[n_4] += (((pos[n_1] - pos[n_3]) % n_ij * iar_j) * strength);
        ergBending += (lambda * util::pow(convex * angleCurrent - mesh->angleEq[n][f_i],2.));
      }
    }
  }

  return;
}

}

}

#endif
