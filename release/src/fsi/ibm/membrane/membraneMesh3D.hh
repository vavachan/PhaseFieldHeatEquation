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

#ifndef MESH_3D_HH
#define MESH_3D_HH

#define MAX_NODE_NEIGHBORS 7 // increase this if mesh file has nodes with more than 7 neighbours

#include <iomanip>
#include <fstream>
#include "membraneMesh3D.h"
#include "utilities/calc.h"

namespace olb {

namespace membrane {

/*********************************************************************************
************* READ GEOMETRY ******************************************************
*********************************************************************************/

/// The geometry is read from a mesh file. The initial positions of the nodes are
/// defined. The memberships of the nodes in the faces is computed.
template<typename T>
void CMesh<T>::readGeometry(const std::string& filenameGeometry) {
  std::string dummyString;
  int dummyInt;
  /// Handle geometry file.
  ifstream fileGeometry;
  fileGeometry.open(filenameGeometry.c_str());
  if(!fileGeometry) {
    throw std::runtime_error("Cannot open geometry file " + filenameGeometry + ".");
  }
  /// Ignore first four lines of file (header data).
  getline(fileGeometry, dummyString);
  getline(fileGeometry, dummyString);
  getline(fileGeometry, dummyString);
  getline(fileGeometry, dummyString);
  /// Read number of nodes and resize all relevant arrays.
  fileGeometry >> numNodes;
  pos.resize(numNodes);
  nodeOrientation.resize(numNodes);
  for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
    nodeToFace[n].resize(numNodes);
    nodeNeighbours[n].resize(numNodes);
  }
  /// Read node positions.
  /// An offset for correctly reading the geometry file is required.
  /// This offset is the index of the first node in the mesh file.
  int offset;
  fileGeometry >> offset;
  fileGeometry >> pos[0].x >> pos[0].y >> pos[0].z;
  for(int n_i = 1; n_i < numNodes; ++n_i) {
    fileGeometry >> dummyInt >> pos[n_i].x >> pos[n_i].y >> pos[n_i].z;
  }
  /// Ignore further header lines.
  getline(fileGeometry, dummyString);
  getline(fileGeometry, dummyString);
  getline(fileGeometry, dummyString);
  /// Read number of faces and resize all relevant arrays.
  fileGeometry >> numFaces;
  area.resize(numFaces);
  angleEq[0].resize(numFaces);
  angleEq[1].resize(numFaces);
  angleEq[2].resize(numFaces);
  faceToNode[0].resize(numFaces);
  faceToNode[1].resize(numFaces);
  faceToNode[2].resize(numFaces);
  faceNeighbours[0].resize(numFaces);
  faceNeighbours[1].resize(numFaces);
  faceNeighbours[2].resize(numFaces);
  b0.resize(numFaces);
  a1.resize(numFaces);
  b1.resize(numFaces);
  a2.resize(numFaces);
  b2.resize(numFaces);
  l1Eq.resize(numFaces);
  l2Eq.resize(numFaces);
  phiEq.resize(numFaces);
  cosEq.resize(numFaces);
  sinEq.resize(numFaces);
  sinInvEq.resize(numFaces);
  /// Create face-to-node look-up table and correct face index by the offset.
  int f_a, f_b, f_c;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    fileGeometry >> dummyInt >> dummyInt >> dummyInt >> dummyInt >> dummyInt >> dummyInt;
    fileGeometry >> f_a >> f_b >> f_c;
    faceToNode[0][f_i] = f_a - offset;
    faceToNode[1][f_i] = f_b - offset;
    faceToNode[2][f_i] = f_c - offset;
  }
  fileGeometry.close();
  clout << "Number of nodes/faces: " << numNodes << "/" << numFaces << std::endl;

  return;
}

/*********************************************************************************
************* FIND FACE NEIGHBOURS ***********************************************
*********************************************************************************/

/// The look-up table for the face neighbours is created. This information is
/// extracted from the face-to-node look-up table. Two faces are neighbours if
/// and only if they share two nodes. Every face has exactly three neighbours
/// (as long as the surface is closed). The standard value for the look-up table
/// is negative. This way, invalid entries can be identified later.
template<typename T>
void CMesh<T>::findFaceNeighbours() {
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    faceNeighbours[0][f_i] = -1;
    faceNeighbours[1][f_i] = -1;
    faceNeighbours[2][f_i] = -1;
  }
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_j = faceToNode[0][f_i];
    const int n_k = faceToNode[1][f_i];
    const int n_l = faceToNode[2][f_i];
    /// Find the other faces with two same nodes.
    for(int f_j = 0; f_j < numFaces; ++f_j) {
      if(f_j != f_i) {
        if(isNodeInFace(n_j, f_j) && isNodeInFace(n_k, f_j)) {
          faceNeighbours[0][f_i] = f_j;
        }
        if(isNodeInFace(n_k, f_j) && isNodeInFace(n_l, f_j)) {
          faceNeighbours[1][f_i] = f_j;
        }
        if(isNodeInFace(n_l, f_j) && isNodeInFace(n_j, f_j)) {
          faceNeighbours[2][f_i] = f_j;
        }
      }
    }
  }
  /// Check whether all next neighbours have been found.
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    if(isNodeInFace(-1, f_i)) {
      throw std::runtime_error("Could not find all face neighbours, first occured for face " + std::to_string(f_i) + ".");
    }
  }

  return;
}

/*********************************************************************************
************* FIND NODE NEIGHBOURS ***********************************************
*********************************************************************************/

/// The node neighbours look-up table is created. This information is extracted
/// from the face-to-node look-up table. Two nodes are neighbours when they are
/// found in the same face. Every node has usually six neighbours, but this number
/// is not fixed. The node-to-face look-up table is also created.
template<typename T>
void CMesh<T>::findNodeNeighbours() {
  /// Create node-to-face look-up table.
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
      nodeToFace[n][n_i] = -1; // standard value for 'not used'
    }
    for(int f_j = 0; f_j < numFaces; ++f_j) {
      if(isNodeInFace(n_i, f_j)) {
        /// find number of facet for given node
        int n = 0;
        while(nodeToFace[n][n_i] != -1) {
          ++n;
        }
        nodeToFace[n][n_i] = f_j;
      }
    }
  }
  /// Initialise node neighbours look-up table.
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
      nodeNeighbours[n][n_i] = -1; // standard value for 'not used'
    }
  }
  /// Fill look-up table with data.
  /// An auxilliary array is temporarily used for this purpose.
  int auxarray[3 * MAX_NODE_NEIGHBORS];
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    for(int m = 0; m < 3 * MAX_NODE_NEIGHBORS; ++m) {
      auxarray[m] = -1; // standard value for 'not used'
    }
    /// Fill auxiliary array with candidates for next neighbours.
    for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
      const int f_j = nodeToFace[n][n_i];
      if(f_j != -1) {
        auxarray[3 * n + 0] = faceToNode[0][f_j];
        auxarray[3 * n + 1] = faceToNode[1][f_j];
        auxarray[3 * n + 2] = faceToNode[2][f_j];
      }
    }
    /// Filter candidates by avoiding doubles.
    /// Neglect all '-1' entries.
    /// The reference node and all neighbours are already found.
    int n = 0;
    for(int m = 0; m < 3 * MAX_NODE_NEIGHBORS; ++m) {
      if(auxarray[m] != -1 && auxarray[m] != n_i && !isInList(auxarray[m], n_i)) {
        nodeNeighbours[n][n_i] = auxarray[m];
        ++n;
      }
    }
  }
  /// Check whether all next neighbours have really been found.
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    if(nodeNeighbours[0][n_i] == -1 || nodeToFace[0][n_i] == -1) {
      throw std::runtime_error("Could not find all node neighbours, first occured for node " + std::to_string(n_i) + ".");
    }
  }

  return;
}

/*********************************************************************************
************* CHECKS WHETHER NODE IS IN LIST *************************************
*********************************************************************************/

/// Checks whether a next neighbour node has already been found.
template<typename T>
bool CMesh<T>::isInList(int m, int n_i) {
  for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
    if(nodeNeighbours[n][n_i] == m) {
      return true;
    }
  }

  return false;
}

/*********************************************************************************
************* CHECK NORMAL VECTORS ***********************************************
*********************************************************************************/

/// The normal vectors have to point outwards. If the direction is not already
/// correct in the mesh file, it is corrected in this method. There are two
/// approaches:
/// 1) Start from one known normal vector orientation, check the remaining vectors
///    successively by comparison of angles.
/// 2) For initially concave shapes, it suffices to compare normals with centre
///    vectors.
/// If it is necessary, two nodes in the face-to-node look-up table are swapped.
/// This inverts the normal vector direction.
template<typename T>
void CMesh<T>::checkNormalVectors(int normalMethod) {
  if(normalMethod != 0 && normalMethod != 1) {
    throw std::runtime_error("No valid normalization method choosen: " + std::to_string(normalMethod));
  }
  /// Method 1: use known normal vector.
  if(normalMethod == 0) {
    /// Use temporary look-up table to see whether normal vector has already been checked.
    /// Start with a known normal vector orientation.
    bool *isChecked;
    isChecked = new bool[numFaces];
    for(int f_i = 0; f_i < numFaces; ++f_i) {
      isChecked[f_i] = false;
    }
    /// Compute centre position.
    p_vec<> cent(0, 0, 0);
    for(int n_i = 0; n_i < numNodes; ++n_i) {
      cent += pos[n_i];
    }
    cent /= numNodes;
    /// Check angle for the first face.
    /// Use first face (0) as reference.
    const int f_0 = 0;
    p_vec<> fcent = (pos[faceToNode[0][f_0]] + pos[faceToNode[1][f_0]] + pos[faceToNode[2][f_0]]) / 3. - cent;
    p_vec<> norm = (pos[faceToNode[0][f_0]] - pos[faceToNode[1][f_0]]) % (pos[faceToNode[2][f_0]] - pos[faceToNode[1][f_0]]);
    fcent.norm_to(1);
    norm.norm_to(1);
    /// Swap nodes if necessary.
    if(fcent * norm < 0) {
      const int temp = faceToNode[0][f_0];
      faceToNode[0][f_0] = faceToNode[2][f_0];
      faceToNode[2][f_0] = temp;
    }
    isChecked[f_0] = true;
    /// Successively check all neighbours.
    while(notAllChecked(isChecked)) {
      for(int f_i = 0; f_i < numFaces; ++f_i) {
        if(isChecked[f_i]) {
          const int f_1 = faceNeighbours[0][f_i];
          const int f_2 = faceNeighbours[1][f_i];
          const int f_3 = faceNeighbours[2][f_i];
          if(!isChecked[f_1]) {
            checkNormals(f_1, f_i);
            isChecked[f_1] = true;
          }

          if(!isChecked[f_2]) {
            checkNormals(f_2, f_i);
            isChecked[f_2] = true;
          }

          if(!isChecked[f_3]) {
            checkNormals(f_3, f_i);
            isChecked[f_3] = true;
          }
        }
      }
    }
    /// Check whether all angles are small.
    /// If an angle is larger than 90°, something is wrong.
    for(int f_i = 0; f_i < numFaces - 1; ++f_i) {
      for(int f_j = f_i + 1; f_j < numFaces; ++f_j) {
        if(areFacesNeighbours(f_i, f_j)) {
          const p_vec<> norm_i = (pos[faceToNode[0][f_i]] - pos[faceToNode[1][f_i]]) % (pos[faceToNode[2][f_i]] - pos[faceToNode[1][f_i]]);
          const p_vec<> norm_j = (pos[faceToNode[0][f_j]] - pos[faceToNode[1][f_j]]) % (pos[faceToNode[2][f_j]] - pos[faceToNode[1][f_j]]);
          if(norm_i * norm_j < 0) {
            throw std::runtime_error("Could not correct all normal vectors.");
          }
        }
      }
    }
  }
  /// Method 2: use initially convex membrane.
  if(normalMethod == 1) {
    /// Compute centre position.
    p_vec<> cent(0, 0, 0);
    for(int n_i = 0; n_i < numNodes; ++n_i) {
      cent += pos[n_i];
    }
    cent /= numNodes;
    /// Compare normal and radial vector.
    for(int f_i = 0; f_i < numFaces; ++f_i) {
      const int n_j = faceToNode[0][f_i];
      const int n_k = faceToNode[1][f_i];
      const int n_l = faceToNode[2][f_i];
      const p_vec<> norm = (pos[n_j] - pos[n_k]) % (pos[n_l] - pos[n_k]);
      const p_vec<> rad = (pos[n_j] + pos[n_k] + pos[n_l]) / 3. - cent;
      /// Check sign and swap sign if necessary.
      if(rad * norm < 0) {
        const int temp = faceToNode[0][f_i];
        faceToNode[0][f_i] = faceToNode[2][f_i];
        faceToNode[2][f_i] = temp;
      }
    }
  }

  return;
}

/*********************************************************************************
************* CHECKS WHETHER ALL FACES HAVE BEEN CHECKED *************************
*********************************************************************************/

/// It is checked whether there are faces left to be checked for their normal
/// vector orientation. Returns true if at least one face has not been checked
/// and false otherwise.
template<typename T>
bool CMesh<T>::notAllChecked(bool *isChecked) {
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    if(!isChecked[f_i]) {
      return true;
    }
  }

  return false;
}

/*********************************************************************************
************* CHECKS ANGLE BETWEEN TWO FACES NORMALS *****************************
*********************************************************************************/

/// The angle between face f_i and reference face f_0 is checked. If the angle is
/// larger than 90°, the orientation of f_i has to be changed. For doing so, the
/// first and third elements in the face-to-node look-up table are swapped,
/// which changes the sign of the normal.
template<typename T>
void CMesh<T>::checkNormals(int f_i, int f_0) {
  const p_vec<> norm_0 = (pos[faceToNode[0][f_0]] - pos[faceToNode[1][f_0]]) % (pos[faceToNode[2][f_0]] - pos[faceToNode[1][f_0]]);
  const p_vec<> norm_i = (pos[faceToNode[0][f_i]] - pos[faceToNode[1][f_i]]) % (pos[faceToNode[2][f_i]] - pos[faceToNode[1][f_i]]);
  if(norm_i * norm_0 < 0) {
    const int temp = faceToNode[0][f_i];
    faceToNode[0][f_i] = faceToNode[2][f_i];
    faceToNode[2][f_i] = temp;
  }

  return;
}

/*********************************************************************************
************* SET EQUILIBRIUM VALUES *********************************************
*********************************************************************************/

/// Starting from the geometry in the mesh file, the relevant equilibrium values
/// are computed for unit radius: face areas, surface, volume, normal angles,
/// shape functions.
template<typename T>
void CMesh<T>::setEquilibriumValues() {
  /// Compute temporary membrane centre.
  p_vec<> centreMesh(0, 0, 0);
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    centreMesh += pos[n_i];
  }
  centreMesh /= numNodes;
  /// Compute equilibrium face areas, surface and volume.
  surface = 0;
  volume = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_j = faceToNode[0][f_i];
    const int n_k = faceToNode[1][f_i];
    const int n_l = faceToNode[2][f_i];
    const p_vec<> edge1 = pos[n_j] - pos[n_k];
    const p_vec<> edge2 = pos[n_l] - pos[n_k];
    const p_vec<> normal = edge1 % edge2;
    const p_vec<> faceCentre = (pos[n_j] + pos[n_k] + pos[n_l]) / 3 - centreMesh;
    area[f_i] = normal.length() / 2;
    surface += area[f_i];
    volume += ((computeNormalEq(f_i) * faceCentre) * area[f_i] / 3);
  }
  /// Compute equilibrium face angles.
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    for(int n = 0; n < 3; ++n) {
      const int f_j = faceNeighbours[n][f_i];
      int n_1, n_2;
      nodesSame(f_i, f_j, n_1, n_2);
      const T sign = (computeNormalEq(f_i) % computeNormalEq(f_j)) * (pos[n_1] - pos[n_2]);
      if(sign > 0) {
        angleEq[n][f_i] = acos(computeNormalEq(f_i) * computeNormalEq(f_j));
      }
      else if(sign < 0) {
        angleEq[n][f_i] = -acos(computeNormalEq(f_i) * computeNormalEq(f_j));
      }
      else {
        angleEq[n][f_i] = 0;
      }
    }
  }
  /// Compute shape functions.
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    const int n_0 = faceToNode[0][f_i];
    const int n_1 = faceToNode[1][f_i];
    const int n_2 = faceToNode[2][f_i];
    l1Eq[f_i] = (pos[n_2] - pos[n_1]).length();
    l2Eq[f_i] = (pos[n_0] - pos[n_1]).length();
    phiEq[f_i] = acos((pos[n_2] - pos[n_1]) * (pos[n_0] - pos[n_1]) / (l1Eq[f_i] * l2Eq[f_i]));
    cosEq[f_i] = cos(phiEq[f_i]);
    sinEq[f_i] = sin(phiEq[f_i]);
    sinInvEq[f_i] = 1 / sinEq[f_i];
    b0[f_i] = l1Eq[f_i] / 2;
    a1[f_i] = -l2Eq[f_i] * sinEq[f_i] / 2;
    b1[f_i] = (l2Eq[f_i] * cosEq[f_i] - l1Eq[f_i]) / 2;
    a2[f_i] = l2Eq[f_i] * sinEq[f_i] / 2;
    b2[f_i] = -l2Eq[f_i] * cosEq[f_i] / 2;
  }

  return;
}

/*********************************************************************************
************* ANALYSE MESH STATISTICS ********************************************
*********************************************************************************/

/// Mesh statistics are computed and reported for a given radius: surface,
/// volume, effective radius, reduced volume, face areas, normal angles, edge
/// angles, minimum and maximum number of node neighbours. This function is used
/// for obtaining mesh statistics only.
template<typename T>
void CMesh<T>::meshStatistics(T radius) {
  /// Compute volume, surface and effective radius.
  const T volumeStat = volume * util::pow(radius,3.);
  const T surfaceStat = surface * util::pow(radius,2.);
  const T radiusEff = sqrt(surface * util::pow(radius,2.) / (4 * M_PI));
  const T volumeRed = 3 * volume * util::pow(radius,3) / (4 * M_PI * util::pow(radiusEff,3.));
  /// Compute minimum, maximum and mean areas.
  T areaMin = 10 * util::pow(radius,2.);
  T areaMax = 0;
  T areaMean = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    areaMean += area[f_i];
    if(area[f_i] < areaMin) {
      areaMin = area[f_i];
    }
    if(area[f_i] > areaMax) {
      areaMax = area[f_i];
    }
  }
  areaMean /= numFaces;
  areaMin *= util::pow(radius,2.);
  areaMax *= util::pow(radius,2.);
  /// Compute area variance.
  T areaVar = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    areaVar += util::pow(areaMean - area[f_i], 2.);
  }
  areaVar = sqrt(areaVar / numFaces) / areaMean;
  areaMean *= util::pow(radius,2.);
  /// Compute minimum, maximum and mean angles.
  T angleMin = M_PI;
  T angleMax = 0;
  T angleMean = 0;
  int m = 0;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    for(int n = 0; n < 3; ++n) {
      const T angleMag = abs(angleEq[n][f_i]);
      angleMean += angleMag;
      if(angleMag < angleMin) {
        angleMin = angleMag;
      }
      if(angleMag > angleMax) {
        angleMax = angleMag;
      }
      ++m;
    }
  }
  angleMean /= m;
  /// Compute angle variance.
  T angleVar = 0;
  for(int f_i = 0; f_i < numFaces - 1; ++f_i) {
    for(int n = 0; n < 3; ++n) {
      angleVar += util::pow(angleMean - angleEq[n][f_i],2.);
    }
  }
  angleVar = sqrt(angleVar / m) / angleMean;
  /// Compute minimum and maximum angles and angle variance.
  const T angleEdgeMean = M_PI / 3.; // mean angle must be 180° / 3
  T angleEdgeMin = M_PI;
  T angleEdgeMax = 0.;
  T angleEdgeVar = 0.;
  for(int f_i = 0; f_i < numFaces; ++f_i) {
    for(int i = 0; i < 3; ++i) {
      const p_vec<> vec1 = pos[faceToNode[(i + 1) % 3][f_i]] - pos[faceToNode[i % 3][f_i]];
      const p_vec<> vec2 = pos[faceToNode[(i + 2) % 3][f_i]] - pos[faceToNode[i % 3][f_i]];
      const T angleMag = abs(acos(vec1 * vec2 / (vec1.length() * vec2.length())));
      if(angleMag < angleEdgeMin) {
        angleEdgeMin = angleMag;
      }
      if(angleMag > angleEdgeMax) {
        angleEdgeMax = angleMag;
      }
      angleEdgeVar += util::pow(angleEdgeMean - angleMag,2);
    }
  }
  angleEdgeVar = sqrt(angleEdgeVar / (3. * numFaces)) / angleEdgeMean;
  /// Compute minimum, maximum and mean edges.
  T lengthMean = 0;
  T lengthMin = 10 * radius;
  T lengthMax = 0;
  T lengthVar = 0;
  int n = 0;
  for(int n_i = 0; n_i < numNodes - 1; ++n_i) {
    for(int n_j = n_i + 1; n_j < numNodes; ++n_j) {
      if(areNodesNeighbours(n_i, n_j)) {
        const T length = (pos[n_i] - pos[n_j]).length();
        lengthMean += length;
        if(length < lengthMin) {
          lengthMin = length;
        }
        if(length > lengthMax) {
          lengthMax = length;
        }
        ++n;
      }
    }
  }
  lengthMean /= n;
  lengthMin *= radius;
  lengthMax *= radius;
  /// Compute edge variance.
  for(int n_i = 0; n_i < numNodes - 1; ++n_i) {
    for(int n_j = n_i + 1; n_j < numNodes; ++n_j) {
      if(areNodesNeighbours(n_i, n_j)) {
        const T length = (pos[n_i] - pos[n_j]).length();
        lengthVar += (util::pow(lengthMean - length,2));
      }
    }
  }
  lengthVar = sqrt(lengthVar / n) / lengthMean;
  lengthMean *= radius;
  /// Number of neighbouring nodes.
  int numMax = 0;
  int numMin = MAX_NODE_NEIGHBORS;
  int numCurr = 0;
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    int l = 0;
    while(nodeToFace[l][n_i] != -1) {
      ++l;
      numCurr = l;
    }
    if(numCurr < numMin) {
      numMin = numCurr;
    }
    if(numCurr > numMax) {
      numMax = numCurr;
    }
  }
  /// Report results.
  clout << "Radius: " << radius << std::endl;
  clout << "Volume: " << volumeStat << std::endl;
  clout << "Surface: " << surfaceStat << std::endl;
  clout << "Effective radius: " << radiusEff << std::endl;
  clout << "Reduced volume: " << volumeRed << std::endl;
  clout << "Number of node neighbours: between " << numMin << " and " << numMax << std::endl;
  clout << setprecision(3);
  clout << "Edge lengths: between " << lengthMin << " and " << lengthMax << ", average = " <<
     lengthMean << ", variance = " << 100 * lengthVar << "%" << std::endl;
  clout << "Face areas: between " << areaMin << " and " << areaMax << ", average = " << areaMean <<
     ", variance = " << 100 * areaVar << "%" << std::endl;
  clout << "Normal-normal angles: between " << angleMin << " and " << angleMax << ", average = " <<
     angleMean << ", variance = " << 100 * angleVar << "%" << std::endl;
  clout << "Edge-edge angles: between " << angleEdgeMin << " and " << angleEdgeMax <<
    ", average = " << angleEdgeMean << ", variance = " << 100 * angleEdgeVar << "%" << std::endl;
  return;
}

/*********************************************************************************
************* SET ORIENTATION NODES **********************************************
*********************************************************************************/

/// Each membrane node gets an index (-1, 0, 1) which identifies it as top,
/// bottom, or neutral node. Top nodes are initially located above the xy-plane
/// (z > 0), bottom nodes initially below the xy-plane (z < 0). Nodes initially
/// located close to the xy-plane (z = 0) are omitted. This information can later
/// be used to identify the top and the bottom of an arbitrarily deformed and
/// rotated membrane.
template<typename T>
void CMesh<T>::initOrientationNodes() {
  numPosNodeOrientation = 0;
  numNegNodeOrientation = 0;
  for(int n_i = 0; n_i < numNodes; ++n_i) {
    if(pos[n_i].z > 0.0001) {
      nodeOrientation[n_i] = 1;
      numPosNodeOrientation++;
    }
    else if(pos[n_i].z < -0.0001) {
      nodeOrientation[n_i] = -1;
      numNegNodeOrientation++;
    }
    else {
      nodeOrientation[n_i] = 0;
    }
  }

  return;
}

/*********************************************************************************
************* COMPUTE EQUILIBRIUM NORMAL VECTOR **********************************
*********************************************************************************/

/// The normal vector of a face in the undeformed mesh is computed.
template<typename T>
inline p_vec<> CMesh<T>::computeNormalEq(int f_i) {
  const int n_j = faceToNode[0][f_i];
  const int n_k = faceToNode[1][f_i];
  const int n_l = faceToNode[2][f_i];
  const p_vec<> normal = (pos[n_j] - pos[n_k]) % (pos[n_l] - pos[n_k]);

  return normal / normal.length();
}

/*********************************************************************************
************* NODE IN ONE FACE BUT NOT IN THE OTHER ******************************
*********************************************************************************/

/// Providing two neighbouring faces, the node belonging to the first face but
/// not to the second face is returned. It is not checked whether both faces are
/// neighbours at all.
template<typename T>
int CMesh<T>::nodeDiff(int f_i, int f_j) {
  const int n_j = faceToNode[0][f_i];
  const int n_k = faceToNode[1][f_i];
  const int n_l = faceToNode[2][f_i];
  if(!isNodeInFace(n_j, f_j)) {
    return n_j;
  }
  else if(!isNodeInFace(n_k, f_j)) {
    return n_k;
  }
  else {
    return n_l;
  }
}

/*********************************************************************************
************* COMMON NODES IN FACES **********************************************
*********************************************************************************/

/// The two common nodes of two neighbouring faces are computed. It is not checked
/// whether both faces are neighbours. The cycling order of nodes is
/// n_j -> n_k -> n_l -> n_j.
template<typename T>
void CMesh<T>::nodesSame(int f_i, int f_j, int& n_1, int& n_2) {
  const int n_j = faceToNode[0][f_i];
  const int n_k = faceToNode[1][f_i];
  const int n_l = faceToNode[2][f_i];
  if(n_j == nodeDiff(f_i, f_j)) {
    n_1 = n_k;
    n_2 = n_l;
  }
  else if(n_k == nodeDiff(f_i, f_j)) {
    n_1 = n_l;
    n_2 = n_j;
  }
  else {
    n_1 = n_j;
    n_2 = n_k;
  }

  return;
}

/*********************************************************************************
************* CHECK WHETHER FACES ARE NEIGHBOURS *********************************
*********************************************************************************/

/// Checks whether two faces are neighbours.
template<typename T>
inline bool CMesh<T>::areFacesNeighbours(int f_i, int f_j) {
  if(f_j == faceNeighbours[0][f_i] || f_j == faceNeighbours[1][f_i] || f_j == faceNeighbours[2][f_i]) {
    return true;
  }
  else {
    return false;
  }
}

/*********************************************************************************
************* CHECK WHETHER NODES ARE NEIGHBOURS *********************************
*********************************************************************************/

/// Checks whether two nodes are neighbours.
template<typename T>
inline bool CMesh<T>::areNodesNeighbours(int n_i, int n_j) {
  for(int n = 0; n < MAX_NODE_NEIGHBORS; ++n) {
    if(n_i == nodeNeighbours[n][n_j]) {
      return true;
    }
  }

  return false;
}

/*********************************************************************************
************* CHECK WHETHER NODE IS IN FACE **************************************
*********************************************************************************/

/// Checks whether a node is member of a face.
template<typename T>
inline bool CMesh<T>::isNodeInFace(int n_i, int f_j) {
  if(n_i == faceToNode[0][f_j] || n_i == faceToNode[1][f_j] || n_i == faceToNode[2][f_j]) {
    return true;
  }
  else {
    return false;
  }
}

}

}

#endif
