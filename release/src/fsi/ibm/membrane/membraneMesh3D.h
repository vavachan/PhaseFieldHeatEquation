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

#ifndef MESH_3D_H
#define MESH_3D_H

#define MAX_NODE_NEIGHBORS 7 // increase this if mesh file has nodes with more than 7 neighbours

#include "utilities/pvecten.h"
#include "io/ostreamManager.h"

namespace olb {

namespace membrane {

template<typename T>
struct CMesh {
  mutable OstreamManager clout;
  int numNodes; // Number of nodes
  int numFaces; // Number of faces
  T surface; // Surface of mesh (for radius 1)
  T volume; // Volume of mesh (for radius 1)
  std::vector<p_vec<>> pos; // Node positions (for radius 1)
  std::vector<T> area; // Face areas (for radius 1)
  std::vector<int> faceToNode[3]; // Face-to-node look-up table
  std::vector<int> nodeToFace[MAX_NODE_NEIGHBORS]; // Node-to-face look-up table
  std::vector<int> faceNeighbours[3]; // Face-neighbour look-up table
  std::vector<int> nodeNeighbours[MAX_NODE_NEIGHBORS]; // Node-neighbour look-up table
  std::vector<T> angleEq[3]; // Equilibrium angles between neighbouring faces
  std::vector<T> b0, a1, b1, a2, b2; // Shape functions
  std::vector<T> l1Eq, l2Eq, phiEq, cosEq, sinEq, sinInvEq; // Reference geometry (for radius 1)
  std::vector<int> nodeOrientation; // Defining "top" and "bottom" of the mesh
  int numPosNodeOrientation; // Number of nodes in "top" half
  int numNegNodeOrientation; // Number of nodes in "bottom" half

  CMesh(std::string meshFilename, T radius, int normalMethod) : clout(std::cout,"CMesh")
  {
    readGeometry(meshFilename);
    findFaceNeighbours();
    findNodeNeighbours();
    checkNormalVectors(normalMethod);
    setEquilibriumValues();
    meshStatistics(radius);
    initOrientationNodes();
  }

  void readGeometry(const std::string&); // Read geometry from mesh file
  void findFaceNeighbours(); // Find face neighbours
  void findNodeNeighbours(); // Find node neighbours
  bool isInList(int, int); // Auxiliary method to check whether node is already in list
  void checkNormalVectors(int); // Check face normal orientations
  bool notAllChecked(bool*); // Checks whether all faces have been checked
  void checkNormals(int, int); // Checks angle between normals
  void setEquilibriumValues(); // Set equilibrium values
  void meshStatistics(T); // Compute mesh statistics
  void initOrientationNodes(); // Nodes are sorted into "top" and "bottom" nodes
  p_vec<> computeNormalEq(int); // Compute equilibrium face normal
  int nodeDiff(int, int); // Finds node belonging to first but not to second face
  void nodesSame(int, int, int&, int&); // Finds common nodes of two given faces
  bool areFacesNeighbours(int, int); // Checks whether faces are neighbours
  bool areNodesNeighbours(int, int); // Checks whether nodes are neighbours
  bool isNodeInFace(int, int); // Checks whether node is in face
};

}

}

#endif
