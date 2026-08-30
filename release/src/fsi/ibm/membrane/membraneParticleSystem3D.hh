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

#ifndef SUPER_MEMBRANE_3D_HH
#define SUPER_MEMBRANE_3D_HH

#include "membraneParticle3D.h"
#include "core/superD.h"
#include "core/platform/platform.h"
#include "descriptor/fields.h"
#include "io/ostreamManager.h"
#include <tinyxml2.h>

namespace olb {

namespace membrane {

template<typename T, typename DESCRIPTOR>
MembraneParticleSystem3D<T,DESCRIPTOR>::MembraneParticleSystem3D() : clout(std::cout,"MembraneParticleSystem3D")
{ }

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::readMeshParameters(std::string filename)
{
  OstreamManager clout(std::cout,"MembraneParticleSystem3D::readMeshParameters");
  clout << "Reading mesh simulation parameters..." << std::endl;
  // Open XML file and check validity.
  tinyxml2::XMLDocument xml_doc;
  bool eResult = xml_doc.LoadFile(filename.c_str());
  if (eResult != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Cannot load or parse XML mesh parameter file.");
  }
  else {
    clout << "Successfully accessed " << filename << "." << std::endl;
  }
  tinyxml2::XMLElement *pElement;
  // Read initialisation parameters.
  clout << "Reading initialisation parameters..." << std::endl;
  pElement = xml_doc.FirstChildElement("init");
  if(pElement == NULL) {
    throw std::runtime_error("Cannot find initialisation parameters.");
  }
  // Read overlap check.
  clout << "Reading overlap check parameters..." << std::endl;
  pElement = xml_doc.FirstChildElement("init")->FirstChildElement("overlapCheck");
  if(pElement == NULL) {
    clout << "WARNING! Cannot find overlap checks; setting all to true." << std::endl;
    _checkInitialParticleOverlap = true;
  }
  else {
    // Read particle-particle overlap check.
    int flagCheckParticleParticleOverlap;
    eResult = pElement->QueryIntAttribute("particleParticle", &flagCheckParticleParticleOverlap);
    if(eResult != tinyxml2::XML_SUCCESS) {
      clout << "WARNING! Cannot find particle-particle overlap check; setting to true." << std::endl;
      _checkInitialParticleOverlap = true;
    }
    else {
      if(flagCheckParticleParticleOverlap > 0) {
        _checkInitialParticleOverlap = true;
      }
      else {
        _checkInitialParticleOverlap = false;
      }
      clout << "Particle-particle overlap check: " << flagCheckParticleParticleOverlap << std::endl;
    }
  }
  // Read normalisation method.
  clout << "Reading normalisation parameters..." << std::endl;
  pElement = xml_doc.FirstChildElement("init")->FirstChildElement("normalisation");
  if(pElement == NULL) {
    clout << "WARNING! Cannot find normalisation parameters; setting to 0." << std::endl;
    _normalMethod = 0;
  }
  else {
    // Read normalisation method.
    eResult = pElement->QueryIntAttribute("method", &_normalMethod);
    if(eResult != tinyxml2::XML_SUCCESS) {
      clout << "WARNING! Cannot find normalisation method; setting to 0." << std::endl;
      _normalMethod = 0;
    }
    else {
      if(_normalMethod == 0 || _normalMethod == 1) {
        clout << "Normalisation method: " << _normalMethod << std::endl;
      }
      else {
        clout << "Normalisation method needs to be 0 or 1; setting to 0." << std::endl;
        _normalMethod = 0;
      }
    }
  }

  // Read mesh information.
  clout << "Reading meshes..." << std::endl;
  unsigned _numMeshes = 0;
  unsigned _numParticles = 0;
  tinyxml2::XMLElement* root = xml_doc.FirstChildElement("mesh");
  // Loop over as many meshes as are found in parameter file.
  for(tinyxml2::XMLElement* e = root; e != NULL; e = e->NextSiblingElement("mesh")) {
    clout << "Found mesh:" << std::endl;
    // Read general mesh properties.
    pElement = e->FirstChildElement("general");
    if(pElement == NULL) {
      clout << "WARNING! No mesh properties provided; ignoring mesh." << std::endl;
      continue;
    }
    // Get number of particles for mesh.
    int numParticlesTemp;
    eResult = pElement->QueryIntAttribute("numParticles", &numParticlesTemp);
    if(eResult != tinyxml2::XML_SUCCESS) {
      clout << "WARNING! Cannot find number of particles; ignoring mesh." << std::endl;
      continue;
    }
    if(numParticlesTemp <= 0) {
      clout << "WARNING! Non-positive number of particles specified; ignoring mesh." << std::endl;
      continue;
    }
    // Read mesh filename.
    const char *filename;
    filename = pElement->Attribute("file");
    if(filename == NULL) {
      clout << "WARNING! No mesh filename found; ignoring mesh." << std::endl;
      continue;
    }
    // Get mesh radius.
    double radiusTemp;
    eResult = pElement->QueryDoubleAttribute("radius", &radiusTemp);
    if(eResult != tinyxml2::XML_SUCCESS) {
      clout << "WARNING! Cannot find mesh radius; ignoring mesh." << std::endl;
      continue;
    }
    if(radiusTemp <= 0) {
      clout << "WARNING! Non-positive radius specified; ignoring mesh." << std::endl;
      continue;
    }
    _numMeshes++;
    _numParticles += numParticlesTemp;
    _meshFilename.push_back(filename);
    _meshRadius.push_back(static_cast<T>(radiusTemp));
    _meshNumParticles.push_back(numParticlesTemp);
    clout << "Mesh file: " << filename << std::endl;
    clout << "Number of particles: " << numParticlesTemp << std::endl;
    clout << "Mesh radius: " << radiusTemp << std::endl;
    // Read physical properties.
    pElement = e->FirstChildElement("physics");
    if(pElement == NULL) {
      clout << "WARNING! No physical properties provided; setting all to zero." << std::endl;
      _meshKV.push_back(0.);
      _meshKSph.push_back(0.);
      _meshKA.push_back(0.);
      _meshKAlpha.push_back(0.);
      _meshKS.push_back(0.);
      _meshKB.push_back(0.);
    }
    else {
      // Volume modulus.
      double kVTemp;
      eResult = pElement->QueryDoubleAttribute("kV", &kVTemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kV specified; setting to 0." << std::endl;
        _meshKV.push_back(T(0.));
      }
      else {
        _meshKV.push_back(static_cast<T>(kVTemp));
        clout << "kV: " << kVTemp << std::endl;
      }
      // Sphere modulus.
      double kSphTemp;
      eResult = pElement->QueryDoubleAttribute("kSph", &kSphTemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kSph specified; setting to 0." << std::endl;
      // Parse parameters from XML file
        _meshKSph.push_back(T(0.));
      }
      else {
        _meshKSph.push_back(static_cast<T>(kSphTemp));
        clout << "kSph: " << kSphTemp << std::endl;
      }
      // Surface modulus.
      double kATemp;
      eResult = pElement->QueryDoubleAttribute("kA", &kATemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kA specified; setting to 0." << std::endl;
        _meshKA.push_back(T(0.));
      }
      else {
        _meshKA.push_back(static_cast<T>(kATemp));
        clout << "kA: " << kATemp << std::endl;
      }
      // Area dilation modulus.
      double kAlphaTemp;
      eResult = pElement->QueryDoubleAttribute("kalpha", &kAlphaTemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kalpha specified; setting to 0." << std::endl;
        _meshKAlpha.push_back(T(0.));
      }
      else {
        _meshKAlpha.push_back(static_cast<T>(kAlphaTemp));
        clout << "kalpha: " << kAlphaTemp << std::endl;
      }
      // Strain modulus.
      double kSTemp;
      eResult = pElement->QueryDoubleAttribute("kS", &kSTemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kS specified; setting to 0." << std::endl;
        _meshKS.push_back(T(0.));
      }
      else {
        _meshKS.push_back(static_cast<T>(kSTemp));
        clout << "kS: " << kSTemp << std::endl;
      }
      // Bending modulus.
      double kBTemp;
      eResult = pElement->QueryDoubleAttribute("kB", &kBTemp);
      if(eResult != tinyxml2::XML_SUCCESS) {
        clout << "WARNING! No kB specified; setting to 0." << std::endl;
        _meshKB.push_back(T(0.));
      }
      else {
        _meshKB.push_back(static_cast<T>(kBTemp));
        clout << "kB: " << kBTemp << std::endl;
      }
    }
  }
  clout << "Found " << _numMeshes << " meshes with " << _numParticles << " particles in total." << std::endl;
}

/// The meshes are initialized
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::initializeMeshesFromXML(std::string filename)
{
  clout << "Initialising Lagrangian meshes from file \"" << filename << "\"." << std::endl;
  readMeshParameters(filename);

  // Create mesh table for particles
  for (unsigned m_i = 0; m_i < _meshFilename.size(); ++m_i) {
    for (unsigned n = 0; n < _meshNumParticles[m_i]; ++n) {
      _meshForParticle.push_back(m_i);
    }
  }

  // Create meshes.
  for(unsigned n = 0; n < _meshFilename.size(); ++n) {
    clout << "Setting up mesh " << n << "..." << std::endl;
    _meshes.emplace_back(std::make_shared<CMesh<T>>(_meshFilename[n], _meshRadius[n], _normalMethod));
  }
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::readPositionParameters(std::string filename)
{
  OstreamManager clout(std::cout,"MembraneParticleSystem3D::readPositionParameters");
  clout << "Reading particle positions..." << std::endl;
  // Open XML file and check validity.
  tinyxml2::XMLDocument xml_doc;
  bool eResult = xml_doc.LoadFile(filename.c_str());
  if(eResult != tinyxml2::XML_SUCCESS) {
    throw std::runtime_error("Cannot load or parse parameter file.Pos");
  }
  else {
    clout << "Successfully accessed " << filename << "." << std::endl;
  }
  // Getting positions from file.
  clout << "Getting positions..." << std::endl;
  unsigned numPositions = 0;
  initParticleStruct particleTemp;
  tinyxml2::XMLElement* root = xml_doc.FirstChildElement("particle");
  for(tinyxml2::XMLElement* e = root; e != NULL; e = e->NextSiblingElement("particle")) {
    eResult = e->QueryDoubleAttribute("X", &particleTemp.centre.x);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get x-position.");
    }
    eResult = e->QueryDoubleAttribute("Y", &particleTemp.centre.y);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get y-position.");
    }
    eResult = e->QueryDoubleAttribute("Z", &particleTemp.centre.z);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get z-position.");
    }
    eResult = e->QueryDoubleAttribute("angle", &particleTemp.angle);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get angle.");
    }
    eResult = e->QueryDoubleAttribute("axisX", &particleTemp.axis.x);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get x-axis.");
    }
    eResult = e->QueryDoubleAttribute("axisY", &particleTemp.axis.y);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get y-axis.");
    }
    eResult = e->QueryDoubleAttribute("axisZ", &particleTemp.axis.z);
    if(eResult != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("Cannot get z-axis.");
    }
    clout << "Found position " << numPositions << ": pos = " << particleTemp.centre <<
      ", angle = " << particleTemp.angle << ", axis = " << particleTemp.axis << std::endl;
    _initParticles.push_back(particleTemp);
    numPositions++;
  }
  if(numPositions < _initParticles.size()) {
    throw std::runtime_error("Number of particle positions (" + std::to_string(numPositions) +
      ") is smaller than number of particles (" + std::to_string(_initParticles.size()) + ").");
  }

  if (_meshes.size() == 0) {
    throw std::runtime_error("No meshes found. Number meshes is equal zero.");
  } else if (_initParticles.size() == 0) {
    throw std::runtime_error("No particles found. Number particles is equal zero.");
  }
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::checkParticleParticleOverlap()
{
  clout << "Check particle-particle overlap." << std::endl;
  for(unsigned c_i = 0; c_i < _initParticles.size(); ++c_i) {
    for(unsigned c_j = 0; c_j < c_i; ++c_j) {
      if (c_i != c_j) {
        const p_vec<> distVec = _initParticles[c_j].centre - _initParticles[c_i].centre;
        if(distVec.length() < (_initParticles[c_i].radius_pos + _initParticles[c_j].radius_pos) + 1) {
          throw std::runtime_error("Particle " + std::to_string(c_i) + " overlaps at least with particle"
              + std::to_string(c_j) + ". \n Change particle positions or switch off overlap check.");
        }
      }
    }
  }
}

/// It is checked if any node is moving too fast, in such a case a warning or
/// error message is printed.
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::checkParticles()
{
  const T speedOfSound = 1./3.; // fluid speed of sound
  for(unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    for(int n_i = 0; n_i < _particle[c_i]->numNodes; ++n_i) {
      const p_vec<> velocity = _particle[c_i]->vel[n_i];
      const p_vec<> position = _particle[c_i]->pos[n_i];
      if(velocity.length() > speedOfSound) {
        clout << "Velocity of node " << n_i << " in particle " << c_i << " is " << velocity << std::endl;
        clout << "  Position:          " << position << std::endl;
        clout << "  Total force:       " << _particle[c_i]->force[n_i] << std::endl;
        clout << "  Strain force:      " << _particle[c_i]->forceStrain[n_i] << std::endl;
        clout << "  Bending force:     " << _particle[c_i]->forceBending[n_i] << std::endl;
        clout << "  Volume force:      " << _particle[c_i]->forceVolume[n_i] << std::endl;
        clout << "  Surface force:     " << _particle[c_i]->forceSurface[n_i] << std::endl;
        clout << "  Interaction force: " << _particle[c_i]->forceInteraction[n_i] << std::endl;
        clout << "The Simulation is aborted; trying to dump data." << std::endl;
      }
    }
  }
}

/// The particle configurations are updated. The order of updates is crucial:
/// 1) Update node positions
/// 2) Update face areas and particle surface
/// 3) Update node weights
/// 4) Update particle centre
/// 5) Update particle volume
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::updateParticles()
{
  // Update particle configuration.
  for(unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    _particle[c_i]->updateNodes();
    _particle[c_i]->updateAreas();
    _particle[c_i]->updateWeights();
    _particle[c_i]->updateCentre();
    _particle[c_i]->updateVolume();
  }
  // Check particle configuration.
  checkParticles();
}

/// The particles are initialized and positioned at the beginning of the simulation.
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::initializeParticlesFromXML(std::string filename)
{
  clout << "Initialising Lagrangian particles from file \"" << filename << "\"." << std::endl;
  readPositionParameters(filename);

  int tempNodeOffset = 0;
  for(unsigned c_i = 0; c_i < _initParticles.size(); ++c_i) {
    _initParticles[c_i].radius_pos = _meshRadius[_meshForParticle[c_i]];
    _initParticles[c_i].nodeOffset = tempNodeOffset;
    _initParticles[c_i].meshType = _meshForParticle[c_i];
    tempNodeOffset += _meshes[_meshForParticle[c_i]]->numNodes;
    clout << "Particle " << c_i << " with mesh " << _meshForParticle[c_i] <<
      " at " << _initParticles[c_i].centre << "." << std::endl;
  }
  if(_checkInitialParticleOverlap) {
    checkParticleParticleOverlap();
  }

  // Initializing particles from initParticleStruct
  for(unsigned c_i = 0; c_i < _initParticles.size(); ++c_i){
    _particle.emplace_back(std::make_shared<MembraneParticle3D<T>>(_meshes[_initParticles[c_i].meshType]));
    _particle[c_i]->setRadius(_initParticles[c_i].radius_pos, true);
    _particle[c_i]->setCentre(_initParticles[c_i].centre);
    _particle[c_i]->nodeOffset = _initParticles[c_i].nodeOffset;
    _particle[c_i]->meshType = _initParticles[c_i].meshType;
    _particle[c_i]->rotate(_initParticles[c_i].angle, _initParticles[c_i].axis.x, _initParticles[c_i].axis.y, _initParticles[c_i].axis.z);
    _particle[c_i]->modVolume = _meshKV[_meshForParticle[c_i]];
    _particle[c_i]->modSphere = _meshKSph[_meshForParticle[c_i]];
    _particle[c_i]->modSurface = _meshKA[_meshForParticle[c_i]];
    _particle[c_i]->modAlpha = _meshKAlpha[_meshForParticle[c_i]];
    _particle[c_i]->modShear = _meshKS[_meshForParticle[c_i]];
    _particle[c_i]->modBending = _meshKB[_meshForParticle[c_i]];
  }

  // Find number of faces and nodes, and number of nodes in largest mesh.
  _numTotalNodes = 0;
  _numTotalFaces = 0;
  for(unsigned n = 0; n < _meshes.size(); ++n) {
    _numTotalNodes += _meshNumParticles[n] * _meshes[n]->numNodes;
    _numTotalFaces += _meshNumParticles[n] * _meshes[n]->numFaces;
  }

  // Update particles.
  updateParticles();
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::initializeData(LoadBalancer<T>& loadBalancer)
{
  // Initialize particle data required for ibm communication
  clout << "Initialize particle data for ibm communication" << std::endl;
  _data = std::make_shared<SuperD<T,DESCRIPTOR>>(loadBalancer);
  _data->getBlock(0).resize({(int)_numTotalNodes,1,1});
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::computeVolumeFraction()
{
  clout << "Computing particle volume fraction..." << std::endl;
  T volumeFraction = 0.0;
  for(unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    volumeFraction += _particle[c_i]->volumeRef;
  }
  clout << "Volume fraction in number cells: " << setprecision(5) << volumeFraction << std::endl;
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::computeForcesInternal()
{
  for (unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    _particle[c_i]->forcesCompute();
  }
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::combineForces()
{
  for (unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    for (int n_i = 0; n_i < _particle[c_i]->numNodes; ++n_i) {
      _particle[c_i]->force[n_i] = _particle[c_i]->forceStrain[n_i]
                                 + _particle[c_i]->forceBending[n_i]
                                 + _particle[c_i]->forceVolume[n_i]
                                 + _particle[c_i]->forceSurface[n_i]
                                 + _particle[c_i]->forceInteraction[n_i];
    }
  }
}

/// For efficient IBM communication, the data is pushed into an efficient
/// data container. Therefore, the data required for IBM-coupling is copied
/// in each step. At intial call, the velocity of the particles is set to zero.
/// Currently this is allowed as the initial particle velocity is always zero.
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::downloadData()
{
  for (unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    for (int n_i = 0; n_i < _particle[c_i]->numNodes; ++n_i) {
      unsigned i = c_i * _particle[c_i]->numNodes + n_i;
      auto view = _data->getBlock(0).get(i);
      auto tmp = view.template getField<fields::membrane::VELOCITY>();
      p_vec<> velD{tmp[0], tmp[1], tmp[2]};
      _particle[c_i]->vel[n_i] = velD;
    }
  }
}

/// For efficient IBM communication, the data is pushed into an efficient
/// data container. Therefore, the data required for IBM-coupling is copied
/// in each step.
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::uploadData()
{
  for (unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    for (int n_i = 0; n_i < _particle[c_i]->numNodes; ++n_i) {
      unsigned i = c_i * _particle[c_i]->numNodes + n_i;
      auto view = _data->getBlock(0).get(i);
      Vector<T,3> nodeCoord {_particle[c_i]->pos[n_i].x,
                             _particle[c_i]->pos[n_i].y,
                             _particle[c_i]->pos[n_i].z};
      view.template setField<fields::PHYS_R>(nodeCoord);
      Vector<T,3> nodeForce {_particle[c_i]->force[n_i].x,
                             _particle[c_i]->force[n_i].y,
                             _particle[c_i]->force[n_i].z};
      view.template setField<fields::membrane::FORCE>(nodeForce);
      Vector<T,3> nodeVel {_particle[c_i]->vel[n_i].x,
                           _particle[c_i]->vel[n_i].y,
                           _particle[c_i]->vel[n_i].z};
      view.template setField<fields::membrane::VELOCITY>(nodeVel);
    }
  }
}

/// Method executing necessary steps in each time step for FEM
template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::updateForcesFromDisplacements()
{
  // Pull data from container containing displacements from ibm communication
  downloadData();

  // Apply displacements of the membrane nodes
  updateParticles();

  // Compute forces via FEM
  computeForcesInternal();

  // Summarize all force contributions
  combineForces();

  // Copy data into container for ibm communication
  uploadData();
}

template<typename T, typename DESCRIPTOR>
void MembraneParticleSystem3D<T,DESCRIPTOR>::updateMomenta()
{
  for(unsigned c_i = 0; c_i < _particle.size(); ++c_i) {
    _particle[c_i]->updateMomenta();
  }
}

template<typename T, typename DESCRIPTOR>
std::shared_ptr<MembraneParticle3D<T>> MembraneParticleSystem3D<T,DESCRIPTOR>::get(int c_i)
{
  return _particle.at(c_i);
}

template<typename T, typename DESCRIPTOR>
SuperD<T,DESCRIPTOR>& MembraneParticleSystem3D<T,DESCRIPTOR>::getData()
{
  return *(_data);
}

template<typename T, typename DESCRIPTOR>
unsigned MembraneParticleSystem3D<T,DESCRIPTOR>::getNumberMeshes()
{
  return _meshes.size();
}

template<typename T, typename DESCRIPTOR>
unsigned MembraneParticleSystem3D<T,DESCRIPTOR>::getNumberParticles()
{
  return _particle.size();
}

template<typename T, typename DESCRIPTOR>
unsigned MembraneParticleSystem3D<T,DESCRIPTOR>::getTotalNumberNodes()
{
  return _numTotalNodes;
}


template<typename T, typename DESCRIPTOR>
unsigned MembraneParticleSystem3D<T,DESCRIPTOR>::getTotalNumberFaces()
{
  return _numTotalFaces;
}

}  // namespace membrane

}  // namespace olb

#endif
