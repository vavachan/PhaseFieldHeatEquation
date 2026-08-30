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

#ifndef IO_H
#define IO_H

//#define GSL  // required for computing the eigenvector of the deformation

#include <fstream>
#include <iomanip>
#include <vector>
#include "membrane_diskIO.h"
#include "utilities/calc.h"

#ifdef GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#endif

namespace olb {

namespace membrane {


/// Particles are written to VTP files.
/// The following data is included:
/// * Node locations and face connectivity
/// * Area deviations
/// * Maximum angles
/// * Particle index
/// * Process rank
template<typename T, typename DESCRIPTOR>
void writeParticleVTK(int iT, std::string directory, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane) {
  /// Check execution condition.
  if (sMembrane.getNumberParticles() == 0) {
    return;
  }
  /// Create and open file.
  std::stringstream outputFilename;
  outputFilename << directory << "/Particles_t" << iT << ".vtp";
  ofstream outputFile;
  outputFile.open(outputFilename.str().c_str());
  /// Write header.

  outputFile << "<?xml version=\"1.0\"?>\n";
  outputFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  outputFile << "<PolyData>\n";
  outputFile << "<Piece NumberOfPoints=\"" << sMembrane.getTotalNumberNodes()
             << "\" NumberOfPolys=\""
             << sMembrane.getTotalNumberFaces()<<"\">\n";

  outputFile << "<PointData>\n";
  outputFile << "</PointData>\n";
  //write area deviation
  outputFile << "<CellData Scalars=\"areaDeviation\">\n";
  outputFile << "<DataArray type=\"Float32\" Name=\"areaDeviation\" format=\"ascii\">\n";
  outputFile.unsetf(ios_base::fixed);
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    for (int f_i = 0; f_i < sMembrane.get(c_i)->numFaces; ++f_i) {
      outputFile << setprecision(3) << sMembrane.get(c_i)->area[f_i] / (sMembrane.get(c_i)->mesh->area[f_i] * util::pow(sMembrane.get(c_i)->radius, 2.)) - 1 << " ";
      if( (f_i+1)%6 == 0){
       outputFile <<"\n";
     }
    }
  }
  outputFile << "\n</DataArray>\n";
  //write max Angle
  outputFile <<"<DataArray type=\"Float32\" Name=\"maxAngle\" format=\"ascii\">\n";
  double maxAngle;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    for (int f_i = 0; f_i < sMembrane.get(c_i)->numFaces; ++f_i) {
      maxAngle = 0;
      for (int i = 0; i < 3; ++i) {
        if(sMembrane.get(c_i)->computeAngleNormals(f_i, sMembrane.get(c_i)->mesh->faceNeighbours[i][f_i]) > maxAngle) {
          maxAngle = sMembrane.get(c_i)->computeAngleNormals(f_i, sMembrane.get(c_i)->mesh->faceNeighbours[i][f_i]);

        }
      }
      outputFile << setprecision(3) << maxAngle * 180 / M_PI << " ";
      if( (f_i+1)%6 == 0){
       outputFile <<"\n";
     }
    }
  }
  outputFile << "\n</DataArray>\n";
  //write particle Index
  outputFile << "<DataArray type=\"Int32\" Name=\"particleIndex\" format=\"ascii\" >\n";
  outputFile.unsetf(ios_base::fixed);
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    for (int f_i = 0; f_i < sMembrane.get(c_i)->numFaces; ++f_i) {
      outputFile << c_i << " ";
      if( (f_i+1)%6 == 0){
       outputFile <<"\n";
     }
    }
  }
  outputFile << "\n</DataArray>\n";
  outputFile << "</CellData>\n";
  outputFile << "<Points>\n";
  outputFile << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    for (int n_i = 0; n_i < sMembrane.get(c_i)->numNodes; ++n_i) {
      outputFile << setprecision(6) << fixed << sMembrane.get(c_i)->pos[n_i] << " ";
      if((n_i+1)%2 == 0){
       outputFile <<"\n";
     }
    }
  }
  outputFile << "</DataArray>\n";
  outputFile << "</Points>\n";
  //write connectivity
  outputFile << "<Polys>\n";
  outputFile << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\"> \n";
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    for (int f_i = 0; f_i < sMembrane.get(c_i)->numFaces; ++f_i) {
      outputFile << sMembrane.get(c_i)->nodeOffset + sMembrane.get(c_i)->mesh->faceToNode[0][f_i] << " "
                 << sMembrane.get(c_i)->nodeOffset + sMembrane.get(c_i)->mesh->faceToNode[1][f_i] << " "
                 << sMembrane.get(c_i)->nodeOffset + sMembrane.get(c_i)->mesh->faceToNode[2][f_i] << " ";
      if( (f_i+1)%2 == 0){
       outputFile <<"\n";
     }
    }
  }
  outputFile << "</DataArray>\n";
  //write offsets. These are simply counting up to the number of faces multiplied by 3
  outputFile << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for(unsigned i= 1; i<sMembrane.getTotalNumberFaces() +1 ; i++ ){
    outputFile << i*3 <<" ";
    if( (i)%6 == 0 ){
      outputFile <<"\n";
    }
  }
  outputFile << "\n</DataArray>\n";
  outputFile << "</Polys>\n";
  outputFile << "</Piece>\n";
  outputFile << "</PolyData>\n";
  outputFile << "</VTKFile>\n";
  outputFile.close();
}

/// Individual particle statistics are written to the disk:
/// 1) centre position
/// 2) centre velocity
/// 3) angular momentum
/// 4) angular velocity
/// 5) total force
/// 6) total torque
/// 7) energy contributions
/// 8) volume and surface area
/// 9) minimum and maximum area deviations
/// 10) minimum and maximum edge lengths
/// 11) maximum normal-normal angle
template<typename T, typename DESCRIPTOR>
void writeParticleData(int iT, std::string directory, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane) {
  /// Check execution condition.
  if(sMembrane.getNumberParticles() == 0) {
    return;
  }
  /// Create files.
  const unsigned numParticles = sMembrane.getNumberParticles();
  std::stringstream *filename;
  filename = new std::stringstream[numParticles];
  int counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    filename[counter] << directory << "/Particle_" << c_i << ".dat";
    counter++;
  }
  /// Write headers.
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    std::stringstream header;
    header << "# individual data for particle " << c_i << " out of " << sMembrane.getNumberParticles() << "\n";
    header << "# column 1: lattice time\n";
    header << "# columns 2-4: centre position\n";
    header << "# columns 5-7: centre velocity\n";
    header << "# columns 8-10: angular momentum\n";
    header << "# columns 11-13: angular velocity\n";
    header << "# columns 14-16: total force\n";
    header << "# columns 17-19: total torque\n";
    header << "# column 20: total energy\n";
    header << "# column 21: strain energy\n";
    header << "# column 22: bending energy\n";
    header << "# column 23: surface energy\n";
    header << "# column 24: volume energy\n";
    header << "# columns 25-33: stresslet\n";
    header << "# column 34: volume\n";
    header << "# column 35: surface\n";
    header << "# column 36: minimum area deviation\n";
    header << "# column 37: maximum area deviation\n";
    header << "# column 38: minimum edge length\n";
    header << "# column 39: maximum edge length\n";
    header << "# column 40: maximum normal-normal angle\n";
    header << "iT x y z v_x v_y v_z L_x L_y L_z omega_x omega_y omega_z F_x F_y F_z T_x T_y T_z E_tot E_S E_B E_A E_V sigma_xx sigma_xy sigma_xz sigma_yx sigma_yy sigma_yz sigma_zx sigma_zy sigma_zz volume surface areaDevMin areaDevMax edgeMin edgeMax angleMax";
    DIO::writeHeader(filename[counter].str(), header.str());
    counter++;
  }
  /// Calculate data on the fly and write data.
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    /// Calculate particle stresslet.
    p_ten2 stresslet;
    stresslet.reset();
    for(int n_i = 0; n_i < sMembrane.get(c_i)->numNodes; ++n_i) {
      const p_vec<> force = sMembrane.get(c_i)->forceStrain[n_i] + sMembrane.get(c_i)->forceBending[n_i] + sMembrane.get(c_i)->forceVolume[n_i] + sMembrane.get(c_i)->forceSurface[n_i];
      const p_vec<> pos = sMembrane.get(c_i)->posOld[n_i];
      stresslet.xx -= (force.x * pos.x);
      stresslet.xy -= (force.x * pos.y);
      stresslet.xz -= (force.x * pos.z);
      stresslet.yx -= (force.y * pos.x);
      stresslet.yy -= (force.y * pos.y);
      stresslet.yz -= (force.y * pos.z);
      stresslet.zx -= (force.z * pos.x);
      stresslet.zy -= (force.z * pos.y);
      stresslet.zz -= (force.z * pos.z);
    }
    /// Calculate area deviations and maximum angle.
    double areaDev;
    double areaDevMax = 0.;
    double areaDevMin = 0.;
    double angle;
    double angleMax = 0.;
    for(int f_i = 0; f_i < sMembrane.get(c_i)->numFaces; ++f_i) {
      /// Find maximum and minimum area deviation.
      areaDev = sMembrane.get(c_i)->area[f_i] / (sMembrane.get(c_i)->mesh->area[f_i] * util::pow(sMembrane.get(c_i)->radius, 2.)) - 1.;
      if(areaDev > areaDevMax) {
        areaDevMax = areaDev;
      }
      if(areaDev < areaDevMin) {
        areaDevMin = areaDev;
      }
      /// Find maximum normal-normal angle.
      for(int i = 0; i < 3; ++i) {
        int f_j = sMembrane.get(c_i)->mesh->faceNeighbours[i][f_i];
        if(f_j > f_i) {
          angle = sMembrane.get(c_i)->computeAngleNormals(f_i, f_j);
          if(angle > angleMax) {
            angleMax = angle;
          }
        }
      }
    }
    /// Calculate minimum and maximum edge lengths.
    double edge;
    double edgeMin = sMembrane.get(c_i)->radius;
    double edgeMax = 0.;
    for(int n_i = 0; n_i < sMembrane.get(c_i)->numNodes; ++n_i) {
      for(int i = 0; i < MAX_NODE_NEIGHBORS; ++i) {
        int n_j = sMembrane.get(c_i)->mesh->nodeNeighbours[i][n_i];
        if(n_j > n_i) {
          edge = (sMembrane.get(c_i)->pos[n_i] - sMembrane.get(c_i)->pos[n_j]).length();
          if(edge > edgeMax) {
            edgeMax = edge;
          }
          if(edge < edgeMin) {
            edgeMin = edge;
          }
        }
      }
    }
    /// Append line to data file.
    std::stringstream line;
    line << iT << " ";
    line << setprecision(15) << sMembrane.get(c_i)->centre.x << " " << sMembrane.get(c_i)->centre.y << " " << sMembrane.get(c_i)->centre.z << " ";
    line << setprecision(4) << sMembrane.get(c_i)->velocity << " ";
    line << sMembrane.get(c_i)->angularMomentum << " ";
    line << sMembrane.get(c_i)->computeAngularVelocity() << " ";
    line << sMembrane.get(c_i)->Force << " ";
    line << sMembrane.get(c_i)->torque << " ";
    line << sMembrane.get(c_i)->ergStrain + sMembrane.get(c_i)->ergBending + sMembrane.get(c_i)->ergSurface + sMembrane.get(c_i)->ergVolume << " ";
    line << sMembrane.get(c_i)->ergStrain << " ";
    line << sMembrane.get(c_i)->ergBending << " ";
    line << sMembrane.get(c_i)->ergSurface << " ";
    line << sMembrane.get(c_i)->ergVolume << " ";
    line << stresslet << " ";
    line << sMembrane.get(c_i)->volume << " ";
    line << sMembrane.get(c_i)->surface << " ";
    line << areaDevMin << " ";
    line << areaDevMax << " ";
    line << edgeMin << " ";
    line << edgeMax << " ";
    line << angleMax * 180. / M_PI << " ";
    DIO::appendLine(filename[counter].str(), line.str());
    counter++;
  }
}



template<typename T, typename DESCRIPTOR>
void createMasterFile(OstreamManager clout, std::string& _name, int iT, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane)
{
  std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                  + createFileName(_name) + "_master.pvd";
  preamblePVD( fullNamePVDmaster, sMembrane);
  closePVD( fullNamePVDmaster, sMembrane);

}

template<typename T, typename DESCRIPTOR>
void preamblePVD( const std::string& fullNamePVD, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}

template<typename T, typename DESCRIPTOR>
void closePVD( const std::string& fullNamePVD, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);

  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, typename DESCRIPTOR>
void dataPVDmaster(OstreamManager clout, int iT,
    const std::string& fullNamePVDMaster, const std::string& namePiece, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane)
{
  std::ofstream fout(fullNamePVDMaster.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD
    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\"\" file=\"" << namePiece << "\"/>\n";
    fout.close();
  }
  else{
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
  closePVD( fullNamePVDMaster, sMembrane);
}





/// The particle inertia tensors and principal orientations are computed
/// using GNU Scientific Library (GSL).
/// Results are written in two directories:
/// 1) for data analysis (Particles)
/// 2) for visualization (VTKParticles)
#ifdef GSL
template<typename T, typename DESCRIPTOR>
void writeParticleInertia(int iT, std::string directoryVTK, std::string directoryData, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane) {
  /// Check execution condition.
  if(sMembrane.getNumberParticles() == 0) {
    return;
  }
  /// Create filenames.
  const unsigned numParticles = sMembrane.getNumberParticles();
  std::stringstream *filename;
  filename = new std::stringstream[numParticles];
  int counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    filename[counter] << directoryData << "/Axes_" << c_i << ".dat";
    counter++;
  }
  /// Write headers.
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    std::stringstream header;
    header << "# orientation data for particle " << c_i << " out of " << sMembrane.getNumberParticles() << "\n";
    header << "# column 1: lattice time\n";
    header << "# columns 2-4: principal axes\n";
    header << "# columns 5-7: eigenvector 0\n";
    header << "# columns 8-10: eigenvector 1\n";
    header << "# columns 11-13: eigenvector 2\n";
    header << "iT a b c eigvec0x eigvec0y eigvec0z eigvec1x eigvec1y eigvec1z eigvec2x eigvec2y eigvec2z";
    DIO::writeHeader(filename[counter].str(), header.str());
    counter++;
  }
  /// Compute inertia tensor.
  p_ten2 *inertiaTensor;
  inertiaTensor = new p_ten2[numParticles];
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    inertiaTensor[counter] = sMembrane.get(c_i)->computeInertiaTensorVolume();
    counter++;
  }
  /// Diagonalise inertia tensors using GSL.
  double *eigval0;
  double *eigval1;
  double *eigval2;
  p_vec<> *eigvec0;
  p_vec<> *eigvec1;
  p_vec<> *eigvec2;
  double *a;
  double *b;
  double *c;
  eigval0 = new double[numParticles];
  eigval1 = new double[numParticles];
  eigval2 = new double[numParticles];
  eigvec0 = new p_vec<>[numParticles];
  eigvec1 = new p_vec<>[numParticles];
  eigvec2 = new p_vec<>[numParticles];
  a = new double[numParticles];
  b = new double[numParticles];
  c = new double[numParticles];
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    /// Translate tensor to plain array for GSL, using tensor symmetry.
    double inertiaTensorGSL[] = {inertiaTensor[counter].xx, inertiaTensor[counter].xy, inertiaTensor[counter].xz,
                                 inertiaTensor[counter].xy, inertiaTensor[counter].yy, inertiaTensor[counter].yz,
                                 inertiaTensor[counter].xz, inertiaTensor[counter].yz, inertiaTensor[counter].zz};
    /// Set up GSL for tensor diagonalisation.
    gsl_matrix_view m = gsl_matrix_view_array(inertiaTensorGSL, 3, 3); // create readable matrix
    gsl_vector *eval = gsl_vector_alloc(3); // create vectors for eigenvalues
    gsl_matrix *evec = gsl_matrix_alloc(3, 3); // create matrix for eigenvectors
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3); // allocate workspace for symmetric tensors
    gsl_eigen_symmv(&m.matrix, eval, evec, w); // calculate eigenvalues and eigenvectors
    gsl_eigen_symmv_free(w); // free memory
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); // sort eigenvalues by ascending magnitude
    /// Pass eigenvalues and eigenvectors to useable variables.
    /// The eigenvalue with the smallest index has the smallest magnitude (eigval0).
    /// The eigenvectors denote the principal directions of the inertia ellipsoid in space.
    eigval0[counter] = gsl_vector_get(eval, 0); // smallest
    eigval1[counter] = gsl_vector_get(eval, 1); // intermediate
    eigval2[counter] = gsl_vector_get(eval, 2); // largest
    eigvec0[counter].set_to(gsl_matrix_get(evec, 0, 0), gsl_matrix_get(evec, 1, 0), gsl_matrix_get(evec, 2, 0));
    eigvec1[counter].set_to(gsl_matrix_get(evec, 0, 1), gsl_matrix_get(evec, 1, 1), gsl_matrix_get(evec, 2, 1));
    eigvec2[counter].set_to(gsl_matrix_get(evec, 0, 2), gsl_matrix_get(evec, 1, 2), gsl_matrix_get(evec, 2, 2));
    /// Free memory allocated by GSL.
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    /// Compute effective radii from eigenvalues, assuming ellipsoidal shapes.
    /// 'a' is largest radius, 'c' is smallest.
    /// The third eigenvector (eigvec2) corresponds to the direction of the largest momentum of inertia.
    a[counter] = sqrt(2.5 * (eigval1[counter] + eigval2[counter] - eigval0[counter]) / sMembrane.get(c_i)->volume);
    b[counter] = sqrt(2.5 * (eigval0[counter] + eigval2[counter] - eigval1[counter]) / sMembrane.get(c_i)->volume);
    c[counter] = sqrt(2.5 * (eigval0[counter] + eigval1[counter] - eigval2[counter]) / sMembrane.get(c_i)->volume);
    /// Correct orientation of eigenvectors relative to absolute particle orientation.
    /// The inertia tensor offers no access to the absolute orientation of the particle,
    /// i.e., eigenvectors are defined only up to their signs.
    /// Using a mesh-fixed quantity, one can correct the orientation in such a way that the sign of
    /// the eigenvectors is defined as well.
    if(eigvec2[counter] * sMembrane.get(c_i)->computeOrientationVector() < 0.) {
      eigvec2[counter] *= -1;
    }
    counter++;
  }
  /// Write data to file.
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    std::stringstream line;
    line << iT << " "
         << a[counter] << " " << b[counter] << " " << c[counter] << " "
         << eigvec0[counter] << " " << eigvec1[counter] << " " << eigvec2[counter];
    DIO::appendLine(filename[counter].str(), line.str());
    counter++;
  }
  /// Write VTK files.
  /// Eigenvectors are scaled by principal radii (a, b, c) in order to
  /// correctly visualise the inertia ellipsoid.
  std::stringstream outputFilename;
  outputFilename << directoryVTK << "/Axes_t" << iT << ".vtk";
  ofstream outputFile;
  outputFile.open(outputFilename.str().c_str());
  /// Write header.
  outputFile << "# vtk DataFile Version 3.0\n";
  outputFile << "Axes\n";
  outputFile << "ASCII\n";
  outputFile << "DATASET POLYDATA\n";
  /// Write Points.
  outputFile << "POINTS " << 4 * numParticles << " float\n";
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    outputFile << sMembrane.get(c_i)->centre << "\n";
    outputFile << sMembrane.get(c_i)->centre + eigvec0[counter] * a[counter] << "\n";
    outputFile << sMembrane.get(c_i)->centre + eigvec1[counter] * b[counter] << "\n";
    outputFile << sMembrane.get(c_i)->centre + eigvec2[counter] * c[counter] << "\n";
    counter++;
  }
  /// Write lines.
  outputFile << "LINES " << 3 * numParticles << " " << numParticles * 9 << "\n";
  counter = 0;
  for (unsigned c_i = 0; c_i < sMembrane.getNumberParticles(); ++c_i) {
    outputFile << "2" << " " << counter * 4 << " " << counter * 4 + 1 << "\n";
    outputFile << "2" << " " << counter * 4 << " " << counter * 4 + 2 << "\n";
    outputFile << "2" << " " << counter * 4 << " " << counter * 4 + 3 << "\n";
    counter++;
  }
  outputFile.close();
}
#endif

/// Call all disk writing functions.
template<typename T, typename DESCRIPTOR>
void writeData(int iT, MembraneParticleSystem3D<T,DESCRIPTOR>& sMembrane) {
  // Create directory.
  DIO::createDirectory("tmp/vtkData/VTKParticles");
  DIO::createDirectory("tmp/particleData");
  DIO::createDirectory("tmp/particleData/Particles");
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  string name="Particles";
  if(rank== 0){
    OstreamManager clout(std::cout, name );
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                  + createFileName(name) + "_master.pvd";

    std::string namePiece =  "VTKparticles/" + name+ "_t"+ std::to_string(iT) + ".vtp";
    if(iT==0){
      createMasterFile(clout, name, iT, sMembrane);
      dataPVDmaster(clout, iT, fullNamePVDmaster, namePiece, sMembrane );
     }
     else{
      dataPVDmaster(clout, iT, fullNamePVDmaster, namePiece, sMembrane );
    }

  }

  writeParticleVTK(iT, "tmp/vtkData/VTKParticles", sMembrane);
  writeParticleData(iT, "tmp/particleData/Particles", sMembrane);
#ifdef GSL
  writeParticleInertia(iT, "tmp/particleData/VTKParticles", "tmp/particleData/Particles", sMembrane);
#endif
}



}

}

#endif
