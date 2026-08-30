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

#ifndef DISK_IO_HH
#define DISK_IO_HH
#include <cstdlib>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include "membrane_diskIO.h"

namespace olb {

namespace membrane {

namespace DIO {

/// Check whether specified directory exists.
/// Rreturns true if it does, false elsewise.
bool directoryExists(std::string directory) {
  bool exists = false;
  DIR *p_dir;
  p_dir = opendir(directory.c_str());
  if(p_dir != NULL) {
    exists = true;
    closedir(p_dir);
  }
  return exists;
}

/// Checks whether specified file exists.
/// Returns true if it does, false elsewise.
bool fileExists(std::string filename) {
  bool exists = false;
  struct stat file_info;
  int int_stat;
  int_stat = stat(filename.c_str(), &file_info);
  if(int_stat == 0) {
    exists = true;
  }
  return exists;
}

/// Specified directory is created if it does not already exist.
/// Only root can create a directory.
void createDirectory(std::string directory) {
  if(!directoryExists(directory.c_str())) {
    std::stringstream message;
    message << "mkdir " << directory.c_str();
    int ret = system(message.str().c_str());
    if (ret != 0) {
      throw std::runtime_error("System call failed with return code: " + std::to_string(ret));
    }
  }
}

/// Specified file is created if it does not already exist.
/// Only root can create a file.
void createFile(std::string filename) {
  if(!fileExists(filename.c_str())) {
    ofstream file(filename.c_str());
    file.close();
  }
}

/// Header line is written to specified file.
/// This is only done if file does not already exist.
void writeHeader(std::string filename, std::string header) {
  if(!fileExists(filename.c_str())) {
    createFile(filename);
    ofstream file(filename.c_str());
    file << header.c_str() << endl;
    file.close();
  }
}

/// Line is appended to specified file.
/// This is only done if file already exists.
void appendLine(std::string filename, std::string line) {
  if(fileExists(filename.c_str())) {
    ofstream file(filename.c_str(), fstream::app);
    file << line.c_str() << endl;
    file.close();
  }
}

}

}

}

#endif
