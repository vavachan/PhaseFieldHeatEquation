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

#ifndef DISK_IO_H
#define DISK_IO_H

#include <sstream>

namespace olb {

namespace membrane {

namespace DIO {
  bool directoryExists(std::string); // Check existence of directory
  bool fileExists(std::string); // Check existence of file
  void createDirectory(std::string); // Create directory
  void createFile(std::string); // Create file
  void writeHeader(std::string, std::string); // Write header to file
  void appendLine(std::string, std::string); // Append line to file
}

}

}

#endif
