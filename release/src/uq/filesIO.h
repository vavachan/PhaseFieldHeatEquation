/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
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

#pragma once

#include <unistd.h>
#include <limits.h>
#include <libgen.h>   // For dirname
#include <sys/stat.h> // For stat, mkdir
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <memory>
#include <algorithm>
#include <numeric> // for std::inner_product

// #include "matrix.h"

// Utility function to check if a directory exists
bool directoryExists(const std::string& path);

// Utility function to create a directory
bool createDirectory(const std::string& path);

// Utility function to delete a directory
bool deleteDirectory(const std::string& path);

// Utility function to check if a file exists
bool fileExists(const std::string& name);

// Utility function to save a 1D vector to a binary file
template <typename T>
void saveVector1D(const std::string& filePath, const std::vector<T>& vec);

// Utility function to read a 1D vector from a binary file
template <typename T>
void readVector1D(const std::string& filePath, std::vector<T>& vec);

// Utility function to save a 3D vector to a binary file
template <typename T>
void saveVector3D(const std::string& filePath, const std::vector<std::vector<std::vector<T>>>& vec);

// Utility function to read a 3D vector from a binary file
template <typename T>
void readVector3D(const std::string& filePath, std::vector<std::vector<std::vector<T>>>& vec);
