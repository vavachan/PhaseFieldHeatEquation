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

#ifndef FILES_IO_H
#define FILES_IO_H

#include "filesIO.h"
#include <cstdlib> // For std::system

bool directoryExists(const std::string& path)
{
  struct stat statbuf;
  return (stat(path.c_str(), &statbuf) == 0 && S_ISDIR(statbuf.st_mode));
}

bool createDirectory(const std::string& path)
{
  return (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0);
}

bool deleteDirectory(const std::string& path)
{
  std::string command = "rm -rf " + path;
  return (std::system(command.c_str()) == 0);
}

bool fileExists(const std::string& name)
{
  std::ifstream f(name.c_str());
  return f.good();
}

template <typename T>
void saveVector1D(const std::string& filePath, const std::vector<T>& vec)
{
  std::ofstream out(filePath, std::ios::binary);
  if (!out.is_open())
    throw std::runtime_error("Cannot open file for writing: " + filePath);

  std::size_t size = vec.size();
  out.write(reinterpret_cast<const char*>(&size), sizeof(size));
  out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(T));
}

template <typename T>
void readVector1D(const std::string& filePath, std::vector<T>& vec)
{
  std::ifstream in(filePath, std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("Cannot open file for reading: " + filePath);

  std::size_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(size));
  vec.resize(size);
  in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(T));
}

template <typename T>
void saveVector3D(const std::string& filePath, const std::vector<std::vector<std::vector<T>>>& vec)
{
  std::ofstream out(filePath, std::ios::binary);
  if (!out.is_open())
    throw std::runtime_error("Cannot open file for writing: " + filePath);

  std::size_t outerSize = vec.size();
  out.write(reinterpret_cast<const char*>(&outerSize), sizeof(outerSize));

  for (const auto& midVec : vec) {
    std::size_t midSize = midVec.size();
    out.write(reinterpret_cast<const char*>(&midSize), sizeof(midSize));

    for (const auto& innerVec : midVec) {
      std::size_t innerSize = innerVec.size();
      out.write(reinterpret_cast<const char*>(&innerSize), sizeof(innerSize));
      out.write(reinterpret_cast<const char*>(innerVec.data()), innerSize * sizeof(T));
    }
  }
}

template <typename T>
void readVector3D(const std::string& filePath, std::vector<std::vector<std::vector<T>>>& vec)
{
  std::ifstream in(filePath, std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("Cannot open file for reading: " + filePath);

  std::size_t outerSize;
  in.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));
  vec.resize(outerSize);

  for (auto& midVec : vec) {
    std::size_t midSize;
    in.read(reinterpret_cast<char*>(&midSize), sizeof(midSize));
    midVec.resize(midSize);

    for (auto& innerVec : midVec) {
      std::size_t innerSize;
      in.read(reinterpret_cast<char*>(&innerSize), sizeof(innerSize));
      innerVec.resize(innerSize);
      in.read(reinterpret_cast<char*>(innerVec.data()), innerSize * sizeof(T));
    }
  }
}

#endif // FILES_IO_hh