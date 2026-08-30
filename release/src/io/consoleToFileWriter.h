/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

#ifndef CONSOLETOFILEWRITER_H
#define CONSOLETOFILEWRITER_H

#include <iostream>
#include <fstream>
#include <streambuf>

namespace olb {

class DoubleBuffer : public std::streambuf {
public:
  DoubleBuffer(std::streambuf* consoleBuf, std::streambuf* fileBuf) : consoleBuffer(consoleBuf), fileBuffer(fileBuf) {}

protected:
  // This function is called for every character written to the stream
  virtual int overflow(int c) override {
    if (c != EOF) {
      // Write the character to both buffers
      if (consoleBuffer->sputc(c) == EOF || fileBuffer->sputc(c) == EOF) {
        return EOF;
      }
    }
    return c;
  }

  // Synchronize both buffers
  virtual int sync() override {
    if (consoleBuffer->pubsync() == 0 && fileBuffer->pubsync() == 0) {
      return 0;
    }
    return -1;
  }

private:
  std::streambuf* consoleBuffer; // Original console buffer
  std::streambuf* fileBuffer;    // File buffer

};

}

#endif
