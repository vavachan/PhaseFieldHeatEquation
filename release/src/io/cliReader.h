/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#ifndef IO_CLI_READER_H_
#define IO_CLI_READER_H_

#include <string>
#include <vector>
#include <charconv>

namespace olb {

/// Very simple CLI argument parser
class CLIreader {
private:
  std::vector<std::string> _tokens;

public:
  CLIreader(int& argc, char** argv) {
    for (int i=1; i < argc; ++i) {
      _tokens.emplace_back(argv[i]);
    }
  }

  /// Returns true iff name is specified
  bool contains(const std::string& name) const {
    return std::find(_tokens.begin(), _tokens.end(), name) != _tokens.end();
  }

  /// Returns value of token after name (i.e. by convention the value assigned to name)
  std::string operator[](const std::string& name) const {
    auto iter = std::find(_tokens.begin(), _tokens.end(), name);
    if (iter != _tokens.end() && ++iter != _tokens.end()){
      return *iter;
    } else {
      return std::string{};
    }
  }

  /// Return value of name as TYPE or fallback if not provided
  template <typename TYPE>
  TYPE getValueOrFallback(const std::string& name, TYPE fallback) const {
    if (contains(name)) {
      const std::string str = operator[](name);
      TYPE value{};
      // Handling required due to inconsistent support for std::from_chars
             if constexpr (std::is_integral_v<TYPE>) {
        auto result = std::from_chars(str.data(),
                                      str.data() + str.size(),
                                      value);
        if (result.ec != std::errc()) {
          return fallback;
        }
        return value;
      } else if constexpr (std::is_floating_point_v<TYPE>) {
        try {
          std::size_t nProcessed{};
          value = static_cast<TYPE>(std::stod(str, &nProcessed));
          if (nProcessed != str.size()) {
            return fallback;
          }
          return value;
        } catch (...) {
          return fallback;
        }
      } else {
        static_assert([](){ return false; }, "Unsupported type");
      }
    } else {
      return fallback;
    }
  }

};

template <>
std::string CLIreader::getValueOrFallback<std::string>(
  const std::string& name, std::string fallback) const {
  if (contains(name)) {
    return operator[](name);
  } else {
    return fallback;
  }
}

}

#endif
