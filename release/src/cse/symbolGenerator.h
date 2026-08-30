/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Shota Ito
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

#ifndef CSE_SYMBOL_GENERATOR_H
#define CSE_SYMBOL_GENERATOR_H

namespace olb {

namespace cse {

/// Generator for placeholder expressions
class SymbolGenerator {
private:
  std::string _symbol;
  unsigned int _count = 0;
  std::vector<std::string> _generated;
public:
  explicit SymbolGenerator(std::string symbol = "x") :
    _symbol(symbol)
  { };

  std::string next() {
    unsigned int tmp = _count;
    ++_count;
    std::string next = _symbol + std::to_string(tmp);
    _generated.push_back(next);
    return next;
  }

  std::string symbol() {
    return _symbol;
  }

  std::vector<std::string> list() {
    return _generated;
  }

  auto find_symbol(std::string symbol) {
    return std::find_if(_generated.begin(), _generated.end(),
                        [&symbol](const std::string& s){
                          return s == symbol;
                        });
  }

  bool contains_symbol(std::string symbol) {
    auto it = find_symbol(symbol);
    return it != _generated.end();
  }

  void remove_symbol_from_list(std::string symbol) {
    auto it = find_symbol(symbol);
    if (it != _generated.end()) {
      _generated.erase(it);
    } else {
      throw std::runtime_error("Attempt to remove non-existing symbol.");
    }
  }
};

/*
 *  Helper functions for string manipulation to provide expressions
 *  in the form SymPy can process it to perform common-subexpression
 *  elimination.
 */

/// Print as symbols
std::string symbol(const Expr& symbol) {
  return "Symbol(\"" + symbol.describe() + "\")";
}

/// Print as symbols
std::string symbol(const std::string& symbol) {
  if (symbol.find("Symbol") != std::string::npos) {
    return symbol;
  }
  return "Symbol(\"" + symbol + "\")";
}

std::string strip_symbol(const std::string& symbol) {
  std::string_view raw{symbol};
  raw = std::string_view(raw.cbegin() + 8,
                         raw.cend() - 2);
  return std::string(raw);
}

std::string extract_c_array_name(const std::string& input) {
  size_t pos = input.rfind('[');
  if (pos != std::string::npos) {
    return input.substr(0, pos);
  }
  return input;
}

/// Print as assignments
std::string assignment(const std::string& lhs, const std::string& rhs) {
  return "Assignment(" + lhs + "," + rhs + ")";
}

/// Print as functions
std::string function(const std::string& body, const std::string& arg) {
  return "Function(" + body + ")(" + arg + ")";
}

/// Unwrap FIELD type from FieldArray
std::string remove_array_from_name(const std::string& arg) {
  std::string_view raw{arg};
  raw = std::string_view(raw.cbegin() + 6,
                         raw.cend() - 1);
  return std::string(raw);
}

}

}

#endif
