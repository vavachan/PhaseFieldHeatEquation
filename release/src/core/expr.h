/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_EXPR_H
#define CORE_EXPR_H

#include <string>
#include <sstream>
#include <memory>
#include <variant>
#include <functional>
#include <numeric>

#include "core/meta.h"

namespace olb {

/// Basic value-substitute enabling extraction of expression trees for code generation
class Expr : public ExprBase {
public:
  enum struct Op : std::size_t {
    Add, Mul, Sub, Div, Sqrt, Abs, Pow, Exp,
    End // keep last
  };

private:
  static std::array<std::size_t, static_cast<std::size_t>(Op::End)> _stats;

  struct Symbol {
    std::string name;

    void describe(std::stringstream& out) const;
    std::size_t size() const;
  };

  struct Constant {
  private:
    double value;

  public:
    Constant();
    Constant(double x);

    void describe(std::stringstream& out) const;
    std::size_t size() const;
  };

  class Binary {
  private:
    std::shared_ptr<Expr> _lhs;
    Op op;
    std::shared_ptr<Expr> _rhs;

  public:
    Binary(Expr lhs, Op op, Expr rhs);

    void describe(std::stringstream& out) const;
    std::size_t size() const;
  };

  class Unary {
  private:
    Op op;
    std::shared_ptr<Expr> _arg;

  public:
    Unary(Op op, Expr arg);

    void describe(std::stringstream& out) const;
    std::size_t size() const;
  };

  std::variant<Constant,Symbol,Binary,Unary> _payload;

public:
  Expr();
  Expr(Expr&& rhs);
  Expr(const Expr& rhs);
  Expr(double v);
  Expr(std::string name);
  Expr(Expr lhs, Op op, Expr rhs);
  Expr(Op op, Expr rhs);

  Expr& operator=(const Expr& rhs);

  Expr& operator+=(Expr rhs);
  Expr& operator-=(Expr rhs);
  Expr& operator*=(Expr rhs);
  Expr& operator/=(Expr rhs);

  /// Writes the serialized expression tree to out
  void describe(std::stringstream& out) const;
  /// Returns the serialized expression tree
  std::string describe() const;

  /// Helper function to quickly check if Expr is a specific symbol
  bool isSymbol(std::string name) const;

  /// Returns expanded size of tree
  std::size_t size() const;

  /// Resets FLOP counter
  static void reset() {
    _stats.fill(0);
  }

  /// Increment FLOP counter
  static void increment(Op op) {
    ++_stats[static_cast<std::size_t>(op)];
  }

  /// Return accumulated FLOPs
  static std::size_t count() {
    return std::accumulate(_stats.cbegin(), _stats.cend(), 0);
  }

  /// Return accumulated FLOPs of type
  static std::size_t count(Op op) {
    return _stats[static_cast<std::size_t>(op)];
  }

};

Expr operator+(Expr lhs, Expr rhs);
Expr operator-(Expr lhs, Expr rhs);
Expr operator*(Expr lhs, Expr rhs);
Expr operator/(Expr lhs, Expr rhs);
Expr operator-(Expr rhs);

Expr operator%(Expr lhs, int rhs);

bool operator==(const Expr& rhs, const Expr& lhs);
bool operator!=(const Expr& rhs, const Expr& lhs);
bool operator>(const Expr& rhs, const Expr& lhs);
bool operator<(const Expr& rhs, const Expr& lhs);
bool operator>=(const Expr& rhs, const Expr& lhs);
bool operator<=(const Expr& rhs, const Expr& lhs);

namespace util {

Expr sqrt(Expr x);
Expr fabs(Expr x);
Expr pow(Expr base, Expr exp);
Expr exp(Expr x);


}

}

#endif
