# This file is part of the OpenLB library
#
# Copyright (C) 2021 Adrian Kummerlaender
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

from sympy import *

from sympy.codegen.ast import CodeBlock, Assignment, Return
from sympy.printing.c import C11CodePrinter
from sympy.codegen.rewriting import ReplaceOptim
from sympy.simplify import cse_main
from sympy.utilities.iterables import numbered_symbols
from sympy.core.function import UndefinedFunction

import itertools

class SymbolGenerator:
    def __init__(self, prefix):
        self.prefix = prefix
        self.symbols = set()

    def __iter__(self):
        self.generator = numbered_symbols(cls=Symbol, prefix=self.prefix)
        return self

    def __next__(self):
        symbol = next(self.generator)
        self.symbols.add(symbol)
        return symbol

class CodeBlockPrinter(C11CodePrinter):
    def __init__(self, subexprs):
        super(CodeBlockPrinter, self).__init__()
        self.known_functions['Abs'] = 'util::fabs'
        self.known_functions['exp'] = 'util::exp'
        self.known_functions['Pow'] = 'util::pow'
        self.subexprs = subexprs

    def _print_Indexed(self, expr):
        if len(expr.indices) == 1:
            return f"{expr.base.name}[{expr.indices[0]}]"
        elif len(expr.indices) == 2:
            return f"{expr.base.name}[{expr.indices[0]}][{expr.indices[1]}]"

    def _print_Float(self, flt):
        return "V{%.15g}" % flt.evalf()

    def _print_Pow(self, expr):
        if expr.exp == -0.5:
            return "V{1} / util::sqrt(%s)" % self.doprint(expr.base)
        elif expr.exp == -1:
            return "V{1} / (%s)" % self.doprint(expr.base)
        elif expr.exp == -2:
            x = self.doprint(expr.base)
            return "V{1} / ((%s)*(%s))" % (x, x)
        elif expr.exp == 0.5:
            return "util::sqrt(%s)" % self.doprint(expr.base)
        elif expr.exp == 2:
            x = self.doprint(expr.base)
            return "((%s)*(%s))" % (x, x)
        else:
            return "util::pow(%s, %s)" % (self.doprint(expr.base), self.doprint(expr.exp))

    def _print_Assignment(self, expr):
        if isinstance(expr.rhs, Symbol):
            if "DECLARE:" in expr.rhs.name:
                return expr.rhs.name.replace("DECLARE:", "")

        if hasattr(expr.rhs, 'func'):
            if isinstance(expr.rhs.func, UndefinedFunction):
                call = expr.rhs.name
                if "defineRho" in call or "defineU" in call:
                    return f"{self.doprint(expr.rhs)}"
                elif "computeRho" in call:
                    return f"V {expr.lhs} = {call}()"
                elif "computeJ" in call or "computeU" in call:
                    decl = Assignment(expr.lhs,Symbol("DECLARE:V "+expr.lhs.name+" [DESCRIPTOR::d]"))
                    return f"{self.doprint(decl)}; {call}({expr.lhs})"
                elif "computeStress" in call:
                    decl = Assignment(expr.lhs,Symbol("DECLARE:V "+expr.lhs.name+" [util::TensorVal<CELL::descriptor_t>::n]"))
                    args = expr.rhs.args + (expr.lhs,)
                    return f"{self.doprint(decl)}; {call}{args}"
                elif "computeEquilibrium" in call:
                    decl = Assignment(expr.lhs,Symbol("DECLARE:V "+expr.lhs.name+" [DESCRIPTOR::q]"))
                    args = expr.rhs.args + (expr.lhs,)
                    return f"{self.doprint(decl)}; {call}{args}"
                elif "FILL_ARRAY_ELEMENT" in call:
                    return f"{self.doprint(Assignment(expr.lhs,expr.rhs.args[1]))}"
                elif "DECLARE_WITH_ARGUMENT:" in call:
                    lhs = expr.rhs.name.replace("DECLARE_WITH_ARGUMENT:", "")
                    return f"{lhs} = {expr.rhs.args[0]}"

        if expr.lhs in self.subexprs:
            return f"auto {self.doprint(expr.lhs)} = {self.doprint(expr.rhs.evalf())};"
        else:
            return f"{self.doprint(expr.lhs)} = {self.doprint(expr.rhs.evalf())};"

    def _print_Function(self, expr):
        return "%s(%s)" % (expr.name, self.stringify(expr.args, ", "))

def code(expr, symbols={ }):
    return CodeBlockPrinter(symbols).doprint(expr)

def cse(block, generator):
    expand_pos_square = ReplaceOptim(
        lambda e: e.is_Pow and e.exp.is_integer and e.exp == 2,
        lambda p: UnevaluatedExpr(Mul(p.base, p.base, evaluate = False)),
    )
    custom_opti = cse_main.basic_optimizations + [
        (expand_pos_square, expand_pos_square)
    ]
    return block.cse(symbols=generator, optimizations=custom_opti, order='none')

