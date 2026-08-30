# This file is part of the OpenLB library
#
# Copyright (C) 2024 Shota Ito, Adrian Kummerlaender
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
from sympy.codegen.ast import CodeBlock, Assignment
from sympy.utilities.iterables import numbered_symbols
from source.cse_utils import SymbolGenerator, CodeBlockPrinter, code, cse
import re

def optimize_dynamics(inputFile):
    # Load expression tree from ouput file
    with open(inputFile, "r") as file:
        load_expressions = file.read()
    results = { }
    exec(load_expressions, globals(), results)

    # Operator information
    collisionO = results['collisionO']
    # Assignments
    cell_assignments = results['cell_assignments']
    return_assignments = results['return_assignments']
    fields_assignments = results['fields_assignments']
    # Symbols
    params_symbols = results['params_symbols']

    # Retrieve dynamics tuple information
    collisionO = collisionO.replace("Expr", "T")
    collisionO = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', collisionO)

    # Combine all assignments
    assignments = cell_assignments
    optional_symbols = params_symbols

    # Current limitations of sympy assignments
    generator = iter(SymbolGenerator('x'))
    cell_backsubstitutions = [ ]
    sub_assignments = [ ]
    for assignment in assignments:
        alias = next(generator)
        cell_backsubstitutions.append((alias, assignment.lhs))
        sub_assignments.append(Assignment(alias, assignment.rhs))

    # Creating aliases and assignments for parameter calls
    optional_assignments = [ ]
    optional_substitutions = [ ]
    for symbol in optional_symbols:
        alias = next(generator)
        optional_assignments.append(Assignment(alias, symbol))
        optional_substitutions.append((symbol, alias))

    # Create code block to perform cse
    block = CodeBlock(*sub_assignments,*fields_assignments,*return_assignments)
    block = block.subs(optional_substitutions)

    # Perform cse
    block = cse(block, generator)

    # These are required due to current CSE limititations
    substitutions = [ ]
    post_collision = [ ]
    simple_post_collision = True
    for assignment in block.args:
        if assignment.lhs in cell_backsubstitutions[0]:
            atoms = assignment.rhs.free_symbols
            for (lhs, rhs) in cell_backsubstitutions:
                if lhs in atoms or rhs in atoms:
                    simple_post_collision = False
    if simple_post_collision:
        substitutions = cell_backsubstitutions
    else:
        for subs in cell_backsubstitutions:
            post_collision.append(Assignment(subs[1], subs[0]))

    # Apply back-substitution of placeholders
    block = block.subs(substitutions)

    # Detect which optional assignments are required to close the symbolic collision
    optional = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)

    # return cse optimized expression
    return '\n'.join([
        "template <typename T, typename... FIELDS>",
        f"""struct CSE<{ collisionO }> {{""",
        "template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>",
        "CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {",
        code(CodeBlock(*optional, *block.args[:-2], *post_collision), symbols=generator.symbols),
        "return { %s, %s };" % (code(block.args[-2].rhs), code(block.args[-1].rhs)),
        "}\n};"
    ])

def optimize_operator(inputFile):
    # Load expression tree from ouput file
    with open(inputFile, "r") as file:
        load_expressions = file.read()
    results = { }
    exec(load_expressions, globals(), results)

    # Operator information
    operatorO = results['operatorO']
    descriptor = results['descriptor']
    isOperatorWithParameter = results['isOperatorWithParameter']
    # Assignments
    cell_assignments = results['cell_assignments']
    fields_assignments = results['fields_assignments']
    function_calls = results['function_calls']
    # Symbols
    params_symbols = results['params_symbols']

    operatorO = operatorO.replace("Expr", "T")
    operatorO = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', operatorO)
    descriptor = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', descriptor)

    # Combine all assignments
    assignments = cell_assignments
    optional_symbols = params_symbols

    # Current limitations of sympy assignments
    generator = iter(SymbolGenerator('x'))
    cell_backsubstitutions = [ ]
    sub_assignments = [ ]
    for assignment in assignments:
        alias = next(generator)
        cell_backsubstitutions.append((alias, assignment.lhs))
        sub_assignments.append(Assignment(alias, assignment.rhs))

    # Creating aliases and assignments for parameter calls
    optional_assignments = [ ]
    optional_substitutions = [ ]
    for symbol in optional_symbols:
        alias = next(generator)
        optional_assignments.append(Assignment(alias, symbol))
        optional_substitutions.append((symbol, alias))

    # Create code block to perform cse
    block = CodeBlock(*function_calls,*sub_assignments,*fields_assignments)
    block = block.subs(optional_substitutions)

    # Perform cse
    block = cse(block, generator)

    # These are required due to current CSE limititations
    substitutions = [ ]
    post_apply = [ ]
    simple_post_apply = True
    for assignment in block.args:
        if assignment.lhs in cell_backsubstitutions[0]:
            atoms = assignment.rhs.free_symbols
            for (lhs, rhs) in cell_backsubstitutions:
                if lhs in atoms or rhs in atoms:
                    simple_post_apply = False
    if simple_post_apply:
        substitutions = cell_backsubstitutions
    else:
        for subs in cell_backsubstitutions:
            post_apply.append(Assignment(subs[1], subs[0]))

    # Apply back-substitution of placeholders
    block = block.subs(substitutions)

    # Detect which optional assignments are required to close the symbolic apply
    optional = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)

    # Remove helper assignments to ensure symbol availability
    args = [ ]
    for assignment in block.args:
        if hasattr(assignment.rhs, 'name'):
            if "EXPOSE_ARRAY_ELEMENT" not in assignment.rhs.name:
                args.append(assignment)
        else:
            args.append(assignment)

    block = CodeBlock(*args)

    # return cse optimized expression
    header = [ ]
    if isOperatorWithParameter:
        header = "\n".join([
            "template <concepts::Cell CELL, concepts::Parameters PARAMETERS>",
            "void apply(CELL& cell, PARAMETERS& parameters) any_platform {"
        ])
    else:
        header = "\n".join([
            "template <concepts::Cell CELL>",
            "void apply(CELL& cell) any_platform {"
        ])

    return '\n'.join([
        "template <typename T, typename... FIELDS>",
        f"""struct CSE_O<{ operatorO },{ descriptor }> {{""",
        header,
        "using V = typename CELL::value_t;",
        "using DESCRIPTOR = typename CELL::descriptor_t;",
        code(CodeBlock(*optional, *block.args, *post_apply), symbols=generator.symbols),
        "}\n};"
    ])
