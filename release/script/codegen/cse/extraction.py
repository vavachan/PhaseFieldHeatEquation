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

from argparse import ArgumentParser
from source.cse_modules import extractDynamicsExpressions, extractOperatorExpressions

# Parse name of dynamics to be optimized
parser = ArgumentParser()
parser.add_argument('target', help='full C++ typename of the target instance')
parser.add_argument('template', help='path and filename to mako template')
parser.add_argument('input', help='path and filename to .txt file containing target info')
parser.add_argument('cpp_output', help='path where the .cpp files are created')
parser.add_argument('out_output', help='path where the .out files are created')
parser.add_argument('type', help='dynamics or operator')
args = parser.parse_args()

# Create .cpp files to extract expression tree
if args.type == 'dynamics':
    extractDynamicsExpressions(args.target, args.template, args.input, args.cpp_output, args.out_output)
elif args.type == 'operator':
    extractOperatorExpressions(args.target, args.template, args.input, args.cpp_output, args.out_output)
else:
    print("Invalid type.")
