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
from source.cse_modules import executeCSE

# Parse name of dynamics to be optimized
parser = ArgumentParser()
parser.add_argument('inputfile', help='filename containing expressions')
parser.add_argument('cse_template', help='filename of mako template to execute cse')
parser.add_argument('outputfile',   help='path and filename of output file')
args = parser.parse_args()

# Generate cse-optimized C++ code from tree and fill template
print('Generate from ' + args.inputfile)
executeCSE(args.inputfile, args.cse_template, args.outputfile)
