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

from mako.template import Template
from pathlib import Path
import subprocess, os

def extractDynamicsExpressions(dynamics_type, template, input_f, output_p, output_o):
    # Fill Mako template of C++ program to retrieve expression tree
    template_extraction = Template(filename=template)
   
    # create subdirectories to prevent confusion by linker
    executable = Path(input_f).stem;

    # specify output for the cpp mako template
    output_cpp = output_p + Path(input_f).stem + ".cpp"
    output_template = output_o + Path(input_f).stem + ".out"

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(template_extraction.render(dynamics=dynamics_type, output=output_template))
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+";", shell=True, cwd=output_p)
    subprocess.call("./"+executable, shell=True, cwd=output_p)

def extractOperatorExpressions(operator_with_descriptor_type, template, input_f, output_p, output_o):
    # Fill Mako template of C++ program to retrieve expression tree
    template_extraction = Template(filename=template)
   
    # create subdirectories to prevent confusion by linker
    executable = Path(input_f).stem;

    # specify output for the cpp mako template
    output_cpp = output_p + Path(input_f).stem + ".cpp"
    output_template = output_o + Path(input_f).stem + ".out"

    # split operator and descriptor
    operator_type, descriptor_type = operator_with_descriptor_type.split(';',1)

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(template_extraction.render(operator=operator_type, descriptor=descriptor_type, output=output_template))
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+";", shell=True, cwd=output_p)
    subprocess.call("./"+executable, shell=True, cwd=output_p)

def executeCSE(expressions, template, result):
    # Generate cse-optimized C++ code from tree and fill template
    template_cse = Template(filename=template)
    with open(result, 'w') as f:
        f.write(template_cse.render(filename=expressions))


