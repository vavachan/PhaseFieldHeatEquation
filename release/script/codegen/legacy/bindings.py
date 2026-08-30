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

import cppyy

cppyy.add_include_path('src')

cppyy.cppdef('#define USING_LEGACY_CODEGEN')
cppyy.cppdef('#define PLATFORM_CPU_SISD')
cppyy.cppdef('#define DISABLE_CSE')
cppyy.cppdef('#define any_platform')

cppyy.include('core/platform/platform.h')

cppyy.include('descriptor/descriptor.h')

cppyy.include('core/expr.cpp')
cppyy.include('core/concepts.h')
cppyy.include('core/fields.h')
cppyy.include('core/data.h')
cppyy.include('core/cellD.h')
cppyy.include('core/fieldParametersD.hh')

cppyy.include('dynamics/interface.h')

cppyy.include('dynamics/equilibrium.h')
cppyy.include('dynamics/collision.h')
cppyy.include('dynamics/collisionLES.h')
cppyy.include('dynamics/collisionMRT.h')
cppyy.include('dynamics/guoZhaoDynamics.h')
cppyy.include('dynamics/freeEnergyDynamics.h')
cppyy.include('dynamics/interactionPotential.h')

from cppyy.gbl import olb

# Using FORCE field as dummy due to instantiation failure on empty variadic pack in cppyy
descriptors = {
    'D2Q5': olb.descriptors.D2Q5[olb.descriptors.FORCE,olb.descriptors.POROSITY,olb.descriptors.VELOCITY],
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.FORCE,olb.descriptors.POROSITY,olb.descriptors.VELOCITY],
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.FORCE,olb.descriptors.POROSITY,olb.descriptors.VELOCITY],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.FORCE,olb.descriptors.POROSITY,olb.descriptors.VELOCITY],
    'D3Q27': olb.descriptors.D3Q27[olb.descriptors.FORCE,olb.descriptors.POROSITY,olb.descriptors.VELOCITY],
}

mrtDescriptors = {
    'D2Q5': olb.descriptors.D2Q5[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.tag.MRT,olb.descriptors.FORCE],
}

rtlbmDescriptors = {
    'D3Q7': olb.descriptors.D3Q7[olb.descriptors.tag.RTLBM],
    'D3Q15': olb.descriptors.D3Q15[olb.descriptors.tag.RTLBM],
    'D3Q27': olb.descriptors.D3Q27[olb.descriptors.tag.RTLBM],
}

freeEnergyDescriptors = {
    'D2Q9': olb.descriptors.D2Q9[olb.descriptors.CHEM_POTENTIAL,olb.descriptors.FORCE],
    'D3Q19': olb.descriptors.D3Q19[olb.descriptors.CHEM_POTENTIAL,olb.descriptors.FORCE],
}
