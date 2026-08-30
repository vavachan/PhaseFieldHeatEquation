# Example build config for OpenLB using GNU C++ and OpenMPI on Apple MacOS
#
# Tested on Apple M1
#
# Also see the guide available at [1].
#
# [1]: https://www.openlb.net/wp-content/uploads/2022/06/olb-tr6.pdf
#
# Usage:
#  - Install GCC, OpenMPI and tinyxml2 using homebrew
#  - Copy this file to OpenLB root as `config.mk`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun ./cavity3d`

CXX := mpic++
CC := gcc-12

CXXFLAGS := -O3 -Wall -march=native -mtune=native
CXXFLAGS += -std=c++20
CXXFLAGS += -I/opt/homebrew/include

LDFLAGS := -L/opt/homebrew/lib

PARALLEL_MODE := MPI
OMPFLAGS := -fopenmp

PLATFORMS := CPU_SISD

USE_EMBEDDED_DEPENDENCIES := OFF
