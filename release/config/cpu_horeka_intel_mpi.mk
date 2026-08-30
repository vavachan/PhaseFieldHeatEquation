# Example build config for OpenLB using Intel C++ and MPI on HoreKa

CXX             := mpicxx
CC              := icx

CXXFLAGS        := -O3 -Wall
CXXFLAGS        += -std=c++20 --gcc-toolchain=/opt/gcc/13/

LDFLAGS += --gcc-toolchain=/opt/gcc/13/ 

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD

USE_EMBEDDED_DEPENDENCIES := ON
