# Example build config for OpenLB using CUDA on single GPU systems
#
# Tested using CUDA 11.4
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Adjust CUDA_ARCH to match your specifc GPU
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `./cavity3d`

CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++20 --forward-unknown-to-host-compiler

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD GPU_CUDA

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
