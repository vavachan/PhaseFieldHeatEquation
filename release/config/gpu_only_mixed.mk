# Example build config for OpenLB using CUDA on single GPU systems with separate compiler for non-CUDA parts
#
# Tested using CUDA 11.4
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Adjust CUDA_ARCH to match your specifc GPU
#  - Adjust CUDA_LDFLAGS s.t. CXX can link CUDA libraries
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `./cavity3d`

CXX             := g++
CC              := gcc

CXXFLAGS        := -O3 -std=c++20

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_CXX        := nvcc
# Adjust to CUDA / driver library path so that CXX can resolve libcuda, libcudart and libcudadevrt
CUDA_LDFLAGS    := -L/run/opengl-driver/lib
# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

FLOATING_POINT_TYPE        := double

USE_EMBEDDED_DEPENDENCIES := ON
