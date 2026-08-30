# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Adapted for the Karolina supercomputer at IT4I
#
# Load modules:
#
# module load make/4.4.1-GCCcore-13.2.0
# module load CUDA/12.4.0
# module load OpenMPI/4.1.6-GCC-12.2.0-CUDA-12.4.0
# module load Boost/1.83.0-GCC-13.2.0

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20

CUDA_ARCH       := 80

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
