# Example build config for OpenLB using CUDA and OpenMPI
#
# Tested using the laptop ROG Zephyrus G14 GA401QE GA401QE which has AMD Ryzen 7 and  NVIDIA GeForce RTX 3050 in WSL2 ubuntu 2204
#
# Usage:
#  - $sudo apt install libopenmpi-dev openmpi-bin or download and install openMPI4.1 from https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html
#  - check package location $dpkg -L libopenmpi-dev | grep mpi
#  - Copy this file to OpenLB root as `config.mk`
#  - If you need, you shoud change CXXFLAGS and LDFLAGS which you check at No2 bullet.
#  - Adjust CUDA_ARCH to match your specifc GPU
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using `mpirun -np 1 ./cavity3d` (All processes share default GPU, not optimal)
#
# CXXFLAGS and LDFLAGS may need to be adjusted depending on the specific MPI installation.
# Compare to `mpicxx --showme:compile` and `mpicxx --showme:link` when in doubt.



CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

MPIFLAGS        := -lmpi_cxx -lmpi

CXXFLAGS += -I/usr/lib/x86_64-linux-gnu/openmpi/include
LDFLAGS += -L/usr/lib/x86_64-linux-gnu/openmpi/lib

PLATFORMS       := CPU_SISD GPU_CUDA

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON
