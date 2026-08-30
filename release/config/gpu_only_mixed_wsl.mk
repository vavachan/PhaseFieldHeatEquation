# Example build config for OpenLB using CUDA and OpenMPI
#
# Tested using the laptop ROG Zephyrus G14 GA401QE GA401QE which has AMD Ryzen 7 and  NVIDIA GeForce RTX 3050 in WSL2 ubuntu 2204
#
# Usage:
#  - $sudo apt install libopenmpi-dev openmpi-bin or download and install openMPI4.1 from https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html
#  - Copy this file to OpenLB root as `config.mk`
#  - check library path $sudo find /usr -name libcuda.so (for defining CUDA_LIB_PATH and if you use wsl, set WSL_LIB_PATH)
#  - check library path $sudo find /usr -name libmpi.so (for defining MPI_PATH)
#  - check package location $dpkg -L libopenmpi-dev | grep mpi
#  - Adjust `CUDA_PATH` to match your CUDA installation
#  - Adjust CUDA_ARCH to match your specifc GPU architecture
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3d`
#  - Run `make`
#  - Start the simulation using `mpirun -np 1 ./cavity3d`

# One advantage of mixed compilation is that `make no-cuda-recompile` is available,
# significantly speeding up the compilation if the set of used dynamics / post
# processors has not been expanded (i.e. if no new GPU kernels need costly
# instantiation).

CUDA_PATH=/usr/local/cuda
MPI_PATH=/usr/lib/x86_64-linux-gnu/openmpi
CUDA_LIB_PATH=/usr/local/cuda/targets/x86_64-linux/lib
WSL_LIB_PATH=/usr/lib/wsl/lib

# Compiler settings
CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3  -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20  -I${CUDA_PATH}/include

PARALLEL_MODE   := MPI
PLATFORMS       := CPU_SISD GPU_CUDA

# CUDA settings
CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20  -I${MPI_PATH}/include -I${CUDA_PATH}/include

CUDA_LDFLAGS    := -L${CUDA_LIB_PATH}  -lcuda -lcudart -lcudadevrt \
                   -L${WSL_LIB_PATH} \
                   -L${MPI_PATH}/lib -lmpi -lmpi_cxx \
                   -Wl,-rpath=${WSL_LIB_PATH} \
                   -Wl,-rpath=${CUDA_LIB_PATH} \
                   -Wl,-rpath=.

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 86

FLOATING_POINT_TYPE := double

USE_EMBEDDED_DEPENDENCIES := ON
