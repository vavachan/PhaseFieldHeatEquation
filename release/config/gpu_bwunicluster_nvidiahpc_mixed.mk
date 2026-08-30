# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Tested using module toolkit/nvidia-hpc-sdk/23.9 in bwUniCluster2.0 in June 2024
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Adjust CUDA_ARCH to match your specifc GPU
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/laminar/cavity3dBenchmark`
#  - Run `make`
#  - Start the simulation using sbatch sbatch testslum.sh
#
#
# One advantage of mixed compilation is that `make no-cuda-recompile` is available,
# significantly speeding up the compilation if the set of used dynamics / post
# processors has not been expanded (i.e. if no new GPU kernels need costly
# instantiation).

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20 --gcc-toolchain=/opt/gcc/12
CXXFLAGS        += -Xcompiler -I/pfs/data5/software_uc2/bwhpc/common/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/cuda/12.2/targets/x86_64-linux/include
LDFLAGS         += --gcc-toolchain=/opt/gcc/12

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

# Compiler to use for CUDA-enabled files
CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20 --compiler-bindir=/opt/gcc/12/bin

# Adjust to enable resolution of libcuda, libcudart, libcudadevrt
CUDA_LDFLAGS   := -L/pfs/data5/software_uc2/bwhpc/common/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/cuda/12.2/targets/x86_64-linux/lib \
		  -L/pfs/data5/software_uc2/bwhpc/common/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/REDIST/cuda/12.2/targets/x86_64-linux/lib/stubs \
		 -lstdc++

# for e.g. RTX 30* (Ampere), see table in `rules.mk` for other options
CUDA_ARCH       := 70

FLOATING_POINT_TYPE := double

USE_EMBEDDED_DEPENDENCIES := ON
