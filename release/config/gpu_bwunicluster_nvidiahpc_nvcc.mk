CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3  --compiler-bindir=/opt/gcc/12/bin -std=c++20
CXXFLAGS        += -I/pfs/data5/software_uc2/bwhpc/common/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/comm_libs/openmpi/openmpi-3.1.5/include


LDFLAGS         := -L/pfs/data5/software_uc2/bwhpc/common/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/comm_libs/openmpi/openmpi-3.1.5/lib \
                 --compiler-bindir=/opt/gcc/12/bin -std=c++20

PARALLEL_MODE   := MPI

MPIFLAGS        := -lmpi_cxx -lmpi

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_ARCH       := 70

FLOATING_POINT_TYPE := double

USE_EMBEDDED_DEPENDENCIES := ON


