# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Adapted for Nvidia HPC SDK 22.3 on the HoreKa supercomputer at KIT
#
# Usage:
#  - Copy this file to OpenLB root as `config.mk`
#  - Load Nvidia HPC toolkit via `module load toolkit/nvidia-hpc-sdk/22.3`
#  - Run `make clean; make`
#  - Switch to example directory, e.g. `examples/turbulence/nozzle3d`
#  - Run `make`
#  - Use `mpirun --map-by ppr:2:socket:pe=19 --bind-to core bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./nozzle3d' for launch
#  - See example SLURM config at the end


CXX             := nvcc
CC              := nvcc

CXXFLAGS        := -O3  --compiler-bindir=/opt/gcc/12/bin -std=c++20
CXXFLAGS        += -I/hkfs/home/software/all/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/comm_libs/openmpi/openmpi-3.1.5/include


LDFLAGS         := -L/hkfs/home/software/all/toolkit/nvidia_hpc_sdk/23.9/Linux_x86_64/23.9/comm_libs/openmpi/openmpi-3.1.5/lib \
		 --compiler-bindir=/opt/gcc/12/bin -std=c++20

PARALLEL_MODE   := MPI

MPIFLAGS        := -lmpi_cxx -lmpi

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_ARCH       := 90

FLOATING_POINT_TYPE := double

USE_EMBEDDED_DEPENDENCIES := ON
# Example SLURM config:
#
# #!/bin/bash
# #SBATCH --account="hk-project-cpe"
# #SBATCH --partition="accelerated"
# #SBATCH --nodes=1
# #SBATCH --tasks-per-node=4
# #SBATCH --cpus-per-task=38
# #SBATCH --gres=gpu:4
# #SBATCH --exclusive
# #SBATCH --time=...
#
# module purge
# module load toolkit/nvidia-hpc-sdk/23.9
# export LD_LIBRARY_PATH=/opt/gcc/12/lib64/:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.


