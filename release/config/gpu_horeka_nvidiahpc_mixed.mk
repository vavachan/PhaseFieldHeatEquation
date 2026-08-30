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

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3
CXXFLAGS        += -std=c++20
CXXFLAGS        += -Xcompiler -I/hkfs/home/software/all/toolkit/nvidia_hpc_sdk/22.3/Linux_x86_64/22.3/cuda/include

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

# Compiler to use for CUDA-enabled files
CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20
CUDA_LDFLAGS    := -L/hkfs/home/software/all/toolkit/nvidia_hpc_sdk/22.3/Linux_x86_64/22.3/cuda/11.6/lib64 -L/hkfs/home/software/all/toolkit/nvidia_hpc_sdk/22.3/Linux_x86_64/22.3/REDIST/cuda/11.6/targets/x86_64-linux/lib/stubs
CUDA_ARCH       := 80

USE_EMBEDDED_DEPENDENCIES := ON

# Example SLURM config:
#
# #!/usr/bin/env bash
# #SBATCH --account="hk-project-cpe"
# #SBATCH --partition="accelerated"
# #SBATCH --nodes=1
# #SBATCH --tasks-per-node=4
# #SBATCH --cpus-per-task=38
# #SBATCH --gres=gpu:4
# #SBATCH --exclusive
# #SBATCH --time=...
#
# module load toolkit/nvidia-hpc-sdk/22.3
#
# export WORKSPACE=...
#
# if [ -d tmp ]; then
#   if [ -L tmp ]; then
#     unlink tmp
#   else
#     rm -rf tmp
#   fi
# fi
# mkdir $WORKSPACE/job_${SLURM_JOB_ID}
# ln -s $WORKSPACE/job_${SLURM_JOB_ID} tmp
#
# mpirun --map-by ppr:2:socket:pe=19 --bind-to core bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./nozzle3d'
