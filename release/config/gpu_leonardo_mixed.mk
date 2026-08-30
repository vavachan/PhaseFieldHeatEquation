# Example build config for OpenLB using mixed compilation of CUDA with OpenMPI
#
# Adapted for the Leonardo supercomputer at CINECA (booster module)
#
# Load modules:
#
# module load gcc/12.2.0
# module load openmpi/4.1.6--gcc--12.2.0
# module load cuda/12.3

CXX             := mpic++
CC              := gcc

CXXFLAGS        := -O3 -Wall -march=native -mtune=native
CXXFLAGS        += -std=c++20

PARALLEL_MODE   := MPI

PLATFORMS       := CPU_SISD GPU_CUDA

CUDA_CXX        := nvcc
CUDA_CXXFLAGS   := -O3 -std=c++20
CUDA_LDFLAGS    := -L/leonardo/prod/opt/compilers/cuda/12.3/none/targets/x86_64-linux/lib/stubs
CUDA_LDFLAGS    += -L/leonardo/prod/opt/compilers/cuda/12.3/none/lib64

CUDA_ARCH       := 80

FLOATING_POINT_TYPE := float

USE_EMBEDDED_DEPENDENCIES := ON

# Example SLURM config:
#
# #!/bin/bash
#
# #SBATCH -J olb-test
# #SBATCH -A TBD
# #SBATCH -p boost_usr_prod
# #SBATCH --qos boost_qos_lprod
# #SBATCH --time 00:30:00
# #SBATCH -o test.log
# #SBATCH --nodes 1
# #SBATCH --ntasks-per-node 4
# #SBATCH --cpus-per-task 8
# #SBATCH --gpus 4
#
# module load gcc/12.2.0
# module load openmpi/4.1.6--gcc--12.2.0
# module load cuda/12.3
#
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
#
# mpirun --map-by ppr:4:node:pe=8 --bind-to core bash -c 'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}; ./cavity3d --steps 1000 --size 500 --no-results'
