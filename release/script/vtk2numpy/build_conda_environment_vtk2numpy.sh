#!/usr/bin/env bash

# Script to install a conda environment whereby vtk2numpy can work.

CONDA_ENV_NAME=env-vtk2numpy

echo "---> Creating the conda environment..."
conda create -n $CONDA_ENV_NAME -y

echo "---> Activating the conda environment..."
source activate $CONDA_ENV_NAME

echo "---> Installing a newer version of python and pip (not necessary but may be safer to avoid inconsistencies)..."
conda install python -y
conda install pip -y

echo "---> Restarting the conda environment..."
conda deactivate
source activate $CONDA_ENV_NAME

echo "---> Installing all the python packages..."
python3 -m pip install matplotlib numpy vtk

echo "---> All done!"
