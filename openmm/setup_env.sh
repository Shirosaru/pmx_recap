#!/bin/bash
set -e

# setup_env.sh
# Usage: ./setup_env.sh
# Installs necessary dependencies for OpenMM Antibody Binding Analysis

echo "Checking for Conda..."
if command -v conda &> /dev/null; then
    echo "Conda found. Creating/Updating environment 'openmm_binding'..."
    conda create -n openmm_binding -c conda-forge openmm pdbfixer parmed mdtraj python=3.9 -y
    echo "Environment 'openmm_binding' created."
    echo "Activate with: conda activate openmm_binding"
else
    echo "Conda not found. Attempting pip install (Ensure Python 3.9+ is active)..."
    if ! command -v python3 &> /dev/null; then
        echo "Error: Python3 not found."
        exit 1
    fi
    
    echo "Installing requirements via pip..."
    pip install openmm pdbfixer parmed mdtraj numpy
    # Note: openmm via pip might be tricky on some linux distros without Conda, 
    # but official wheels exist for valid python versions.
    echo "Installation complete."
fi

echo "Setup Done. You can now run 'python compare_binding.py'."
