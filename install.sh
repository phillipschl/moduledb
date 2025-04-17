#!/bin/bash
# Installation script for NRPS C-A-T Module Database Generation Pipeline

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Miniconda or Anaconda first."
    exit 1
fi

# Create conda environment
echo "Creating conda environment..."
conda env create -f environment.yml

# Activate environment
echo "To activate the environment, run:"
echo "    conda activate moduledb"

echo "Installation complete!"
echo "Usage: python main.py --uniref90 /path/to/uniref90.fasta --config config.yaml" 