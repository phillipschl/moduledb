#!/bin/bash

# Set environment variables
export PFAM_DB=/media/data/moduledb/Pfam-A.hmm

# Load the conda environment
eval "$(conda shell.bash hook)"
conda activate moduledb

# Run the pipeline with pfam path
python main.py --uniref90 /media/data/moduledb/uniref90.fasta --log_level DEBUG

# Exit with the pipeline's exit code
exit $? 